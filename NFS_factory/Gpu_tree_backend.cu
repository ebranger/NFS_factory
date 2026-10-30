//The tree computations on the GPU (TREE_BACKEND_GPU, in the GPU build only). Each depth of the product and remainder trees is
//done with batched kernels for all its nodes: the NTT products of int_ntt with one plan per depth (tree_level_plan), or schoolbook
//products for small nodes, exactly as Ntt_tree_backend does on the CPU, so that every node gets the same number. The forward
//transforms of the nodes from the product tree are kept on the GPU for the middle products of the remainder tree.
//The root and the leaves are done on the CPU, by the code of Cpu_tree_backend: the root with its large products on the GPU, and
//the leaves from the fractions of their parents, which are copied back.

#include "Tree_backend.h"
#include "ntt_gpu.h"
#include <stdio.h>
#include <stdlib.h>

//The products of the root on the GPU, one at a time. Small products are left to MPIR.
class Gpu_tree_products : public Tree_products
{
public:
	Gpu_tree_products()
	{
		in_limbs = 0;
		out_limbs = 0;
		transform_words = 0;
		limbs[0] = limbs[1] = NULL;
		sizes = NULL;
		X = Y = NULL;
		CUDA_CHECK(cudaMalloc(&sizes, 4 * sizeof(uint32_t)));
	}

	~Gpu_tree_products()
	{
		cudaFree(limbs[0]);
		cudaFree(limbs[1]);
		cudaFree(sizes);
		cudaFree(X);
		cudaFree(Y);
	}

	void Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b)
	{
		mp_size_t a_limbs = abs(a->_mp_size), b_limbs = abs(b->_mp_size);
		Ntt_plan plan;
		if (a_limbs + b_limbs < 512 || a_limbs < 2 || b_limbs < 2 || !ntt_plan_product(plan, mpz_sizeinbase(a, 2), mpz_sizeinbase(b, 2)))
		{
			mpz_mul(result, a, b);
			return;
		}
		const bool negative = (a->_mp_size < 0) != (b->_mp_size < 0);
		Product(result, a, b, plan, false);
		if (negative)
		{
			mpz_neg(result, result);
		}
	}

	mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L)
	{
		Ntt_plan plan;
		if (!ntt_plan_cyclic(plan, min_L, mpz_sizeinbase(y, 2), mpz_sizeinbase(x, 2)))
		{
			return cpu.Product_mod(result, y, x, min_L);
		}
		if (mpz_sgn(y) == 0 || mpz_sgn(x) == 0)
		{
			mpz_set_ui(result, 0);
			return plan.L();
		}
		Product(result, y, x, plan, true);
		return plan.L();
	}

	void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision)
	{
		cpu.Middle_product_pair(result0, result1, y, x0, x1, shift0, shift1, precision); //Not used by the GPU tree.
	}

	Gpu_ntt gpu;

private:
	//result = a * b, exactly or modulo 2^L - 1, with the plan.
	void Product(mpz_ptr result, mpz_srcptr a, mpz_srcptr b, const Ntt_plan& plan, bool cyclic)
	{
		size_t a_limbs = mpz_size(a), b_limbs = mpz_size(b);
		size_t stride = (a_limbs > b_limbs) ? a_limbs : b_limbs;
		size_t out_stride = Gpu_ntt::Combine_limbs(plan);
		Reserve(2 * stride, out_stride, (size_t)plan.primes * 2 * plan.Length());
		CUDA_CHECK(cudaMemcpy(limbs[0], a->_mp_d, a_limbs * sizeof(uint64_t), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(limbs[0] + stride, b->_mp_d, b_limbs * sizeof(uint64_t), cudaMemcpyHostToDevice));
		uint32_t size[2] = { (uint32_t)a_limbs, (uint32_t)b_limbs };
		CUDA_CHECK(cudaMemcpy(sizes, size, sizeof(size), cudaMemcpyHostToDevice));
		Gpu_level in = { limbs[0], stride, sizes, 2 };
		Gpu_level out = { limbs[1], out_stride, sizes + 2, 1 };
		gpu.Forward(plan, in, X, 0);
		gpu.Pointwise_adjacent(plan, X, Y, 1, 0);
		gpu.Combine(plan, Y, 1, cyclic, out, 0);
		uint32_t n;
		CUDA_CHECK(cudaMemcpy(&n, sizes + 2, sizeof(uint32_t), cudaMemcpyDeviceToHost));
		if (result->_mp_alloc < (int)n)
		{
			mpz_realloc2(result, n * GMP_LIMB_BITS);
		}
		CUDA_CHECK(cudaMemcpy(result->_mp_d, limbs[1], n * sizeof(uint64_t), cudaMemcpyDeviceToHost));
		result->_mp_size = (int)n;
	}

	void Reserve(size_t in, size_t out, size_t words)
	{
		if (in > in_limbs)
		{
			cudaFree(limbs[0]);
			CUDA_CHECK(cudaMalloc(&limbs[0], in * sizeof(uint64_t)));
			in_limbs = in;
		}
		if (out > out_limbs)
		{
			cudaFree(limbs[1]);
			CUDA_CHECK(cudaMalloc(&limbs[1], out * sizeof(uint64_t)));
			out_limbs = out;
		}
		if (words > transform_words)
		{
			cudaFree(X);
			cudaFree(Y);
			CUDA_CHECK(cudaMalloc(&X, words * sizeof(uint32_t)));
			CUDA_CHECK(cudaMalloc(&Y, words / 2 * sizeof(uint32_t)));
			transform_words = words;
		}
	}

	Ntt_product cpu; //for products beyond the longest transform
	uint64_t* limbs[2];
	uint32_t* sizes;
	uint32_t *X, *Y;
	size_t in_limbs, out_limbs, transform_words;
};

//Device memory that grows as needed and is kept.
template <class T> static void reserve(T*& p, size_t& capacity, size_t count)
{
	if (count > capacity)
	{
		cudaFree(p);
		CUDA_CHECK(cudaMalloc(&p, count * sizeof(T)));
		capacity = count;
	}
}

class Gpu_tree_backend : public Cpu_tree_backend
{
public:
	Gpu_tree_backend(int levels, int max_bits_in_input);
	~Gpu_tree_backend();
	mp_bitcnt_t Build_product_tree();
	void Remainder_tree();

private:
	//The nodes of a depth: the product tree, and then the fractions of the remainder tree.
	struct Depth
	{
		Gpu_level level;
		size_t limb_capacity, size_capacity;
		uint32_t* bits; //of the product tree nodes
		size_t bits_capacity;
		uint32_t* transform; //of the nodes, kept from the product tree where the depth uses the NTT
		size_t transform_capacity;
		Ntt_plan plan;
		bool uses_ntt;
	};

	void Set_stride(Depth& d, size_t count, size_t stride);
	void Upload(int depth);
	void Download(int depth);
	void Download_bits(int depth);
	void Shifts(int depth);

	Gpu_ntt* gpu;
	std::vector<Depth> depth_data;
	uint32_t* parent_transform;
	size_t parent_transform_capacity;
	uint32_t* pair_transform;
	size_t pair_transform_capacity;
	Gpu_level scratch_level; //cyclic products, and exact truncated products
	size_t scratch_limb_capacity, scratch_size_capacity;
	uint32_t *shift, *precision;
	size_t shift_capacity, precision_capacity;
	std::vector<uint64_t> host_limbs;
	std::vector<uint32_t> host_sizes, host_words;
};

Gpu_tree_backend::Gpu_tree_backend(int levels, int max_bits_in_input) : Cpu_tree_backend(new Gpu_tree_products, levels, max_bits_in_input)
{
	gpu = &static_cast<Gpu_tree_products*>(products)->gpu;
	depth_data.resize(levels);
	for (int d = 0; d < levels; d++)
	{
		Depth& x = depth_data[d];
		x.level.limbs = NULL;
		x.level.size = NULL;
		x.level.stride = 0;
		x.level.count = (size_t)1 << d;
		x.limb_capacity = x.size_capacity = x.bits_capacity = x.transform_capacity = 0;
		x.bits = NULL;
		x.transform = NULL;
		x.uses_ntt = false;
		reserve(x.level.size, x.size_capacity, x.level.count);
		reserve(x.bits, x.bits_capacity, x.level.count);
	}
	parent_transform = pair_transform = NULL;
	parent_transform_capacity = pair_transform_capacity = 0;
	scratch_level.limbs = NULL;
	scratch_level.size = NULL;
	scratch_limb_capacity = scratch_size_capacity = 0;
	shift = precision = NULL;
	shift_capacity = precision_capacity = 0;
}

Gpu_tree_backend::~Gpu_tree_backend()
{
	for (Depth& x : depth_data)
	{
		cudaFree(x.level.limbs);
		cudaFree(x.level.size);
		cudaFree(x.bits);
		cudaFree(x.transform);
	}
	cudaFree(parent_transform);
	cudaFree(pair_transform);
	cudaFree(scratch_level.limbs);
	cudaFree(scratch_level.size);
	cudaFree(shift);
	cudaFree(precision);
}

void Gpu_tree_backend::Set_stride(Depth& d, size_t count, size_t stride)
{
	d.level.stride = stride;
	d.level.count = count;
	reserve(d.level.limbs, d.limb_capacity, count * stride);
}

//Copy the nodes of a depth from the CPU to the GPU.
void Gpu_tree_backend::Upload(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	const size_t count = d.level.count, stride = d.level.stride;
	host_limbs.assign(count * stride, 0);
	host_sizes.resize(count);
	for (size_t r = 0; r < count; r++)
	{
		mpz_srcptr a = node[first + r].get_mpz_t();
		size_t n = mpz_size(a);
		if (n > stride)
		{
			fprintf(stderr, "Gpu_tree_backend: a node of %zu limbs does not fit in %zu\n", n, stride);
			exit(3);
		}
		memcpy(host_limbs.data() + r * stride, a->_mp_d, n * sizeof(uint64_t));
		host_sizes[r] = (uint32_t)n;
	}
	CUDA_CHECK(cudaMemcpy(d.level.limbs, host_limbs.data(), host_limbs.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d.level.size, host_sizes.data(), count * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

//Copy the numbers of a depth from the GPU to the CPU.
void Gpu_tree_backend::Download(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	const size_t count = d.level.count, stride = d.level.stride;
	host_limbs.resize(count * stride);
	host_sizes.resize(count);
	CUDA_CHECK(cudaMemcpy(host_limbs.data(), d.level.limbs, host_limbs.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(host_sizes.data(), d.level.size, count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	for (size_t r = 0; r < count; r++)
	{
		mpz_ptr a = node[first + r].get_mpz_t();
		size_t n = host_sizes[r];
		if (a->_mp_alloc < (int)n)
		{
			mpz_realloc2(a, n * GMP_LIMB_BITS);
		}
		memcpy(a->_mp_d, host_limbs.data() + r * stride, n * sizeof(uint64_t));
		a->_mp_size = (int)n;
	}
}

void Gpu_tree_backend::Download_bits(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	host_words.resize(d.level.count);
	gpu->Bits(d.level, d.bits, 0);
	CUDA_CHECK(cudaMemcpy(host_words.data(), d.bits, d.level.count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	for (size_t r = 0; r < d.level.count; r++)
	{
		node_bits[first + r] = host_words[r];
	}
}

mp_bitcnt_t Gpu_tree_backend::Build_product_tree()
{
	//The leaves, with room for their fractions too.
	const int leaf_depth = levels - 1;
	size_t max_leaf = 1;
	for (int i = first_leaf; i < first_leaf + num_leaves; i++)
	{
		node_bits[i] = mpz_sizeinbase(node[i].get_mpz_t(), 2);
		size_t n = mpz_size(node[i].get_mpz_t());
		max_leaf = (n > max_leaf) ? n : max_leaf;
	}
	Set_stride(depth_data[leaf_depth], num_leaves, max_leaf + 1);
	Upload(leaf_depth);

	for (int depth = levels - 2; depth >= 0; depth--)
	{
		Depth& children = depth_data[depth + 1];
		Depth& parents = depth_data[depth];
		const size_t count = parents.level.count;
		children.uses_ntt = tree_level_plan(children.plan, &node_bits[(1 << (depth + 1)) - 1], children.level.count, guard);
		if (children.uses_ntt)
		{
			const Ntt_plan& p = children.plan;
			const size_t T = p.Length();
			//The transforms of the children are kept for the remainder tree, except those of the leaves.
			uint32_t* X;
			if (depth + 1 <= levels - 2)
			{
				reserve(children.transform, children.transform_capacity, (size_t)p.primes * 2 * count * T);
				X = children.transform;
			}
			else
			{
				reserve(pair_transform, pair_transform_capacity, (size_t)p.primes * 2 * count * T);
				X = pair_transform;
			}
			reserve(parent_transform, parent_transform_capacity, (size_t)p.primes * count * T);
			Set_stride(parents, count, Gpu_ntt::Combine_limbs(p));
			gpu->Forward(p, children.level, X, 0);
			gpu->Pointwise_adjacent(p, X, parent_transform, count, 0);
			gpu->Combine(p, parent_transform, count, false, parents.level, 0);
		}
		else
		{
			Set_stride(parents, count, 2 * children.level.stride);
			gpu->School_products(children.level, parents.level, 0);
		}
		Download_bits(depth);
	}
	Download(0);
	return node_bits[0];
}

//The shifts and precisions of the middle products of the nodes of a depth: precision = parent bits + guard, and
//shift = parent bits - node bits.
void Gpu_tree_backend::Shifts(int depth)
{
	const int first = (1 << depth) - 1;
	const size_t count = (size_t)1 << depth;
	host_words.resize(2 * count);
	for (size_t r = 0; r < count; r++)
	{
		const int i = first + (int)r, parent = (i - 1) / 2;
		host_words[r] = (uint32_t)(node_bits[parent] - node_bits[i]);
		host_words[count + r] = (uint32_t)(node_bits[parent] + guard);
	}
	reserve(shift, shift_capacity, count);
	reserve(precision, precision_capacity, count);
	CUDA_CHECK(cudaMemcpy(shift, host_words.data(), count * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(precision, host_words.data() + count, count * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Gpu_tree_backend::Remainder_tree()
{
	if (levels < 3)
	{
		return; //The root is the parent of the leaves, and is on the CPU.
	}
	Upload(0);
	for (int depth = 1; depth <= levels - 2; depth++)
	{
		Depth& nodes = depth_data[depth];
		Depth& parents = depth_data[depth - 1];
		const size_t count = nodes.level.count;
		Shifts(depth);
		if (nodes.uses_ntt)
		{
			const Ntt_plan& p = nodes.plan;
			const size_t T = p.Length();
			reserve(parent_transform, parent_transform_capacity, (size_t)p.primes * (count / 2) * T);
			reserve(pair_transform, pair_transform_capacity, (size_t)p.primes * count * T);
			const size_t combine_limbs = Gpu_ntt::Combine_limbs(p);
			reserve(scratch_level.limbs, scratch_limb_capacity, count * combine_limbs);
			reserve(scratch_level.size, scratch_size_capacity, count);
			scratch_level.stride = combine_limbs;
			scratch_level.count = count;
			gpu->Forward(p, parents.level, parent_transform, 0);
			gpu->Pointwise_swap(p, parent_transform, nodes.transform, pair_transform, count / 2, 0);
			gpu->Combine(p, pair_transform, count, true, scratch_level, 0);
			gpu->Extract_middle(scratch_level, shift, precision, parents.level.size, 2, nodes.level, 0);
		}
		else
		{
			//The exact truncated products need the nodes themselves, so they are written to the scratch level, which then
			//takes the place of the nodes.
			reserve(scratch_level.limbs, scratch_limb_capacity, count * nodes.level.stride);
			reserve(scratch_level.size, scratch_size_capacity, count);
			scratch_level.stride = nodes.level.stride;
			scratch_level.count = count;
			gpu->School_middle(parents.level, nodes.level, shift, precision, scratch_level, 0);
			std::swap(nodes.level.limbs, scratch_level.limbs);
			std::swap(nodes.level.size, scratch_level.size);
			std::swap(nodes.limb_capacity, scratch_limb_capacity);
			std::swap(nodes.size_capacity, scratch_size_capacity);
		}
	}
	Download(levels - 2);
}

Tree_backend* create_gpu_tree_backend(int levels, int max_bits_in_input)
{
	return new Gpu_tree_backend(levels, max_bits_in_input);
}
