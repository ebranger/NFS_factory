//a class for big-int polynomials and some operations using such polynomials.
//Originally intended for a lattice-siever, so most of this functionality is not needed in NFS_factory.

#include "Polynomial.h"
#include <iostream>

using namespace std;

#define POLYNOMIAL_MAX_DEGREE 10

uint8 Polynomial::get_degree()
{
	return degree;
}

Polynomial::Polynomial()
{
	//create a new polynomial class, initialize coefficients and set to 0.
	degree = 0;
	coefficients = new mpz_class*[2*POLYNOMIAL_MAX_DEGREE + 1];

	for (int i = 0; i <= 2*POLYNOMIAL_MAX_DEGREE; i++)
	{
		coefficients[i] = new mpz_class;
		*coefficients[i] = 0;
	}
}

Polynomial::Polynomial(Polynomial* input)
{
	//Copy the degree and all coefficients
	degree = input->degree;
	coefficients = new mpz_class*[2*POLYNOMIAL_MAX_DEGREE + 1];

	for (int i = 0; i <= 2 * POLYNOMIAL_MAX_DEGREE; i++)
	{
		coefficients[i] = new mpz_class;
		*coefficients[i] = 0;
	}

	for(int i = 0; i <= degree; i++)
	{
		*coefficients[i] = *input->coefficients[i];
	}
}

Polynomial::~Polynomial(void)
{
	for(int i = 0; i <= 2*POLYNOMIAL_MAX_DEGREE; i++)
	{
		delete coefficients[i];
	}
	delete[] coefficients;
}

void Polynomial::copy(Polynomial* input)
{
	//copy an entire polynoial.
	degree = input->degree;

	for (int i = 2*POLYNOMIAL_MAX_DEGREE; i > degree; i--)
	{
		*coefficients[i] = 0;
	}

	for(int i = 0; i <= degree; i++)
	{
		*coefficients[i] = *input->coefficients[i];
	}
}

void Polynomial::evaluate_polynomial(mpz_class* result, mpz_class* value)
{
	//Evaluate P(x) using Horners method
	*result = 0;

	for(int i = degree; i >=0; i--) 
	{
		*result = *result * *value + *coefficients[i];	
	}
}

void Polynomial::evaluate_polynomial(mpz_class* result, uint32 value)
{
	//Evaluate P(x) using Horners method
	*result = 0;

	for(int i = degree; i >=0; i--) 
	{
		*result = *result * value + *coefficients[i];	
	}
}

void Polynomial::set_coeff(uint8 coeff, mpz_class* value)
{
	//Set one coefficient in the polynomial
	if (coeff > POLYNOMIAL_MAX_DEGREE)
	{
		return;
	}

	if (coeff > degree)
	{
		degree = coeff;
	}

	*coefficients[coeff] = *value;

	if(*coefficients[coeff] == 0 && coeff == degree)
	{
		update_degree();
	}
}

void Polynomial::set_coeff(uint8 coeff, int32 value)
{
	//Set one coefficient in the polynomial
	if (coeff > POLYNOMIAL_MAX_DEGREE)
	{
		return;
	}

	if (coeff > degree)
	{
		degree = coeff;
	}

	*coefficients[coeff] = value;

	if(*coefficients[coeff] == 0 && coeff == degree)
	{
		update_degree();
	}
}

void Polynomial::clear()
{
	//Set all coefficients to 0.
	for (int i = 0; i <= 2*POLYNOMIAL_MAX_DEGREE; i++)
	{
		*coefficients[i] = 0;
	}
	degree = 0;
}

void Polynomial::evaluate_polynomial_modulo(mpz_class* result, mpz_class* value, mpz_class* module)
{
	//evaluate P(x) mod module using Horners method
	*result = 0;

	for(int i = degree; i >=0; i--) 
	{
		*result = (*result * *value + *coefficients[i]) % *module;	
	}
}

void Polynomial::evaluate_polynomial_modulo(mpz_class* result, uint32 value, uint32 module)
{
	//evaluate P(x) mod module using Horners method
	*result = 0;

	for(int i = degree; i >=0; i--) 
	{
		*result = (*result * value + *coefficients[i]) % module;	
	}
}

void Polynomial::evaluate_polynomial_homogenous(mpz_class* result, int32 a, int32 b)
{
	//Evaluate b^deg*f((a/b)^deg). 

	if (degree == 0)
	{
		//should not happen, but better check anyways...
		*result = 0;
		return;
	}

	*result = a* *coefficients[degree] + b* *coefficients[degree-1];

	mpz_class bpower;
	bpower = b;

	for (int i = degree - 2; i >= 0; i--)
	{
		bpower = bpower * b;
		*result = *result * a + *coefficients[i] * bpower;
	}

}

void Polynomial::evaluate_polynomial_homogenous(mpz_class* result, int64 a, int64 b)
{
	//Evaluate b^deg*f((a/b)^deg). 

	if (degree == 0)
	{
		//should not happen, but better check anyways...
		*result = 0;
		return;
	}

	*result = a* *coefficients[degree] + b* *coefficients[degree - 1];

	mpz_class bpower;
	bpower = b;

	for (int i = degree - 2; i >= 0; i--)
	{
		bpower = bpower * b;
			*result = *result * a + *coefficients[i] * bpower;

	}
}

bool Polynomial::poly_is_binomial_SNFS_deg4()
{
	//check if the polynomial is of the form x^4-k, if so there is a fast SNFS relation value calculation routine to use. 
	if (degree == 4 && *coefficients[1] == 0 && *coefficients[2] == 0 && *coefficients[3] == 0)
	{
		return true; 
	}
	else
	{
		return false;
	}
}

void Polynomial::evaluate_binomial_polynomial_homogeneous_SNFS_deg4(mpz_class* result, int64 a, int64 b, mpz_class* temp)
{
	//Evaluate b^deg*f((a/b)^deg). 
	//fast version for deg-4 polynomial with x, x^2 and x^3 = 0, which happends to be many SNFS numbers i am interested in...
	//also, use external already defined temp number to store information, to not have to declare a new scratchpad memory.

		*temp = a;
		*temp = *temp * *temp;
		*temp = *temp * *temp;
		*result = *coefficients[degree] * *temp;

		*temp = b;
		*temp = *temp * *temp;
		*temp = *temp * *temp;
		*result = *result + *coefficients[0] * *temp;
		return;

}

void Polynomial::debug_print_poly()
{
	cout << "Printing polynomial" << endl;
	for(int i = degree; i >=0; i--)
	{
		cout << "x^" << i << " " << coefficients[i]->get_str() << endl;
	}
}

void Polynomial::update_degree()
{
	//Update the degree if it has changed.
	int i = 2*POLYNOMIAL_MAX_DEGREE;

	while (i >= 0 && *coefficients[i] == 0)
	{
		i--;
	}

	if (i <= 0)
	{
		i = 0;
	}

	if (i > POLYNOMIAL_MAX_DEGREE)
	{
		cout << "A polynomial of too high degree was found!" << endl;
		return;
	}
	degree = i;
}

void Polynomial::mult_modpoly(Polynomial* a, Polynomial* modpoly, uint32 p)
{
	//Multiply this polynomial with a, take the polynomial mod modpoly and take all coefficients mod p.
	Polynomial* result = new Polynomial;

	for(int i = 0; i <= a->degree; i++)
	{
		for(int j = 0; j <= this->degree; j++)
		{
			*result->coefficients[i+j] = (*result->coefficients[i+j] + *a->coefficients[i] * *this->coefficients[j]) % p;
		}
	}

	result->degree = a->degree + this->degree;
	this->copy(result);
	delete result;

	this->reduce_modpoly(modpoly, p);
}

void Polynomial::reduce_modp(uint32 p)
{
	//reduce all coefficients mod p
	mpz_class* prime = new mpz_class;
	*prime = p;
	for(int i = 0; i <= degree; i++)
	{
		mpz_mod(coefficients[i]->get_mpz_t(), coefficients[i]->get_mpz_t(), prime->get_mpz_t());
	}
	delete prime;
}

void Polynomial::reduce_modp(mpz_class* prime)
{
	//reduce all coefficients mod p
	for(int i = 0; i <= degree; i++)
	{
		mpz_mod(coefficients[i]->get_mpz_t(), coefficients[i]->get_mpz_t(), prime->get_mpz_t());
	}
}

void Polynomial::make_monic(uint32 p)
{
	//make a polynomial monic mod p

	if (*coefficients[degree] % p == 0)
	{
		cout << "Cannot make polynomial monic!" << endl;
		return;
	}

	mpz_class* temp = new mpz_class;
	mpz_class* prime = new mpz_class;

	*prime = p;
	mpz_invert(temp->get_mpz_t(), coefficients[degree]->get_mpz_t(), prime->get_mpz_t());

	if (*temp == 0)
	{
		cout << "Inversion failed in poly_make_moic" << endl;
		return;
	}

	for(int i = 0; i <= degree; i++)
	{
		*coefficients[i] = *coefficients[i] * *temp;
		mpz_mod(coefficients[i]->get_mpz_t(), coefficients[i]->get_mpz_t(), prime->get_mpz_t());
	}

	delete temp;
	delete prime;
}

void Polynomial::make_monic(mpz_class* prime)
{
	//make a polynomial monic mod p

	mpz_class* temp = new mpz_class;
	mpz_mod(temp->get_mpz_t(), coefficients[degree]->get_mpz_t(), prime->get_mpz_t());

	if (*temp == 0)
	{
		cout << "Cannot make polynomial monic!" << endl;
		delete temp;
		return;
	}

	mpz_invert(temp->get_mpz_t(), coefficients[degree]->get_mpz_t(), prime->get_mpz_t());
	if (*temp == 0)
	{
		cout << "Inversion failed in poly_make_moic" << endl;
		delete temp;
		return;
	}
	for(int i = 0; i <= degree; i++)
	{
		*coefficients[i] = *coefficients[i] * *temp;
		mpz_mod(coefficients[i]->get_mpz_t(), coefficients[i]->get_mpz_t(), prime->get_mpz_t());
	}
	delete temp;
}

void Polynomial::reduce_modpoly(Polynomial* modulus, uint32 p)
{
	//reduce this polynomial mod modulus, take coefficients mod p.
	if (modulus->degree == 0) 
	{
		clear();
		return;
	}

	int32 i;
	mpz_class* top_value = new mpz_class;
	Polynomial* mod = new Polynomial;
	mpz_class* c = new mpz_class;
	uint32 j;

	mod->copy(modulus);
	mod->make_monic(p);

	while (degree >= mod->degree) 
	{
		*top_value = *coefficients[degree];

		*coefficients[degree] = 0;
		for (i = mod->degree-1; i >= 0; i--) 
		{
			*c = *top_value * *mod->coefficients[i];
			mpz_mod_ui(c->get_mpz_t(), c->get_mpz_t(), p);
			j = degree - (mod->degree - i);
			*c = *coefficients[j] - *c;
			mpz_mod_ui(coefficients[j]->get_mpz_t(), c->get_mpz_t(), p);
		}
		this->update_degree();
	}
	delete top_value;
	delete mod;
	delete c;

	this->update_degree();
}

void Polynomial::gcd(Polynomial* a, Polynomial* b, uint32 p)
{
	//gcd of two polynomials.

	Polynomial* h = new Polynomial;
	Polynomial* r = new Polynomial;
	Polynomial* result = new Polynomial;
	Polynomial* test = new Polynomial;

	if (a->degree > b->degree) 
	{
		result->copy(a);
		h->copy(b);
	}
	else 
	{
		h->copy(a);
		result->copy(b);
	}

	while ((h->degree > 0) || (*h->coefficients[h->degree] != 0)) 
	{
		r->copy(result);
		r->reduce_modpoly(h,p);
		result->copy(h);
		h->copy(r);
	}

	this->copy(result);

	if (this->degree == 0)
	{
		*this->coefficients[0] = 1;
	}

	delete h;
	delete r;
	delete result;
	delete test;
	update_degree();
}

void Polynomial::square_poly_modpoly(Polynomial* modpoly, uint32 p)
{
	//mutliply a polynomial with itself modulo modpoly
	this->mult_modpoly(this, modpoly, p);
} 

void Polynomial::poly_exponentiate(uint32 power, Polynomial* modpoly, uint32 p)
{
	//calculate poly^power mod modpoly

	Polynomial* powerpol = new Polynomial;
	Polynomial* result = new Polynomial;

	powerpol->copy(this);
	result->set_coeff(0,1);
	uint32 temppower = power;

	while(temppower > 0) 
	{
		if(temppower % 2 == 1) 
		{
			result->mult_modpoly(powerpol, modpoly, p);
			temppower--;
		}
		powerpol->square_poly_modpoly(modpoly,p);
		temppower = temppower / 2;
	}

	this->copy(result);
	delete result;
	delete powerpol;
} 

void Polynomial::divide_out_root(uint32 root, uint32 p)
{
	//divide polynomial by (x - root)

	Polynomial* result = new Polynomial;

	for(int i = degree; i > 0; i--)
	{
		result->set_coeff(i-1, coefficients[i]->get_ui());
		*coefficients[i-1] = *coefficients[i-1] + *coefficients[i] * root;
		mpz_mod_ui(coefficients[i-1]->get_mpz_t(),coefficients[i-1]->get_mpz_t(),p);
		*coefficients[i] = 0;
	}

	copy(result);
	delete result;
}


void Polynomial::get_polyroots_modp(uint32* numroots, uint32* rootlist, uint32 p)
{
	//find all root of the polynomial modulo p. return the number of roots and a list with the roots

	uint32 currentroot =  0;
	*numroots = 0;

	if(p <= 25)
	{ 
		mpz_class* mpztemp = new mpz_class;
		//Handle very small primes by brute force. Limit 25 is arbitrary but small.

		for(int i = 0; i < p; i++)
		{
			this->evaluate_polynomial_modulo(mpztemp,i,p);

			if(*mpztemp == 0)
			{
				*numroots = *numroots + 1;
				rootlist[currentroot] = i;
				currentroot++;
			}
		}
		delete mpztemp;
		return;
	}

	//treat root x=0 separately...
	mpz_class* mpztemp = new mpz_class;
	mpz_mod_ui(mpztemp->get_mpz_t(), this->coefficients[0]->get_mpz_t(), p);
	if(*mpztemp == 0)
	{
		//x=0 is a root
		*numroots = *numroots + 1;
		rootlist[0] = 0;
		currentroot = 1;
	}

	Polynomial* allroots = new Polynomial;
	Polynomial* temp = new Polynomial;

	allroots->set_coeff(1,1);
	allroots->poly_exponentiate(p-1,this,p);

	*allroots->coefficients[0] = (*allroots->coefficients[0] - 1) ;
	mpz_mod_ui(allroots->coefficients[0]->get_mpz_t(), allroots->coefficients[0]->get_mpz_t(), p);
	//temp is now x^(p-1)-1

	allroots->gcd(this, allroots,p);

	if(allroots->degree == 0)
	{
		delete mpztemp;
		delete allroots;
		delete temp;
		return;
	}
	else
	{
		*numroots = *numroots + allroots->degree;
	}

	if(allroots->degree == 1)
	{
		//only one more root in [1, p-1] so we find it directly.
		allroots->make_monic(p);
		mpz_mod_ui(mpztemp->get_mpz_t(),allroots->coefficients[0]->get_mpz_t(),p);
		rootlist[currentroot] = (p - mpz_get_ui(mpztemp->get_mpz_t()));

		delete temp;
		delete mpztemp;
		delete allroots;
		return;
	}

	//more than one root in [1, p-1] so we need to find them... Cantor-Zassenhaus to the rescue!
	//allroots contain gcd((x^(p-1) - 1), polynomial) and is used to track the roots found

	uint32 tries = 0; // not sure if needed, but just to make sure that the loop stops at some point
	int32 root;
	allroots->make_monic(p);

	while (allroots->degree > 0 && tries < p) //tries < p should guarantee finding the root, since it will essentially try all roots... early abort may be considered.
	{
		tries++;
		temp->clear();
		temp->set_coeff(0,tries);
		temp->set_coeff(1,1);

		temp->poly_exponentiate((p-1)/2,allroots,p);
		*temp->coefficients[0] = *temp->coefficients[0] - 1;
		mpz_mod_ui(temp->coefficients[0]->get_mpz_t(), temp->coefficients[0]->get_mpz_t(), p);
		//temp is now x^((p-1)/2) - 1 mod allroots.

		temp->gcd(temp, this, p);

		if (temp->degree == 1)
		{
			//one root isolated
			temp->make_monic(p);
			root = p - temp->coefficients[0]->get_ui();

			if(root < 0)
			{
				root += p;
			}

			rootlist[currentroot] = root;
			allroots->divide_out_root(rootlist[currentroot],p);
			allroots->make_monic(p);
			currentroot++;

		}

		//if only one root remains
		if(allroots->degree == 1)
		{
			rootlist[currentroot] = p - allroots->coefficients[0]->get_ui();
			allroots->divide_out_root(rootlist[currentroot],p);
		}
	}

	if(tries == p)
	{
		cout << "Polynomial rootfinding failed to find roots, this is a problem... " << endl;
	}

	delete allroots;
	delete temp;
	delete mpztemp;
	return;
} 

