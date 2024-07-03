#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <cstdlib>
#include <stdio.h>
#include <string.h>

#include "dataTools.h"
using namespace NTL;
using namespace std;

/**********************************************/

typedef unsigned long long int usfixn64;
int verbose = 0;
double totalElapsed = 0;

/**********************************************/

void
inline
hrule (int len = 25)
{
	for (int i = 0; i < len; i++)
		printf ("=");
	printf ("\n");
}

/**********************************************/

void
computeBigPrime (ZZ& bigPrime, ZZ radix, usfixn64 coefficientSize)
{

	power (bigPrime, radix, coefficientSize);
	add (bigPrime, bigPrime, 1);

}

/**********************************************/
//convert to big field element
void
toBigFieldElement (usfixn64 *bigFieldElement, ZZ bigInt, ZZ prime, ZZ radix,
									 usfixn64 coefficientSize)
{
	rem (bigInt, bigInt, prime);
	ZZ m;
	for (int i = 0; i < coefficientSize; i++)
	{
		rem (m, bigInt, radix);
//		cout<<"m " <<m<<endl;
		bigFieldElement[i] = to_ulong (m);
//		cout<<bigFieldElement[i]<<endl;
		div (bigInt, bigInt, radix);
	}
}

/**********************************************/
void
toBigInt (ZZ& bigInt, usfixn64* bigFieldElement, ZZ radix,
					usfixn64 coefficientSize)
{
	bigInt = 0;
	ZZ t;
	char numStr[64];

	for (int i = 0; i < coefficientSize; i++)
	{
		sprintf (numStr, "%llu", bigFieldElement[i]);
		conv (t, numStr);
//		t = conv<ZZ>(bigFieldElement[i]);
		mul (t, t, power (radix, i));
		add (bigInt, bigInt, t);
	}

}
/**********************************************/
void
print_big_field_element (usfixn64* bigFieldElement, usfixn64 coefficientSize)
{
	for (int i = 0; i < coefficientSize; i++)
		cout << "[" << i << "] = " << bigFieldElement[i] << endl;
}
/**********************************************/

int
main (int argc, char** argv)
{
	int K = 16; //2*coefficientSize
	int e = 1; //from input
	int coefficientSize = K / 2;
	verbose = 0;
	char helpMessage[1024] =
			"\n"
					"\t-----------------------------------------------------------------------"
					"\n"
					"\t -Precomputing powers of omega for vectors of K^e coefficients"
					"\n"
					"\t used in big-prime-field-FFT p=r^{K/2}+1"
					"\n"
					"\t 1:K,  2:e,  3: radix_high_pow,  4: radix_low_pow,  5:path_omega,  6:path_pow_omega,   7:-v "
					"\n"
					"\n"
					"\t-----------------------------------------------------------------------"
					"\n";

	int radix_lower_pow, radix_higher_pow;
	radix_higher_pow = 63;
	radix_lower_pow = 34;

	char omega_path[1024] = "omega";
	char powers_omega_path[1024] = "pow_omega";
	//N=K^e = inputVectorSize
	//1:K
	//2:e
	//3: radix_higher_power
	//4: radix_higher_lower
	//5: path_to_omega
	//6: path_to_computed_powers_of_omega

	if (argc == 1)
	{
		printf ("%s\n", helpMessage);
		return -1;
	}
	if (argc > 1)
	{
		if (strcmp (argv[1], "-h") == 0)
		{
			printf ("%s\n", helpMessage);
			return -1;
		}
		else
		{
			K = atoi (argv[1]);
			coefficientSize = K / 2;
		}
	}
	if (argc > 2)
	{
		e = atoi (argv[2]);
	}
	if (argc > 3)
	{
		radix_higher_pow = atoi (argv[3]);
	}
	if (argc > 4)
	{
		radix_lower_pow = atoi (argv[4]);
	}
	if (argc > 5)
	{
		sprintf (omega_path, argv[5]);
	}
	if (argc > 6)
	{
		sprintf (powers_omega_path, argv[6]);
	}
	if (argc > 7)
	{
		if (strcmp (argv[7], "-v"))
			verbose = 1;
		if (strcmp (argv[7], "-verbose"))
			verbose = 1;
	}
//	r = (1UL << radix_higher_pow) + (1 << radix_lower_pow);
	ZZ r, radix_high, radix_low;
	power (radix_high, 2, radix_higher_pow);
	power (radix_low, 2, radix_lower_pow);
	add (r, radix_high, radix_low);

	ZZ p (0);
	ZZ omega (0);
	ZZ m (0);

	printf ("\n... K=%d, e=%d, Radix = 2^%d + 2^%d (%lu), coefficientSize=%d\n",
					K, e, radix_higher_pow, radix_lower_pow, to_ulong (r),
					coefficientSize);
	computeBigPrime (p, r, coefficientSize);

	usfixn64 *field_element = (usfixn64*) malloc (
			coefficientSize * sizeof(usfixn64));

	readVectorFromFile (field_element, 1, coefficientSize, omega_path);

	toBigInt (omega, field_element, r, coefficientSize);

	printVectorToFile (field_element, 1, coefficientSize, powers_omega_path, 0);
	m = omega;
	for (int i = 2; i < pow (K, e - 1) + 1; i++)
	{
		mul (m, m, omega);
		rem (m, m, p);
		toBigFieldElement (field_element, m, p, r, coefficientSize);
		printVectorToFile (field_element, 1, coefficientSize, powers_omega_path, 1);
	}

	printf ("\n... Precomputed powers of omega are written to %s\n\n", powers_omega_path);
	return 0;
}
