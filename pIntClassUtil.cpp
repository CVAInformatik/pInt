
/*
Copyright  © 2024 Claus Vind - Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/

#include <iostream>
//#define PERF
#ifdef PERF
//#include "winnt.h"
#include "windows.h"
#include "profileapi.h"
#endif

#include "PrimeTable.h"
#include "pIntClass.h"
#include "pIntClassRandom.h"
#include "pIntClassUtil.h"


pIntClass& exponentiation(int  a, int exp)
{
    pIntClass* res = new pIntClass();

    unsigned int  mask = exp >> 1;
    if (exp == 0)
    {
        *res = 1;
    }
    else if (mask == 0) // exponent == 1
    {
        *res = a;
    }
    else {
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        mask++;  // mask is now a single 1 at the postion of the leftmost 1 in exp
        mask = mask >> 1;
        *res = a;
        do {
            *res *= *res;
            if (exp & mask)
                *res *= a;
            mask = mask >> 1;
        } while (mask);
    }
    return *res;
}


pIntClass& RemQuotient(const pIntClass& A, const pIntClass& M, pIntClass* Quotient)
{
	pIntClass  _Quotient;
	pIntClass* Rem = new pIntClass();
	int counter = 0;

	if (M.IsZero()) {
		std::cout << "divison by zero" << std::endl;
		return *Rem;
	}
	else if ((A.value.size() == 1) && (M.value.size() == 1)) {
		/* small numbers both less than RMOD */
		if (Quotient != NULL)
		{
			Quotient->value.clear();
			Quotient->value.push_back(A.value[0] / M.value[0]);
		}

		Rem->value.push_back( A.value[0] % M.value[0]);
        Rem->normalize(Rem->value);
		return *Rem;
	}
	else {
        pIntClass  _dividend(A);
        pIntClass  _divisor(M);
        int     reciprocal = pIntClass::MODULUS / (2 + _divisor.value.back());
        int     shift = (int)M.value.size();
		pIntClass  _reciprocal(reciprocal);

        *Rem = reciprocal;
		*Rem *= _dividend;

		if (!Rem->IsZero())  
            for (int i = 0; i < shift;i++)	Rem->DivModulus();

		while (1)
		{
			_Quotient = *Rem;
			*Rem *= _divisor;
            Rem->ChSignBit();
			*Rem += _dividend;

			if (_divisor.IsBiggerNummerically(*Rem)) 	
                break;

			*Rem *= _reciprocal;
			for (int i = 0; i < shift;i++) {
				Rem->DivModulus();
			}
			if (Rem->IsZero())
				*Rem = 1;
			*Rem += _Quotient;
		}
		/* we are done */
		if (Rem->IsNeg())
		{
			*Rem += _divisor;
			_Quotient += -1;
		}
		if (Quotient)  *Quotient = _Quotient;
		return *Rem;
	}
}

int Jacobi(const pIntClass& a, const pIntClass& b)
{

	pIntClass A(a);
	pIntClass M(b);
	int ResSign = 1;
	if (M.IsZero())
		return 1;
	else {
		A =RemQuotient(A, M, NULL);
		if (A.IsZero()) {  //A|M 
			return 0;
		}
		else {
			while (!A.IsZero()) {
				while (!A.IsZero() && ((A[0] & 1) == 0))
				{
					A >>= 1;
					switch (M[0] & 0x7)
					{
					case 3: 
                    case 5:   
                        ResSign = -ResSign;
						break;
					default:
						break;
					}
				}
				pIntClass temp;
				temp = A;
				A = M;
				M = temp;
				if (!A.IsZero() && !M.IsZero() && (3 == (A[0] & 0x3)) && (3 == (M[0] & 0x3))) 
                    ResSign = -ResSign;

                A = RemQuotient(A, M, NULL);
			}
			if ((M == 1))
				return ResSign;
			else
				return 0;
		}
	}
};


void testTonelliShanks()
{
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif

    pIntClass  P("26959946667150639794667015087019630673557916260026308143510066298881");
    pIntClass  A("18958286285566608000408668544493926415504680968679321075787234672564");
    pIntClass  Res;
#define P224 1
#if P224
    // NIST P-224 
    for (int i = 0; i < 5; i++) {
        pIntClass t("2021");
        pIntClass  A("18958286285566608000408668544493926415504680968679321075787234672564");
        t = t * t * t;
        A += t;
        A += -1*(3 * 2021);
#ifdef PERF
        QueryPerformanceFrequency(&Frequency);
        QueryPerformanceCounter(&StartingTime);
#endif
        TonelliShanks(A, P, Res);
#ifdef PERF
        QueryPerformanceCounter(&EndingTime);
        ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
        ElapsedMicroseconds.QuadPart *= 1000000;
        ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
        std::cout << "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif

        std::cout << std::endl << " A:   " << A.ToString() << std::endl;
        std::cout << " P:   " << P.ToString() << std::endl;
        std::cout << " Res: " << Res.ToString() << std::endl;
        Res = Res * Res;
        std::cout << " Res * Res  : " << Res.ToString()<< std::endl;
        Res = RemQuotient(Res, P, NULL);
        std::cout << " Res * Res  mod P: " << Res.ToString() << std::endl;
    }
    return;
#else
#endif

}

bool  TonelliShanks(const pIntClass& n, const pIntClass& p, pIntClass& res) {

	pIntClass* Res = new pIntClass();
	if (Jacobi(n, p) != 1) {
		*Res = 0;
		res = *Res;
		return false;
	}

	// Factor out powers of 2 from p - 1
	pIntClass q = p; q -= 1;
	int s = 0;
	while (((q[0] & 1) == 0)) {
		if (q.IsZero()) {
            res = 0;
			std::cout << "not a square " << std::endl;
			return false;
		}
		q >>= 1;
		s++;
	}

	// Find a non-square z such as ( z | p ) = -1
	pIntClass z = 2;
	pIntClassRandom Rands(pIntClass::MODULUS);

	Rands.SetSeed(p);
	do {
		z = Rands.Rand();
	} while (Jacobi(z, p) != -1);

	pIntClass c = modpow(z, q, p);
	//std::cout << " z: " << z.ToString() << std::endl;
	//std::cout << " q: " << q.ToString() << std::endl;
	pIntClass p1 = p;
	//std::cout << " p: " << p1.ToString() << std::endl;
	//std::cout << " c: " << c.ToString() << std::endl;
	pIntClass  t = modpow(n, q, p);
	pIntClass  m = s;
	pIntClass  m_1 = m; m_1;
	pIntClass  result = modpow(n, (q + 1) >>= 1, p);

	//std::cout << " c: " << c.ToString() << std::endl;
	//std::cout << " t: " << t.ToString() << std::endl;
	//std::cout << " m: " << m.ToString() << std::endl;
	//std::cout << " R: " << result.ToString() << std::endl;

loop:
	if (t.IsZero()) {
		res = 0;
		std::cout << "Root is 0" << std::endl;
		return false;
	}
	if (t.IsOne()) {
        res = result;
		std::cout << "Root is " << result.ToString() << std::endl;
		return true;
	}
	int  i = 0;
	pIntClass t1 = t;
	do {
		i++;
		t1 = modmult(t1, t1, p);
	} while (!t1.IsOne());



	pIntClass temp1 = m;
	temp1 -= 1;
	temp1 -= i;

	pIntClass temp2;

	temp2 = modpow(2, temp1, p);

	pIntClass  b = modpow(c, temp2, p);
	c = modmult(b, b, p);
	t = modmult(t, c, p);
	m = i;
	result = modmult(result, b, p);
	//std::cout << " c: " << c.ToString() << std::endl;
	//std::cout << " t: " << t.ToString() << std::endl;
	//std::cout << " R: " << result.ToString() << std::endl;

	goto loop;
	return  true;
}




pIntClass modmult(const pIntClass &_a, const pIntClass &_b, const pIntClass &mod) {  // Compute a*b % mod
    pIntClass result;
    pIntClass a = _a;
    pIntClass b = _b;
    while (!b.IsZero()) {

        if ((b[0] & 1) == 1) {
            result += a;
            result = RemQuotient(result, mod, NULL);
        }
        a = a + a;
        a = RemQuotient(a, mod, NULL);
        b >>= 1;
    }
    return result;
}

pIntClass modpow(const pIntClass& _a, const pIntClass& _b, const pIntClass& mod) {  // Compute a^b % mod
    pIntClass result;
    pIntClass  a = _a;
    pIntClass  b = _b;
    pIntClass t;
    ++(result);

    while (!b.IsZero()) {
        if ((b[0] & 1) == 1) {
            t = modmult(result, a, mod);
            result = t;
        }
        a = modmult(a, a, mod);
        b >>= 1;
    }

    return result;
}

void testMR1(int &npcount, int &pcount, int width,  pIntClass& p)
{
    pIntClass prime = p;

    std::string s = prime.ToString();
    std::string prefix = "";

    while (prefix.length() + s.length() < width) prefix = prefix + " ";

    if (MillerRabin(prime, 30)) {
        std::cout << "n " << prefix << s << " is probably prime  " <<   ++pcount << std::endl;
    }
    else {
        std::cout << "n " << prefix << s << " is not prime  " <<  ++npcount + pcount <<  std::endl;
    };
}


void testModMult()
{
    pIntClass mod("2147483647");
    pIntClass result("4026531840");
    std::cout << "mod : " << mod.ToString() << std::endl;
    std::cout << "result: " <<result.ToString() << std::endl;
     result = RemQuotient(result, mod, NULL);
     std::cout << "result: " << result.ToString() << std::endl;


}



void testMR()
{
    pIntClass prime("26959946667150639794667015087019630673557916260026308143510066298881");

    pIntClass prime1("5127821565631733");


    pIntClass prime2("2147483647");
#define WIDTH 75
    int np = 0, p = 0;
    testMR1(np, p, WIDTH, prime);
    testMR1(np, p, WIDTH, prime1);
    testMR1(np, p, WIDTH, prime2);

    pIntClassRandom Rands(pIntClass::MODULUS);

    pIntClass r = prime;
    r -= 3;
    Rands.SetSeed(r);

#define COUNT 1000

    np = 0;
    p = 0;
    for (int i = 0; i < COUNT; i++) {
#if 1
        pIntClass pc = Rands.Rand();
        if ((pc[0] & 1) == 0) pc++;
 //       std::cout << "mod : " << p.ToString() << std::endl;
        testMR1(np, p, WIDTH, pc);
#else
        testModMult();
#endif

    }
}

bool MillerRabin(const pIntClass& number, int witnesses)
{
    pIntClass m = number;

    if (!m.IsZero() && ((m[0] & 1) == 0 )) {
        std::cout << "argument must be odd " << std::endl;
        return false;
    }
    else {

        pIntClassRandom Rands(pIntClass::MODULUS);

        pIntClass r = number;
        r -= 3;

        Rands.SetSeed(r);

        pIntClass d = m;  
        
        d -= 1;
        
        int s = 0;

        while (!d.IsZero() && ((d[0] & 1) == 0)) {
            d >>= 1;
            s++;
        }

        for (size_t ix = 0; ix < witnesses; ix++)
        {
            pIntClass a = Rands.Rand();
            a += 2;

            //std::cout << "a: " << a.ToString() << std::endl;

            pIntClass x = modpow(a, d, m);

            for (int i = 0; i < s; i++) {
               pIntClass y = modmult(x, x, m);
               pIntClass t = m;
               t -= x;
               if (y.IsOne() && !x.IsOne() && !t.IsOne() ) {
                   std::cout << "mr fail at " << ix << " ";
                    return false;  
                }
                x = y;
            }
            if (!x.IsOne()) {
                std::cout << "mr fail at " << ix << " ";
                return false;
            }
        }
    }
    return true;
}


bool  CheckedTonelliShanks(const pIntClass& n, const pIntClass& p, pIntClass& res) {

    if (MillerRabin(p, 30))
        return TonelliShanks(n, p, res);
    return
        false;
}


