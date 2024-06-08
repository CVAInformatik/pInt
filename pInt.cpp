

/*
Copyright  © 2024 Claus Vind - Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


// pInt.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include "pIntClass.h"
#include "pIntClassUtil.h"


#ifdef PERF
#include "windows.h"
#include "profileapi.h"
#endif


void test13QuotientReminder()
{


    pIntClass num1("777777777777777");
    pIntClass num("5555555555555");



#define LIMIT1 100000
#define STEP  177395769

    for (int i = 0; i < LIMIT1; i++)
    {
        pIntClass t1 = num1;
        pIntClass t = num;
        pIntClass step(STEP);
        pIntClass Quotient;
        pIntClass Remainder;
        step *= i;
        t += step;
        Remainder =RemQuotient(t1, t, &Quotient);
        t = num;
        step = STEP;
        step *= i;
        t += step;
        std::cout << "Quotient: " << Quotient.ToString() << std::endl;
        t *= Quotient;
        t += Remainder;
        if (num1 != t) {
            std::cout << "num * Quotient + Remainder: " << t.ToString() << std::endl;
            std::cout << "Num1                      : " << num1.ToString() << std::endl;
        }

    }
}


void test15Jacobi()
{

    for (int j = 3; j < 60; j += 2)
        for (int i = 1; i < 31; i++) {
            std::cout << "Jacobi ( ";
            pIntClass cj(j);
            pIntClass ci(i);
            std::cout << ci.ToString() << "/ ";
            std::cout << cj.ToString() << " ) = ";            
            std::cout << Jacobi( ci, cj) << std::endl;
        }

}

void testSubtraction();
void testAddition();
void testMultiplication();
void testTonelliShanks();
void testMR();
void testModMult();

void testpIntClass2(){

    pIntClass  a("4999999999999999999999999999");
    pIntClass  b("5000000000000000000000000000");

    pIntClass a1(a);
    pIntClass b1(b);

    std::cout << "1 a1: " << a1.ToString() << std::endl;
    std::cout << "2 b1: " << b1.ToString() << std::endl;

    a1 -= b1;

    std::cout << "3 a1: " << a1.ToString() << std::endl;
    std::cout << "4 b1: " << b1.ToString() << std::endl;

    a1 += b1;

    std::cout << "5 a1: " << a1.ToString() << std::endl;
    std::cout << "6 b1: " << b1.ToString() << std::endl;

    --b1;
    a1 -= b1;

    std::cout << "7 a1: " << a1.ToString() << std::endl;
    std::cout << "8 b1: " << b1.ToString() << std::endl;

    a1 = a;
    b1 = b;

    a1 = b + b;
    b1 = b1 + b1;
    std::cout << "9  a1: " << a1.ToString() << std::endl;
    std::cout << "10 b1: " << b1.ToString() << std::endl;

    a1 = pIntClass("4999999999");
    b1 = b;

    a1 = a1 + a1;
    b1 = a1 + b1;
    std::cout << "11 a1: " << a1.ToString() << std::endl;
    std::cout << "12 b1: " << b1.ToString() << std::endl;

    a1 = pIntClass("4999999999");
    b1 = b;
    std::cout << "13 a1: " << a1.ToString() << std::endl;
    std::cout << "14 b1: " << b1.ToString() << std::endl;

    a1 = a1 + a1;
    std::cout << "15 a1: " << a1.ToString() << std::endl;
    b1 = a;
    std::cout << "16 b1: " << b1.ToString() << std::endl;
    b1 += a1;
    std::cout << "17 b1: " << b1.ToString() << std::endl;
    b1 += b1;
    std::cout << "18 b1: " << b1.ToString() << std::endl;
    b1 += 7;
    std::cout << "19 b1: " << b1.ToString() << std::endl;
    b1 += 10;
    std::cout << "20 b1: " << b1.ToString() << std::endl;

    a1 = pIntClass("4999999999");
    b1 = b;
    std::cout << "21 a1: " << a1.ToString() << std::endl;
    std::cout << "22 b1: " << b1.ToString() << std::endl;

    a1 = a1 + a1;
    std::cout << "23 a1: " << a1.ToString() << std::endl;
    b1 = a + a1 ;
    std::cout << "24 b1: " << b1.ToString() << std::endl;
    b1 = b1 + b1 +  7 + 10;
    std::cout << "25 b1: " << b1.ToString() << std::endl;

}

void testExponentiation()
{
    pIntClass res;

    res = exponentiation(2, 86243);
    --res;
    std::cout << "res: " << res.ToString() << std::endl;

    res = exponentiation(2, 1257787);
    --res;
    std::cout << "res: " << res.ToString() << std::endl;

    res = exponentiation(2, 3021377);
    --res;
    std::cout << "res: " << res.ToString() << std::endl;

    res = exponentiation(2, 20996011);
    --res;
    std::cout << "res: " << res.ToString() << std::endl;

}

#define TESTLENGTH 100
void testLeaks()
{
    pIntClass* A = new pIntClass();

    for (int i = 0; i < TESTLENGTH; i++)
    {
        delete A;
        A = new pIntClass();
        *A = 100;
        std::string s;
        s = A->ToString();
        ++(*A);
        std::cout << "A " << A->ToString() << std::endl;
    }
    for (int i = 0; i < TESTLENGTH; i++)
    {
        delete A;
        A = new pIntClass();
        *A = (*A) + 7;
        *A *= 100;
        *A = (*A) * 100;
        std::string s;
        s = A->ToString();
        std::cout << "A " << A->ToString() << std::endl;
    }

    delete A;
}




int main()
{


    testLeaks();
    testMR();
    testExponentiation();
    testSubtraction();
    testAddition();
    testMultiplication();
    testModMult();
    testTonelliShanks();



}
