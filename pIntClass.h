#pragma once

/*
Copyright  © 2024 Claus Vind - Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


/*
        pIntClass, 
                 an experiment in constructing a portable large int using modern C++

*/

#include <vector>
#include <string>


#if OS_WINDOWS   // windows
//#define PERF
typedef unsigned __int64  u64;
typedef __int64  s64;
#else   //linux
typedef  uint64_t  u64;
typedef  int64_t  s64;
#endif



/*
*     The pIntClass represents multiprecision integers as normalized std::vector<int> of digits 
*     in the range [-999.999.999....999.999.999]. 
*     This leaves a 1 bit headroom in both ends, which makes for convenient add/subtract 
*     operations, e.g. no need for looking at sizes of operands and such.
*     
* 
*     A pIntClass number is normalized if
* 
*      value.size() == 0        -> meaning the number has the value 0
* 
*      all non-zero entries in value have the same sign and are in the open
*                interval ] -1000000000,...., 1000000000[
*
*      and the most significant digit is non-zero.
*
*      that means that the sign of a non-zero pIntClass2 number is the sign of the
*          value.back(), the last non-zero entry).
*
*	   Addition and Subtraction may temporarily make a value not normalized,
*	   normalize() will bring it back in order.
* 
*/


/* undefine if you don't want to include Schönhage-Strassen multiplication */
#define SSLIMIT 220


class pIntClass {

public:

    pIntClass();
    pIntClass(const pIntClass &x);
    pIntClass(const std::string  &x);
    pIntClass(const int x);

    ~pIntClass();

    static const int ReservationSize = 8;
    static const int MODULUS = 1000000000;

    std::string& ToString();

    pIntClass& operator+=(const pIntClass& rhs);
    pIntClass& operator*=(const pIntClass& rhs);

    pIntClass& operator-=(const pIntClass& rhs);
    pIntClass& operator=(const pIntClass& rhs);

    pIntClass& operator+=(const int rhs);
    pIntClass& operator*=(const int rhs);

    pIntClass& operator-=(const int rhs);
    pIntClass& operator=(const int  rhs);




    pIntClass& operator<<=(const unsigned int shift);
    pIntClass& operator>>=(const unsigned int shift);

    pIntClass& operator++();
    pIntClass operator++(int dummy);
    pIntClass& operator--();
    pIntClass operator--(int dummy);

    int operator[](int index); // not really an array
    int Size() {   return (int) value.size();   };
    int Sign();   // 0 for 0, 1 for positive, -1 for negative
    int ChSignBit(); // returns the new value of the sign
    inline bool IsZero() const { return value.size() == 0; } // == 0
    inline bool IsPos() const { return value.size() == 0 || ( value.back() > 0); }//  =< 0
    inline bool IsNeg() const { return value.size() != 0 && (value.back() < 0); } // < 0
    inline bool IsOne() const {  return (value.size() == 1) && (value[0] ==1 ); }
    inline bool IsMinusOne() const { return (value.size() == 1) && ( value[0] == -1); }

    bool operator!=(const int b);
    bool operator!=(const pIntClass& b);
    bool IsBiggerNummerically(const pIntClass& b);

    friend pIntClass& operator+(pIntClass &a, const pIntClass &b);
    friend pIntClass& operator*(pIntClass& a, const pIntClass& b);

    friend bool operator<(const pIntClass& a, const pIntClass& b);
    friend bool operator==(const pIntClass& a, const pIntClass& b);
    friend bool operator!=(const pIntClass& a, const pIntClass& b);

    friend pIntClass& operator+(const int a, const pIntClass &b);
    friend pIntClass& operator-(const int a, const pIntClass& b);
    friend pIntClass& operator*(const int a, const pIntClass &b);


    friend int Jacobi(const pIntClass& a, const pIntClass& b);
    friend pIntClass& RemQuotient(const pIntClass& A, const pIntClass& M, pIntClass *Quotient);

    friend class pIntClassRandom;

private:

    void mul10();
    void normalize(std::vector<int> &value);
    void Scale(int scale);
    void DivModulus();

    /* schoolbook multiplication */
    pIntClass& SchoolbookMultiplication(const pIntClass& rhs);
    //pIntClass& KaratsubaMultiplication(const pIntClass& rhs);


#ifdef SSLIMIT
    pIntClass& SchoenhageStrassenMultiplication(const pIntClass& rhs);

    void LoadFFT(const std::vector<int> &A, double* Buffer);
    void Carry(s64 size, double* Buffer);
#endif

    std::vector<int> value;
};


