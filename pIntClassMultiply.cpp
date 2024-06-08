
/*
Copyright  © 2024 Claus Vind - Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/

#include <iostream>
#ifndef OS_WINDOWS
#include <math.h>
#endif

#include "pIntClass.h"

#ifdef SSLIMIT
#include "PrimeFactorDFT.h"
#endif

void pIntClass::mul10()
{
	int carry = 0;
	s64 t = 0;

	for (int i = 0; i < value.size(); i++)
	{
		t = (s64)value[i];
		t *= 10;
		t += carry;
		carry = 0;
		if (t >= MODULUS) {
			carry = (int)(t / MODULUS);
			t = t % MODULUS;
		}
		value[i] = (int)t;
	}
	if (carry) value.push_back(carry);
}


pIntClass pIntClass::operator*=(const pIntClass& rhs)
{
	size_t sz = value.size() + rhs.value.size();
#ifdef SSLIMIT
	if(sz >= SSLIMIT)
		return SchoenhageStrassenMultiplication(rhs);
	else
#endif
		return SchoolbookMultiplication(rhs);

}

pIntClass pIntClass::SchoolbookMultiplication(const pIntClass& rhs)
{
	std::vector<s64> buffer;
	std::vector<int> acc;
	if (value.size() && rhs.value.size()) {
		int mysign =     (value.back() >= 0) ? 1 : -1;
		int rhsign = (rhs.value.back() >= 0) ? 1 : -1;

		buffer.reserve(value.size() + rhs.value.size() + 1);
		acc.reserve(value.size() + rhs.value.size() + 1);

		for (int i1 = 0; i1 < value.size() + rhs.value.size() + 1; i1++) {
			acc.push_back(0);
			buffer.push_back(0);
		}


		for (int j = 0; j < rhs.value.size(); j++) {
			s64 multiplier = (s64)( (rhsign > 0) ? rhs.value[j] : -rhs.value[j]);
			for (int ix = 0; ix < value.size(); ix++)  
				buffer[ix] = ((s64)( (mysign > 0) ? value[ix] : -value[ix])) * multiplier;
			s64 carry = 0;
			for (int cx = 0; cx < value.size(); cx++) {
				buffer[cx] = buffer[cx] + acc[j+cx]+  carry;
				carry = buffer[cx] / MODULUS;
				acc[j + cx] = buffer[cx] % MODULUS;
				buffer[cx] = 0;
			}
			if (carry) acc[j + value.size()] = (int) carry;
		}

		value.clear();
		for (int vx = 0; vx < acc.size(); vx++) value.push_back(acc[vx]);
		normalize(value);
		if (( (mysign < 0) && (rhsign > 0)) ||	((mysign > 0) && (rhsign < 0)))
			for (int ix = 0; ix < value.size(); ix++) value[ix] = -value[ix];

		buffer.clear();
		acc.clear();
	}
	else
		value.clear();

	return *this;
}




void pIntClass::Scale(int scale)
{
	if (value.size() == 0) 	return; //scaling 0 -> 0
	else if (scale == 0) { value.clear();	return; }  //scaling with 0 -> 0
	else if (scale == 1) return; // scaling with 1 -> no change
	else if (scale == -1) { 
		for(size_t ix = 0; ix < value.size(); ix++ )
			value[ix] = -value[ix] ;	
		return; 
	} //scaling with -1 -> change sign
	else {
		s64 carry = 0;
		s64 t = 0;
		s64 s = scale;
		int sign = (value.back()>= 0) ? 1 : -1 ;  // save sign

		if (s < 0) {
			s = -1 * s;
			sign = -sign ;
		}

		for (int i = 0; i < value.size(); i++)
		{
			t = (s64)value[i];
			t *= s;
			t += carry;
			carry = 0;
			if (t >= MODULUS) {
				carry = (int)(t / MODULUS);
				t = t % MODULUS;
			}
			value[i] = (int)t;
		}

		while (carry) {
			value.push_back(carry % MODULUS);
			carry = carry / MODULUS;
		}
		if(sign == -1)
			for (size_t ix = 0; ix < value.size(); ix++)
				value[ix] = -value[ix];

		normalize(value);
		return ;
	}
}



pIntClass& pIntClass::operator<<=(const unsigned int shift)
{
	unsigned int _shift = shift;
	if (value.size()) {
		int sign = (value.back()) >= 0 ? 1 : -1;

		if (sign == -1)
			for (size_t i = 0; i < value.size(); i++)
				value[i] = -value[i];

		while (_shift--) {
			int carry = 0;
			for (size_t i = 0; i < value.size(); i++)
			{
				value[i] = value[i] << 1;
				value[i] += carry;
				if (value[i] >= MODULUS) {
					value[i] -= MODULUS;
					carry = 1;
				}
				else
					carry = 0;
			}
			if (carry) value.push_back(1);
		}

		if (sign == -1)
			for (size_t i = 0; i < value.size(); i++)
				value[i] = -value[i];
	}
	return *this;
}

pIntClass& pIntClass::operator>>=(const unsigned int shift)
{
	unsigned int _shift = shift;
	if (_shift && value.size()) {
		int sign = (value.back()) >= 0 ? 1 : -1;

		if (sign == -1)
			for (size_t i = 0; i < value.size(); i++)
				value[i] = -value[i];

		while (_shift--) {
			int carry = 0;
			for (size_t i = value.size(); i > 0; i--)
			{
				if (carry)  value[i - 1] += MODULUS;
				carry = (value[i - 1] & 1);
				value[i - 1] = value[i - 1] >> 1;
			}
			if (value.back() == 0) value.pop_back();
		}

		if (sign == -1)
			for (size_t i = 0; i < value.size(); i++)
				value[i] = -value[i];
		normalize(value);
	}
	return *this;
}

#ifdef SSLIMIT

#define OVERALLOCATION 2

//#define PREALLOC
#ifdef PREALLOC
#define PREALLOCE 2000000
Data* real1 = new Data[PREALLOCE];
Data* imag1 = new Data[PREALLOCE];
Data* real3 = new Data[PREALLOCE];
Data* imag3 = new Data[PREALLOCE];
#endif


pIntClass pIntClass::SchoenhageStrassenMultiplication(const  pIntClass& rhs)
{

	PrimeFactorDFT pf;
	

	int MySign = value.back() >= 0 ? 1 : -1;
	int rhsSign = rhs.value.back() >= 0 ? 1 : -1;

	factorSeq  factors;

	u64 min_sz = value.size() +  rhs.value.size();

	pf.FastCalcFactors((unsigned int) (3 /** 2*/ * min_sz), factors);

	pf.SetFactors(factors);
	//std::cout << "FFT length " << pf.Status() << std::endl;

	if (pf.Status() > 0) {
		/* this should be less dynamic, but right now it is OK */
#ifndef PREALLOC
		Data* real1 = new Data[pf.Status() + OVERALLOCATION];
		Data* imag1 = new Data[pf.Status() + OVERALLOCATION];
		Data* real3 = new Data[pf.Status() + OVERALLOCATION];
		Data* imag3 = new Data[pf.Status() + OVERALLOCATION];
#endif
		for (s64 i = 0; i < pf.Status() + OVERALLOCATION; i++) {
			real1[i] = 0;
			imag1[i] = 0;
			real3[i] = 0;
			imag3[i] = 0;
		}

		LoadFFT(value, real1);
		LoadFFT(rhs.value, imag1);

		pf.forwardFFT(real1, imag1);

		real3[0] = real1[0] * imag1[0];
		s64 size = pf.Status();
		for (s64 i = 1; i < size; i++)
		{
			Data X01Real = real1[i];
			Data X01Imag = imag1[i];
			Data X02Real = real1[size - i];
			Data X02Imag = imag1[size - i];
			Data X1Real = (X01Real + X02Real) / 2.0;
			Data X1Imag = (X01Imag - X02Imag) / 2.0;
			Data X2Imag = -1 * (X01Real - X02Real) / 2;
			Data X2Real = (X01Imag + X02Imag) / 2;
			Data X3Real = X1Real * X2Real - X1Imag * X2Imag;
			Data X3Imag = X1Real * X2Imag + X1Imag * X2Real;
			real3[i] = X3Real;
			imag3[i] = X3Imag;
		}
		for (s64 i = 0; i < size; i++)
		{
			real1[i] = real3[i];
			imag1[i] = imag3[i];
		}

		pf.ScaledInverseFFT(real1, imag1);

		Carry(size, real1);

		/* convert back to radix 10^9 from radix 10^3 double*/
		value.clear();
		for (int i = 0; i < pf.Status(); i += 3)
		{
			int t0 = (int)real1[i];
			int t1 = 1000 * (int)real1[i + 1];  // these values are 0 when we reach
			int t2 = 1000 * 1000 * (int)real1[i + 2];// the end of buffer, due to
			value.push_back(t0 + t1 + t2);//OVERALLOCATION
		}
#ifndef PREALLOC
		delete[] real1;
		delete[] imag1;
		delete[] real3;
		delete[] imag3;
#endif

		while (value.size() && (value.back() == 0)) value.pop_back();
		if((MySign * rhsSign) < 0) /* negative result*/
			for( size_t ix = 0 ; ix < value.size(); ix++) 
				value[ix] =  -value[ix];
		
		return *this;
	}
	return *this;
}

#define RMOD3 1000

void pIntClass::LoadFFT(const std::vector<int>& A, double* Buffer)
{
	int carry = 0;
	int FFTIndex = 0;
	for (int ix = 0; ix < A.size(); ix++) {
		int temp = (A.back() >= 0) ? A[ix] : -A[ix];   // temp is  radix 10^9 
		for (int i = 0; i < 3; i++) {
			int tmp = (temp % RMOD3);// tmp is  radix 10^3
			temp = temp / RMOD3;
			tmp = tmp + carry; carry = 0;
			if (tmp > (RMOD3 / 2) - 1) {  // tmp is 'balanced' radix 10^3 
				tmp = tmp - RMOD3;
				carry++;
			}
			Buffer[FFTIndex++] = (double)tmp;
		}
	}
	if (carry)
		Buffer[FFTIndex] = 1.0;
}

void pIntClass::Carry(s64 size, double* Buffer)
{
	s64 carry = 0;
	s64 tmp = 0;

	/* conversion from 'balanced' Radix 10^3 double  to unbalanced Radix 10^3 double */
	for (s64 i = 0; i < size; i++)
	{
		tmp = (s64)std::round(Buffer[i]);

		tmp += carry; carry = 0;
		while (tmp < (-1 * (RMOD3 / 2))) { tmp += RMOD3; carry--; }
		while (tmp > ((RMOD3 / 2) - 1)) { tmp -= RMOD3; carry++; }
		Buffer[i] = (double)tmp;
	}
	carry = 0;
	/* carry we are still in Radix 10^3 double*/
	for (s64 i = 0; i < size; i++)
	{
		tmp = (s64)std::round(Buffer[i]);
		tmp += carry;  carry = 0;
		if (tmp < 0) { tmp += RMOD3; carry--; }
		Buffer[i] = double(tmp);
	}

}

#endif

pIntClass& pIntClass::operator*=(const int rhs)
{
	pIntClass temp(rhs);
	*this *= temp;
	return *this;
}


pIntClass operator*(const pIntClass& a, const pIntClass& b)
{
	pIntClass temp;// = new pIntClass();
	temp = a;
	temp *= b;
	return temp;
}

bool operator<(const pIntClass& a, const pIntClass& b)
{
	if (a.IsPos() && b.IsPos()) { // a >= 0 and b >= 0
		if (a.IsZero() &&  b.IsZero())  return false;
		if (a.IsZero() && !b.IsZero())  return true;
		if (!a.IsZero() && b.IsZero())  return false;
		if (!a.IsZero() && !b.IsZero())
		{ // same sign !
			if (a.value.size() < b.value.size()) return true;
			if (a.value.size() > b.value.size()) return false;
		  // same size !
			for (size_t i = a.value.size(); i > 0; i--)
				if (a.value[i - 1] == b.value[i - 1]) continue;
				else if (a.value[i - 1] > b.value[i - 1]) return false;
				else return true;
			return false;
		}

	}
	else if (a.IsNeg() && b.IsPos()) return true;
	else if (a.IsPos() && b.IsNeg()) return false;
	else // a < 0 amd b < 0 
	{
		if (a.value.size() < b.value.size()) return false;
		if (a.value.size() > b.value.size()) return true;
		// same size !
		for (size_t i = a.value.size(); i > 0; i--)
			if (a.value[i - 1] == b.value[i - 1]) continue;
			else if (a.value[i - 1] < b.value[i - 1]) return false;
			else return  true;
		return false;
	}
	return false;
}

bool operator==(const pIntClass& a, const pIntClass& b)
{
	if (a.value.size() != b.value.size()) return false;
	if (a.value.size() == 0) return true;
	for (size_t i = 0; i < a.value.size(); i++)
		if (a.value[i] != b.value[i]) return false;
	return true;
}

bool operator!=(const pIntClass& a, const pIntClass& b) {	return !operator==(a, b);  }


pIntClass operator*(const int a, const pIntClass& b)
{
	pIntClass temp;// = new pIntClass();
	temp = b;
	(temp) *= a;
	return temp;
}



void testMultiplication()
{

	pIntClass a("2628461924971");

	pIntClass b;

	std::string s;
	b = 7 * a;
	s = b.ToString();
	std::cout << " b = 7 * a " << s << std::endl;

	s.clear();
#if 0
	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = a * 7;
	std::cout << " b = a * 7 " << b.ToString() << std::endl;


	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = 10 * 7;
	std::cout << " b = 10 * 7 " << b.ToString() << std::endl;

	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = a * 7;
	b.ChSignBit();
	std::cout << " b = -(a * 7) " << b.ToString() << std::endl;

#endif

}
