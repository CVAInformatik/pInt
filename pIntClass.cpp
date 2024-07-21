
/*
Copyright  © 2024 Claus Vind - Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/

#include <iostream>
#include "pIntClass.h"




/******************************************************************************

       Below is (some of)  the implementation of pIntClass

****************************************************************************/


pIntClass::pIntClass()
{
	value.clear();
	value.reserve(ReservationSize);
}

pIntClass::pIntClass(const pIntClass& x)
{
	value.clear();
	value.reserve(ReservationSize);
	for (int i = 0; i < x.value.size(); i++) value.push_back(x.value[i]);
}

pIntClass::pIntClass(const int x)
{
	int _x = x;
	value.clear();
	value.reserve(ReservationSize);
	if (x < 0) _x = -_x;
	for (; _x; _x = _x / MODULUS) {
			value.push_back(_x % MODULUS);
	}
	if (x < 0)
		for (size_t ix = 0; ix < value.size();ix++)
			value[ix] = -value[ix];
	normalize(value);
}


/*
*     A pIntClass number is normalized if
*      value.size() == 0        -> meaning the number has the value 0
*      all non-zero entries in value[] have the same sign and is in the open
*                interval ] -modulus,...., modulus[
*      (it means that the sign of a non-zero pIntClass2 number is the sign of
*          value.back(), the most significant 'digit' non-zero entry).
*
*	   Addition and Subtraktion may temporarily make a value not normalized,
*	   normalize will bring it back in order.
*/
void pIntClass::normalize(std::vector<int>& val) {
	// drop leading zeroes
	while (val.size() && (val.back() == 0)) val.pop_back();
	if (val.size()) {
		/* low to high negative carry */
		int carry = 0;
		for (size_t i = 0; i < val.size(); i++)
		{
			val[i] += carry;
			if (val[i] >= MODULUS) {
				carry = 1;
				val[i] -= MODULUS;
			}
			else
				carry = 0;
		}
		if (carry) val.push_back(carry);
		carry = 0;
		/* low to high positive carry */
		for (size_t i = 0; i < val.size(); i++)
		{
			val[i] += carry;
			if (val[i] <= -MODULUS) {
				carry = -1;
				val[i] += MODULUS;
			}
			else
				carry = 0;
		}
		if (carry) val.push_back(carry);
		/* all entries in value now is in ]-modulus...modulus[ */
		/* now we need to make sure they have the same sign     */
		carry = 0;
		if (val.back() < 0) { /* negative numberss*/
			for (size_t i = 0; i < val.size(); i++) {
				val[i] += carry;
				carry = 0;
				if (val[i] > 0) {
					val[i] -= MODULUS;
					carry = 1;
				}
			}
			if (carry) val.push_back(carry);
		}
		else if (val.back() > 0) { /* positive numbers */
			for (size_t i = 0; i < val.size(); i++) {
				val[i] += carry;
				carry = 0;
				if (val[i] < 0) {
					val[i] += MODULUS;
					carry = -1;
				}
			}
			if (carry) val.push_back(carry);
		}
		while (val.back() == 0) {
			val.pop_back();
		}
	}
}



pIntClass::~pIntClass() 
{
	value.clear();
}



pIntClass& pIntClass::operator=(const pIntClass& rhs)
{
	if ( this != &rhs)
	{
		value.clear();
		for (int i = 0; i < rhs.value.size(); i++) value.push_back(rhs.value[i]);
	}
	return *this;
}


pIntClass pIntClass::operator=(const int  rhs)
{
	value.clear();
	s64 tmp = 0;
	if (rhs < 0) {
		tmp = -((s64) rhs);
	}
	else
		tmp = rhs;
	do {
		if(rhs < 0)
			value.push_back(-(tmp % MODULUS));
		else
			value.push_back(tmp % MODULUS);
		tmp = tmp / MODULUS;
	} while (tmp > 0);
	normalize(value);
	return *this;
}



bool pIntClass::operator!=( const pIntClass& b) {	return !( *this == b);  }

bool pIntClass::operator!=(const int b) { return !(*this == b); }



bool pIntClass::IsBiggerNummerically(const pIntClass& b)
{
	if (IsZero() && b.IsZero()) return false;
	if (IsZero() && !b.IsZero()) return false;
	if (!IsZero() && b.IsZero()) return true;
	if (value.size() > b.value.size()) return true;
	if (value.size() < b.value.size()) return false;
	/* this and b are the same size, but may differ in signs */
	for (size_t i = value.size(); i > 0; i--)
	{
		int a  =  ( value.back() >= 0) ?   value[i-1] :   -value[i-1];
		int _b = (b.value.back() >= 0) ? b.value[i-1] : -b.value[i-1];
		if (a > _b) return true ;
		else if (a < _b) return false;
	}
	return false;
}

void pIntClass::DivModulus() {
	switch (value.size())
	{
	case 0: 		break;
	default:
		for (int ix = 0; value.size() && (ix < (value.size() - 1)); ix++)
		{
			value[ix] = value[ix + 1];
		}
		value.pop_back();
		break;
	}
}


int pIntClass::operator[](int index)
{
	if (index < value.size())
		return value[index];  // we strip the sign !
	else 
		return 0;  // elements of value[] are all positive !
}

int pIntClass::Sign()
{
	if (value.size() == 0)  return 0;// 0 is positive by definition
	else  return value.back() > 0 ? 1 : -1;
}

int pIntClass::ChSignBit()
{
	if (value.size() == 0) 	return 0;
	else
		for (size_t i = 0; i < value.size(); i++) value[i] = -value[i];
	return (value.back() >= 0) ? 1 : -1;
}

pIntClass& pIntClass::operator++()
{
	if (value.size() == 0) value.push_back(1);
	else {
		value[0]++;
		normalize(value);
	}
	return *this;
}


pIntClass pIntClass::operator++(int)
{
	pIntClass res = *this;
	++*this;
	return res;
}

pIntClass& pIntClass::operator--()
{

	if (value.size() == 0) value.push_back(-1);
	else {
		value[0]--;
		normalize(value);
	}
	return *this;
}

pIntClass pIntClass::operator--(int)
{
	pIntClass res = *this;
	--*this;
	return res;
}
