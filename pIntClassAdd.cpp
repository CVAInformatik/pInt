
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

pIntClass& pIntClass::operator+=(const pIntClass& rhs)
{
	size_t sz = std::min(value.size(), rhs.value.size());

	size_t i = 0;
	for (; i < sz; i++)	value[i] += rhs.value[i];

	while (i < rhs.value.size())  value.push_back(rhs.value[i++]);

	normalize(value);
	return *this;
}

pIntClass& pIntClass::operator-=(const pIntClass& rhs)
{
	size_t sz = std::min(value.size(), rhs.value.size());

	size_t i = 0;
	for (; i < sz; i++)	value[i] -= rhs.value[i];

	while (i < rhs.value.size())  value.push_back(-rhs.value[i++]);

	normalize(value);

	return *this;
}



pIntClass& pIntClass::operator-=(const int rhs)
{
	pIntClass temp;
	temp = rhs;
	temp.ChSignBit();
	*this += temp;
	return *this;
}


pIntClass& pIntClass::operator+=(const int rhs)
{
	pIntClass temp(rhs);
	*this += temp;
	return *this;

}


pIntClass& operator+(pIntClass& a, const pIntClass& b)
{
	pIntClass* t = new pIntClass();
	*t = a;
	*t += b;
	return *t;
}

pIntClass& operator-(pIntClass& a, const pIntClass& b)
{
	pIntClass* t = new pIntClass();
	*t = a;
	*t -= b;
	return *t;
}


pIntClass& operator+(const int a, const pIntClass& b)
{
	pIntClass* t = new pIntClass();
	*t = a;
	*t += b;
	return *t;
}

pIntClass& operator-(const int a, const pIntClass& b)
{
	pIntClass* t = new pIntClass();
	*t = a;
	*t -= b;
	return *t;
}

void testAddition()
{

	pIntClass a("2628461924971");

	pIntClass b;

	b = 7 + a;
	std::cout << " b = 7 + a " << b.ToString() << std::endl;

	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = a + 7;
	std::cout << " b = a + 7 " << b.ToString() << std::endl;


	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = 10 + 7;
	std::cout << " b = 10 + 7 " << b.ToString() << std::endl;

	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = a + 7;
	b.ChSignBit();
	std::cout << " b = -(a + 7) " << b.ToString() << std::endl;
}

void testSubtraction()
{

	pIntClass a("2628461924971");

	pIntClass b;

	b = 7 - a;
	std::cout << " b = 7 - a " << b.ToString() << std::endl;

	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = a - 7;
	std::cout << " b = a - 7 " << b.ToString() << std::endl;


	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = 10 - 7;
	std::cout << " b = 10 - 7 " << b.ToString() << std::endl;

	b = pIntClass("99999999999999999");

	std::cout << " b " << b.ToString() << std::endl;

	b = a - 7;
	b.ChSignBit();
	std::cout << " b = -(a - 7) " << b.ToString() << std::endl;
}

