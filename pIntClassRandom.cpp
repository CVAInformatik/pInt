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
#include "pIntClassRandom.h"

#define DIGITS 9 

pIntClassRandom::pIntClassRandom(unsigned int _Modulus) 
{
	Modulus = _Modulus;
	dist = new std::uniform_int_distribution<unsigned int>(0, Modulus - 1);
}

pIntClassRandom::~pIntClassRandom() {
	delete dist;
}

void pIntClassRandom::SetSeed(const pIntClass& seed) {
	sz = DIGITS * ((int) seed.value.size() - 1);
	s = seed.value.size() ? seed.value.back() : 0;
	if (s == 0) {
		std::cout << "Zero Argument " << std::endl;
		return;
	}
	int t = s;
	while (t > 0) { t = t / 10; sz++; } //count the numbers of digits in the MSInt


}

pIntClass& pIntClassRandom::Rand() {

	pIntClass* res = new pIntClass();

	Rand(*res);

	return *res;

}

void pIntClassRandom::Rand(pIntClass& res) {

	int sztemp = sz;
	while (sztemp > 0) {
		// approx 90% of all numbers smaller than sz is  in range sz ... sz/10
		// approx 99% of all numbers smaller than sz is  in range sz ... sz/100
		// approx 99.9% of all numbers smaller than sz is  in range sz ... sz/1000.....
		if (_Rand(Modulus) >= (Modulus / 10)) break;
		sztemp--;
	}
	int size = std::max(1, sztemp);

	res.value.clear();
	for (; size >= DIGITS; size = size - DIGITS)
		res.value.push_back(_Rand(Modulus)); //  _Rand() returns an integer in the range 0..RMOD-1


	int i = _Rand(Modulus);

	while (DIGITS > size) {
		i = i / 10;
		size++;
	}
	i = i % s;

	res.value.push_back(i);

	while (res.value.size() && (res.value.back() == 0)) res.value.pop_back();
}

unsigned int  pIntClassRandom::_Rand(unsigned int UpperBound)
{
	unsigned int ix2 = dist->operator()(rd);
	return ix2 % UpperBound;
}

