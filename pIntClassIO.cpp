
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
#include <algorithm>
#endif
#include "pIntClass.h"


pIntClass::pIntClass(const std::string& x)
{
	value.clear();
	value.reserve(ReservationSize);

	int index1 = 0;
	int sign = 1; //default positive

	while (index1 < x.length() && isspace(x[index1])) index1++;

	if (x[index1] == '-') {
		sign = -1;
		index1++;
	}

	int index2 = index1;

	while (index2 < x.length() && isdigit(x[index2])) index2++;

	if (index2 == index1)
		std::cerr << "unknown format: " << x << std::endl;
	// we assume we have something number like 
	value.push_back(0);
	while (index1 < x.length() && isdigit(x[index1])) {
		mul10();
		value[0] += x[index1] - '0';
		index1++;
	}
	normalize(value);
	for (int i = 0; i < value.size(); i++) value[i] *= sign;

}

std::string pIntClass::ToString()
{
#define FORMATSTRING "%09d"
#define DIGITS 9

	if (value.size() == 0) {
		std::string* s = new std::string();
		s->append("0");
		return *s;
	}
	else {
		char buffer[12];
		std::string s;// = new std::string();
		int sign = value.back();

		for (int i = 0; i < value.size(); i++) {
			sprintf(buffer, FORMATSTRING, (sign<0 )? -value[i] : value[i]);
			char* c1 = buffer;
			char* c2 = buffer + DIGITS - 1;
			while (c1 < c2)
			{
				char t = *c1;
				*c1 = *c2;
				*c2 = t;
				c1++; c2--;
			}
			s.append(buffer);
		}
		while (s.size() && (s.back() == '0')) s.pop_back();
		if (s.size() == 0) 
			s.append("0");
		else 
		    if (sign < 0)  s.append("-");
		std::reverse(s.begin(), s.end());
		return s;
	}
}
