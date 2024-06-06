
CC = g++
CFLAGS = -g 
CPPFLAGS =  -O1  

%.o  :  %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@


clean:
	rm *.o


PrimeFactorDFT.o : PrimeFactorDFT.cpp  PrimeFactorDFT.h


pInt.o : pInt.cpp
PrimeFactorDFT.o :   PrimeFactorDFT.cpp 
pIntClass.o  :       pIntClass.cpp pIntClass.h
pIntClassAdd.o  :    pIntClassAdd.cpp pIntClass.h
pIntClassIO.o  :     pIntClassIO.cpp pIntClass.h
pIntClassMultiply.o :   pIntClassMultiply.cpp pIntClass.h
pIntClassRandom.o  : pIntClassRandom.cpp pIntClassRandom.h pIntClass.h
pIntClassUtil.o :    pIntClassUtil.cpp pIntClassUtil.h pIntClass.h

pInt :  pInt.o PrimeFactorDFT.o pIntClass.o pIntClassAdd.o pIntClassIO.o pIntClassMultiply.o pIntClassRandom.o pIntClassUtil.o


