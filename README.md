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
*      and the most significant (last) digit is non-zero.
*
*      that means that the sign of a non-zero pIntClass number is the sign of the
*          value.back(), the last non-zero entry).
*
*	   Addition and Subtraction may temporarily make a value not normalized,
*	   normalize() will bring it back in order.
* 
*/
  it comppiles in a Visual Studio 2022 x64 Console App project, define OS_WINDOWS and _CRT_SECURE_NO_WARNINGS for the entire project, compile for 'release' for speed.
  
  it compiles and run in WSL/Ubuntu with the attached makefile, the MRtest is killed after 300-400 iterations,   I suspect a leak, but otherwise it appears to work OK
  
ADDED: Looks like I have found most of the leaks, it now runs the full test to end without being killed...
