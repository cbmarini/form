extern "C" {
	#include "form3.h"
}

#include <iostream>
using namespace std;
// #include<ginac/ginac.h>

int EvaluateLin(PHEAD WORD *term, WORD level, WORD par) {
	DUMMYUSE(par);
	printf("Hello from the GiNaC wrapper\n");
	return(Generator(BHEAD term, level));
}
