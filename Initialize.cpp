#include <iostream>
#include <math.h>
#include "dlfcn.h"
#include "mex.h"

#include <N_DEV_ADMSbsimcmg.h>

using namespace std;

void Initialize();

typedef void (*Initializeance)();

void mexFunction(
		 int          nlhs,
		 mxArray      *[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
  Initialize();
}
