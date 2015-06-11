#include <iostream>
#include <math.h>
#include "dlfcn.h"
#include "mex.h"

using namespace std;

void BSIMCMG(double Vd, double Vg, double Vs, double Ve);

extern void _main();

/****************************/
class MyData {

public:
  void display();
  void set_data(double v1, double v2);
  MyData(double v1 = 0, double v2 = 0);
  ~MyData() { }
private:
  double val1, val2;
};

MyData::MyData(double v1, double v2)
{
  val1 = v1;
  val2 = v2;
}

void MyData::display()
{
#ifdef _WIN32
	mexPrintf("Value1 = %g\n", val1);
	mexPrintf("Value2 = %g\n\n", val2);
#else
  cout << "Value1 = " << val1 << "\n";
  cout << "Value2 = " << val2 << "\n\n";
#endif
}

void MyData::set_data(double v1, double v2) { val1 = v1; val2 = v2; }

/*********************/

static
void mexcpp(
	    double num1,
	    double num2
	    )
{
#ifdef _WIN32
	mexPrintf("\nThe initialized data in object:\n");
#else
  cout << "\nThe initialized data in object:\n";
#endif
  MyData *d = new MyData; // Create a  MyData object
  d->display();           // It should be initialized to
                          // zeros
  d->set_data(num1,num2); // Set data members to incoming
                          // values
#ifdef _WIN32
  mexPrintf("After setting the object's data to your input:\n");
#else
  cout << "After setting the object's data to your input:\n";
#endif
  d->display();           // Make sure the set_data() worked
  delete(d);
  flush(cout);
  return;
}

void mexFunction(
		 int          nlhs,
		 mxArray      *[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
  double      *vin1, *vin2;

  /* Check for proper number of arguments */

  if (nrhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "MEXCPP requires two input arguments.");
  } else if (nlhs >= 1) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "MEXCPP requires no output argument.");
  }
 
  void* lib_handle = dlopen("/mnt/hgfs/kminhsu/Desktop/Xyce-Build/Xyce-6.2-build/src/DeviceModelPKG/test/libtest_dyn.so", RTLD_LAZY); 
  cout << "error:" << dlerror() << endl;
  cout << lib_handle << endl; 

  BSIMCMG(1.0, 0.3, 0.0, 0.0 );

  vin1 = (double *) mxGetPr(prhs[0]);
  vin2 = (double *) mxGetPr(prhs[1]);

  mexcpp(*vin1, *vin2);
  return;
}
