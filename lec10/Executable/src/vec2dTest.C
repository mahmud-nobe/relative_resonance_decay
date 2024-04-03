#include <iostream>
#include "TROOT.h" //always needed
#include "TStyle.h" //usually needed
#include "TH2.h"
#include "../../Lib/vec2d.h"

using namespace std;

void vec2dTest()
{
  // Load our library into
  // the ROOT system to be
  // usable! Not needed if compiled!
  //gSystem->Load("../Lib/libVec2d.so");
  
  vec2d *v=new vec2d(2,2);
  v->Print();
  cout<<v->Length()<<endl;       
    
  delete v;
}

int main()
{
  vec2dTest();
}
