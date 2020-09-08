#include <iostream>
#include <fstream>
#include <cstdlib>

#include "utilities.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  
  try
  {
    if(argc <2)
    {
      throw("# No parameter file is specified.");
    }
  }
  catch(const char *s)
  {
    cout<<s<<endl;
    exit(1);
  }
  
  cali = new Cali("data/ngc5548_cont.txt", "data/ngc5548_line.txt");
  
  cali->mcmc();
  cali->get_best_params();
  cali->align_with_error();
  cali->output();
  cali->recon();

  return EXIT_SUCCESS;
}