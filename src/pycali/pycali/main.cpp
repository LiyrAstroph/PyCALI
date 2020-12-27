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

  Config cfg(argv[1]);
  cfg.print_cfg();

  Cali cali(cfg);
  cali.mcmc();
  cali.get_best_params();
  cali.output();
  cali.recon();

  return EXIT_SUCCESS;
}