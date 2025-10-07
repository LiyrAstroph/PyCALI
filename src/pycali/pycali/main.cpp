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
  //cali.mcmc();
  cali.get_best_params();
  cali.output();
  cali.recon();
  
  cout<<"================================"<<endl;
  cout<<"One can execute 'python plot_for_cali.py param.txt' to plot the result,\n"
        "which output a file 'PyCALI_results.pdf'. "<<endl;

  return EXIT_SUCCESS;
}