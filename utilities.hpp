#ifndef _UTILITIES_HPP

#define _UTILITIES_HPP

#include <string>
#include <vector>

#include <dnestvars.h>

#define WhiteSpace " \t\v\r"

enum PRIOR_TYPE {GAUSSIAN=1, UNIFORM=2};
enum PAR_FIX {NOFIXED=false, FIXED=true};

using namespace std;

double prob_cali(const void *model);
void from_prior_cali(void *model);
void print_particle_cali(FILE *fp, const void *model);
double perturb_cali(void *model);

/* 
 * Data class for light curves.
 */
class Data 
{
  public:
    Data();
    Data(const string& fname);
    ~Data();
    void load(const string& fname);
    void normalize();
    void sort_data();
    /* variables */
    double norm;
    vector<double> flux_org, error_org;
    vector<double> time, flux, error;
    vector<int> code;
    vector<int> index;
    vector<int> num_code;
    vector<string> code_list;
};

class Cali
{
  public:
    Cali();
    Cali(const string& fcont, const string& fline="");
    ~Cali();
    void mcmc();
    void align(double *model);
    void get_best_params();
    

    Data cont, line;
    size_t size_max;
    size_t ncode;

    int num_params;
    int num_params_var;
    double **par_range_model;
    bool *par_fix;
    double *par_fix_val;
    int *par_prior_model;
    double **par_prior_gaussian;
    double *best_params, *best_params_std;

    double *workspace;

    DNestFptrSet *fptrset;
};

extern Cali *cali;
#endif