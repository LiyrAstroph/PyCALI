#ifndef _UTILITIES_HPP

#define _UTILITIES_HPP

#include <string>
#include <vector>

#include "./cdnest/dnestvars.h"

#define WhiteSpace " \t\v\r"

enum PRIOR_TYPE {GAUSSIAN=1, UNIFORM=2, LOG=3};
enum PAR_FIX {NOFIXED=false, FIXED=true};

using namespace std;

double prob_cali(const void *model, const void *arg);
void from_prior_cali(void *model, const void *arg);
void print_particle_cali(FILE *fp, const void *model, const void *arg);
double perturb_cali(void *model, const void *arg);

class Config;
class DataLC;
class Data;
class Cali;

class Config
{
  public:
    Config();
    Config(const string& fname);
    ~Config();
    void load(const string& fname);
    void setup(const string& fcont, const string& fline="", 
             int nmcmc=10000, double pdiff=0.1, 
             double scale_range_low= 0.5, double scale_range_up=1.5,
             double shift_range_low= -1.0, double shift_range_up= 1.0,
             double syserr_range_low = 0.0, double syserr_range_up = 0.1,
             double errscale_range_low=0.1, double errscale_range_up = 2.0,
             double sigma_range_low= 1.0e-4, double sigma_range_up=1.0,
             double tau_range_low = 1.0, double tau_range_up = 1.0e4,
             bool fixed_scale = false, bool fixed_shift = false,
             bool fixed_syserr=true, bool fixed_error_scale=true);
    string get_param_filename();
    void print_cfg();
    vector<double> test(const vector<double>& range = {0.5, 1.5});

    string fname;
    size_t nmcmc;
    double pdiff;
    double scale_range_up, scale_range_low;
    double shift_range_up, shift_range_low;
    double syserr_range_up, syserr_range_low;
    double errscale_range_up, errscale_range_low;
    double sigma_range_up, sigma_range_low;
    double tau_range_up, tau_range_low;
    char *fcont, *fline;

    bool fixed_scale, fixed_shift;
    bool fixed_syserr, fixed_error_scale;
};

class DataLC
{
  public:
    DataLC();
    DataLC(size_t n);
    ~DataLC();
    void resize(size_t n);

    vector<double> time, flux, error;
};

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
    void check_code(Data& data);
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
    Cali(Config& cfg);
    ~Cali();
    void mcmc();
    void align(double *model);
    void align_with_error();
    void get_best_params();
    void output();
    void recon();
    void set_covar_Umat_cont(double sigma, double tau, double *USmat);
    void set_covar_Umat_line(double sigma, double tau, double *USmat);
    double get_norm_cont();
    double get_norm_line();

    string fcont, fline;
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
    double *best_params, *best_params_std, *best_params_covar;

    double *Larr_data;
    double *workspace;

    size_t nmcmc;
    double pdiff;
    /* reconstruction */
    DataLC cont_recon, line_recon;
    size_t size_recon_max;

    DNestFptrSet *fptrset;
};

#endif