#ifndef _UTILITIES_HPP

#define _UTILITIES_HPP

#include <string>
#include <vector>
#include <list>

#include "../cdnest/dnestvars.h"

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
    void setup(const string& fcont, const list<string>& fline=list<string>({}), 
             int nmcmc=10000, double ptol=0.1, 
             double scale_range_low= 0.5, double scale_range_up=1.5,
             double shift_range_low= -1.0, double shift_range_up= 1.0,
             double syserr_range_low = 0.0, double syserr_range_up = 0.1,
             double errscale_range_low=0.1, double errscale_range_up = 2.0,
             double sigma_range_low= 1.0e-4, double sigma_range_up=1.0,
             double tau_range_low = 1.0, double tau_range_up = 1.0e4,
             bool fixed_scale = false, bool fixed_shift = false,
             bool fixed_syserr=true, bool fixed_error_scale=true,
             const vector<int>& fixed_codes=vector<int>({}),
             const vector<int>& fixed_scalecodes=vector<int>({}),
             bool flag_norm=true);
    string get_param_filename();
    void check_directory();
    void print_cfg();
    void parse_fline_str(const string& fline_str);
    void parse_fixed_codes_str(const string& fixed_codes_str);
    void parse_fixed_scalecodes_str(const string& fixed_scalecodes_str);
    vector<double> test(const vector<double>& range = vector<double>({0.5, 1.5}));

    string fname;
    size_t nmcmc;
    double ptol;
    double scale_range_up, scale_range_low;
    double shift_range_up, shift_range_low;
    double syserr_range_up, syserr_range_low;
    double errscale_range_up, errscale_range_low;
    double sigma_range_up, sigma_range_low;
    double tau_range_up, tau_range_low;
    char *fcont;
    list<string> fline;

    bool fixed_scale, fixed_shift;
    bool fixed_syserr, fixed_error_scale;
    bool flag_norm;
    vector<int> fixed_codes;  /* fix some specific codes */
    vector<int> fixed_scalecodes; /* fix scale of some specific codes */
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
    Data(const string& fname, bool flag_norm);
    ~Data();
    void load(const string& fname, bool flag_norm);
    void normalize();
    void sort_data();
    void check_code(Data& data);
    /* variables */
    double norm;
    vector<double> flux_org, error_org;
    vector<double> time, flux, error;
    vector<int> syserrflag;
    vector<int> syserrflag_list;
    vector<int> code;
    vector<int> index;
    vector<int> num_code;
    vector<int> num_flag;
    vector<double>mean_code;
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
    void align_cont(double *model);
    void align_line(double *model, int il);
    void align_with_error();
    void get_best_params();
    void output();
    void recon();
    void set_covar_Umat_cont(double sigma, double tau, double *USmat);
    void set_covar_Umat_line(double sigma, double tau, double *USmat, unsigned int il);
    double get_norm_cont();
    double get_norm_line(unsigned int il);
    void check_directory();
    void check_fixed_codes(Config& cfg);
    void check_fixed_scalecodes(Config& cfg);

    string fcont;
    list<string> fline;
    Data cont;
    list<Data> lines;
    size_t size_max;
    size_t ncode;
    size_t nsyserr_flag;

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
    double ptol;
    /* reconstruction */
    DataLC cont_recon;
    list<DataLC> lines_recon;
    size_t size_recon_max;

    DNestFptrSet *fptrset;

    int stat_type;
};

#endif