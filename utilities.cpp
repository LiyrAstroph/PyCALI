
#include <iostream>
#include <iomanip>
#include <fstream> 
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <cblas.h>
#include <float.h>

#include "utilities.hpp"
#include "mathfun.h"

/*=====================================================*/
Config::Config()
{
  nmcmc = 10000;
  ptol = 0.1;
  scale_range_low = 0.5;
  scale_range_up = 1.5;
  shift_range_low = -1.0;
  shift_range_up = 1.0;
  syserr_range_low = 0.0;
  syserr_range_up = 0.1;
  errscale_range_low = 0.1;
  errscale_range_up = 2.0;
  sigma_range_low = 1.0e-4;
  sigma_range_up = 1.0;
  tau_range_low = 1.0;
  tau_range_up = 1.0e4;

  fixed_scale = false;
  fixed_shift = false;
  fixed_syserr= true;
  fixed_error_scale = true;

  fcont = new char [256];
  fline = new char [256];
  strcpy(fcont, "\0");
  strcpy(fline, "\0");
}
Config::Config(const string& fname)
      :fname(fname)
{
  nmcmc = 10000;
  ptol = 0.1;
  scale_range_low = 0.5;
  scale_range_up = 1.5;
  shift_range_low = -1.0;
  shift_range_up = 1.0;
  syserr_range_low = 0.0;
  syserr_range_up = 0.1;
  errscale_range_low = 0.1;
  errscale_range_up = 2.0;
  sigma_range_low = 1.0e-4;
  sigma_range_up = 1.0;
  tau_range_low = 1.0;
  tau_range_up = 1.0e4;

  fixed_scale = false;
  fixed_shift = false;
  fixed_syserr= true;
  fixed_error_scale = true;
  
  fcont = new char [256];
  fline = new char [256];
  strcpy(fcont, "\0");
  strcpy(fline, "\0");

  load(fname);
}
Config::~Config()
{
  delete[] fcont;
  delete[] fline;
}
void Config::load(const string& fname)
{
  ifstream fin;

  fin.open(fname);
  if(fin.fail())
  {
    cout<<"cannot open file "<<fname<<endl;
    exit(-1);
  }

  #define MAXTAGS 30
  #define DOUBLE 1
  #define STRING 2
  #define INT 3

  int i, j, nt;
  char str[256], buf1[256], buf2[256], buf3[256];
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];

  nt = 0;
  strcpy(tag[nt], "FileCont");
  addr[nt] = fcont;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileLine");
  addr[nt] = fline;
  id[nt++] = STRING;

  strcpy(tag[nt], "NMcmc");
  addr[nt] = &nmcmc;
  id[nt++] = INT;

  strcpy(tag[nt], "PTol");
  addr[nt] = &ptol;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ScaleRangeLow");
  addr[nt] = &scale_range_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ScaleRangeUp");
  addr[nt] = &scale_range_up;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ShiftRangeLow");
  addr[nt] = &shift_range_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ShiftRangeUp");
  addr[nt] = &shift_range_up;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SyserrRangeLow");
  addr[nt] = &syserr_range_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SyserrRangeUp");
  addr[nt] = &syserr_range_up;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ErrscaleRangeLow");
  addr[nt] = &errscale_range_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ErrscaleRangeUp");
  addr[nt] = &errscale_range_up;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SigmaRangeLow");
  addr[nt] = &sigma_range_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SigmaRangeUp");
  addr[nt] = &sigma_range_up;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "TauRangeLow");
  addr[nt] = &tau_range_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "TauRangeUp");
  addr[nt] = &tau_range_up;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FixedScale");
  addr[nt] = &fixed_scale;
  id[nt++] = INT;

  strcpy(tag[nt], "FixedShift");
  addr[nt] = &fixed_shift;
  id[nt++] = INT;

  strcpy(tag[nt], "FixedSyserr");
  addr[nt] = &fixed_syserr;
  id[nt++] = INT;

  strcpy(tag[nt], "FixedErrorScale");
  addr[nt] = &fixed_error_scale;
  id[nt++] = INT;

  while(!fin.eof())
  {
    sprintf(str,"empty");

    fin.getline(str, 256);
    if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
      continue;
    if(buf1[0]=='#')
      continue;
    for(i=0, j=-1; i<nt; i++)
      if(strcmp(buf1, tag[i]) == 0)
      {
        j = i;
        tag[i][0] = 0;
        //printf("%s %s\n", buf1, buf2);
        break;
      }
    if(j >=0)
    {
      switch(id[j])
      {
        case DOUBLE:
          *((double *) addr[j]) = atof(buf2);
          break;
        case STRING:
          strcpy((char *)addr[j], buf2);
          break;
        case INT:
          *((int *)addr[j]) = (int) atof(buf2);
          break;
      }
    }
    else
    {
      fprintf(stderr, "# Error in file %s: Tag '%s' is not allowed or multiple defined.\n", 
                    fname.c_str(), buf1);
      exit(0);
    }
  }

  if(strlen(fcont) == 0)
  {
    cout<<"fcont is empty."<<endl;
    exit(-1);
  }

  if(scale_range_low >= scale_range_up)
  {
    cout<<"Incorrect settings in ScaleRangeLow and ScaleRangeUp."<<endl;
    exit(-1);
  }

  if(shift_range_low >= shift_range_up)
  {
    cout<<"Incorrect settings in ShiftRangeLow and ShiftRangeUp."<<endl;
    exit(-1);
  }

  if(syserr_range_low >= syserr_range_up)
  {
    cout<<"Incorrect settings in SyserrRangeLow and SyserrRangeUp."<<endl;
    exit(-1);
  }

  if(errscale_range_low >= errscale_range_up)
  {
    cout<<"Incorrect settings in ErrscaleRangeLow and ErrscaleRangeUp."<<endl;
    exit(-1);
  }

  if(sigma_range_low >= sigma_range_up)
  {
    cout<<"Incorrect settings in SigmaRangeLow and SigmaRangeUp."<<endl;
    exit(-1);
  }

  if(tau_range_low >= tau_range_up)
  {
    cout<<"Incorrect settings in TauRangeLow and TauRangeUp."<<endl;
    exit(-1);
  }

  if(fixed_scale && fixed_shift)
  {
    cout<<"Better not to fix both Scale and Shift parameters."<<endl;
    exit(-1);
  }
  fin.close();
}
void Config::setup(const string& fcont_in, const string& fline_in, 
             int nmcmc_in, double ptol_in, 
             double scale_range_low_in, double scale_range_up_in,
             double shift_range_low_in, double shift_range_up_in,
             double syserr_range_low_in, double syserr_range_up_in,
             double errscale_range_low_in, double errscale_range_up_in,
             double sigma_range_low_in, double sigma_range_up_in,
             double tau_range_low_in, double tau_range_up_in,
             bool fixed_scale_in, bool fixed_shift_in,
             bool fixed_syserr_in, bool fixed_error_scale_in)
{
  strcpy(fcont, fcont_in.c_str());
  strcpy(fline, fline_in.c_str());
  nmcmc = nmcmc_in;
  ptol = ptol_in;
  scale_range_low = scale_range_low_in;
  scale_range_up = scale_range_up_in;
  shift_range_low = shift_range_low_in;
  shift_range_up = shift_range_up_in;
  syserr_range_low = syserr_range_low_in;
  syserr_range_up = syserr_range_up_in;
  errscale_range_low = errscale_range_low_in;
  errscale_range_up = errscale_range_up_in;
  sigma_range_low = sigma_range_low_in;
  sigma_range_up = sigma_range_up_in;
  tau_range_low = tau_range_low_in;
  tau_range_up = tau_range_up_in;

  fixed_scale = fixed_scale_in;
  fixed_shift = fixed_shift_in;
  fixed_syserr = fixed_syserr_in;
  fixed_error_scale = fixed_error_scale_in;

  if(strlen(fcont) == 0)
  {
    cout<<"fcont is empty."<<endl;
    exit(-1);
  }

  if(scale_range_low >= scale_range_up)
  {
    cout<<"Incorrect settings in ScaleRangeLow and ScaleRangeUp."<<endl;
    exit(-1);
  }

  if(shift_range_low >= shift_range_up)
  {
    cout<<"Incorrect settings in ShiftRangeLow and ShiftRangeUp."<<endl;
    exit(-1);
  }

  if(sigma_range_low >= sigma_range_up)
  {
    cout<<"Incorrect settings in SigmaRangeLow and SigmaRangeUp."<<endl;
    exit(-1);
  }

  if(tau_range_low >= tau_range_up)
  {
    cout<<"Incorrect settings in TauRangeLow and TauRangeUp."<<endl;
    exit(-1);
  }

  if(fixed_scale && fixed_shift)
  {
    cout<<"Better not to fix both Scale and Shift parameters."<<endl;
    exit(-1);
  }

  fname.clear();
}

string Config::get_param_filename()
{
  return fname;
}

void Config::print_cfg()
{
  cout<<setw(20)<<"fname: "<<fname<<endl;
  cout<<setw(20)<<"fcont: "<<fcont<<endl;
  cout<<setw(20)<<"fline: "<<fline<<endl;
  cout<<setw(20)<<"nmcmc: "<<nmcmc<<endl;
  cout<<setw(20)<<"scale_range_low: "<<scale_range_low<<endl;
  cout<<setw(20)<<"scale_range_up: "<<scale_range_up<<endl;
  cout<<setw(20)<<"shift_range_low: "<<shift_range_low<<endl;
  cout<<setw(20)<<"shift_range_up: "<<shift_range_up<<endl;
  cout<<setw(20)<<"syserr_range_low: "<<syserr_range_low<<endl;
  cout<<setw(20)<<"syserr_range_up: "<<syserr_range_up<<endl;
  cout<<setw(20)<<"errscale_range_low: "<<errscale_range_low<<endl;
  cout<<setw(20)<<"errscale_range_up: "<<errscale_range_up<<endl;
  cout<<setw(20)<<"sigma_range_low: "<<sigma_range_low<<endl;
  cout<<setw(20)<<"sigma_range_up: "<<sigma_range_up<<endl;
  cout<<setw(20)<<"tau_range_low: "<<tau_range_low<<endl;
  cout<<setw(20)<<"tau_range_up: "<<tau_range_up<<endl;
  cout<<setw(20)<<"fixed_scale: "<<fixed_scale<<endl;
  cout<<setw(20)<<"fixed_shift: "<<fixed_shift<<endl;
  cout<<setw(20)<<"fixed_syserr: "<<fixed_syserr<<endl;
  cout<<setw(20)<<"fixed_error_scale: "<<fixed_error_scale<<endl;

  ofstream fout;
  fout.open("data/param_input");
  fout<<setw(20)<<left<<"fname"<<" = "<<fname<<endl;
  fout<<setw(20)<<left<<"fcont"<<" = "<<fcont<<endl;
  fout<<setw(20)<<left<<"fline"<<" = "<<fline<<endl;
  fout<<setw(20)<<left<<"nmcmc"<<" = "<<nmcmc<<endl;
  fout<<setw(20)<<left<<"scale_range_low"<<" = "<<scale_range_low<<endl;
  fout<<setw(20)<<left<<"scale_range_up"<<" = "<<scale_range_up<<endl;
  fout<<setw(20)<<left<<"shift_range_low"<<" = "<<shift_range_low<<endl;
  fout<<setw(20)<<left<<"shift_range_up"<<" = "<<shift_range_up<<endl;
  fout<<setw(20)<<left<<"syserr_range_low"<<" = "<<syserr_range_low<<endl;
  fout<<setw(20)<<left<<"syserr_range_up"<<" = "<<syserr_range_up<<endl;
  fout<<setw(20)<<left<<"errscale_range_low"<<" = "<<errscale_range_low<<endl;
  fout<<setw(20)<<left<<"errscale_range_up"<<" = "<<errscale_range_up<<endl;
  fout<<setw(20)<<left<<"sigma_range_low"<<" = "<<sigma_range_low<<endl;
  fout<<setw(20)<<left<<"sigma_range_up"<<" = "<<sigma_range_up<<endl;
  fout<<setw(20)<<left<<"tau_range_low"<<" = "<<tau_range_low<<endl;
  fout<<setw(20)<<left<<"tau_range_up"<<" = "<<tau_range_up<<endl;
  fout<<setw(20)<<left<<"fixed_scale"<<" = "<<fixed_scale<<endl;
  fout<<setw(20)<<left<<"fixed_shift"<<" = "<<fixed_shift<<endl;
  fout<<setw(20)<<left<<"fixed_syserr"<<" = "<<fixed_syserr<<endl;
  fout<<setw(20)<<left<<"fixed_error_scale"<<" = "<<fixed_error_scale<<endl;
  fout.close();
}

vector<double> Config::test(const vector<double>& range)
{
  return vector<double>(range);
}
/*=====================================================*/
DataLC::DataLC()
{

}
DataLC::DataLC(size_t n)
{
  time.resize(n);
  flux.resize(n);
  error.resize(n);
}
DataLC::~DataLC()
{
  time.clear();
  flux.clear();
  error.clear();
}
void DataLC::resize(size_t n)
{
  time.resize(n);
  flux.resize(n);
  error.resize(n);
}
/*=====================================================*/
Data::Data()
{
}

Data::Data(const string& fname)
{
  load(fname);
}

Data::~Data()
{
  flux_org.clear();
  error_org.clear();
  time.clear();
  flux.clear();
  error.clear();
  code.clear();

  num_code.clear();
  mean_code.clear();
  code_list.clear();
}

void Data::load(const string& fname)
{
  /* first clear all vectors */
  flux_org.clear();
  error_org.clear();
  time.clear();
  flux.clear();
  error.clear();
  code.clear();

  num_code.clear();
  mean_code.clear();
  code_list.clear();

  /* now read data */
  fstream fin;
  string line;
  stringstream ss;

  fin.open(fname);
  if(fin.fail())
  {
    cout<<"cannot open file "<<fname<<endl;
    exit(-1);
  }

  int j, num;
  int idx, idx_str;
  double t, f, e, mean;
  string cstr;
  
  idx = 0;
  while(1)
  {
    getline(fin, line);
    if(fin.fail())
      break;
    idx_str = line.find_first_not_of(WhiteSpace);
    if(idx_str>=line.size())
      break;
    if(line[0]!='#')
    {
      cout<<"Incorrect line format in "<<fname<<endl;
      exit(-1);
    }

    line.erase(0, 1);  /* skip the first charactor '#' */
    idx_str = line.find_first_not_of(WhiteSpace);
    line.erase(0, idx_str);  /* skip the left white spaces */

    idx_str = line.find_first_of(WhiteSpace); /* extract the code string */
    cstr = line.substr(0, idx_str);
    code_list.push_back(cstr);

    line.erase(0, idx_str); /* extract the number */
    num = stoi(line);
    num_code.push_back(num);
    cout<<code_list[idx]<<"   "<<num_code[idx]<<endl;
    mean = 0.0;
    for(j=0; j<num_code[idx]; j++)
    {
      getline(fin, line);
      if(fin.fail())
      {
        cout<<"# Wrong in reading "<<fname<<endl;
        exit(-1);
      }

      ss.str(line);
      ss>>t>>f>>e;
      if(ss.fail())
      {
        cout<<"# Wrong in reading "<<fname<<endl;
        exit(-1);
      }
        
      time.push_back(t);
      flux_org.push_back(f);
      error_org.push_back(e);
      code.push_back(idx);
      ss.clear();

      mean += f;
    }
    mean /= num;
    mean_code.push_back(mean);
    idx++;
  }
  if(idx==0)
  {
    cout<<"Error: an empty file "<<fname<<endl;
    exit(-1);
  }
  cout<<time.size()<<" points, "<<code_list.size()<<" codes."<<endl;
  fin.close();

  normalize();
  sort_data();
  return;
}

void Data::normalize()
{
  int i;
  norm = mean_code[0]; /* the first dataset is the reference */
  for(i=0; i<flux_org.size(); i++)
  {
    flux_org[i] /= mean_code[code[i]];
    error_org[i] /= mean_code[code[i]];
  }
  return;
}

void Data::sort_data()
{
  vector<int> code_tmp;
  vector<double> time_tmp;

  index.resize(time.size());
  iota(index.begin(), index.end(), 0);
  stable_sort(index.begin(), index.end(), [&](size_t i1, size_t i2) {return time[i1] < time[i2];});
  
  time_tmp = time;
  flux = flux_org;
  error = error_org;
  code_tmp = code;
  int i;
  for(i=0; i<index.size(); i++)
  {
    time[i] = time_tmp[index[i]];
    flux_org[i] = flux[index[i]];
    error_org[i] = error[index[i]];
    code[i] = code_tmp[index[i]];
  }
  code_tmp.clear();
  time_tmp.clear();
}

void Data::check_code(Data& data)
{
  int i;
  if(code_list.size() != data.code_list.size())
  {
    cout<<"# Numbers of codes do not match."<<endl;
    exit(-1);
  }

  for(i=0; i<code_list.size(); i++)
  {
    if(code_list[i] != data.code_list[i])
    {
      cout<<"# Codes do not match or in different orders."<<endl;
      exit(-1);
    }
  }
}
/*=====================================================*/
/* class for calibration */
Cali::Cali()
{
  num_params = 0;
  par_range_model = NULL;
  par_fix = NULL;
  par_fix_val = NULL;

  par_prior_model = NULL;
  par_prior_gaussian = NULL;

  best_params = NULL;
  best_params_std = NULL;
  best_params_covar = NULL;

  workspace = NULL;
  Larr_data = NULL;
}

Cali::Cali(Config& cfg)
     :fcont(cfg.fcont), fline(cfg.fline), cont(cfg.fcont),
      nmcmc(cfg.nmcmc), ptol(cfg.ptol)
{
  int i, j;
  bool isfixed;

  num_params_var = 2;
  size_max = cont.time.size();
  ncode = cont.code_list.size();
  if(!fline.empty())
  {
    line.load(fline);
    size_max = fmax(size_max, line.time.size());
    ncode = fmax(ncode, line.code_list.size());
    num_params_var += 2;
    
    /* check whether cont and line have the same codes */
    cont.check_code(line);

    /* special treatment with line data, line and cont should be scaled with a same factor */
    for(i=0; i<line.time.size(); i++)
    {
      line.flux_org[i] *=  (line.mean_code[line.code[i]]/line.mean_code[0]) *  (cont.mean_code[0]/cont.mean_code[line.code[i]]);
      line.error_org[i] *= (line.mean_code[line.code[i]]/line.mean_code[0]) * (cont.mean_code[0]/cont.mean_code[line.code[i]]);
    }
    for(i=0; i<ncode; i++)
    {
      line.mean_code[i] *= (line.mean_code[i]/line.mean_code[0]) * (cont.mean_code[0]/cont.mean_code[i]);
    }
  }
  /* variability, scale, shift, syserr, error scale */
  num_params = num_params_var + ncode*2 + ncode + ncode;
  if(!fline.empty())
  {
    /* syserr and error scale of line */
    num_params += ncode + ncode;
  }
  par_range_model = new double * [num_params];
  for(i=0; i<num_params; i++)
  {
    par_range_model[i] = new double [2];
  }
  par_fix = new bool [num_params];
  par_fix_val = new double [num_params];
  par_prior_model = new int [num_params];
  par_prior_gaussian = new double * [num_params];
  for(i=0; i<num_params; i++)
  {
    par_prior_gaussian[i] = new double [2];
  } 
  best_params = new double[num_params];
  best_params_std = new double[num_params];
  best_params_covar = new double[num_params*num_params];

  /* set parameter ranges */
  double tau_up, tau_low;
  tau_up = fmin(cfg.tau_range_up, (cont.time[cont.time.size()-1]-cont.time[0]));
  tau_low = fmax(cfg.tau_range_low,(cont.time[cont.time.size()-1]-cont.time[0])/cont.time.size());
  i=0; /* sigma */
  par_range_model[i][0] = log(cfg.sigma_range_low);
  par_range_model[i][1] = log(cfg.sigma_range_up);
  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;
  
  i+=1; /* tau */
  par_range_model[i][0] = log(tau_low);
  par_range_model[i][1] = log(tau_up);
  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  if(!fline.empty())
  {
    tau_up = fmin(cfg.tau_range_up, (line.time[line.time.size()-1]-line.time[0]));
    tau_low = fmin(cfg.tau_range_low, (line.time[line.time.size()-1]-line.time[0])/line.time.size());
    i+=1; /* sigma */
    par_range_model[i][0] = log(cfg.sigma_range_low);
    par_range_model[i][1] = log(cfg.sigma_range_up);
    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  
    i+=1; /* tau */
    par_range_model[i][0] = log(tau_low);
    par_range_model[i][1] = log(tau_up);
    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  
  /* scale */
  for(j=0; j<ncode; j++)
  {
    i+=1;
    par_range_model[i][0] = cfg.scale_range_low;
    par_range_model[i][1] = cfg.scale_range_up;
    par_prior_model[i] = LOG;
  }
  /* shift */
  for(j=0; j<ncode; j++)
  {
    i+=1;
    par_range_model[i][0] = cfg.shift_range_low;
    par_range_model[i][1] = cfg.shift_range_up;
    par_prior_model[i] = UNIFORM;
  }
  /* syserr of continuum */
  for(j=0; j<ncode; j++)
  {
    i+=1;
    par_range_model[i][0] = cfg.syserr_range_low;
    par_range_model[i][1] = cfg.syserr_range_up;
    par_prior_model[i] = UNIFORM;
  }
  /* error scale of continuum */
  for(j=0; j<ncode; j++)
  {
    i+=1;
    par_range_model[i][0] = cfg.errscale_range_low;
    par_range_model[i][1] = cfg.errscale_range_up;
    par_prior_model[i] = LOG;
  }

  if(!fline.empty())
  {
    /* syserr of line */
    for(j=0; j<ncode; j++)
    {
      i+=1;
      par_range_model[i][0] = cfg.syserr_range_low;
      par_range_model[i][1] = cfg.syserr_range_up;
      par_prior_model[i] = UNIFORM;
    }
    /* error scale of line */
    for(j=0; j<ncode; j++)
    {
      i+=1;
      par_range_model[i][0] = cfg.errscale_range_low;
      par_range_model[i][1] = cfg.errscale_range_up;
      par_prior_model[i] = LOG;
    }
  }

  for(i=0; i<num_params; i++)
  {
    par_fix[i] = NOFIXED;
    par_fix_val[i] = -DBL_MAX;
  }
  par_fix[num_params_var] = FIXED;
  par_fix_val[num_params_var] = 1.0;
  par_fix[num_params_var+ncode] = FIXED;
  par_fix_val[num_params_var+ncode] = 0.0;

  if(cfg.fixed_scale)
  {
    for(i=1; i<ncode; i++)
    {
      par_fix[num_params_var+i] = FIXED;
      par_fix_val[num_params_var+i] = 1.0;
    }
  }

  if(cfg.fixed_shift)
  {
    for(i=1; i<ncode; i++)
    {
      par_fix[num_params_var+i+ncode] = FIXED;
      par_fix_val[num_params_var+i+ncode] = 0.0;
    }
  }

  if(cfg.fixed_syserr)
  {
    for(i=0; i<ncode; i++)
    {
      par_fix[num_params_var+2*ncode+i] = FIXED;
      par_fix_val[num_params_var+2*ncode+i] = 0.0;
    }

    if(!fline.empty())
    {
      for(i=0; i<ncode; i++)
      {
        par_fix[num_params_var+4*ncode+i] = FIXED;
        par_fix_val[num_params_var+4*ncode+i] = 0.0;
      }
    }
  }

  if(cfg.fixed_error_scale)
  {
    for(i=0; i<ncode; i++)
    {
      par_fix[num_params_var+3*ncode+i] = FIXED;
      par_fix_val[num_params_var+3*ncode+i] = 1.0;
    }

    if(!fline.empty())
    {
      for(i=0; i<ncode; i++)
      {
        par_fix[num_params_var+5*ncode+i] = FIXED;
        par_fix_val[num_params_var+5*ncode+i] = 1.0;
      }
    }
  }

  workspace = new double[10*size_max];
  Larr_data = new double[size_max];
  for(i=0; i<size_max; i++)
    Larr_data[i] = 1.0;
  
  fptrset = dnest_malloc_fptrset();
  fptrset->from_prior = from_prior_cali;
  fptrset->perturb = perturb_cali;
  fptrset->print_particle = print_particle_cali;
  fptrset->log_likelihoods_cal = prob_cali;


  /* reconstruction */
  double t1, t2, tspan;
  tspan = cont.time[cont.time.size()-1] - cont.time[0];
  t1 = cont.time[0] - 0.05*tspan;
  t2 = cont.time[cont.time.size()-1] + 0.05*tspan;
  cont_recon.resize(cont.time.size()*2);
  size_recon_max = cont_recon.time.size();
  for(i=0; i<cont_recon.time.size(); i++)
  {
    cont_recon.time[i] = t1 + (t2 - t1)/(cont_recon.time.size()-1.0) * i;
  }
  if(!fline.empty())
  {
    line_recon.resize(line.time.size()*2);
    tspan = line.time[line.time.size()-1] - line.time[0];
    t1 = line.time[0] - 0.05*tspan;
    t2 = line.time[line.time.size()-1] + 0.05*tspan;
    for(i=0; i<line_recon.time.size(); i++)
    {
      line_recon.time[i] = t1 + (t2 - t1)/(line_recon.time.size()-1.0) * i;
    }
    size_recon_max = fmax(size_recon_max, line_recon.time.size());
  }
}

Cali::~Cali()
{
  int i;

  delete[] par_fix;
  delete[] par_fix_val;
  delete[] par_prior_model;
  for(i=0; i<num_params; i++)
  {
    delete[] par_range_model[i];
    delete[] par_prior_gaussian[i];
  }
  delete[] par_range_model;
  delete[] par_prior_gaussian;
  delete[] best_params;
  delete[] best_params_std;
  delete[] best_params_covar;

  delete[] workspace;
  delete[] Larr_data;
  dnest_free_fptrset(fptrset);
}

void Cali::align(double *model)
{
  int i, idx;
  double *ps_scale = model+num_params_var;
  double *es_shift = ps_scale + ncode;
  double *syserr = es_shift + ncode;
  double *error_scale = syserr + ncode;
  for(i=0; i<cont.time.size(); i++)
  {
    idx = cont.code[i];
    cont.flux[i] = cont.flux_org[i] * ps_scale[idx] - es_shift[idx];
    /* note that this error does not include errors of scale and shift */
    cont.error[i] = sqrt(cont.error_org[i]*cont.error_org[i]*error_scale[idx]*error_scale[idx] 
                    + syserr[idx]*syserr[idx]) * ps_scale[idx];
  }

  if(!fline.empty())
  {
    syserr = error_scale + ncode;
    error_scale = syserr + ncode;
    for(i=0; i<line.time.size(); i++)
    {
      idx = line.code[i];
      line.flux[i] = line.flux_org[i] * ps_scale[idx];
      line.error[i] = sqrt(line.error_org[i]*line.error_org[i]*error_scale[idx]*error_scale[idx] 
                      + syserr[idx]*syserr[idx] ) * ps_scale[idx];
    }
  }
}

void Cali::align_with_error()
{
  int i, idx;
  double *ps_scale = best_params+num_params_var;
  double *es_shift = ps_scale + ncode;
  double *syserr = es_shift + ncode;
  double *error_scale = syserr + ncode;
  double *ps_scale_err = best_params_std + num_params_var;
  double *es_shift_err = ps_scale_err + ncode;
  for(i=0; i<cont.time.size(); i++)
  {
    idx = cont.code[i];
    cont.flux[i] = cont.flux_org[i] * ps_scale[idx] - es_shift[idx];
    cont.error[i] = sqrt((cont.error_org[i]*cont.error_org[i]*error_scale[idx]*error_scale[idx] +syserr[idx]*syserr[idx]) 
                         *ps_scale[idx]*ps_scale[idx]
                        +pow(cont.flux_org[i]*ps_scale_err[idx], 2.0)
                        +pow(es_shift_err[idx], 2.0)
                        -2.0*cont.flux_org[i]*best_params_covar[(num_params_var+idx)*num_params + (num_params_var+idx+ncode)]
                        );
  }

  if(!fline.empty())
  {
    syserr = error_scale + ncode;
    error_scale = syserr + ncode;
    for(i=0; i<line.time.size(); i++)
    {
      idx = line.code[i];
      line.flux[i] = line.flux_org[i] * ps_scale[idx];
      line.error[i] = sqrt((line.error_org[i]*line.error_org[i]*error_scale[idx]*error_scale[idx] + syserr[idx]*syserr[idx]) 
                          *ps_scale[idx]*ps_scale[idx]
                          +pow(line.flux_org[i] * ps_scale_err[idx], 2.0)
                          );
    }
  }
}

void Cali::mcmc()
{
  int i, argc=0;
  char **argv;
  double logz_con;
  char dnest_options_file[256];

  argv = new char * [9];
  for(i=0; i<9; i++)
  {
    argv[i] = new char [256];
  }
  
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], "./");
  strcat(argv[argc++], "/data/restart_dnest.txt");

  strcpy(dnest_options_file, "OPTIONS");
  logz_con = dnest(argc, argv, fptrset, num_params, "data/", nmcmc, ptol, (void *)this);

  for(i=0; i<9; i++)
  {
    delete[] argv[i];
  }
  delete[] argv;
}

void Cali::get_best_params()
{
  int i, j, num_ps;
  FILE *fp;
  char posterior_sample_file[256];
  double *post_model, *posterior_sample;
  double *pm, *pmstd;

  strcpy(posterior_sample_file, "data/posterior_sample.txt");

  /* open file for posterior sample */
  fp = fopen(posterior_sample_file, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
    exit(0);
  }

  /* read number of points in posterior sample */
  if(fscanf(fp, "# %d", &num_ps) < 1)
  {
    fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
    exit(0);
  }
  printf("# Number of points in posterior sample: %d\n", num_ps);

  if(num_ps < 500)
  {
    cout<<"########################################################\n"
          "# Too few effective posterior samples.\n"
          "# Try to increse nmcmc, or decrease ptol,\n"
          "# or set a more appropriate range for scale and shift.\n"
          "########################################################"<<endl;
    exit(-1);
  }

  post_model = new double[num_params*sizeof(double)];
  posterior_sample = new double[num_ps * num_params*sizeof(double)];
  
  for(i=0; i<num_ps; i++)
  {
    for(j=0; j<num_params; j++)
    {
      if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
      {
        fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
        exit(0);
      }
    }
    fscanf(fp, "\n");

    memcpy(posterior_sample+i*num_params, post_model, num_params*sizeof(double));

  }
  fclose(fp);

  /* calcaulte mean and standard deviation of posterior samples. */
  pm = (double *)best_params;
  pmstd = (double *)best_params_std;
  for(j=0; j<num_params; j++)
  {
    pm[j] = pmstd[j] = 0.0;
  }
  for(i=0; i<num_ps; i++)
  {
    for(j =0; j<num_params; j++)
      pm[j] += *((double *)posterior_sample + i*num_params + j );
  }

  for(j=0; j<num_params; j++)
    pm[j] /= num_ps;

  for(i=0; i<num_ps; i++)
  {
    for(j=0; j<num_params; j++)
      pmstd[j] += pow( *((double *)posterior_sample + i*num_params + j ) - pm[j], 2.0 );
  }

  for(j=0; j<num_params; j++)
  {
    if(num_ps > 1)
      pmstd[j] = sqrt(pmstd[j]/(num_ps-1.0));
    else
      pmstd[j] = 0.0;
  }  

  /* find out the largest likelihood */
  ifstream fin;
  vector<double> prob(num_ps);
  string str;
  double prob_max;
  int ip_max;
  fin.open("data/posterior_sample_info.txt");
  if(!fin.good())
  {
    cout<<"Error: Cannot open file data/posterior_sample_info.txt.\n"<<endl;
    exit(0);
  }
  getline(fin, str);
  prob_max = -DBL_MAX;
  for(i=0; i<num_ps; i++)
  {
    fin>>prob[i];
    if(prob_max < prob[i])
    {
      prob_max = prob[i];
      ip_max = i;
    }
  }
  fin.close();
  cout<<"Lmax:"<<prob_max<<" at "<<ip_max<<"th posterior sample."<<endl;
  printf("The params with the maximum likelihood:\n");
  for(j = 0; j<num_params; j++)
    printf("%d %f\n", j, *((double *)posterior_sample + ip_max*num_params + j)); 

  /* calculate covariance */
  double covar;
  int k;
  for(i=0; i<num_params; i++)
  {
    for(j=0; j<i; j++)
    {
      covar = 0.0;
      for(k=0; k<num_ps; k++)
      {
        covar += (*((double *)posterior_sample + k*num_params + i )) * (*((double *)posterior_sample + k*num_params + j ));
      }
      best_params_covar[i*num_params+j] = best_params_covar[j*num_params+i] = covar/num_ps - best_params[i]*best_params[j];
    }
    best_params_covar[i*num_params+i] = best_params_std[i] * best_params_std[i];
  }

  /* directly calculate flux and error */
  int stat_type = 0;
  double error_scale_shift;
  DataLC cont_output(cont.time.size());
  if(stat_type == 0)  /* mediate values */
  {
    double *flux, *error;
    flux = new double [cont.time.size()*num_ps];
    error = new double [cont.time.size()*num_ps];

    for(i=0; i<num_ps; i++)
    {
      align((double *)posterior_sample + i*num_params);
      for(j=0; j<cont.time.size(); j++)
      {
        flux[j * num_ps + i] = cont.flux[j];
        error[j * num_ps + i] = cont.error[j];
      }
    }
    /* sort, and use the media value */
    for(j=0; j<cont.time.size(); j++)
    {
      /* ascending order */
      qsort(flux+j*num_ps, num_ps, sizeof(double), compare);
      cont_output.flux[j] = flux[j*num_ps + (int)(0.5*num_ps)];
      qsort(error+j*num_ps, num_ps, sizeof(double), compare);
      cont_output.error[j] = error[j*num_ps + (int)(0.5*num_ps)];
      /* include error of scale and shift */
      error_scale_shift = 0.5*( (cont_output.flux[j] - flux[j*num_ps + (int)(0.1585*num_ps)])
                               +(flux[j*num_ps + (int)(0.8415*num_ps)] - cont_output.flux[j]) );
      cont_output.error[j] = sqrt(cont_output.error[j]*cont_output.error[j] + error_scale_shift*error_scale_shift);
    }
  
    cont.flux = cont_output.flux;
    cont.error = cont_output.error;
    delete[] flux;
    delete[] error;

    if(!fline.empty())
    {
      DataLC line_output(line.time.size());
      flux = new double [line.time.size()*num_ps];
      error = new double [line.time.size()*num_ps];
      for(i=0; i<num_ps; i++)
      {
        align((double *)posterior_sample + i*num_params);
        for(j=0; j<line.time.size(); j++)
        {
          flux[j * num_ps + i] = line.flux[j];
          error[j * num_ps + i] = line.error[j];
        }
      }
      /* sort, and use the media value */
      for(j=0; j<line.time.size(); j++)
      {
        qsort(flux+j*num_ps, num_ps, sizeof(double), compare);
        line_output.flux[j] = flux[j*num_ps + (int)(0.5*num_ps)];
        qsort(error+j*num_ps, num_ps, sizeof(double), compare);
        line_output.error[j] = error[j*num_ps + (int)(0.5*num_ps)];
        /* include error of scale and shift */
        error_scale_shift = 0.5*( (line_output.flux[j] - flux[j*num_ps + (int)(0.1585*num_ps)])
                               +(flux[j*num_ps + (int)(0.8415*num_ps)] - line_output.flux[j]) );
        line_output.error[j] = sqrt(line_output.error[j]*line_output.error[j] + error_scale_shift*error_scale_shift);
      }
    
      line.flux = line_output.flux;
      line.error = line_output.error;
      delete[] flux;
      delete[] error;
    }
  }
  else if(stat_type == 1)  /* mean values */
  {
    vector<double> flux_rms(cont.time.size());
    for(j=0; j<cont.time.size(); j++)
    {
      cont_output.flux[j] = 0.0;
      cont_output.error[j] = 0.0;

      flux_rms[j] = 0.0;
    }

    for(i=0; i<num_ps; i++)
    {
      align((double *)posterior_sample + i*num_params);
      for(j=0; j<cont.time.size(); j++)
      {
        cont_output.flux[j] += cont.flux[j];
        cont_output.error[j] += cont.error[j];

        flux_rms[j] += cont.flux[j]*cont.flux[j];
      }
    }
    for(j=0; j<cont.time.size(); j++)
    {
      cont_output.flux[j] /= num_ps;
      cont_output.error[j] /= num_ps;
      
      flux_rms[j] /= num_ps;
      /* include error of scale */
      error_scale_shift = fmax(0.0, flux_rms[j] - cont_output.flux[j]*cont_output.flux[j]);
      cont_output.error[j] = sqrt(cont_output.error[j]*cont_output.error[j] + error_scale_shift);
    }
  
    cont.flux = cont_output.flux;
    cont.error = cont_output.error;

    if(!fline.empty())
    {
      DataLC line_output(line.time.size());
      flux_rms.resize(line.time.size());

      for(j=0; j<line.time.size(); j++)
      {
        line_output.flux[j] = 0.0;
        line_output.error[j] = 0.0;

        flux_rms[j] = 0.0;
      }

      for(i=0; i<num_ps; i++)
      {
        align((double *)posterior_sample + i*num_params);
        for(j=0; j<line.time.size(); j++)
        {
          line_output.flux[j] += line.flux[j];
          line_output.error[j] += line.error[j];

          flux_rms[j] += line.flux[j]*line.flux[j];
        }
      }
      for(j=0; j<line.time.size(); j++)
      {
        line_output.flux[j] /= num_ps;
        line_output.error[j] /= num_ps;

        flux_rms[j] /= num_ps;
        /* include error of scale */
        error_scale_shift = fmax(0.0, flux_rms[j] - line_output.flux[j]*line_output.flux[j]);
        line_output.error[j] = sqrt(line_output.error[j]*line_output.error[j] + error_scale_shift);
      }
    
      line.flux = line_output.flux;
      line.error = line_output.error;
    }
  }
  else /* error propagate */
  {
    align_with_error();
  }

  delete[] post_model;
  delete[] posterior_sample;
}

void Cali::output()
{
  int i;
  ofstream fout;
  fout.open(fcont+"_cali");
  for(i=0; i<cont.time.size(); i++)
  {
    fout<<fixed<<cont.time[i]
        <<scientific<<" "<<cont.flux[i]*cont.norm<<"  "<<cont.error[i]*cont.norm<<"  "<<cont.code_list[cont.code[i]]<<endl;
  }
  fout.close();

  /* output indices that sort cont */
  fout.open("data/cont_sort_index.txt");
  for(i=0; i<cont.time.size(); i++)
  {
    fout<<cont.index[i]<<endl;
  }
  fout.close();
  
  if(!fline.empty())
  {
    fout.open(fline+"_cali");
    for(i=0; i<line.time.size(); i++)
    {
      fout<<fixed<<line.time[i]
          <<scientific<<" "<<line.flux[i]*line.norm<<"  "<<line.error[i]*line.norm<<"  "<<line.code_list[line.code[i]]<<endl;
    }
    fout.close();

    /* output indices that sort line */
    fout.open("data/line_sort_index.txt");
    for(i=0; i<line.time.size(); i++)
    {
      fout<<line.index[i]<<endl;
    }
    fout.close();
  }

  fout.open("data/factor.txt");
  fout<<"Code \t Scale  \t Error  \t Shift  \t Error    \t     Cov"<<endl;
  for(i=0; i<ncode;i++)
  {
    fout<<scientific
        <<cont.code_list[i]<<"\t"<<best_params[i+num_params_var]<<"\t"<<best_params_std[i+num_params_var]
        <<"\t"<<best_params[i+ncode+num_params_var]<<"\t"<<best_params_std[i+ncode+num_params_var]
        <<"\t"<<best_params_covar[(i+num_params_var)*num_params+(i+num_params_var+ncode)]<<endl;
  }
  fout.close();

  fout.open("data/pyCALI_output.txt");
  fout<<"# mean of continuum:"<<endl;
  for(i=0; i<ncode; i++)
  {
    fout<<cont.code_list[i]<<"\t"<<cont.mean_code[i]<<"\t"<<cont.num_code[i]<<endl;
  }
  if(!fline.empty())
  {
    fout<<"# mean of line:"<<endl;
    for(i=0; i<ncode; i++)
    {
      fout<<line.code_list[i]<<"\t"<<line.mean_code[i]<<"\t"<<line.num_code[i]<<endl;
    }
  }
  fout.close();
}

void Cali::recon()
{
  double *Lbuf, *ybuf, *y, *Cq, *yq, *W, *D, *phi;
  double *USmat, *PEmat1, *PEmat2;
  double syserr;

  double *pm = (double *)best_params;
  double sigma, sigma2, tau;
  int i, info, nq;
  int nd_cont = cont.time.size(), nd_cont_recon = cont_recon.time.size();

  syserr = 0.0;
  tau = exp(pm[1]);
  sigma = exp(pm[0]) * sqrt(tau);
  sigma2 = sigma*sigma;
  
  nq = 1;
  Lbuf = workspace;
  ybuf = Lbuf + nd_cont*nq; 
  y = ybuf + nd_cont;
  Cq = y + nd_cont;
  yq = Cq + nq*nq;

  W = new double [size_recon_max];
  D = new double [size_recon_max];
  phi = new double [size_recon_max];
  USmat = new double [size_recon_max * size_max];
  PEmat1 = new double [size_recon_max * size_max];
  PEmat2 = new double [size_recon_max * size_recon_max];

  compute_semiseparable_drw(cont.time.data(), nd_cont, sigma2, 1.0/tau, cont.error.data(), syserr, W, D, phi);
  // Cq^-1 = L^TxC^-1xL
  multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_cont, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_cont);

  // L^TxC^-1xy
  multiply_matvec_semiseparable_drw(cont.flux.data(), W, D, phi, nd_cont, sigma2, ybuf);
  multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, nd_cont);

  // (hat q) = Cqx(L^TxC^-1xy)
  inverse_pomat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);
  for(i=0; i<nq; i++)
    yq[i] = ybuf[i];
  
  // y = yc - Lxq
  multiply_matvec_MN(Larr_data, nd_cont, nq, yq, ybuf);
  for(i=0; i<nd_cont; i++)
  {
    y[i] = cont.flux[i] - ybuf[i];
  }
  
  set_covar_Umat_cont(sigma, tau, USmat);
  // (hat s) = SxC^-1xy
  multiply_matvec_semiseparable_drw(y, W, D, phi, nd_cont, sigma2, ybuf);
  multiply_matvec_MN(USmat, nd_cont_recon, nd_cont, ybuf, cont_recon.flux.data());

  // SxC^-1xS^T
  multiply_mat_transposeB_semiseparable_drw(USmat, W, D, phi, nd_cont, nd_cont_recon, sigma2, PEmat1);
  multiply_mat_MN(USmat, PEmat1, PEmat2, nd_cont_recon, nd_cont_recon, nd_cont);

  for(i=0; i<nd_cont_recon; i++)
  {
    cont_recon.error[i] = sqrt(sigma2 + syserr*syserr - PEmat2[i*nd_cont_recon + i]);
  }

  for(i=0; i<nd_cont_recon; i++)
  {
    cont_recon.flux[i] += yq[0];
  }

  ofstream fout;
  fout.open("data/cont_recon.txt");
  for(i=0; i<nd_cont_recon; i++)
  {
    fout<<scientific
        <<cont_recon.time[i]<<"   "<<cont_recon.flux[i]*cont.norm<<"  "<<cont_recon.error[i]*cont.norm<<endl;
  }
  fout.close();

  if(!fline.empty())
  {
    int nd_line = line.time.size(), nd_line_recon = line_recon.time.size();
    
    syserr = 0.0;
    tau = exp(pm[3]);
    sigma = exp(pm[2]) * sqrt(tau);
    sigma2 = sigma*sigma;
    
    nq = 1;
    Lbuf = workspace;
    ybuf = Lbuf + nd_line*nq; 
    y = ybuf + nd_line;
    Cq = y + nd_line;
    yq = Cq + nq*nq;

    compute_semiseparable_drw(line.time.data(), nd_line, sigma2, 1.0/tau, line.error.data(), syserr, W, D, phi);
    // Cq^-1 = L^TxC^-1xL
    multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_line, nq, sigma2, Lbuf);
    multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_line);

    // L^TxC^-1xy
    multiply_matvec_semiseparable_drw(line.flux.data(), W, D, phi, nd_line, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, nd_line);

    // (hat q) = Cqx(L^TxC^-1xy)
    inverse_pomat(Cq, nq, &info);
    multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);
    for(i=0; i<nq; i++)
      yq[i] = ybuf[i];
  
    // y = yc - Lxq
    multiply_matvec_MN(Larr_data, nd_line, nq, yq, ybuf);
    for(i=0; i<nd_line; i++)
    {
      y[i] = line.flux[i] - ybuf[i];
    }
  
    set_covar_Umat_line(sigma, tau, USmat);
    // (hat s) = SxC^-1xy
    multiply_matvec_semiseparable_drw(y, W, D, phi, nd_line, sigma2, ybuf);
    multiply_matvec_MN(USmat, nd_line_recon, nd_line, ybuf, line_recon.flux.data());

    // SxC^-1xS^T
    multiply_mat_transposeB_semiseparable_drw(USmat, W, D, phi, nd_line, nd_line_recon, sigma2, PEmat1);
    multiply_mat_MN(USmat, PEmat1, PEmat2, nd_line_recon, nd_line_recon, nd_line);

    for(i=0; i<nd_line_recon; i++)
    {
      line_recon.error[i] = sqrt(sigma2 + syserr*syserr - PEmat2[i*nd_line_recon + i]);
    }

    for(i=0; i<nd_line_recon; i++)
    {
      line_recon.flux[i] += yq[0];
    }

    ofstream fout;
    fout.open("data/line_recon.txt");
    for(i=0; i<nd_line_recon; i++)
    {
      fout<<scientific
          <<line_recon.time[i]<<"   "<<line_recon.flux[i]*line.norm<<"  "<<line_recon.error[i]*line.norm<<endl;
    }
    fout.close();
  }

  delete[] D;
  delete[] W;
  delete[] phi;
  delete[] USmat;
  delete[] PEmat1;
  delete[] PEmat2;
}
void Cali::set_covar_Umat_cont(double sigma, double tau, double *USmat)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<cont_recon.time.size(); i++)
  {
    t1 = cont_recon.time[i];
    for(j=0; j<cont.time.size(); j++)
    {
      t2 = cont.time[j];
      USmat[i*cont.time.size()+j] = sigma*sigma * exp (- fabs(t1-t2) / tau );
    }
  }
  return;
}
void Cali::set_covar_Umat_line(double sigma, double tau, double *USmat)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<line_recon.time.size(); i++)
  {
    t1 = line_recon.time[i];
    for(j=0; j<line.time.size(); j++)
    {
      t2 = line.time[j];
      USmat[i*line.time.size()+j] = sigma*sigma * exp (- fabs(t1-t2) / tau );
    }
  }
  return;
}
double Cali::get_norm_cont()
{
  return cont.norm;
}
double Cali::get_norm_line()
{
  return line.norm;
}
/*=============================================================*/
double prob_cali(const void *model, const void *arg)
{
  Cali *cali = (Cali *)arg;
  double prob, prob1=0.0, prob2=0.0, lambda, ave_con, lndet, sigma, sigma2, tau;
  double lndet_n;
  double *ybuf, *W, *D, *phi, *Cq, *Lbuf, *yq;
  int i, nq;
  double *workspace = cali->workspace;
  double *pm = (double *)model;
  double *ps_scale = pm + cali->num_params_var;
  Data &cont = cali->cont;
  int nd_cont = cont.time.size();
  double *Larr_data = cali->Larr_data;

  nq = 1;
  Lbuf = workspace;
  ybuf = Lbuf + nd_cont*nq;
  W = ybuf + nd_cont;
  D = W + nd_cont;
  phi = D + nd_cont;
  Cq = phi + nd_cont;
  yq = Cq + nq*nq;

  tau = exp(pm[1]);
  sigma = exp(pm[0]) * sqrt(tau);
  sigma2 = sigma*sigma;
  
  cali->align(pm);

  compute_semiseparable_drw(cont.time.data(), nd_cont, sigma2, 1.0/tau, cont.error.data(), 0.0, W, D, phi);
  lndet = 0.0;
  for(i=0; i<nd_cont; i++)
    lndet += log(D[i]);

  /* calculate L^T*C^-1*L */
  multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_cont, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_cont);

  /* calculate L^T*C^-1*y */
  multiply_matvec_semiseparable_drw(cont.flux.data(), W, D, phi, nd_cont, sigma2, ybuf);
  multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, nd_cont);
  
  lambda = Cq[0];
  ave_con = yq[0]/Cq[0];

/* get the probability */
  for(i=0;i<nd_cont;i++)
  {
    ybuf[i] = cont.flux[i] - ave_con;
  }
  multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_cont, sigma2, Lbuf);
  prob1 = -0.5 * cblas_ddot(nd_cont, ybuf, 1, Lbuf, 1);
  
  lndet_n = 0.0;
  for(i=0; i<cont.num_code.size(); i++)
  {
    lndet_n += 2.0*log(ps_scale[i]) * cont.num_code[i];
  }
  prob1 = prob1 - 0.5*lndet - 0.5*log(lambda) + 0.5 * lndet_n;
  
  if(!cali->fline.empty())
  {
    Data& line = cali->line;
    int nd_line = line.time.size();

    Lbuf = workspace;
    ybuf = Lbuf + nd_line*nq;
    W = ybuf + nd_line;
    D = W + nd_line;
    phi = D + nd_line;
    Cq = phi + nd_line;
    yq = Cq + nq*nq;

    tau = exp(pm[3]);
    sigma = exp(pm[2])*sqrt(tau);
    sigma2 = sigma*sigma;

    compute_semiseparable_drw(line.time.data(), nd_line, sigma2, 1.0/tau, line.error.data(), 0.0, W, D, phi);
    lndet = 0.0;
    for(i=0; i<nd_line; i++)
      lndet += log(D[i]);

    /* calculate L^T*C^-1*L */
    multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_line, nq, sigma2, Lbuf);
    multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_line);

    /* calculate L^T*C^-1*y */
    multiply_matvec_semiseparable_drw(line.flux.data(), W, D, phi, nd_line, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, nd_line);
  
    lambda = Cq[0];
    ave_con = yq[0]/Cq[0];

    for(i=0;i<nd_line;i++)
    {
      ybuf[i] = line.flux[i] - ave_con;
    }
    multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_line, sigma2, Lbuf);
    prob2 = -0.5 * cblas_ddot(nd_line, ybuf, 1, Lbuf, 1);

    lndet_n =  0.0;
    for(i=0; i<line.num_code.size(); i++)
    {
      lndet_n += 2.0*log(ps_scale[i]) * line.num_code[i];
    }

    prob2 = prob2 - 0.5*lndet - 0.5*log(lambda) + 0.5 * lndet_n;
  }
  
  prob = prob1 + prob2;
  
  return prob;
}
void from_prior_cali(void *model, const void *arg)
{
  int i;
  double *pm = (double *)model;
  
  Cali *cali = (Cali *)arg;

  for(i=0; i<cali->num_params; i++)
  {
    if(cali->par_prior_model[i] == GAUSSIAN )
    {
      pm[i] = dnest_randn() * cali->par_prior_gaussian[i][1] + cali->par_prior_gaussian[i][0];
      dnest_wrap(&pm[i], cali->par_range_model[i][0], cali->par_range_model[i][1]);
    }
    else if(cali->par_prior_model[i] == LOG)
    {
      pm[i] = log(cali->par_range_model[i][0]) + dnest_rand()*(log(cali->par_range_model[i][1]) - log(cali->par_range_model[i][0]));
      pm[i] = exp(pm[i]);
    }
    else 
    {
      pm[i] = cali->par_range_model[i][0] + dnest_rand()*(cali->par_range_model[i][1] - cali->par_range_model[i][0]);
    }
  }

  for(i=0; i<cali->num_params; i++)
  {
    if(cali->par_fix[i] == FIXED)
      pm[i] = cali->par_fix_val[i];
  }
}
void print_particle_cali(FILE *fp, const void *model, const void *arg)
{
  int i;
  double *pm = (double *)model;

  Cali *cali = (Cali *)arg;

  for(i=0; i<cali->num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
}
double perturb_cali(void *model, const void *arg)
{
  double *pm = (double *)model;
  double logH = 0.0, width;
  int which;
  
  Cali *cali = (Cali *)arg;

  /* sample variability parameters more frequently */
  do
  {
    which = dnest_rand_int(cali->num_params);
  }while(cali->par_fix[which] == FIXED);
 
  width = ( cali->par_range_model[which][1] - cali->par_range_model[which][0] );

  if(cali->par_prior_model[which] == UNIFORM)
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), cali->par_range_model[which][0], cali->par_range_model[which][1]);
  }
  else if(cali->par_prior_model[which] == LOG)
  {
    logH -= (-log(pm[which]));
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), cali->par_range_model[which][0], cali->par_range_model[which][1]);
    logH += (-log(pm[which]));
  }
  else
  {
    logH -= (-0.5*pow((pm[which] - cali->par_prior_gaussian[which][0])/cali->par_prior_gaussian[which][1], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], cali->par_range_model[which][0], cali->par_range_model[which][1]);
    logH += (-0.5*pow((pm[which] - cali->par_prior_gaussian[which][0])/cali->par_prior_gaussian[which][1], 2.0) );
  }
  
  return logH;
}