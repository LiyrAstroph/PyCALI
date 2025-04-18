
#include <iostream>
#include <iomanip>
#include <fstream> 
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <float.h>
#include <math.h>
#include <sys/stat.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_filter.h>

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

  fixed_codes.clear();
  fixed_scalecodes.clear();

  fcont = new char [256];
  strcpy(fcont, "\0");
  fline.clear();
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

  flag_norm = true;

  fixed_codes.clear();
  fixed_scalecodes.clear();
  
  fcont = new char [256];
  strcpy(fcont, "\0");
  fline.clear();

  load(fname);
}
Config::~Config()
{
  fline.clear();
  delete[] fcont;
}
void Config::load(const string& fname)
{
  ifstream fin;
  char fbuf[256], fbuf_codes[256], fbuf_scalecodes[256];
  /* empty bufs */
  fbuf_codes[0]='\0';
  fbuf_scalecodes[0]='\0';

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
  #define BOOL 4

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
  addr[nt] = fbuf;
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
  id[nt++] = BOOL;

  strcpy(tag[nt], "FixedShift");
  addr[nt] = &fixed_shift;
  id[nt++] = BOOL;

  strcpy(tag[nt], "FixedSyserr");
  addr[nt] = &fixed_syserr;
  id[nt++] = BOOL;

  strcpy(tag[nt], "FixedErrorScale");
  addr[nt] = &fixed_error_scale;
  id[nt++] = BOOL;

  strcpy(tag[nt], "FixedCodes");
  addr[nt] = fbuf_codes;
  id[nt++] = STRING;

  strcpy(tag[nt], "FixedScaleCodes");
  addr[nt] = fbuf_scalecodes;
  id[nt++] = STRING;

  strcpy(tag[nt], "FlagNorm");
  addr[nt] = &flag_norm;
  id[nt++] = BOOL;

  // default values 
  strcpy(fbuf,"\0");

  while(!fin.eof())
  {
    sprintf(str,"empty");

    fin.getline(str, 256);
    if(str[0]=='#')
      continue;
    if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
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
        case BOOL:
          *((bool *)addr[j]) = (bool) atof(buf2);
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

  /* parse fline string */
  parse_fline_str(fbuf);

  parse_fixed_codes_str(fbuf_codes);

  parse_fixed_scalecodes_str(fbuf_scalecodes);
}

void Config::parse_fline_str(const string& fline_str)
{
  fline.clear();
  if(fline_str.empty())
    return;

  string strbuf(fline_str);
  char buf1[256];
  char *pstr, *pchr;
  
  pstr = (char *)strbuf.c_str();
  pchr = strchr(pstr, ',');
  if(pchr == NULL)
  {
    fline.push_back(pstr);
  }
  else 
  {
    strncpy(buf1, pstr, pchr-pstr);
    buf1[pchr-pstr] = '\0';
    fline.push_back(buf1);
    pstr = pchr + 1;
    while(1)
    {
      pchr = strchr(pstr, ',');
      if(pchr == NULL)
      {
        fline.push_back(pstr);
        break;
      }
      else 
      {
        strncpy(buf1, pstr, pchr-pstr);
        buf1[pchr-pstr] = '\0';
        fline.push_back(buf1);
        pstr = pchr+1;
      }     
    }
  }
  cout<<fline.size()<<" line(s) input."<<endl;
  return;
}

void Config::parse_fixed_codes_str(const string &fixed_codes_str)
{
  fixed_codes.clear();
  if(fixed_codes_str.empty())
    return;
  
  string strbuf(fixed_codes_str);
  char buf1[256];
  char *pstr, *pchr;
  
  pstr = (char *)strbuf.c_str();
  pchr = strchr(pstr, ',');
  if(pchr == NULL)
  {
    fixed_codes.push_back(strtol(pstr, NULL, 10));
  }
  else 
  {
    strncpy(buf1, pstr, pchr-pstr);
    buf1[pchr-pstr] = '\0';
    fixed_codes.push_back(strtol(buf1, NULL, 10));
    pstr = pchr + 1;
    while(1)
    {
      pchr = strchr(pstr, ',');
      if(pchr == NULL)
      {
        fixed_codes.push_back(strtol(pstr, NULL, 10));
        break;
      }
      else 
      {
        strncpy(buf1, pstr, pchr-pstr);
        buf1[pchr-pstr] = '\0';
        fixed_codes.push_back(strtol(buf1, NULL, 10));
        pstr = pchr+1;
      }     
    }
  }
  return;
}

void Config::parse_fixed_scalecodes_str(const string &fixed_scalecodes_str)
{
  fixed_scalecodes.clear();
  if(fixed_scalecodes_str.empty())
    return;
  
  string strbuf(fixed_scalecodes_str);
  char buf1[256];
  char *pstr, *pchr;
  
  pstr = (char *)strbuf.c_str();
  pchr = strchr(pstr, ',');
  if(pchr == NULL)
  {
    fixed_scalecodes.push_back(strtol(pstr, NULL, 10));
  }
  else 
  {
    strncpy(buf1, pstr, pchr-pstr);
    buf1[pchr-pstr] = '\0';
    fixed_scalecodes.push_back(strtol(buf1, NULL, 10));
    pstr = pchr + 1;
    while(1)
    {
      pchr = strchr(pstr, ',');
      if(pchr == NULL)
      {
        fixed_scalecodes.push_back(strtol(pstr, NULL, 10));
        break;
      }
      else 
      {
        strncpy(buf1, pstr, pchr-pstr);
        buf1[pchr-pstr] = '\0';
        fixed_scalecodes.push_back(strtol(buf1, NULL, 10));
        pstr = pchr+1;
      }     
    }
  }
  return;
}

void Config::setup(const string& fcont_in, const list<string>& fline_in, 
             int nmcmc_in, double ptol_in, 
             double scale_range_low_in, double scale_range_up_in,
             double shift_range_low_in, double shift_range_up_in,
             double syserr_range_low_in, double syserr_range_up_in,
             double errscale_range_low_in, double errscale_range_up_in,
             double sigma_range_low_in, double sigma_range_up_in,
             double tau_range_low_in, double tau_range_up_in,
             bool fixed_scale_in, bool fixed_shift_in,
             bool fixed_syserr_in, bool fixed_error_scale_in,
             const vector<int>& fixed_codes_in,
             const vector<int>& fixed_scalecodes_in,
             bool flag_norm_in)
{
  strcpy(fcont, fcont_in.c_str());
  fline=fline_in;
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

  fixed_codes = fixed_codes_in;
  fixed_scalecodes = fixed_scalecodes_in;

  flag_norm = flag_norm_in;

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

void Config::check_directory()
{
  /* check if ./data exists
   * if not, create it;
   * if exists, check if it is a directory;
   * if not, throw an error.*/
  struct stat st;
  int status;
  status = stat("./data", &st);
  if(status != 0)
  {
    cout<<"================================"<<endl
        <<"Directory './data' not exist! PyCALI create it."<<endl;
    status = mkdir("./data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(status!=0)
    {
      cout<<"Cannot create './data'"<<endl
          <<"================================"<<endl;
    }
  }
  else
  {
    if(!S_ISDIR(st.st_mode))
    {
      cout<<"================================"<<endl
          <<"'./data' is not a direcotry!"<<endl
          <<"================================"<<endl;
      exit(-1);
    }
  }
  return;
}

void Config::print_cfg()
{
  list<string>::iterator it;
  vector<int>::iterator ic;
  
  check_directory();

  cout<<"=======Input parameters========="<<endl;
  cout<<setw(20)<<"fname: "<<fname<<endl;
  cout<<setw(20)<<"fcont: "<<fcont<<endl;
  if(fline.empty())
  {
    cout<<setw(20)<<"fline: "<<endl;
  }
  else 
  {
    cout<<setw(20)<<"fline: ";
    it = fline.begin();
    cout<<*it;
    for(++it; it != fline.end(); ++it)
      cout<<","<<*it;
    cout<<endl;
  }
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
  if(fixed_codes.empty())
  {
    cout<<setw(20)<<"fixed_codes: "<<endl;
  }
  else
  {
    cout<<setw(20)<<"fixed_codes: ";
    ic = fixed_codes.begin();
    cout<<*ic;
    for(++ic; ic != fixed_codes.end(); ++ic)
      cout<<","<<*ic;
    cout<<endl;

  }
  if(fixed_scalecodes.empty())
  {
    cout<<setw(20)<<"fixed_scalecodes: "<<endl;
  }
  else
  {
    cout<<setw(20)<<"fixed_scalecodes: ";
    ic = fixed_scalecodes.begin();
    cout<<*ic;
    for(++ic; ic != fixed_scalecodes.end(); ++ic)
      cout<<","<<*ic;
    cout<<endl;

  }
  cout<<setw(20)<<"flag_norm: "<<flag_norm<<endl;
  cout<<"================================"<<endl;

  ofstream fout;
  fout.open("data/param_input");
  fout<<setw(20)<<left<<"fname"<<" = "<<fname<<endl;
  fout<<setw(20)<<left<<"fcont"<<" = "<<fcont<<endl;
  if(fline.empty())
  {
    fout<<setw(20)<<left<<"fline"<<" = "<<endl;
  }
  else 
  {
    fout<<setw(20)<<left<<"fline"<<" = ";
    it = fline.begin();
    fout<<*it;
    for(++it; it!=fline.end(); ++it)
      fout<<","<<*it;
    fout<<endl;
  }
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
  if(fixed_codes.empty())
  {
    fout<<setw(20)<<"fixed_codes"<<" = "<<endl;
  }
  else
  {
    fout<<setw(20)<<"fixed_codes"<<" = ";
    ic = fixed_codes.begin();
    fout<<*ic;
    for(++ic; ic != fixed_codes.end(); ++ic)
      fout<<","<<*ic;
    fout<<endl;

  }
  if(fixed_scalecodes.empty())
  {
    fout<<setw(20)<<"fixed_scalecodes"<<" = "<<endl;
  }
  else
  {
    fout<<setw(20)<<"fixed_scalecodes"<<" = ";
    ic = fixed_scalecodes.begin();
    fout<<*ic;
    for(++ic; ic != fixed_scalecodes.end(); ++ic)
      fout<<","<<*ic;
    fout<<endl;

  }
  fout<<setw(20)<<left<<"flag_norm"<<" = "<<flag_norm<<endl;
  fout.flush();
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
  load(fname, true);
}

Data::Data(const string& fname, bool flag_norm)
{
  load(fname, flag_norm);
}

Data::~Data()
{
  flux_org.clear();
  error_org.clear();
  time.clear();
  flux.clear();
  error.clear();
  code.clear();
  
  syserrflag.clear();
  syserrflag_list.clear();

  num_code.clear();
  num_flag.clear();
  mean_code.clear();
  code_list.clear();
}

void Data::load(const string& fname, bool flag_norm=true)
{
  /* first clear all vectors */
  flux_org.clear();
  error_org.clear();
  time.clear();
  flux.clear();
  error.clear();
  code.clear();

  syserrflag.clear();
  syserrflag_list.clear();

  num_code.clear();
  num_flag.clear();
  mean_code.clear();
  code_list.clear();

  /* now read data */
  fstream fin;
  string line;
  stringstream ss;
  ss.clear();

  fin.open(fname);
  if(fin.fail())
  {
    cout<<"cannot open file "<<fname<<endl;
    exit(-1);
  }
  cout<<fname<<endl;

  int i, j, num;
  int idx, idx_str;
  double t, f, e, mean;
  int flag; /* flag */
  bool ifflag, setflag;
  vector<int>::iterator it;
  vector<int> flag_list; /* recode flags of each code */
  string cstr;
  
  idx = 0;
  setflag = false; /* not set flag at begining */ 
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
    cout<<idx<<"  "<<code_list[idx]<<"   "<<num_code[idx]<<endl;
    mean = 0.0;
    flag_list.clear();
    for(j=0; j<num_code[idx]; j++)
    {
      getline(fin, line);
      if(fin.fail())
      {
        cout<<"# Wrong in reading "<<fname<<endl;
        exit(-1);
      }
      
      if(setflag == false) /* set flag at first read */
      {
        ss.str(line);
        ss>>t>>f>>e>>flag;

        if(ss.fail())
        {
          ss.clear();
          ss.str(line);
          ss>>t>>f>>e;
          if(ss.fail())
          {
            cout<<"# Wrong in reading "<<fname<<endl;
            exit(-1);
          }
          else 
          {
            ifflag = false; /* no flag */
          }
        }
        else 
        {
          ifflag = true;
        }
        setflag=true;
      }
      
      /* now do read */
      ss.clear();
      ss.str(line);
      if(ifflag==true)
      {
        ss>>t>>f>>e>>flag;
      }
      else 
      {
        ss>>t>>f>>e;
        flag = 1; /* all point use the same flag */
      }
      if(ss.fail())
      {
        cout<<"# Wrong in reading "<<fname<<endl;
        exit(-1);
      }
      
      /* check if the flag exists */
      // it = find(syserrflag_list.begin(), syserrflag_list.end(), flag);
      // if(it == syserrflag_list.end())  /* not exist */
      // {
      //   syserrflag.push_back(syserrflag_list.size());
      //   syserrflag_list.push_back(flag);
      // }
      // else /* exist, set the index of the flag */
      // {
      //   syserrflag.push_back(it-syserrflag_list.begin());
      // }
      it = find(flag_list.begin(), flag_list.end(), flag);
      if(it == flag_list.end())  /* not exist */
      {
        syserrflag.push_back(flag_list.size());
        flag_list.push_back(flag);
      }
      else /* exist, set the index of the flag */
      {
        syserrflag.push_back(it-flag_list.begin());
      }

      time.push_back(t);
      flux_org.push_back(f);
      error_org.push_back(e);
      code.push_back(idx);
      ss.clear();

      mean += f;
    }
    
    num_flag.push_back(flag_list.size());
    /* store the flags */
    if(syserrflag_list.size() < flag_list.size())
    {
      for(i=syserrflag_list.size(); i<flag_list.size(); i++)
      {
        syserrflag_list.push_back(i);
      }
    }

    /* cope with the case of zero point */
    if(num_code[idx]>0)
      mean /= num;
    else
      mean = 1.0;
    
    if(mean <= 0.0)
    {
      cout<<"Error in code "<<j<<": negative mean flux!"<<endl;
      exit(-1);
    }
    mean_code.push_back(mean);
    idx++;
  }
  if(idx==0)
  {
    cout<<"Error: an empty file "<<fname<<endl;
    exit(-1);
  }
  
  /* check if the errors are positive */
  for(j=0; j<time.size(); j++)
  {
    if(error_org[j]<=0.0)
    {
      cout<<"Zero or negative error in code "<<code_list[code[j]]<<" of "<<fname<<"!"<<endl;
      exit(-1);
    }
  }
  cout<<"  "<<time.size()<<" points, "<<code_list.size()<<" codes, "<<syserrflag_list.size()<<" flags."<<endl;
  cout<<"================================"<<endl;
  fin.close();
  
  /* if not do normalization, mean_code set to that of 0th code */
  if(!flag_norm)
  {
    for(j=1; j<mean_code.size(); j++)
    {
      mean_code[j] = mean_code[0];
    }
  }

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
  vector<int> syserrflag_tmp;

  index.resize(time.size());
  iota(index.begin(), index.end(), 0);
  stable_sort(index.begin(), index.end(), [&](size_t i1, size_t i2) {return time[i1] < time[i2];});
  
  time_tmp = time;
  flux = flux_org;
  error = error_org;
  code_tmp = code;
  syserrflag_tmp = syserrflag;

  int i;
  for(i=0; i<index.size(); i++)
  {
    time[i] = time_tmp[index[i]];
    flux_org[i] = flux[index[i]];
    error_org[i] = error[index[i]];
    code[i] = code_tmp[index[i]];
    syserrflag[i] = syserrflag_tmp[index[i]];
  }
  code_tmp.clear();
  time_tmp.clear();
  syserrflag_tmp.clear();
}

void Data::check_code(Data& data)
{
  int i;
  if(code_list.size() != data.code_list.size())
  {
    cout<<"# Numbers of codes do not match.\n Note that a dataset is permitted to have none points."<<endl;
    exit(-1);
  }

  for(i=0; i<code_list.size(); i++)
  {
    if(code_list[i] != data.code_list[i])
    {
      cout<<"# Codes do not match or in different orders.\n Note that a dataset is permitted to have none points."<<endl;
      exit(-1);
    }
  }
}
/*=====================================================*/
/* class for calibration */
Cali::Cali()
{
  check_directory();

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
     :fcont(cfg.fcont), fline(cfg.fline), cont(cfg.fcont, cfg.flag_norm),
      nmcmc(cfg.nmcmc), ptol(cfg.ptol)
{
  int i, j, m, k;
  bool isfixed;
  
  check_directory();

  num_params_var = 2;
  size_max = cont.time.size();
  ncode = cont.code_list.size();
  nsyserr_flag = cont.syserrflag_list.size();
  if(!fline.empty())
  {
    Data line;
    list<string>::iterator it; 
    for(it=fline.begin(); it!=fline.end(); ++it)
    {
      line.load(*it, cfg.flag_norm);
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
        line.mean_code[i] /= (line.mean_code[i]/line.mean_code[0]) * (cont.mean_code[0]/cont.mean_code[i]);
      }
      lines.push_back(line);
      
      /* use the largest number of flags */
      if(nsyserr_flag < line.syserrflag_list.size())
        nsyserr_flag = line.syserrflag_list.size();
    }
  }

  check_fixed_codes(cfg);
  check_fixed_scalecodes(cfg);
  
  /* variability, scale, shift, syserr, error scale */
  num_params = num_params_var + ncode*2 + (ncode + ncode)*nsyserr_flag;
  if(!fline.empty())
  {
    /* syserr and error scale of line */
    num_params += (ncode + ncode)*nsyserr_flag*lines.size();
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
  /* check tau_range */
  if(cfg.tau_range_low >= (cont.time[cont.time.size()-1]-cont.time[0]))
  {
    cfg.tau_range_low = (cont.time[cont.time.size()-1]-cont.time[0]) / cont.time.size();
    cout<<"# warning: tau_range_low is too large, reduce it to be time length of "<<cfg.fcont
        <<"="<<cfg.tau_range_low<<endl;
    //exit(-1);
  }
  if(cfg.tau_range_up <= (cont.time[cont.time.size()-1]-cont.time[0])/cont.time.size())
  {
    cfg.tau_range_up = (cont.time[cont.time.size()-1]-cont.time[0]);
    cout<<"# warning: tau_range_up is too small, increase it to be mean sampling interval of "<<cfg.fcont
        <<"="<<cfg.tau_range_up<<endl;
    //exit(-1);
  }
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
    list<Data>::iterator it;
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *it;
      if(cfg.tau_range_low >= (line.time[line.time.size()-1]-line.time[0]))
      {
        cfg.tau_range_low = (line.time[line.time.size()-1]-line.time[0])/line.time.size();
        cout<<"# warning: tau_range_low is too large, reduce it to be time length of "<<cfg.fcont
            <<"="<<cfg.tau_range_low<<endl;
        //exit(-1);
      }
      if(cfg.tau_range_up <= (line.time[line.time.size()-1]-line.time[0])/line.time.size())
      {
        cfg.tau_range_up = (line.time[line.time.size()-1]-line.time[0]);
        cout<<" warning: tau_range_up is too small, increase it to be mean sampling interval of "<<cfg.fcont
            <<"="<<cfg.tau_range_up<<endl;
        //exit(-1);
      }
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
    for(k=0; k<nsyserr_flag; k++)
    {
      i+=1;
      if(cont.num_code[j] > 5) /* for a large number, use a uniform prior */
      {
        par_range_model[i][0] = cfg.syserr_range_low;
        par_range_model[i][1] = cfg.syserr_range_up;
        par_prior_model[i] = UNIFORM;
      }
      else /* for a small number, use the Gaussian prior to preset a strong constraint */
      {
        par_range_model[i][0] = cfg.syserr_range_low;
        par_range_model[i][1] = cfg.syserr_range_up;
        par_prior_model[i] = GAUSSIAN;
        par_prior_gaussian[i][0] =  cfg.syserr_range_low;
        par_prior_gaussian[i][1] = (cfg.syserr_range_up - cfg.syserr_range_low)/3.0;
      }
    }
  }
  /* error scale of continuum */
  for(j=0; j<ncode*nsyserr_flag; j++)
  {
    i+=1;
    par_range_model[i][0] = cfg.errscale_range_low;
    par_range_model[i][1] = cfg.errscale_range_up;
    par_prior_model[i] = LOG;
  }

  if(!fline.empty())
  {
    list<Data>::iterator it;
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      /* syserr of line */
      for(j=0; j<ncode; j++)
      {
        for(k=0; k<nsyserr_flag; k++)
        {
          i+=1;
          if(line.num_code[j] > 5) /* for a large number, use a uniform prior */
          {
            par_range_model[i][0] = cfg.syserr_range_low;
            par_range_model[i][1] = cfg.syserr_range_up;
            par_prior_model[i] = UNIFORM;
          }
          else  /* for a small number, use the Gaussian prior to preset a strong constraint */
          {
            par_range_model[i][0] = cfg.syserr_range_low;
            par_range_model[i][1] = cfg.syserr_range_up;
            par_prior_model[i] = GAUSSIAN;
            par_prior_gaussian[i][0] =  cfg.syserr_range_low;
            par_prior_gaussian[i][1] = (cfg.syserr_range_up - cfg.syserr_range_low)/3.0;
          }
        }
      }
      /* error scale of line */
      for(j=0; j<ncode*nsyserr_flag; j++)
      {
        i+=1;
        par_range_model[i][0] = cfg.errscale_range_low;
        par_range_model[i][1] = cfg.errscale_range_up;
        par_prior_model[i] = LOG;
      }
    }
  }
  
  for(i=0; i<num_params; i++)
  {
    par_fix[i] = NOFIXED;
    par_fix_val[i] = -DBL_MAX;
  }
  /* scale and shift of the 1st code are always fixed */
  par_fix[num_params_var] = FIXED;
  par_fix_val[num_params_var] = 1.0;
  par_fix[num_params_var+ncode] = FIXED;
  par_fix_val[num_params_var+ncode] = 0.0;
  
  bool flag;
  /* for empty continuum dataset, shifts are fixed */
  for(i=0; i<ncode; i++)
  {
    if(cont.num_code[i] == 0)
    {
      par_fix[num_params_var    +i+ncode] = FIXED;
      par_fix_val[num_params_var+i+ncode] = 0.0;

      if(fline.empty()) /* if no lines, scale also fixed */
      {
        par_fix[num_params_var    +i] = FIXED;
        par_fix_val[num_params_var+i] = 1.0;
      }
      else 
      {
        list<Data>::iterator it;
        it = lines.begin();
        flag = true;
        for(j=0; j<fline.size(); j++)
        {
          Data& line = *(it);
          if(line.num_code[i] != 0)
          {
            flag = false;
            break;
          }
          it++;
        }
        /* for i-th code, all lines have no points, scale fixed. */
        if(flag == true)
        {
          par_fix[num_params_var    +i] = FIXED;
          par_fix_val[num_params_var+i] = 1.0;
        }
      }
    }
  }

  if(cfg.fixed_scale)
  {
    for(i=1; i<ncode; i++)
    {
      par_fix[num_params_var+i] = FIXED;
      par_fix_val[num_params_var+i] = 1.0;
    }
  }
  for(i=0; i<cfg.fixed_scalecodes.size();i++) /* fix scale of specific codes */
  {
    m = cfg.fixed_scalecodes[i];
    par_fix[num_params_var+m] = FIXED;
    /* take into account different normalization of codes */
    par_fix_val[num_params_var+m] = 1.0;
  }
  for(i=0; i<cfg.fixed_codes.size();i++) /* fix the specific codes */
  {
    m = cfg.fixed_codes[i];
    par_fix[num_params_var+m] = FIXED;
    /* take into account different normalization of codes, otherwise, the final flux will change 
     * because all the parameters will be fixed.
     */
    par_fix_val[num_params_var+m] = cont.mean_code[m]/cont.norm;
  }

  if(cfg.fixed_shift)
  {
    for(i=1; i<ncode; i++)
    {
      par_fix[num_params_var+i+ncode] = FIXED;
      par_fix_val[num_params_var+i+ncode] = 0.0;
    }
  }
  for(i=0; i<cfg.fixed_codes.size();i++) /* fix the specific codes */
  {
    m = cfg.fixed_codes[i];
    par_fix[num_params_var+ncode+m] = FIXED;
    par_fix_val[num_params_var+ncode+m] = 0.0;
  }

  /* if the number point of continuum <= 2 and scale is not fixed, 
     and number point of line is nonzero, fix shift 
   */
  for(i=1; i<ncode; i++)
  {
    flag = true;
    if(!fline.empty())
    {
      list<Data>::iterator it;
      it = lines.begin();
      for(j=0; j<fline.size(); j++)
      {
        Data& line = *(it);
        if(line.num_code[i] > 0)
        {
          flag = false;
          break;
        }
        it++;
      }
    }
    if(cont.num_code[i] <= 2 && flag == true)
    {
      if(!cfg.fixed_scale)
      {
        par_fix[num_params_var    +i+ncode] = FIXED;
        par_fix_val[num_params_var+i+ncode] = 0.0;
      }
    }
  }

  if(cfg.fixed_syserr) /* if syserr is fixed */
  {
    for(i=0; i<ncode*nsyserr_flag; i++)
    {
      par_fix[num_params_var    +2*ncode+i] = FIXED;
      par_fix_val[num_params_var+2*ncode+i] = 0.0;
    }

    if(!fline.empty())
    {
      for(j=0; j<lines.size(); j++)
      {
        for(i=0; i<ncode*nsyserr_flag; i++)
        {
          par_fix[num_params_var    +2*ncode+(2+2*j)*nsyserr_flag*ncode+i] = FIXED;
          par_fix_val[num_params_var+2*ncode+(2+2*j)*nsyserr_flag*ncode+i] = 0.0;
        }
      }
    }
  }
  else  /* if syserr is not fixed */
  {
    /* syserrs of flags not in the flags of a code is fixed  */
    for(i=0; i<ncode; i++)
    {
      for(k=cont.num_flag[i]; k<nsyserr_flag; k++)
      {
        par_fix[num_params_var    +2*ncode+i*nsyserr_flag+k] = FIXED;
        par_fix_val[num_params_var+2*ncode+i*nsyserr_flag+k] = 0.0;
      }
    }
    
    /* now fix the specific codes */
    for(i=0; i<cfg.fixed_codes.size();i++) 
    {
      m = cfg.fixed_codes[i];
      for(k=0; k<nsyserr_flag; k++)
      {
        par_fix[num_params_var    +2*ncode+m*nsyserr_flag+k] = FIXED;
        par_fix_val[num_params_var+2*ncode+m*nsyserr_flag+k] = 0.0;
      }
    }

    if(!fline.empty())
    {
      list<Data>::iterator it;
      it = lines.begin();
      for(j=0; j<lines.size(); j++)
      {
        /* syserrs of flags not in the flags of a code is fixed  */
        Data& line = *(it);
        for(i=0; i<ncode; i++)
        {
          for(k=line.num_flag[i]; k<nsyserr_flag; k++)
          {
            par_fix[num_params_var    +2*ncode+(2+2*j)*nsyserr_flag*ncode+i*nsyserr_flag+k] = FIXED;
            par_fix_val[num_params_var+2*ncode+(2+2*j)*nsyserr_flag*ncode+i*nsyserr_flag+k] = 0.0;
          }
        }
        ++it;
        
        /* fix the specific codes */
        for(i=0; i<cfg.fixed_codes.size();i++) 
        {
          m = cfg.fixed_codes[i];
          for(k=0; k<nsyserr_flag; k++)
          {
            par_fix[num_params_var    +2*ncode+(2+2*j)*nsyserr_flag*ncode+m*nsyserr_flag+k] = FIXED;
            par_fix_val[num_params_var+2*ncode+(2+2*j)*nsyserr_flag*ncode+m*nsyserr_flag+k] = 0.0;
          }
        }
      }
    }
  }

  if(cfg.fixed_error_scale) /* if error scale is fixed */
  {
    for(i=0; i<ncode*nsyserr_flag; i++)
    {
      par_fix[num_params_var    +2*ncode+ncode*nsyserr_flag+i] = FIXED;
      par_fix_val[num_params_var+2*ncode+ncode*nsyserr_flag+i] = 1.0;
    }

    if(!fline.empty())
    {
      for(j=0; j<lines.size(); j++)
      {
        for(i=0; i<ncode*nsyserr_flag; i++)
        {
          /* cont: 2xnxns; j-th line: cont + 2xjxnxns + nxns = (3 + 2xj)xnxns */
          par_fix[num_params_var    +2*ncode+(3+2*j)*ncode*nsyserr_flag+i] = FIXED;
          par_fix_val[num_params_var+2*ncode+(3+2*j)*ncode*nsyserr_flag+i] = 1.0;
        }
      }
    }
  }
  else /* if error scale is not fixed */
  {
    for(i=0; i<ncode; i++)
    {
      /* error scale of flags not in the flags of a code is fixed  */
      for(k=cont.num_flag[i]; k<nsyserr_flag; k++)
      {
        par_fix[num_params_var    +2*ncode+ncode*nsyserr_flag+i*nsyserr_flag+k] = FIXED;
        par_fix_val[num_params_var+2*ncode+ncode*nsyserr_flag+i*nsyserr_flag+k] = 1.0;
      }
    }
    
    /* fix the specific codes */
    for(i=0; i<cfg.fixed_codes.size();i++) 
    {
      m = cfg.fixed_codes[i];
      for(k=0; k<nsyserr_flag; k++)
      {
        par_fix[num_params_var    +2*ncode+ncode*nsyserr_flag+m*nsyserr_flag+k] = FIXED;
        par_fix_val[num_params_var+2*ncode+ncode*nsyserr_flag+m*nsyserr_flag+k] = 1.0;
      }
    }

    if(!fline.empty())
    {
      list<Data>::iterator it;
      it = lines.begin();
      for(j=0; j<lines.size(); j++)
      {
        /* error scale of flags not in the flags of a code is fixed  */
        Data& line = *(it);
        for(i=0; i<ncode; i++)
        {
          for(k=line.num_flag[i]; k<nsyserr_flag; k++)
          {
            par_fix[num_params_var    +2*ncode+(3+2*j)*ncode*nsyserr_flag+i*nsyserr_flag+k] = FIXED;
            par_fix_val[num_params_var+2*ncode+(3+2*j)*ncode*nsyserr_flag+i*nsyserr_flag+k] = 1.0;
          }
        }
        ++it;
        
        /* fix the specific codes */
        for(i=0; i<cfg.fixed_codes.size();i++) 
        {
          m = cfg.fixed_codes[i];
          for(k=0; k<nsyserr_flag; k++)
          {
            par_fix[num_params_var    +2*ncode+(3+2*j)*ncode*nsyserr_flag+m*nsyserr_flag+k] = FIXED;
            par_fix_val[num_params_var+2*ncode+(3+2*j)*ncode*nsyserr_flag+m*nsyserr_flag+k] = 1.0;
          }
        }
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
  cont_recon.resize(fmax(cont.time.size()*2, 100));
  size_recon_max = cont_recon.time.size();
  for(i=0; i<cont_recon.time.size(); i++)
  {
    cont_recon.time[i] = t1 + (t2 - t1)/(cont_recon.time.size()-1.0) * i;
  }
  if(!fline.empty())
  {
    list<Data>::iterator it;
    DataLC line_recon;
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      line_recon.resize(fmax(line.time.size()*2, 100));
      tspan = line.time[line.time.size()-1] - line.time[0];
      t1 = line.time[0] - 0.05*tspan;
      t2 = line.time[line.time.size()-1] + 0.05*tspan;
      for(i=0; i<line_recon.time.size(); i++)
      {
        line_recon.time[i] = t1 + (t2 - t1)/(line_recon.time.size()-1.0) * i;
      }
      size_recon_max = fmax(size_recon_max, line_recon.time.size());

      lines_recon.push_back(line_recon);
    }
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

  lines.clear();
  lines_recon.clear();
}

void Cali::check_directory()
{
  /* check if ./data exists
   * if not, create it;
   * if exists, check if it is a directory;
   * if not, throw an error.*/
  struct stat st;
  int status;
  status = stat("./data", &st);
  if(status != 0)
  {
    cout<<"================================"<<endl
        <<"Directory './data' not exist! PyCALI create it."<<endl;
    status = mkdir("./data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(status!=0)
    {
      cout<<"Cannot create './data'"<<endl
          <<"================================"<<endl;
    }
  }
  else
  {
    if(!S_ISDIR(st.st_mode))
    {
      cout<<"================================"<<endl
          <<"'./data' is not a direcotry!"<<endl
          <<"================================"<<endl;
      exit(-1);
    }
  }
  return;
}

void Cali::check_fixed_codes(Config& cfg)
{
  int i;
  for(i=0; i<cfg.fixed_codes.size(); i++)
  {
    if(cfg.fixed_codes[i]<0 || cfg.fixed_codes[i]>=cont.code_list.size())
    {
      cout<<"Incorrect fixed code number "<<cfg.fixed_codes[i]<<endl;
      exit(-1);
    }
  }
  return;
}

void Cali::check_fixed_scalecodes(Config& cfg)
{
  int i;
  for(i=0; i<cfg.fixed_scalecodes.size(); i++)
  {
    if(cfg.fixed_scalecodes[i]<0 || cfg.fixed_scalecodes[i]>=cont.code_list.size())
    {
      cout<<"Incorrect fixed code number "<<cfg.fixed_scalecodes[i]<<endl;
      exit(-1);
    }
  }
  return;
}

void Cali::align(double *model)
{
  int i, idx, idx_err;
  double *ps_scale = model+num_params_var;
  double *es_shift = ps_scale + ncode;
  double *syserr = es_shift + ncode;
  double *error_scale = syserr + ncode*nsyserr_flag;
  for(i=0; i<cont.time.size(); i++)
  {
    idx = cont.code[i];
    cont.flux[i] = cont.flux_org[i] * ps_scale[idx] - es_shift[idx];
    /* note that this error does not include errors of scale and shift */
    idx_err = idx*nsyserr_flag + cont.syserrflag[i];
    cont.error[i] = sqrt(cont.error_org[i]*cont.error_org[i]*error_scale[idx_err]*error_scale[idx_err] 
                    + syserr[idx_err]*syserr[idx_err]) * ps_scale[idx];
  }

  if(!fline.empty())
  {
    list<Data>::iterator it;
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      syserr = error_scale + ncode*nsyserr_flag;
      error_scale = syserr + ncode*nsyserr_flag;
      for(i=0; i<line.time.size(); i++)
      {
        idx = line.code[i];
        line.flux[i] = line.flux_org[i] * ps_scale[idx];

        idx_err = idx*nsyserr_flag + line.syserrflag[i];
        line.error[i] = sqrt(line.error_org[i]*line.error_org[i]*error_scale[idx_err]*error_scale[idx_err] 
                        + syserr[idx_err]*syserr[idx_err] ) * ps_scale[idx];
      }
    }
  }
}

void Cali::align_cont(double *model)
{
  int i, idx, idx_err;
  double *ps_scale = model+num_params_var;
  double *es_shift = ps_scale + ncode;
  double *syserr = es_shift + ncode;
  double *error_scale = syserr + ncode*nsyserr_flag;
  for(i=0; i<cont.time.size(); i++)
  {
    idx = cont.code[i];
    cont.flux[i] = cont.flux_org[i] * ps_scale[idx] - es_shift[idx];
    /* note that this error does not include errors of scale and shift */
    idx_err = idx*nsyserr_flag + cont.syserrflag[i];
    cont.error[i] = sqrt(cont.error_org[i]*cont.error_org[i]*error_scale[idx_err]*error_scale[idx_err] 
                    + syserr[idx_err]*syserr[idx_err]) * ps_scale[idx];
  }
  return;
}

/* il-th line, counting from 0 */
void Cali::align_line(double *model, int il)
{
  int i, idx, idx_err;
  double *ps_scale = model+num_params_var;
  double *es_shift = ps_scale + ncode;
  double *syserr = es_shift + ncode;                 /* continuum's */
  double *error_scale = syserr + ncode*nsyserr_flag; /* continuum's */
  
  if(!fline.empty())
  {
    list<Data>::iterator it = lines.begin();
    syserr = error_scale + ncode*nsyserr_flag;   /* line's begin */
    error_scale = syserr + ncode*nsyserr_flag;   /* line's begin */
    for(i=0; i<il; i++)
    {
      ++it;
      syserr = error_scale + ncode*nsyserr_flag;
      error_scale = syserr + ncode*nsyserr_flag;
    }
    Data& line = *(it);
    for(i=0; i<line.time.size(); i++)
    {
      idx = line.code[i];
      line.flux[i] = line.flux_org[i] * ps_scale[idx];

      idx_err = idx*nsyserr_flag + line.syserrflag[i];
      line.error[i] = sqrt(line.error_org[i]*line.error_org[i]*error_scale[idx_err]*error_scale[idx_err] 
                      + syserr[idx_err]*syserr[idx_err] ) * ps_scale[idx];
    }
  }
  return;
}

void Cali::align_with_error()
{
  int i, idx, idx_err;
  double *ps_scale = best_params+num_params_var;
  double *es_shift = ps_scale + ncode;
  double *syserr = es_shift + ncode;
  double *error_scale = syserr + ncode*nsyserr_flag;
  double *ps_scale_err = best_params_std + num_params_var;
  double *es_shift_err = ps_scale_err + ncode;
  for(i=0; i<cont.time.size(); i++)
  {
    idx = cont.code[i];
    cont.flux[i] = cont.flux_org[i] * ps_scale[idx] - es_shift[idx];

    idx_err = idx*nsyserr_flag + cont.syserrflag[i];
    cont.error[i] = sqrt((cont.error_org[i]*cont.error_org[i]*error_scale[idx_err]*error_scale[idx_err] +syserr[idx_err]*syserr[idx_err]) 
                         *ps_scale[idx]*ps_scale[idx]
                        +pow(cont.flux_org[i]*ps_scale_err[idx], 2.0)
                        +pow(es_shift_err[idx], 2.0)
                        -2.0*cont.flux_org[i]*best_params_covar[(num_params_var+idx)*num_params + (num_params_var+idx+ncode)]
                        );
  }

  if(!fline.empty())
  {
    list<Data>::iterator it;
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      syserr = error_scale + ncode*nsyserr_flag;
      error_scale = syserr + ncode*nsyserr_flag;
      for(i=0; i<line.time.size(); i++)
      {
        idx = line.code[i];
        line.flux[i] = line.flux_org[i] * ps_scale[idx];

        idx_err = idx*nsyserr_flag + line.syserrflag[i];
        line.error[i] = sqrt((line.error_org[i]*line.error_org[i]*error_scale[idx_err]*error_scale[idx_err] + syserr[idx_err]*syserr[idx_err]) 
                            *ps_scale[idx]*ps_scale[idx]
                            +pow(line.flux_org[i] * ps_scale_err[idx], 2.0)
                            );
      }
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
  //strcpy(argv[argc++], "-l");  //level-dependent sampling

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
  int i, j, k, num_ps;
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
  cout<<"================================"<<endl;
  cout<<"Lmax:"<<prob_max<<" at "<<ip_max<<"th posterior sample."<<endl;
  printf("The params with the maximum likelihood:\n");
  printf("scale and shift\n");
  for(i=0; i<ncode; i++)
  {
    cout<<scientific
        <<cont.code_list[i]<<"\t"<<*((double *)posterior_sample + ip_max*num_params + num_params_var + i)
                           <<"\t"<<*((double *)posterior_sample + ip_max*num_params + num_params_var + i + ncode)<<endl;
  }
  printf("syserr and error scale of continuum\n");
  for(i=0; i<ncode; i++)
  {
    for(k=0; k<nsyserr_flag; k++)
    {
      cout<<scientific<<cont.code_list[i]
          <<"\t"<<k<<"\t"<<*((double *)posterior_sample + ip_max*num_params + num_params_var + k + i*nsyserr_flag + 2*ncode)
          <<"\t"<<*((double *)posterior_sample + ip_max*num_params + num_params_var + k + i*nsyserr_flag + ncode*nsyserr_flag + 2*ncode)
          <<endl;
    }
  }

  if(!fline.empty())
  {
    for(j=0; j<fline.size(); j++)
    {
      printf("syserr and error scale of line %d\n", j);
      for(i=0; i<ncode; i++)
      {
        for(k=0; k<nsyserr_flag; k++)
        {
          cout<<scientific<<cont.code_list[i]
            <<"\t"<<k<<"\t"<<*((double *)posterior_sample + ip_max*num_params + num_params_var + k + i*nsyserr_flag + (2+2*j)*ncode*nsyserr_flag + 2*ncode)
            <<"\t"<<*((double *)posterior_sample + ip_max*num_params + num_params_var + k + i*nsyserr_flag + (3+2*j)*ncode*nsyserr_flag + 2*ncode)
            <<endl;
        }
      }
    }
  }

  /* calculate covariance */
  double covar;
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
  stat_type = 2;  /* 0: median, 1: mean, 2, error peak */
  double error_scale_shift;
  DataLC cont_output(cont.time.size());
  if(stat_type == 0)  /* median values */
  {
    double *flux, *error;
    flux = new double [cont.time.size()*num_ps];
    error = new double [cont.time.size()*num_ps];

    for(i=0; i<num_ps; i++)
    {
      align_cont((double *)posterior_sample + i*num_params);
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
      list<Data>::iterator it;
      int il=0;
      for(it=lines.begin(); it!=lines.end(); ++it)
      {
        Data& line = *(it);
        DataLC line_output(line.time.size());
        flux = new double [line.time.size()*num_ps];
        error = new double [line.time.size()*num_ps];
        for(i=0; i<num_ps; i++)
        {
          align_line((double *)posterior_sample + i*num_params, il);
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
        il++;
      }
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
      align_cont((double *)posterior_sample + i*num_params);
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
      list<Data>::iterator it;
      int il=0;
      for(it=lines.begin(); it!=lines.end(); ++it)
      {
        Data& line = *(it);
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
          align_line((double *)posterior_sample + i*num_params, il);
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
        il++;
      }
    }
  }
  else if(stat_type == 2) /* use best likelihood-maximize values */
  {
    /* Gaussian filter */
    size_t nhist = 50; /* number of bins for histograms */
    size_t ngauss = 50;
    double alpha  = (ngauss-1)/2.0/2.0;
    gsl_filter_gaussian_workspace *gauss_p = gsl_filter_gaussian_alloc(ngauss);
    gsl_vector *hist_in = gsl_vector_alloc(nhist);
    gsl_vector *hist_out = gsl_vector_alloc(nhist);

    double *flux, *error;
    double error_min, error_max;
    gsl_histogram * he = gsl_histogram_alloc (nhist);

    flux = new double [cont.time.size()*num_ps];
    error = new double [cont.time.size()*num_ps];
    
    /* calculate errors of scale and shift */
    for(i=0; i<num_ps; i++)
    {
      align_cont((double *)posterior_sample + i*num_params);
      for(j=0; j<cont.time.size(); j++)
      {
        flux[j * num_ps + i] = cont.flux[j];
        error[j * num_ps + i] = cont.error[j];
      }
    }
    for(j=0; j<cont.time.size(); j++)
    {
      error_min =  DBL_MAX;
      error_max = -DBL_MAX;
      for(i=0; i<num_ps; i++)
      {
        error_min = fmin(error_min, error[j * num_ps + i]);
        error_max = fmax(error_max, error[j * num_ps + i]);
      }

      /* error is a single value, no need to do histogram */      
      if(error_min == error_max)
      {
        cont_output.error[j] = error_min;
      }
      else 
      {
        gsl_histogram_reset(he);
        gsl_histogram_set_ranges_uniform(he, error_min, error_max);
        for(i=0; i<num_ps; i++)
        {
          gsl_histogram_increment(he, error[j * num_ps + i]);
        }
        
        /* gaussian smooth */
        memcpy(hist_in->data, he->bin, nhist*sizeof(double));
        gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 0, hist_in, hist_out, gauss_p);
        memcpy(he->bin, hist_out->data, nhist*sizeof(double));
        
        /* locate the peak */
        cont_output.error[j] = 0.5*(he->range[gsl_histogram_max_bin(he)] + he->range[gsl_histogram_max_bin(he)+1]);
      }

      /* ascending order */
      qsort(flux+j*num_ps, num_ps, sizeof(double), compare);
      cont_output.flux[j] = flux[j*num_ps + (int)(0.5*num_ps)];
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
      list<Data>::iterator it;
      int il=0;
      for(it=lines.begin(); it!=lines.end(); ++it)
      {
        Data& line = *(it);
        DataLC line_output(line.time.size());

        flux = new double [line.time.size()*num_ps];
        error = new double [line.time.size()*num_ps];
        for(i=0; i<num_ps; i++)
        {
          align_line((double *)posterior_sample + i*num_params, il);
          for(j=0; j<line.time.size(); j++)
          {
            flux[j * num_ps + i] = line.flux[j];
            error[j * num_ps + i] = line.error[j];
          }
        }
  
        /* sort, and use the media value */
        for(j=0; j<line.time.size(); j++)
        {
          error_min =  DBL_MAX;
          error_max = -DBL_MAX;
          for(i=0; i<num_ps; i++)
          {
            error_min = fmin(error_min, error[j * num_ps + i]);
            error_max = fmax(error_max, error[j * num_ps + i]);
          }
          
          /* error is a single value, no need to do histogram */
          if(error_min == error_max)
          {
            line_output.error[j] = error_min;
          }
          else
          {
            gsl_histogram_reset(he);
            gsl_histogram_set_ranges_uniform(he, error_min, error_max);
            for(i=0; i<num_ps; i++)
            {
              gsl_histogram_increment(he, error[j * num_ps + i]);
            }
            /* gaussian smooth */
            memcpy(hist_in->data, he->bin, nhist*sizeof(double));
            gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 0, hist_in, hist_out, gauss_p);
            memcpy(he->bin, hist_out->data, nhist*sizeof(double));

            /* locate the peak */
            line_output.error[j] = 0.5*(he->range[gsl_histogram_max_bin(he)] + he->range[gsl_histogram_max_bin(he)+1]);
          }
          
          /* ascending order */
          qsort(flux+j*num_ps, num_ps, sizeof(double), compare);
          line_output.flux[j] = flux[j*num_ps + (int)(0.5*num_ps)];
          /* include error of scale and shift */
          error_scale_shift = 0.5*( (line_output.flux[j] - flux[j*num_ps + (int)(0.1585*num_ps)])
                                 +(flux[j*num_ps + (int)(0.8415*num_ps)] - line_output.flux[j]) );
          line_output.error[j] = sqrt(line_output.error[j]*line_output.error[j] + error_scale_shift*error_scale_shift);
        }

        line.flux = line_output.flux;
        line.error = line_output.error;
        delete[] flux;
        delete[] error;
        il++;
      }
    }
    gsl_histogram_free(he);
    gsl_vector_free(hist_in);
    gsl_vector_free(hist_out);
    gsl_filter_gaussian_free(gauss_p);
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
  int i, j, k;
  ofstream fout;
  fout.open(fcont+"_cali");
  for(i=0; i<cont.time.size(); i++)
  {
    fout<<fixed<<cont.time[i]
        <<scientific<<" "<<cont.flux[i]*cont.norm<<"  "<<cont.error[i]*cont.norm<<"  "<<cont.code_list[cont.code[i]]<<endl;
  }
  fout.flush();
  fout.close();

  /* output indices that sort cont */
  fout.open(fcont+"_sort");
  for(i=0; i<cont.time.size(); i++)
  {
    fout<<cont.index[i]<<endl;
  }
  fout.close();

  if(!fline.empty())
  {
    unsigned int il;
    list<Data>::iterator it;
    list<string>::iterator ifl;
    il = 0;
    ifl = fline.begin();
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      fout.open(*(ifl)+"_cali");
      for(i=0; i<line.time.size(); i++)
      {
        fout<<fixed<<line.time[i]
            <<scientific<<" "<<line.flux[i]*line.norm<<"  "<<line.error[i]*line.norm<<"  "<<line.code_list[line.code[i]]<<endl;
      }
      fout.close();
  
      /* output indices that sort line */
      fout.open(*(ifl)+"_sort");
      for(i=0; i<line.time.size(); i++)
      {
        fout<<line.index[i]<<endl;
      }
      fout.close();
      ifl++;
    }
  }

  fout.open("data/factor.txt");
  fout<<"# Note: factors apply to light curves after being nomoralized by their respective means."<<endl
      <<"#       Means are output to the file PyCALI_output.txt."<<endl;
  fout<<"Code \t Scale  \t Error  \t Shift  \t Error    \t     Cov"<<endl;
  for(i=0; i<ncode;i++)
  {
    fout<<scientific
        <<cont.code_list[i]<<"\t"<<best_params[i+num_params_var]<<"\t"<<best_params_std[i+num_params_var]
        <<"\t"<<best_params[i+ncode+num_params_var]<<"\t"<<best_params_std[i+ncode+num_params_var]
        <<"\t"<<best_params_covar[(i+num_params_var)*num_params+(i+num_params_var+ncode)];
        //<<"\t"<<best_params[i+2*ncode+num_params_var]
        //<<"\t"<<best_params[i+3*ncode+num_params_var];
    /*if(!fline.empty())
    {
      for(j=0; j<fline.size(); j++)
      {
        fout<<"\t"<<best_params[i+(4+2*j)*ncode+num_params_var]
            <<"\t"<<best_params[i+(5+2*j)*ncode+num_params_var];
      }
    }*/
    fout<<endl;
  }

  fout<<endl;
  fout<<"Code \tflag \tSyserr_Cont    \t Err_Scale";
  if(!fline.empty())
  {
    for(j=0; j<fline.size(); j++)
    {
      fout<<"\tSyserr_Line   \t Err_Scale";
    }
  }
  fout<<endl;
  for(i=0; i<ncode; i++)
  {
    for(k=0; k<nsyserr_flag; k++)
    {
      fout<<scientific
          <<cont.code_list[i]<<"\t"<<k<<"\t"
          <<best_params[num_params_var + 2*ncode +                      i*nsyserr_flag + k]<<"\t"
          <<best_params[num_params_var + 2*ncode + ncode*nsyserr_flag + i*nsyserr_flag + k];
      
      if(!fline.empty())
      {
        for(j=0; j<fline.size(); j++)
        {
          fout<<"\t"<<best_params[num_params_var + 2*ncode + (2+2*j)*ncode*nsyserr_flag + i*nsyserr_flag + k]
              <<"\t"<<best_params[num_params_var + 2*ncode + (3+2*j)*ncode*nsyserr_flag + i*nsyserr_flag + k];
        }
      }
      fout<<endl;
    }
  }

  fout.close();

  fout.open("data/PyCALI_output.txt");
  fout<<"# normalization factor of continuum: "<<fcont<<endl;
  for(i=0; i<ncode; i++)
  {
    fout<<i<<"\t"<<cont.code_list[i]<<"\t"<<cont.mean_code[i]<<"\t"<<cont.num_code[i]<<endl;
  }
  if(!fline.empty())
  {
    list<Data>::iterator it;
    list<string>::iterator ifl;
    ifl = fline.begin();
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      fout<<"# normalization factor of line: "<<*(ifl)<<endl;
      for(i=0; i<ncode; i++)
      {
        fout<<i<<"\t"<<line.code_list[i]<<"\t"<<line.mean_code[i]<<"\t"<<line.num_code[i]<<endl;
      }
      ifl++;
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
  // Cq^-1 = L^TxC^-1xL, Lbuf = C^-1xL
  multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_cont, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_cont);

  // L^TxC^-1xy
  multiply_matvec_MN_transposeA(Lbuf, nd_cont, nq, cont.flux.data(), yq);

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
  /* PEmat1 = C^-1 x S^T, dimension: nd*nrec */
  multiply_mat_transposeB_semiseparable_drw(USmat, W, D, phi, nd_cont, nd_cont_recon, sigma2, PEmat1);

  // (hat s) = SxC^-1xy
  multiply_matvec_MN_transposeA(PEmat1, nd_cont, nd_cont_recon, y, cont_recon.flux.data());

  // SxC^-1xS^T
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
  int np;
  fout.open(fcont+"_recon");

  /* determine the best precision */
  np = 1 + ceil(-log10((cont_recon.time[1] - cont_recon.time[0])));
  if(np<6)np=6;

  for(i=0; i<nd_cont_recon; i++)
  {
    fout<<fixed<<setprecision(np)
        <<cont_recon.time[i]<<"   "
        <<scientific<<cont_recon.flux[i]*cont.norm<<"  "<<cont_recon.error[i]*cont.norm<<endl;
  }
  fout.flush();
  fout.close();

  if(!fline.empty())
  {
    int il=0;
    list<Data>::iterator it;
    list<DataLC>::iterator itr;
    list<string>::iterator ifl;
    itr = lines_recon.begin();
    ifl = fline.begin();
    for(it=lines.begin(); it!=lines.end(); ++it)
    {
      Data& line = *(it);
      DataLC& line_recon = *(itr);
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
      // Cq^-1 = L^TxC^-1xL, Lbuf = C^-1xL
      multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_line, nq, sigma2, Lbuf);
      multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_line);
  
      // L^TxC^-1xy
      multiply_matvec_MN_transposeA(Lbuf, nd_line, nq, line.flux.data(), yq);
  
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
    
      set_covar_Umat_line(sigma, tau, USmat, il);
      /* PEmat1 = C^-1 x S^T */
      multiply_mat_transposeB_semiseparable_drw(USmat, W, D, phi, nd_line, nd_line_recon, sigma2, PEmat1);

      // (hat s) = SxC^-1xy
      multiply_matvec_MN_transposeA(PEmat1, nd_line, nd_line_recon, y, line_recon.flux.data());
  
      // SxC^-1xS^T
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
      fout.open(*(ifl)+"_recon");

      /* determine the best precision */
      np = 1 + ceil(-log10((line_recon.time[1] - line_recon.time[0])));
      if(np<6)np=6;

      for(i=0; i<nd_line_recon; i++)
      {
        fout<<fixed<<setprecision(np)
            <<line_recon.time[i]<<"   "
            <<scientific<<line_recon.flux[i]*line.norm<<"  "<<line_recon.error[i]*line.norm<<endl;
      }
      fout.flush();
      fout.close();

      itr++;
      ifl++;
      il++;
    }
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
  double t1, t2, sigma2;
  int i, j;
  
  sigma2 = sigma*sigma;
  for(i=0; i<cont_recon.time.size(); i++)
  {
    t1 = cont_recon.time[i];
    for(j=0; j<cont.time.size(); j++)
    {
      t2 = cont.time[j];
      USmat[i*cont.time.size()+j] = sigma2 * exp (- fabs(t1-t2) / tau );
    }
  }
  return;
}
void Cali::set_covar_Umat_line(double sigma, double tau, double *USmat, unsigned int il)
{
  double t1, t2, sigma2;
  int i, j;
  
  list<Data>::iterator it = lines.begin();
  list<DataLC>::iterator itr = lines_recon.begin();
  advance(it, il);
  advance(itr, il);

  Data& line = *(it);
  DataLC& line_recon = *(itr);
  
  sigma2 = sigma*sigma;
  for(i=0; i<line_recon.time.size(); i++)
  {
    t1 = line_recon.time[i];
    for(j=0; j<line.time.size(); j++)
    {
      t2 = line.time[j];
      USmat[i*line.time.size()+j] = sigma2 * exp (- fabs(t1-t2) / tau );
    }
  }
  return;
}
double Cali::get_norm_cont()
{
  return cont.norm;
}
double Cali::get_norm_line(unsigned int il)
{
  list<Data>::iterator it = lines.begin();
  advance(it, il);
  Data& line = *(it);
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

  /* calculate L^T*C^-1*L, note Lbuf = C^-1*L */
  multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_cont, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_cont);

  /* calculate L^T*C^-1*y, note Lbuf = C^-1 x L */
  multiply_matvec_MN_transposeA(Lbuf, nd_cont, nq, cont.flux.data(), yq);
  
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
    unsigned int il;
    list<Data>::iterator it;
    il = 0;
    for(it=cali->lines.begin(); it!=cali->lines.end(); ++it)
    {
      Data& line = *(it);
      int nd_line = line.time.size();

      Lbuf = workspace;
      ybuf = Lbuf + nd_line*nq;
      W = ybuf + nd_line;
      D = W + nd_line;
      phi = D + nd_line;
      Cq = phi + nd_line;
      yq = Cq + nq*nq;
  
      tau = exp(pm[3+il*2]);
      sigma = exp(pm[2+il*2])*sqrt(tau);
      sigma2 = sigma*sigma;
  
      compute_semiseparable_drw(line.time.data(), nd_line, sigma2, 1.0/tau, line.error.data(), 0.0, W, D, phi);
      lndet = 0.0;
      for(i=0; i<nd_line; i++)
        lndet += log(D[i]);
  
      /* calculate L^T*C^-1*L */
      multiply_mat_semiseparable_drw(Larr_data, W, D, phi, nd_line, nq, sigma2, Lbuf);
      multiply_mat_MN_transposeA(Larr_data, Lbuf, Cq, nq, nq, nd_line);
  
      /* calculate L^T*C^-1*y, note Lbuf = C^-1 x L */
      multiply_matvec_MN_transposeA(Lbuf, nd_line, nq, line.flux.data(), yq);
    
      lambda = Cq[0];
      ave_con = yq[0]/Cq[0];
  
      for(i=0;i<nd_line;i++)
      {
        ybuf[i] = line.flux[i] - ave_con;
      }
      multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_line, sigma2, Lbuf);
      prob2 += -0.5 * cblas_ddot(nd_line, ybuf, 1, Lbuf, 1);
  
      lndet_n =  0.0;
      for(i=0; i<line.num_code.size(); i++)
      {
        lndet_n += 2.0*log(ps_scale[i]) * line.num_code[i];
      }
  
      prob2 += - 0.5*lndet - 0.5*log(lambda) + 0.5 * lndet_n;

      il++;
    }    
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
  double logH = 0.0, width, move;
  int which;
  
  Cali *cali = (Cali *)arg;

  /* sample variability parameters more frequently */
  do
  {
    which = dnest_rand_int(cali->num_params);
  }while(cali->par_fix[which] == FIXED);

  width = ( cali->par_range_model[which][1] - cali->par_range_model[which][0] );
  
  move = pm[which];
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
  move = pm[which] - move;

  /* scale (phi) and shift (G) are degenerated
   * phi - G approx constant  
   */
  if(which >= cali->num_params_var && which < cali->num_params_var + cali->ncode && cali->par_fix[which+cali->ncode] == NOFIXED)
  {
    pm[which+cali->ncode] += move + dnest_randh() * (dnest_randn()*0.1);  /* random scatter with a std of 0.1*/
    dnest_wrap(&(pm[which+cali->ncode]), cali->par_range_model[which+cali->ncode][0], cali->par_range_model[which+cali->ncode][1]);
  }
  else if (which >= cali->num_params_var + cali->ncode && which < cali->num_params_var + 2*cali->ncode 
           && cali->par_fix[which-cali->ncode] == NOFIXED)
  {
    logH -= (-log(pm[which-cali->ncode]));
    pm[which-cali->ncode] += move + dnest_randh() * (dnest_randn()*0.1); /* random scatter with a std of 0.1*/
    dnest_wrap(&(pm[which-cali->ncode]), cali->par_range_model[which-cali->ncode][0], cali->par_range_model[which-cali->ncode][1]);
    logH += (-log(pm[which-cali->ncode]));
  }
  return logH;
}