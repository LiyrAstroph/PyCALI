/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */
#ifndef _DNESTVARS_H
#define _DNESTVARS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#define DNEST_MAJOR_VERSION 0  // Dec 2, 2018
#define DNEST_MINOR_VERSION 1
#define DNEST_PATCH_VERSION 0

#define STR_MAX_LENGTH (100)
#define BUF_MAX_LENGTH (200)
#define LEVEL_NUM_MAX (1000)

/* output files */
extern FILE *fsample, *fsample_info;

/* random number generator */
extern const gsl_rng_type * dnest_gsl_T;
extern gsl_rng * dnest_gsl_r;

typedef struct 
{
  double value;
  double tiebreaker;
}LikelihoodType;

typedef struct
{
  LikelihoodType log_likelihood;
  double log_X;
  unsigned long long int visits, exceeds;
  unsigned long long int accepts, tries;
}Level;

// struct for options
typedef struct
{
  unsigned int num_particles;
  unsigned int new_level_interval;
  unsigned int save_interval;
  unsigned int thread_steps;
  unsigned int max_num_levels;
  double lambda, beta;
  unsigned int max_num_saves;
  double max_pdiff;

  char sample_file[STR_MAX_LENGTH];
  char sample_info_file[STR_MAX_LENGTH];
  char levels_file[STR_MAX_LENGTH];
  char sampler_state_file[STR_MAX_LENGTH];
  char posterior_sample_file[STR_MAX_LENGTH];
  char posterior_sample_info_file[STR_MAX_LENGTH];
  char limits_file[STR_MAX_LENGTH];
}Options;
extern Options options;
extern char options_file[STR_MAX_LENGTH];

typedef struct
{
  void (*from_prior)(void *model);
  double (*log_likelihoods_cal)(const void *model);
  double (*log_likelihoods_cal_initial)(const void *model);
  double (*log_likelihoods_cal_restart)(const void *model);
  double (*perturb)(void *model);
  void (*print_particle)(FILE *fp, const void *model);
  void (*read_particle)(FILE *fp, void *model);
  void (*restart_action)(int iflag);
  void (*accept_action)();
  void (*kill_action)(int i, int i_copy);
}DNestFptrSet;

extern void *particles;
extern int dnest_size_of_modeltype;
extern int particle_offset_size, particle_offset_double;

// sampler
extern bool save_to_disk;
extern unsigned int num_threads;
extern double compression;
extern unsigned int regularisation;

extern void *particles;
extern LikelihoodType *log_likelihoods;
extern unsigned int *level_assignments;

// number account of unaccepted times
extern unsigned int *account_unaccepts;

extern int size_levels;  
extern Level *levels;
extern unsigned int count_saves, num_saves, num_saves_restart;
extern unsigned long long int count_mcmc_steps;
extern LikelihoodType *above;
extern unsigned int size_above;

extern int dnest_flag_restart, dnest_flag_postprc, dnest_flag_sample_info, dnest_flag_limits;
extern double dnest_post_temp;
extern char file_restart[STR_MAX_LENGTH], file_save_restart[STR_MAX_LENGTH];

extern double post_logz;
extern int dnest_num_params;
extern char dnest_sample_postfix[STR_MAX_LENGTH], dnest_sample_tag[STR_MAX_LENGTH], dnest_sample_dir[STR_MAX_LENGTH];

//the limits of parameters for each level;
extern double *limits, *copies_of_limits;

extern int dnest_which_particle_update; // which particle to be updated
extern int dnest_which_level_update;    // which level to be updated;
extern int *dnest_perturb_accept;
extern int dnest_root;
//***********************************************
/*                  functions                  */
double mod(double y, double x);
void dnest_wrap(double *x, double min, double max);
void wrap_limit(double *x, double min, double max);
int mod_int(int y, int x);
int dnest_cmp(const void *pa, const void *pb);

void options_load(int max_num_saves, double pdiff);
void setup(int argc, char** argv, DNestFptrSet *fptrset, int num_params, char *sample_dir, int max_num_saves, double pdiff);
void finalise();

double dnest(int argc, char **argv, DNestFptrSet *fptrset,  int num_params, char *sample_dir, int max_num_saves, double pdff);
void dnest_run();
void dnest_mcmc_run();
void update_particle(unsigned int which);
void update_level_assignment(unsigned int which);
double log_push(unsigned int which_level);
bool enough_levels(Level *l, int size_l);
void do_bookkeeping();
void save_levels();
void save_particle();
void save_limits();
void kill_lagging_particles();
void renormalise_visits();
void recalculate_log_X();
double dnest_randh();
double dnest_rand();
double dnest_randn();
int dnest_rand_int(int size);
void dnest_postprocess(double temperature, int max_num_saves, double pdiff);
void postprocess(double temperature);
void initialize_output_file();
void close_output_file();
void dnest_save_restart();
void dnest_restart();
void dnest_restart_action(int iflag);
void dnest_accept_action();
void dnest_kill_action(int i, int i_copy);
void dnest_print_particle(FILE *fp, const void *model);
void dnest_read_particle(FILE *fp, void *model);
int dnest_get_size_levels();
int dnest_get_which_level_update();
int dnest_get_which_particle_update();
void dnest_get_posterior_sample_file(char *fname);
int dnest_check_version(char *verion_str);
unsigned int dnest_get_which_num_saves();
unsigned int dnest_get_count_saves();
unsigned long long int dnest_get_count_mcmc_steps();
void dnest_check_fptrset(DNestFptrSet *fptrset);
DNestFptrSet * dnest_malloc_fptrset();
void dnest_free_fptrset(DNestFptrSet * fptrset);
/*=====================================================*/
// users responsible for following functions
void (*print_particle)(FILE *fp, const void *model);
void (*read_particle)(FILE *fp, void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*log_likelihoods_cal_restart)(const void *model);
double (*perturb)(void *model);
void (*restart_action)(int iflag);
void (*accept_action)();
void (*kill_action)(int i, int i_copy);
/*=====================================================*/

#ifdef __cplusplus
}
#endif

#endif
