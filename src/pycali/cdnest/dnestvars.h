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

#include "progress-bar.h"

#define DNEST_MAJOR_VERSION 0  // Dec 2, 2018
#define DNEST_MINOR_VERSION 1
#define DNEST_PATCH_VERSION 0

#define STR_MAX_LENGTH (100)
#define BUF_MAX_LENGTH (200)
#define LEVEL_NUM_MAX (2000)

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
  double max_ptol;

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
  void (*from_prior)(void *model, const void *arg);
  double (*log_likelihoods_cal)(const void *model, const void *arg);
  double (*log_likelihoods_cal_initial)(const void *model, const void *arg);
  double (*log_likelihoods_cal_restart)(const void *model, const void *arg);
  double (*perturb)(void *model, const void *arg);
  void (*print_particle)(FILE *fp, const void *model, const void *arg);
  void (*read_particle)(FILE *fp, void *model);
  void (*restart_action)(int iflag);
  void (*accept_action)();
  void (*kill_action)(int i, int i_copy);
}DNestFptrSet;

extern void *particles;
extern int dnest_size_of_modeltype;
extern int particle_offset_double;

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

extern void *dnest_arg;
//***********************************************
/*                  functions                  */
extern double mod(double y, double x);
extern void dnest_wrap(double *x, double min, double max);
extern void wrap_limit(double *x, double min, double max);
extern int mod_int(int y, int x);
extern int dnest_cmp(const void *pa, const void *pb);

extern void options_load(int max_num_saves, double ptol);
extern void setup(int argc, char** argv, DNestFptrSet *fptrset, int num_params, char *sample_dir, int max_num_saves, double ptol);
extern void finalise();

extern double dnest(int argc, char **argv, DNestFptrSet *fptrset,  int num_params, char *sample_dir, 
             int max_num_saves, double pdff, const void *arg);
extern void dnest_run();
extern void dnest_mcmc_run();
extern void update_particle(unsigned int which);
extern void update_level_assignment(unsigned int which);
extern double log_push(unsigned int which_level);
extern bool enough_levels(Level *l, int size_l);
extern void do_bookkeeping();
extern void save_levels();
extern void save_particle();
extern void save_limits();
extern void kill_lagging_particles();
extern void renormalise_visits();
extern void recalculate_log_X();
extern double dnest_randh();
extern double dnest_rand();
extern double dnest_randn();
extern int dnest_rand_int(int size);
extern void dnest_postprocess(double temperature, int max_num_saves, double ptol);
extern void postprocess(double temperature);
extern void initialize_output_file();
extern void close_output_file();
extern void dnest_save_restart();
extern void dnest_restart();
extern void dnest_restart_action(int iflag);
extern void dnest_accept_action();
extern void dnest_kill_action(int i, int i_copy);
extern void dnest_print_particle(FILE *fp, const void *model, const void *arg);
extern void dnest_read_particle(FILE *fp, void *model);
extern int dnest_get_size_levels();
extern int dnest_get_which_level_update();
extern int dnest_get_which_particle_update();
extern void dnest_get_posterior_sample_file(char *fname);
extern int dnest_check_version(char *verion_str);
extern unsigned int dnest_get_which_num_saves();
extern unsigned int dnest_get_count_saves();
extern unsigned long long int dnest_get_count_mcmc_steps();
extern void dnest_check_fptrset(DNestFptrSet *fptrset);
extern DNestFptrSet * dnest_malloc_fptrset();
extern void dnest_free_fptrset(DNestFptrSet * fptrset);
/*=====================================================*/
// users responsible for following functions
extern void (*print_particle)(FILE *fp, const void *model, const void *arg);
extern void (*read_particle)(FILE *fp, void *model);
extern void (*from_prior)(void *model, const void *arg);
extern double (*log_likelihoods_cal)(const void *model, const void *arg);
extern double (*log_likelihoods_cal_initial)(const void *model, const void *arg);
extern double (*log_likelihoods_cal_restart)(const void *model, const void *arg);
extern double (*perturb)(void *model, const void *arg);
extern void (*restart_action)(int iflag);
extern void (*accept_action)();
extern void (*kill_action)(int i, int i_copy);
/*=====================================================*/

extern ProgressBar *pb;

#ifdef __cplusplus
}
#endif

#endif
