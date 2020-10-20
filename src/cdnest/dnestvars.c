/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "dnestvars.h"
 
/* output files */
FILE *fsample, *fsample_info;

/* random number generator */
const gsl_rng_type * dnest_gsl_T;
gsl_rng * dnest_gsl_r;

Options options;
char options_file[STR_MAX_LENGTH];

// sampler
bool save_to_disk;
double compression;
unsigned int regularisation;

void *particles;
int dnest_size_of_modeltype;
int particle_offset_size, particle_offset_double;

LikelihoodType *log_likelihoods;
unsigned int *level_assignments;

// number account of unaccepted times
unsigned int *account_unaccepts;

int size_levels;  
Level *levels;
unsigned int count_saves, num_saves, num_saves_restart;
int dnest_which_particle_update; // which particle to be updated
int dnest_which_level_update;    // which level to be updated;
unsigned long long int count_mcmc_steps;
LikelihoodType *above;
unsigned int size_above;

double post_logz;
int dnest_num_params;
char dnest_sample_postfix[STR_MAX_LENGTH], dnest_sample_tag[STR_MAX_LENGTH], dnest_sample_dir[STR_MAX_LENGTH];

double *limits, *copies_of_limits;

int *dnest_perturb_accept;
int dnest_root;

int dnest_flag_restart=0, dnest_flag_postprc=0, dnest_flag_sample_info=0, dnest_flag_limits=0;
double dnest_post_temp=1.0;
char file_restart[STR_MAX_LENGTH], file_save_restart[STR_MAX_LENGTH];

void *dnest_arg;