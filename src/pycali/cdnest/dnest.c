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
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "dnestvars.h"

double dnest(int argc, char** argv, DNestFptrSet *fptrset, int num_params, 
             char *sample_dir, int max_num_saves, double ptol, const void *arg)
{
  int opt;
  extern int optind, opterr, optopt;
  extern char *optarg;
  extern int getopt(int argc, char *const *argv, const char *options);
  
  dnest_arg = (void *)arg;
  
  dnest_check_fptrset(fptrset);

  // cope with argv
  
  dnest_post_temp = 1.0;
  dnest_flag_restart = 0;
  dnest_flag_postprc = 0;
  dnest_flag_sample_info = 0;
  dnest_flag_limits = 0;

  strcpy(file_save_restart, "restart_dnest.txt");
  strcpy(dnest_sample_postfix, "\0");
  strcpy(dnest_sample_tag, "\0");
    
  opterr = 0;
  optind = 0;

  /* MAC getopt and GNU  getopt seem not compatible */
#if defined(__APPLE__) && defined(__MACH__)
  while( (opt = getopt(argc-1, argv+1, "r:s:pt:clx:g:")) != -1)
#else
  while( (opt = getopt(argc, argv, "r:s:pt:clx:g:")) != -1)
#endif
  {
    switch(opt)
    {
      case 'r':
        dnest_flag_restart = 1;
        strcpy(file_restart, optarg);
        printf("# CDnest restarts.\n");
        break;
      case 's':
        strcpy(file_save_restart, optarg);
        //printf("# Dnest sets restart file %s.\n", file_save_restart);
        break;
      case 'p':
        dnest_flag_postprc = 1;
        dnest_post_temp = 1.0;
        printf("# CDnest does postprocess.\n");
        break;
      case 't':
        dnest_post_temp = atof(optarg);
        printf("# CDnest sets a temperature %f.\n", dnest_post_temp);
        if(dnest_post_temp == 0.0)
        {
          printf("# CDnest incorrect option -t %s.\n", optarg);
          exit(0);
        }
        if(dnest_post_temp < 1.0)
        {
          printf("# CDnest temperature should >= 1.0\n");
          exit(0);
        }
        break;
      case 'c':
        dnest_flag_sample_info = 1;
        printf("# CDnest recalculates sample information.\n");
        break;
      case 'l':
        dnest_flag_limits = 1;
        printf("# CDnest level-dependent sampling.\n");
        break;
      case 'x':
        strcpy(dnest_sample_postfix, optarg);
        printf("# CDnest sets sample postfix %s.\n", dnest_sample_postfix);
        break;
      case 'g':
        strcpy(dnest_sample_tag, optarg);
        printf("# CDnest sets sample tag %s.\n", dnest_sample_tag);
        break;
      case '?':
        printf("# CDnest incorrect option -%c %s.\n", optopt, optarg);
        exit(0);
        break;
      default:
        break;
    }
  }
  
  setup(argc, argv, fptrset, num_params, sample_dir, max_num_saves, ptol);

  if(dnest_flag_postprc == 1)
  {
    dnest_postprocess(dnest_post_temp, max_num_saves, ptol);
    finalise();
    return post_logz;
  }

  if(dnest_flag_sample_info == 1)
  {
    dnest_postprocess(dnest_post_temp, max_num_saves, ptol);
    finalise();
    return post_logz;
  }

  if(dnest_flag_restart==1)
    dnest_restart();

  initialize_output_file();
  dnest_run();
  close_output_file();

  dnest_postprocess(dnest_post_temp, max_num_saves, ptol);

  finalise();
  
  return post_logz;
}

// postprocess, calculate evidence, generate posterior sample.
void dnest_postprocess(double temperature, int max_num_saves, double ptol)
{
  options_load(max_num_saves, ptol);
  postprocess(temperature);
}

void dnest_run()
{
  int i, j, k, size_all_above_incr;
  Level *pl, *levels_orig;
  int *buf_size_above, *buf_displs;
  double *plimits;
  extern int fileno(FILE *stream);
  
  printf("# Start diffusive nested sampling.\n");

  while(true)
  {
    //check for termination
    if(options.max_num_saves !=0 &&
        count_saves != 0 && (count_saves%options.max_num_saves == 0))
      break;

    dnest_mcmc_run();

    count_mcmc_steps += options.thread_steps;
    
    if(dnest_flag_limits == 1)
    {
      // limits of smaller levels should be larger than those of higher levels
      for(j=size_levels-2; j >= 0; j--)
        for(k=0; k<particle_offset_double; k++)
        {
          limits[ j * particle_offset_double *2 + k*2 ] = fmin( limits[ j * particle_offset_double *2 + k*2 ],
                  limits[ (j+1) * particle_offset_double *2 + k*2 ] );
          limits[ j * particle_offset_double *2 + k*2 + 1] = fmax( limits[ j * particle_offset_double *2 + k*2 +1 ],
                  limits[ (j+1) * particle_offset_double *2 + k*2 + 1 ] );
        }
    }

    do_bookkeeping();

    if(count_mcmc_steps >= (count_saves + 1)*options.save_interval)
    {
      save_particle();

      // save levels, limits, sync samples when running a number of steps
      if( count_saves % num_saves == 0 )
      {
        save_levels();
        if(dnest_flag_limits == 1)
          save_limits();
        fflush(fsample_info);
        fsync(fileno(fsample_info));
        fflush(fsample);
        fsync(fileno(fsample));
        printf("# Save levels, limits, and sync samples at N= %d.\n", count_saves);
      }

      //if( count_saves % num_saves_restart == 0 )
      //{
      //  dnest_save_restart();
      //}
    }
  }
  
  //dnest_save_restart();

  //save levels
  save_levels();
  if(dnest_flag_limits == 1)
    save_limits();

  /* output state of sampler */
  FILE *fp;
  fp = fopen(options.sampler_state_file, "w");
  fprintf(fp, "%d %d\n", size_levels, count_saves);
  fclose(fp);
}

void do_bookkeeping()
{
  int i;
  //bool created_level = false;

  if(!enough_levels(levels, size_levels) && size_above >= options.new_level_interval)
  {
    // in descending order 
    qsort(above, size_above, sizeof(LikelihoodType), dnest_cmp);
    int index = (int)( (1.0/compression) * size_above);

    Level level_tmp = {above[index], 0.0, 0, 0, 0, 0};
    levels[size_levels] = level_tmp;
    size_levels++;
    
    printf("# Creating level %d with log likelihood = %e.\n", 
               size_levels-1, levels[size_levels-1].log_likelihood.value);

    // clear out the last index records
    for(i=index; i<size_above; i++)
    {
      above[i].value = 0.0;
      above[i].tiebreaker = 0.0;
    }
    size_above = index;

    if(enough_levels(levels, size_levels))
    {
      renormalise_visits();
      options.max_num_levels = size_levels;
      printf("# Done creating levles.\n");
    }
    else
    {
      kill_lagging_particles();
    }
  }
  recalculate_log_X();
}

void recalculate_log_X()
{
  int i;

  levels[0].log_X = 0.0;
  for(i=1; i<size_levels; i++)
  {
    levels[i].log_X = levels[i-1].log_X 
    + log( (double)( (levels[i-1].exceeds + 1.0/compression * regularisation)
                    /(levels[i-1].visits + regularisation)  ) );
  }
}

void renormalise_visits()
{
  size_t i;

  for(i=0; i<size_levels; i++)
  {
    if(levels[i].tries >= regularisation)
    {
      levels[i].accepts = ((double)(levels[i].accepts+1) / (double)(levels[i].tries+1)) * regularisation;
      levels[i].tries = regularisation;
    }

    if(levels[i].visits >= regularisation)
    {
      levels[i].exceeds = ( (double) (levels[i].exceeds+1) / (double)(levels[i].visits + 1) ) * regularisation;
      levels[i].visits = regularisation;
    }
  }
}

void kill_lagging_particles()
{
  static unsigned int deletions = 0;

  bool *good;
  good = (bool *)malloc(options.num_particles * sizeof(bool));

  double max_log_push = -DBL_MAX;

  double kill_probability = 0.0;
  unsigned int num_bad = 0;
  size_t i;

  for(i=0; i<options.num_particles; i++)good[i] = true;

  for(i=0; i<options.num_particles; i++)
  {
    if( log_push(level_assignments[i]) > max_log_push)
      max_log_push = log_push(level_assignments[i]);

    kill_probability = pow(1.0 - 1.0/(1.0 + exp(-log_push(level_assignments[i]) - 4.0)), 3);
    if(gsl_rng_uniform(dnest_gsl_r) <= kill_probability)
    {
      good[i] = false;
      ++num_bad;
    }
  }

  if(num_bad < options.num_particles)
  {
    for(i=0; i< options.num_particles; i++)
    {
      if(!good[i])
      {
        int i_copy;
        do
        {
          i_copy = gsl_rng_uniform_int(dnest_gsl_r, options.num_particles);
        }while(!good[i_copy] || gsl_rng_uniform(dnest_gsl_r) >= exp(log_push(level_assignments[i_copy]) - max_log_push));

        /* sizeof(char *) == sizeof(void *) */
        memcpy((void *)((double *)particles+i*particle_offset_double), 
               (void *)((double *)particles + i_copy*particle_offset_double), 
               dnest_size_of_modeltype);
        log_likelihoods[i] = log_likelihoods[i_copy];
        level_assignments[i] = level_assignments[i_copy];
         
        kill_action(i, i_copy);

        deletions++;

        printf("# Replacing lagging particle.\n");
        printf("# This has happened %d times.\n", deletions);
      }
    }
  }
  else
    printf("# Warning: all particles lagging!.\n");

  free(good);
}

/* save levels */
void save_levels()
{
  if(!save_to_disk)
    return;
  
  int i;
  FILE *fp;

  fp = fopen(options.levels_file, "w");
  fprintf(fp, "# log_X, log_likelihood, tiebreaker, accepts, tries, exceeds, visits\n");
  for(i=0; i<size_levels; i++)
  {
    fprintf(fp, "%14.12g %14.12g %f %llu %llu %llu %llu\n", levels[i].log_X, levels[i].log_likelihood.value, 
      levels[i].log_likelihood.tiebreaker, levels[i].accepts,
      levels[i].tries, levels[i].exceeds, levels[i].visits);
  }
  fclose(fp);

  /* update state of sampler */
  fp = fopen(options.sampler_state_file, "w");
  fprintf(fp, "%d %d\n", size_levels, count_saves);
  fclose(fp);
}

void save_limits()
{
  int i, j;
  FILE *fp;

  fp = fopen(options.limits_file, "w");
  for(i=0; i<size_levels; i++)
  {
    fprintf(fp, "%d  ", i);
    for(j=0; j<particle_offset_double; j++)
      fprintf(fp, "%f  %f  ", limits[i*2*particle_offset_double+j*2], limits[i*2*particle_offset_double+j*2+1]);

    fprintf(fp, "\n");
  }
  fclose(fp);
}

/* save particle */
void save_particle()
{
  count_saves++;

  if(!save_to_disk)
    return;
  
  int whichparticle, whichtask;
  void *particle_message;
  
  if(count_saves%10 == 0)
    printf("#[%.1f%%] Saving sample N= %d.\n", 100.0*count_saves/options.max_num_saves, count_saves);
    
  whichparticle =  gsl_rng_uniform_int(dnest_gsl_r,options.num_particles);

  print_particle(fsample, (void *)((double *)particles + whichparticle * particle_offset_double), dnest_arg);

  fprintf(fsample_info, "%d %e %f %d\n", level_assignments[whichparticle], 
        log_likelihoods[whichparticle].value,
        log_likelihoods[whichparticle].tiebreaker,
        whichparticle);
}

void dnest_mcmc_run()
{
  unsigned int which;
  unsigned int i;
  
  for(i = 0; i<options.thread_steps; i++)
  {

    /* randomly select out one particle to update */
    which = gsl_rng_uniform_int(dnest_gsl_r, options.num_particles);

    dnest_which_particle_update = which;

    //if(count_mcmc_steps >= 10000)printf("FFFF\n");
    //printf("%d\n", which);
    //printf("%f %f %f\n", particles[which].param[0], particles[which].param[1], particles[which].param[2]);
    //printf("level:%d\n", level_assignments[which]);
    //printf("%e\n", log_likelihoods[which].value);

    if(gsl_rng_uniform(dnest_gsl_r) <= 0.5)
    {
      update_particle(which);
      update_level_assignment(which);
    }
    else
    {
      update_level_assignment(which);
      update_particle(which);
    }
        
    if( !enough_levels(levels, size_levels)  && levels[size_levels-1].log_likelihood.value <= log_likelihoods[which].value)
    {
      above[size_above] = log_likelihoods[which];
      size_above++;
    }
  }
}


void update_particle(unsigned int which)
{
  void *particle = (void *)((double *)particles+ which*particle_offset_double);
  LikelihoodType *logl = &(log_likelihoods[which]);
  
  Level *level = &(levels[level_assignments[which]]);

  void *proposal = (void *)malloc(dnest_size_of_modeltype);
  LikelihoodType logl_proposal;
  double log_H;

  memcpy(proposal, particle, dnest_size_of_modeltype);
  dnest_which_level_update = level_assignments[which];
  
  log_H = perturb(proposal, dnest_arg);
  
  logl_proposal.value = log_likelihoods_cal(proposal, dnest_arg);
  logl_proposal.tiebreaker =  (*logl).tiebreaker + gsl_rng_uniform(dnest_gsl_r);
  dnest_wrap(&logl_proposal.tiebreaker, 0.0, 1.0);
  
  if(log_H > 0.0)
    log_H = 0.0;

  dnest_perturb_accept[which] = 0;
  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_H) && level->log_likelihood.value <= logl_proposal.value)
  {
    memcpy(particle, proposal, dnest_size_of_modeltype);
    memcpy(logl, &logl_proposal, sizeof(LikelihoodType));
    level->accepts++;

    dnest_perturb_accept[which] = 1;
    accept_action();
    account_unaccepts[which] = 0; /* reset the number of unaccepted perturb */
  }
  else 
  {
    account_unaccepts[which] += 1; /* number of unaccepted perturb */
  }
  level->tries++;
  
  unsigned int current_level = level_assignments[which];
  for(; current_level < size_levels-1; ++current_level)
  {
    levels[current_level].visits++;
    if(levels[current_level+1].log_likelihood.value <= log_likelihoods[which].value)
      levels[current_level].exceeds++;
    else
      break; // exit the loop if it does not satify higher levels
  }
  free(proposal);
}

void update_level_assignment(unsigned int which)
{
  int i;

  int proposal = level_assignments[which] 
                 + (int)( pow(10.0, 2*gsl_rng_uniform(dnest_gsl_r))*gsl_ran_ugaussian(dnest_gsl_r));

  if(proposal == level_assignments[which])
    proposal =  ((gsl_rng_uniform(dnest_gsl_r) < 0.5)?(proposal-1):(proposal+1));

  proposal=mod_int(proposal, size_levels);

  double log_A = -levels[proposal].log_X + levels[level_assignments[which]].log_X;

  log_A += log_push(proposal) - log_push(level_assignments[which]);
  
  // enforce uniform exploration if levels are enough
  if(enough_levels(levels, size_levels))
    log_A += options.beta*log( (double)(levels[level_assignments[which]].tries +1)/ (levels[proposal].tries +1) );

  if(log_A > 0.0)
    log_A = 0.0;

  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_A) && levels[proposal].log_likelihood.value <= log_likelihoods[which].value)
  {
    level_assignments[which] = proposal;

// update the limits of the level
    if(dnest_flag_limits == 1)
    {
      double *particle = (double *) ((double *)particles+ which*particle_offset_double);
      for(i=0; i<particle_offset_double; i++)
      {
        limits[proposal * 2 * particle_offset_double +  i*2] = 
            fmin(limits[proposal * 2* particle_offset_double +  i*2], particle[i]);
        limits[proposal * 2 * particle_offset_double +  i*2+1] = 
            fmax(limits[proposal * 2 * particle_offset_double +  i*2+1], particle[i]);
      }
    }

  }

}

double log_push(unsigned int which_level)
{
  if(which_level > size_levels)
  {
    printf("level overflow %d %d.\n", which_level, size_levels);
    exit(0);
  }
  if(enough_levels(levels, size_levels))
    return 0.0;

  int i = which_level - (size_levels - 1);
  return ((double)i)/options.lambda;
}

bool enough_levels(Level *l, int size_l)
{
  int i;

  if(options.max_num_levels == 0)
  {
    if(size_l >= LEVEL_NUM_MAX)
    {
      printf("Warning:size of levels approches the limit %d!\nbetter to increase compression!\n", LEVEL_NUM_MAX);
      printf("The default value is exp(1.0)=2.72; the present value is %f.", compression);
      exit(0);
      return true;
    }

    if(size_l < 10)
      return false;

    int num_levels_to_check = 20;
    if(size_l > 80)
      num_levels_to_check = (int)(sqrt(20) * sqrt(0.25*size_l));

    int k = size_l - 1, kc = 0;
    double tot = 0.0;
    double max = -DBL_MAX;
    double diff;

    for(i= 0; i<num_levels_to_check; i++)
    {
      diff = l[k].log_likelihood.value - l[k-1].log_likelihood.value;
      tot += diff;
      if(diff > max)
        max = diff;

      k--;
      kc++;
      if( k < 1 )
        break;
    }
    if(tot/kc < options.max_ptol && max < options.max_ptol*1.1)
      return true;
    else
      return false;
  }
  return (size_l >= options.max_num_levels);
}

void initialize_output_file()
{
  if(dnest_flag_restart !=1)
    fsample = fopen(options.sample_file, "w");
  else
    fsample = fopen(options.sample_file, "a");
  
  if(fsample==NULL)
  {
    fprintf(stderr, "# Cannot open file sample.txt.\n");
    exit(0);
  }
  if(dnest_flag_restart != 1)
    fprintf(fsample, "# \n");

  if(dnest_flag_restart != 1)
    fsample_info = fopen(options.sample_info_file, "w");
  else
    fsample_info = fopen(options.sample_info_file, "a");

  if(fsample_info==NULL)
  {
    fprintf(stderr, "# Cannot open file %s.\n", options.sample_info_file);
    exit(0);
  }
  if(dnest_flag_restart != 1)
    fprintf(fsample_info, "# level assignment, log likelihood, tiebreaker, ID.\n");
}

void close_output_file()
{
  fclose(fsample);
  fclose(fsample_info);
}

void setup(int argc, char** argv, DNestFptrSet *fptrset, int num_params, char *sample_dir, int max_num_saves, double ptol)
{
  int i, j;

  // root task.
  dnest_root = 0;

  // setup function pointers
  from_prior = fptrset->from_prior;
  log_likelihoods_cal = fptrset->log_likelihoods_cal;
  log_likelihoods_cal_initial = fptrset->log_likelihoods_cal_initial;
  log_likelihoods_cal_restart = fptrset->log_likelihoods_cal_restart;
  perturb = fptrset->perturb;
  print_particle = fptrset->print_particle;
  read_particle = fptrset->read_particle;
  restart_action = fptrset->restart_action;
  accept_action = fptrset->accept_action;
  kill_action = fptrset->kill_action;
  strcpy(dnest_sample_dir, sample_dir);

  // random number generator
  dnest_gsl_T = (gsl_rng_type *) gsl_rng_default;
  dnest_gsl_r = gsl_rng_alloc (dnest_gsl_T);
#ifndef Debug
  gsl_rng_set(dnest_gsl_r, time(NULL));
#else
  gsl_rng_set(dnest_gsl_r, 9999);
  printf("# debugging, dnest random seed %d\n", 9999);
#endif  
  
  dnest_num_params = num_params;
  dnest_size_of_modeltype = dnest_num_params * sizeof(double);

  // read options
  options_load(max_num_saves, ptol);

  //dnest_post_temp = 1.0;
  compression = exp(1.0);
  regularisation = options.new_level_interval*sqrt(options.lambda);
  save_to_disk = true;

  // particles
  //particle_offset_size = dnest_size_of_modeltype/sizeof(void);
  particle_offset_double = dnest_size_of_modeltype/sizeof(double);
  particles = (void *)malloc(options.num_particles*dnest_size_of_modeltype);
  
  // initialise sampler
  above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));

  log_likelihoods = (LikelihoodType *)malloc(2*options.num_particles * sizeof(LikelihoodType));
  level_assignments = (unsigned int*)malloc(options.num_particles * sizeof(unsigned int));

  account_unaccepts = (unsigned int *)malloc(options.num_particles * sizeof(unsigned int));
  for(i=0; i<options.num_particles; i++)
  {
    account_unaccepts[i] = 0;
  }

  if(options.max_num_levels != 0)
  {
    levels = (Level *)malloc(options.max_num_levels * sizeof(Level));
    if(dnest_flag_limits == 1)
    {
      limits = (double *)malloc(options.max_num_levels * particle_offset_double * 2 * sizeof(double));
      if(limits == NULL)
      {
        printf("Cannot allocate memory for limits.\n"
               "This usually happens when both the numbers of parameters and levels are extremely large.\n"
               "Please do not switch on '-l' option in argv passed to cdnest.\n");
        exit(EXIT_FAILURE);
      }
      for(i=0; i<options.max_num_levels; i++)
      {
        for(j=0; j<particle_offset_double; j++)
        {
          limits[i*2*particle_offset_double+ j*2] = DBL_MAX;
          limits[i*2*particle_offset_double + j*2 + 1] = -DBL_MAX;
        }
      }
    }
  }
  else
  {
    levels = (Level *)malloc(LEVEL_NUM_MAX * sizeof(Level));

    if(dnest_flag_limits == 1)
    {
      limits = (double *)malloc(LEVEL_NUM_MAX * particle_offset_double * 2 * sizeof(double));
      for(i=0; i<LEVEL_NUM_MAX; i++)
      {
        for(j=0; j<particle_offset_double; j++)
        {
          limits[i*2*particle_offset_double + j*2] = DBL_MAX;
          limits[i*2*particle_offset_double + j*2 + 1] = -DBL_MAX;
        }
      }
    }
  }
  
  dnest_perturb_accept = (int *)malloc(options.num_particles * sizeof(int));
  for(i=0; i<options.num_particles; i++)
  {
    dnest_perturb_accept[i] = 0;
  }

  count_mcmc_steps = 0;
  count_saves = 0;
  num_saves = (int)fmax(0.02*options.max_num_saves, 1.0);
  num_saves_restart = (int)fmax(0.2 * options.max_num_saves, 1.0);

// first level
  size_levels = 0;
  size_above = 0;
  LikelihoodType like_tmp = {-DBL_MAX, gsl_rng_uniform(dnest_gsl_r)};
  Level level_tmp = {like_tmp, 0.0, 0, 0, 0, 0};
  levels[size_levels] = level_tmp;
  size_levels++;
  
  for(i=0; i<options.num_particles; i++)
  {
    dnest_which_particle_update = i;
    dnest_which_level_update = 0;
    from_prior((void *)((double *)particles+i*particle_offset_double), dnest_arg);
    log_likelihoods[i].value = log_likelihoods_cal_initial((void *)((double *)particles+i*particle_offset_double), dnest_arg);
    log_likelihoods[i].tiebreaker = dnest_rand();
    level_assignments[i] = 0;
  }

  /*ModelType proposal;
  printf("%f %f %f \n", particles[0].param[0], particles[0].param[1], particles[0].param[2] );
  proposal = particles[0];
  printf("%f %f %f \n", proposal.param[0], proposal.param[1], proposal.param[2] );*/
}

void finalise()
{
  free(particles);
  free(above);
  free(log_likelihoods);
  free(level_assignments);
  free(levels);

  free(account_unaccepts);

  if(dnest_flag_limits == 1)
    free(limits);

  gsl_rng_free(dnest_gsl_r);

  free(dnest_perturb_accept);

  printf("# Finalizing dnest.\n");
}


void options_load(int max_num_saves, double ptol)
{
  //sscanf(buf, "%d", &options.num_particles);
  options.num_particles = 2;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%d", &options.new_level_interval);
  options.new_level_interval = options.num_particles * dnest_num_params*10;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%d", &options.save_interval);
  options.save_interval = options.new_level_interval;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%d", &options.thread_steps);
  options.thread_steps = options.new_level_interval;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%d", &options.max_num_levels);
  options.max_num_levels = 0;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%lf", &options.lambda);
  options.lambda = 10.0;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%lf", &options.beta);
  options.beta = 100.0;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%d", &options.max_num_saves);
  options.max_num_saves = max_num_saves;

  options.max_ptol = ptol;

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.sample_file);
  strcpy(options.sample_file, dnest_sample_dir);
  strcat(options.sample_file,"/sample");
  strcat(options.sample_file, dnest_sample_tag);
  strcat(options.sample_file, ".txt");
  strcat(options.sample_file, dnest_sample_postfix);
  
  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.sample_info_file);
  strcpy(options.sample_info_file, dnest_sample_dir);
  strcat(options.sample_info_file,"/sample_info");
  strcat(options.sample_info_file, dnest_sample_tag);
  strcat(options.sample_info_file, ".txt");
  strcat(options.sample_info_file, dnest_sample_postfix);
  
  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.levels_file);
  strcpy(options.levels_file, dnest_sample_dir);
  strcat(options.levels_file,"/levels");
  strcat(options.levels_file, dnest_sample_tag);
  strcat(options.levels_file, ".txt");
  strcat(options.levels_file, dnest_sample_postfix);
  
  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.sampler_state_file);
  strcpy(options.sampler_state_file, dnest_sample_dir);
  strcat(options.sampler_state_file,"/sampler_state");
  strcat(options.sampler_state_file, dnest_sample_tag);
  strcat(options.sampler_state_file, ".txt");
  strcat(options.sampler_state_file, dnest_sample_postfix);
  
  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.posterior_sample_file);
  strcpy(options.posterior_sample_file, dnest_sample_dir);
  strcat(options.posterior_sample_file,"/posterior_sample");
  strcat(options.posterior_sample_file, dnest_sample_tag);
  strcat(options.posterior_sample_file, ".txt");
  strcat(options.posterior_sample_file, dnest_sample_postfix);

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.posterior_sample_info_file);
  strcpy(options.posterior_sample_info_file, dnest_sample_dir);
  strcat(options.posterior_sample_info_file,"/posterior_sample_info");
  strcat(options.posterior_sample_info_file, dnest_sample_tag);
  strcat(options.posterior_sample_info_file, ".txt");
  strcat(options.posterior_sample_info_file, dnest_sample_postfix);

  //fgets(buf, BUF_MAX_LENGTH, fp);
  //sscanf(buf, "%s", options.limits_file);
  strcpy(options.limits_file, dnest_sample_dir);
  strcat(options.limits_file,"/limits");
  strcat(options.limits_file, dnest_sample_tag);
  strcat(options.limits_file, ".txt");
  strcat(options.limits_file, dnest_sample_postfix);

  // check options.
  
  if(options.new_level_interval < options.thread_steps)
  {
    printf("# incorrect options:\n");
    printf("# new level interval should be equal to or larger than"); 
    printf("  totaltask * thread step.\n");
    exit(0);
  }

/*  strcpy(options.sample_file, "sample.txt");
  strcpy(options.sample_info_file, "sample_info.txt");
  strcpy(options.levels_file, "levels.txt");
  strcpy(options.sampler_state_file, "sampler_state.txt");*/
}


double mod(double y, double x)
{
  if(x > 0.0)
  {
    return (y/x - floor(y/x))*x;
  }
  else if(x == 0.0)
  {
    return 0.0;
  }
  else
  {
    printf("Warning in mod(double, double) %e\n", x);
    exit(0);
  }
  
}

void dnest_wrap(double *x, double min, double max)
{
  *x = mod(*x - min, max - min) + min;
}

void wrap_limit(double *x, double min, double max)
{

  *x = fmax(fmin(*x, max), min);
}

int mod_int(int y, int x)
{
  if(y >= 0)
    return y - (y/x)*x;
  else
    return (x-1) - mod_int(-y-1, x);
}

double dnest_randh()
{
  return pow(10.0, 1.5 - 3.0*fabs(gsl_ran_tdist(dnest_gsl_r, 2))) * gsl_ran_ugaussian(dnest_gsl_r);
}

double dnest_rand()
{
  return gsl_rng_uniform(dnest_gsl_r);
}

int dnest_rand_int(int size)
{
  return gsl_rng_uniform_int(dnest_gsl_r, size);
}

double dnest_randn()
{
  return gsl_ran_ugaussian(dnest_gsl_r);
}

int dnest_cmp(const void *pa, const void *pb)
{
  LikelihoodType *a = (LikelihoodType *)pa;
  LikelihoodType *b = (LikelihoodType *)pb;

  // in decesending order
  if(a->value > b->value)
    return false;
  if( a->value == b->value && a->tiebreaker > b->tiebreaker)
    return false;
  
  return true;
}


int dnest_get_size_levels()
{
  return size_levels;
}

int dnest_get_which_level_update()
{
  return dnest_which_level_update;
}

int dnest_get_which_particle_update()
{
  return dnest_which_particle_update;
}

unsigned int dnest_get_which_num_saves()
{
  return num_saves;
}
unsigned int dnest_get_count_saves()
{
  return count_saves;
}

unsigned long long int dnest_get_count_mcmc_steps()
{
  return count_mcmc_steps;
}

void dnest_get_posterior_sample_file(char *fname)
{
  strcpy(fname, options.posterior_sample_file);
  return;
}
/* 
 * version check
 * 
 *  1: greater
 *  0: equal
 * -1: lower
 */
int dnest_check_version(char *version_str)
{
  int major, minor, patch;

  sscanf(version_str, "%d.%d.%d", &major, &minor, &patch);
  
  if(major > DNEST_MAJOR_VERSION)
    return 1;
  if(major < DNEST_MAJOR_VERSION)
    return -1;

  if(minor > DNEST_MINOR_VERSION)
    return 1;
  if(minor < DNEST_MINOR_VERSION)
    return -1;

  if(patch > DNEST_PATCH_VERSION)
    return 1;
  if(patch > DNEST_PATCH_VERSION)
    return -1;

  return 0;
}

void dnest_check_fptrset(DNestFptrSet *fptrset)
{
  if(fptrset->from_prior == NULL)
  {
    printf("\"from_prior\" function is not defined.\n");
    exit(0);
  }

  if(fptrset->print_particle == NULL)
  {
    //printf("\"print_particle\" function is not defined. \
    //  \nSet to be default function in dnest.\n");
    fptrset->print_particle = dnest_print_particle;
  }

  if(fptrset->read_particle == NULL)
  {
    //printf("\"read_particle\" function is not defined. \
    //  \nSet to be default function in dnest.\n");
    fptrset->read_particle = dnest_read_particle;
  }

  if(fptrset->log_likelihoods_cal == NULL)
  {
    printf("\"log_likelihoods_cal\" function is not defined.\n");
    exit(0);
  }

  if(fptrset->log_likelihoods_cal_initial == NULL)
  {
    //printf("\"log_likelihoods_cal_initial\" function is not defined. \
    //  \nSet to the same as \"log_likelihoods_cal\" function.\n");
    fptrset->log_likelihoods_cal_initial = fptrset->log_likelihoods_cal;
  }

  if(fptrset->log_likelihoods_cal_restart == NULL)
  {
    //printf("\"log_likelihoods_cal_restart\" function is not defined. \
    //  \nSet to the same as \"log_likelihoods_cal\" function.\n");
    fptrset->log_likelihoods_cal_restart = fptrset->log_likelihoods_cal;
  }

  if(fptrset->perturb == NULL)
  {
    printf("\"perturb\" function is not defined.\n");
    exit(0);
  }

  if(fptrset->restart_action == NULL)
  {
    //printf("\"restart_action\" function is not defined.\
    //  \nSet to the default function in dnest.\n");
    fptrset->restart_action = dnest_restart_action;
  }

  if(fptrset->accept_action == NULL)
  {
    //printf("\"accept_action\" function is not defined.\
    //  \nSet to the default function in dnest.\n");
    fptrset->accept_action = dnest_accept_action;
  }

  if(fptrset->kill_action == NULL)
  {
    //printf("\"kill_action\" function is not defined.\
    //  \nSet to the default function in dnest.\n");
    fptrset->kill_action = dnest_kill_action;
  }

  return;
}

DNestFptrSet * dnest_malloc_fptrset()
{
  DNestFptrSet * fptrset;
  fptrset = (DNestFptrSet *)malloc(sizeof(DNestFptrSet));

  fptrset->from_prior = NULL;
  fptrset->log_likelihoods_cal = NULL;
  fptrset->log_likelihoods_cal_initial = NULL;
  fptrset->log_likelihoods_cal_restart = NULL;
  fptrset->perturb = NULL;
  fptrset->print_particle = NULL;
  fptrset->read_particle = NULL;
  fptrset->restart_action = NULL;
  fptrset->accept_action = NULL;
  fptrset->kill_action = NULL;
  return fptrset;
}

void dnest_free_fptrset(DNestFptrSet * fptrset)
{
  free(fptrset);
  return;
}


/*!
 *  Save sampler state for later restart. 
 */
void dnest_save_restart()
{
  FILE *fp;
  int i, j;
  void *particles_all;
  LikelihoodType *log_likelihoods_all;
  unsigned int *level_assignments_all;
  char str[200];

  
  sprintf(str, "%s_%d", file_save_restart, count_saves);
  fp = fopen(str, "wb");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s. \n", file_save_restart);
    exit(0);
  }

  
  printf("# Save restart data to file %s.\n", str);

    //fprintf(fp, "%d %d\n", count_saves, count_mcmc_steps);
    //fprintf(fp, "%d\n", size_levels_combine);

  fwrite(&count_saves, sizeof(int), 1, fp);
  fwrite(&count_mcmc_steps, sizeof(int), 1, fp);
  fwrite(&size_levels, sizeof(int), 1, fp);

  for(i=0; i<size_levels; i++)
  {
    //fprintf(fp, "%14.12g %14.12g %f %llu %llu %llu %llu\n", levels_combine[i].log_X, levels_combine[i].log_likelihood.value, 
    //  levels_combine[i].log_likelihood.tiebreaker, levels_combine[i].accepts,
    //  levels_combine[i].tries, levels[i].exceeds, levels_combine[i].visits);

    fwrite(&levels[i], sizeof(Level), 1, fp);
  }

  for(i=0; i<options.num_particles; i++)
  {
    //fprintf(fp, "%d %e %f\n", level_assignments_all[j*options.num_particles + i], 
    //  log_likelihoods_all[j*options.num_particles + i].value,
    //  log_likelihoods_all[j*options.num_particles + i].tiebreaker);  
    fwrite(&level_assignments[j*options.num_particles + i], sizeof(int), 1, fp);
    fwrite(&log_likelihoods[j*options.num_particles + i], sizeof(LikelihoodType), 1, fp);    
  }
    
    
  if(dnest_flag_limits == 1)
  {
    for(i=0; i<size_levels; i++)
    {
      //fprintf(fp, "%d  ", i);
      for(j=0; j<particle_offset_double; j++)
      {
        //fprintf(fp, "%f  %f  ", limits[i*2*particle_offset_double+j*2], limits[i*2*particle_offset_double+j*2+1]);
        fwrite(&limits[i*2*particle_offset_double+j*2], sizeof(double), 2, fp);
      }

        //fprintf(fp, "\n");
    }
  }
    

  for(i=0; i<options.num_particles; i++)
  {
    //print_particle(fp, particles_all + (j * options.num_particles + i) * particle_offset_size);
    fwrite((void*)((double *)particles + (j * options.num_particles + i) * particle_offset_double), dnest_size_of_modeltype, 1, fp);
  } 
    
  fclose(fp);

  restart_action(0);
}

void dnest_restart()
{
  FILE *fp;
  int i, j;
  void *particles_all;
  unsigned int *level_assignments_all;
  LikelihoodType *log_likelihoods_all;
  void *particle;

  fp = fopen(file_restart, "rb");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s. \n", file_restart);
    exit(0);
  }

  printf("# Reading %s\n", file_restart);

  fread(&count_saves, sizeof(int), 1, fp);
  fread(&count_mcmc_steps, sizeof(int), 1, fp);
  fread(&size_levels, sizeof(int), 1, fp);
    
  /* consider that the newly input max_num_levels may be different from the save one */
  if(options.max_num_levels != 0)
  {
    if(size_levels > options.max_num_levels)
    {
      printf("# input max_num_levels %d smaller than the one in restart data %d.\n", options.max_num_levels, size_levels);
      size_levels = options.max_num_levels;
    }   
  }
  else  /* not input max_num_levels, directly use the saved size of levels */
  {
    if(size_levels > LEVEL_NUM_MAX)
    {
      printf("# the saved size of levels %d exceeds LEVEL_NUM_MAX %d. \n", size_levels, LEVEL_NUM_MAX);
      exit(EXIT_FAILURE);
    }
  }
  // read levels
  for(i=0; i<size_levels; i++)
  {     
    if(i<size_levels) // not read all the levels
      fread(&levels[i], sizeof(Level), 1, fp);
    else
      fseek(fp, sizeof(Level), SEEK_CUR); /* offset the file point */
  }
    memcpy(levels, levels, size_levels * sizeof(Level));

    // read level assignment
  for(i=0; i<options.num_particles; i++)
  {
    fread(&level_assignments[j*options.num_particles + i], sizeof(int), 1, fp);
    fread(&log_likelihoods[j*options.num_particles + i], sizeof(LikelihoodType), 1, fp); 
    
    /* reset the level assignment that exceeds the present maximum level */
    if(level_assignments[j*options.num_particles + i] > size_levels -1)
    {
      level_assignments[j*options.num_particles + i] = size_levels - 1;
    }
  }

  // read limits
  if(dnest_flag_limits == 1)
  {
    for(i=0; i<size_levels; i++)
    {
      if(i < size_levels)
      {
        for(j=0; j<particle_offset_double; j++)
        {
          fread(&limits[i*2*particle_offset_double+j*2], sizeof(double), 2, fp);      
        }
      }
      else 
      {
        fseek(fp, sizeof(double) * 2 * particle_offset_double, SEEK_CUR); /* offset the file point */
      }
    }
  }

  // read particles
  for(i=0; i<options.num_particles; i++)
  {
    particle = (void *)((double *)particles + (j * options.num_particles + i) * particle_offset_double);
    fread(particle, dnest_size_of_modeltype, 1, fp);
  }

  fclose(fp);

  if(count_saves > options.max_num_saves)
  {
    printf("# Number of samples already larger than the input number, exit!\n");
    exit(0);
  }
  
  num_saves = (int)fmax(0.02*(options.max_num_saves-count_saves), 1.0); /* reset num_saves */
  num_saves_restart = (int)fmax(0.2 * (options.max_num_saves-count_saves), 1.0); /* reset num_saves_restart */

  restart_action(1);

  for(i=0; i<options.num_particles; i++)
  {
    dnest_which_particle_update = i;
    dnest_which_level_update = level_assignments[i];
    //printf("%d %d %f\n", thistask, i, log_likelihoods[i].value);
    log_likelihoods[i].value = log_likelihoods_cal_restart((void *)((double *)particles+i*particle_offset_double), dnest_arg);
    //printf("%d %d %f\n", thistask, i, log_likelihoods[i].value);
    //due to randomness, the original level assignment may be incorrect. re-asign the level
    while(log_likelihoods[i].value < levels[level_assignments[i]].log_likelihood.value)
    {
      printf("# level assignment decrease %d %f %f %d.\n", i, log_likelihoods[i].value, 
        levels[level_assignments[i]].log_likelihood.value, level_assignments[i]);
      level_assignments[i]--;
    }
    
  }
  return;
}

void dnest_print_particle(FILE *fp, const void *model, const void *arg)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<dnest_num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

void dnest_read_particle(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < dnest_num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file %s.\n", options.sample_file);
      exit(0);
    }
  }
  return;
}

void dnest_restart_action(int iflag)
{
  return;
}

void dnest_accept_action()
{
  return;
}

void dnest_kill_action(int i, int i_copy)
{
  return;
}



