#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#include "utils.h"

#include "benchmarks/jacobi1d.h"
#include "benchmarks/jacobi2d.h"
#include "benchmarks/jacobi1d_vslope.h"

/* Problem size (in space) between 2 ^ MIN_POW and 2 ^ MAX_POW */
#define MIN_POW 3
#define MIN_ITER_POW 4
#define DEFAULT_RANGE 9
#define DEFAULT_NRUNS 20
#define DEFAULT_RANGE_ITER 3
/*
*1-Dimension problem size :
* L3 4096kB -> 512k = (1 << 9)k < (1<<19)  doubles (64 bits)  1 << 19
* L2 256kB -> 32k = (1<<5)k > (1 << 15) doubles               1 << 15
* L1 32 k -> 4k = (1<<2)k > (1<<12)                           1 << 12
*/
static const int Pbsize_1d =  8 * 4096 * KB;
static const int Num_iters_1d = 1 << 5;
/*
* 2-Dimension problem size :
* L3 4096kB -> 512k = (1 << 9)k < (1<<19)  doubles (64 bits)  1 << 9
* L2 256kB -> 32k = (1<<5)k > (1 << 15) doubles               1 << 7
* L1 32 k -> 4k = (1<<2)k > (1<<12)                           1 << 6
*/
static const int Pbsize_2d = 1 << 8;
static const int Num_iters_2d = 1 << 4;

/* Options */
static int debug_use_defaults = 0;
static int hdiam_flag = 0;
static int test_with_long_flag = 0;
static int verbose_flag = 0;

static const struct option longopts[] = {
  {"brief",         no_argument,            &verbose_flag,        0},
  {"dimt",          required_argument,      NULL,               't'},
  {"dimx",          required_argument,      NULL,               'x'},
  {"dimy",          required_argument,      NULL,               'y'},
  {"hdiam",         no_argument,            &hdiam_flag,          1},
  {"hdmask",       required_argument,      NULL,               'M'},
  {"help",          no_argument,            NULL,               'h'},
  {"iters-range",   required_argument,      NULL,               'i'},
  {"mask",          required_argument,      NULL,               'm'},
  {"nruns",         required_argument,      NULL,               'n'},
  {"range",         required_argument,      NULL,               'r'},
  {"long",          no_argument,            &test_with_long_flag, 1},
  {"use_defaults",  no_argument,            &debug_use_defaults,  1},
  {"verbose",       no_argument,            &verbose_flag,        1},
  {0, 0, 0, 0}
};

static const char * opts_msg[] = {
  "\tno verbose.",
  "\t\tnumber of iterations.",
  "\t\tfirst space dimension size.",
  "\t\tsecond space dimension size.",
  "\trun benchmarks with half-diamonds",
  "\tspecify mask for half-diamond benchmarks",
  "\t\tprint options and other information.",
  "\t= value with --hdiam, iteration space dimensions ranges \n\
    \t\tfrom 2 ^ 3 to 2 ^(3 + value)",
  "\t\t= value, mask as specified below",
  "\t= value. Number of tests / benchmark.",
  "\t= value with --hdiam, space dimension ranges from 2^5 to 2^(5+value)",
  "\t\t run tests with longs for corectness checking",
  "\t use default values when debugging (small values)",
  "\t",
  "/0"
};

struct args_dimt get2dargs(int, int, int, struct benchspec2d);
struct args_dimt get1dargs(int, int, struct benchspec);
struct args_dimt get1dargs_l(int, int, struct benchspec1d_l);
void print_opts();
double test1d(int, int, int, struct benchspec);
double test2d(int, int, int, int, struct benchspec2d);
double test1d_l(int, int, int, struct benchspec1d_l);
void test_suite_hdiam1d(int, int, int, char *, struct benchspec *,FILE *);

void usage(int, int, char **, struct benchspec *, struct benchspec2d *);

int
main(int argc, char ** argv)
{
  int bs;
  /* Hostname for saving records */
  char hostname[512];
  hostname[511] = '\0';
  gethostname(hostname, 511);
  printf("Host : %s\n", hostname);

  /* 1-D benchmarks and misc benchmarks */
  struct benchspec benchmarks[] = {
    {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap, check_default,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive, check_default,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test, check_default,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_SEQUENTIAL (reference)", djbi1d_sequential, check_default,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_HALF_DIAMONDS", djbi1d_half_diamonds_test, check_low_iter,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_HDIAM(GROUPED TILES)", djbi1d_hdiam_grouped_test, check_low_iter,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_HDIAM (TASKS)", djbi1d_hdiam_tasked_test, check_low_iter,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_FROM_PLUTO", djbi1d_from_pluto, check_default,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_HDIAM_(VAR. SLOPE)", djbi1d_hdiam_vslope_t, check_low_iter,
      Pbsize_1d, Num_iters_1d, 1}
  };

  /* 2-D benchmarks */
  struct  benchspec2d benchmarks2d [] = {
    {"JACOBI2D_HALF_DIAMONDS", djbi2d_half_diamonds, check2d_default,
      Pbsize_2d, Pbsize_2d, Num_iters_2d},
    {"JACOBI2D_SEQ ", djbi2d_seq, check2d_default,
      Pbsize_2d, Pbsize_2d, Num_iters_2d},
  };

  /* 1-D benchmarks with long integers */
  struct benchspec1d_l benchmarks_l [] ={
    {"JACOBI1D_SWAP_SEQ (LONG)", ljbi1d_sequential, check_default,
      Pbsize_1d, Num_iters_1d, 1},
    {"JACOBI1D_HALF_DIAMONDS (LONG)", ljbi1d_half_diamonds_test, check_low_iter,
      Pbsize_1d, Num_iters_1d, 1}
  };

  /* Only half-diamonds in 1d */
  struct benchspec hdiam_benchmarks [] =
  {
    {"JACOBI1D_HALF_DIAMONDS", djbi1d_half_diamonds_test, check_low_iter,
      0, 0, 1},
    {"JACOBI1D_HDIAM_(VAR. SLOPE)", djbi1d_hdiam_vslope_t, check_low_iter,
      0, 0, 1},
    {"JACOBI1D_HDIAM(GROUPED TILES)", djbi1d_hdiam_grouped_test, check_low_iter,
      0, 0, 1},
    {"JACOBI1D_HDIAM (USING TASKS)", djbi1d_hdiam_tasked_test, check_low_iter,
      0, 0, 1}
  };

  int nbench = sizeof(benchmarks) / sizeof(struct benchspec);
  int nbench2d = sizeof(benchmarks2d) / sizeof(struct benchspec2d);
  int nbench_hd = sizeof(hdiam_benchmarks) / sizeof(struct benchspec);

  int opt = -1, option_index = 0;
  int dimx = -1, dimy = -1, dimt = -1;
  int range = -1, range_iters = -1;
  int maskl = 0, hdmaskl = 0;
  int nruns = 0;
  char *benchmask, *hdmask;

  hdmask = "111101";

  while ((opt =
    getopt_long(argc, argv, "hi:m:M:r:t:vx:y:", longopts, &option_index))
    != -1) {
      switch (opt) {
        case '0':
        case 'h':
          print_opts();
          usage(nbench, nbench2d, argv, benchmarks, benchmarks2d);
          break;
        case 'i':
          range_iters = atoi(optarg);
          break;
        case 'm':
          benchmask = optarg;
          if ((maskl = strlen(benchmask)) > nbench + nbench2d) {
            printf("Error : not a valid mask ! Your mask must be %i bits long\n",
            nbench + nbench2d);
            return -1;
          }
          break;
        case 'M':
          hdmask = optarg;
          if ((hdmaskl = strlen(hdmask)) > nbench_hd) {
            printf("Error : the hdiams mask you specified is too long."
              "Must be maximum %i bits long.\n", nbench_hd);
          }
        case 'n':
          nruns = atoi(optarg);
        case 'r':
          range = atoi(optarg);
          break;
        case 't':
          dimt = 2 << atoi(optarg);
          break;
        case 'v':
          verbose_flag++;
          break;
        case 'x':
          dimx = (2 << atoi(optarg)) * KB;
          break;
        case 'y':
          dimy = (2 << atoi(optarg)) * KB;
          break;
        case '?':
          break;
      }
  }

  if(range_iters < 0){
    range_iters = DEFAULT_RANGE_ITER;
  }
  if(range < 0){
    range = DEFAULT_RANGE;
  }
/* If debugging can set default parameters */
#ifdef DEBUG
  if(debug_use_defaults >= 0){
    dimt = DEBUG_ITER;
    dimx = DEBUG_SIZE;
    dimy = DEBUG_SIZE;
  }
#endif
/* -------------------------------------- */
/* If hdiam has been set execute the corresponding tests */
  if (hdiam_flag) {

    /* test sute for half diamonds versions (basic/ grouped tiles/ omp tasks)*/
    int num_benchs = sizeof(hdiam_benchmarks) / sizeof(struct benchspec) ;

    /* Store results in a csv file */
    FILE * csv_file;
    char filename[1024];
    strcpy(filename, "./data/jacobi1D_hdiams_at_");
    strcat(filename, hostname);
    strcat(filename, ".csv");
    csv_file = fopen(filename,"w");

    if (csv_file == NULL) {
      why_fopen(errno);
      fprintf(stderr, "Failed to open %s . Aborting ...\n", filename);
      return -1;
    }
    if(hdmask == NULL){
      hdmask = "11101";
    }
    test_suite_hdiam1d(num_benchs, range, range_iters, hdmask, hdiam_benchmarks,
      csv_file);
    fclose(csv_file);
    return 0;
  }

  if(test_with_long_flag > 0){
    int nbench_l = sizeof(benchmarks_l) / sizeof(struct benchspec1d_l);
    for (bs = 0; bs < nbench_l; bs++) {
      test1d_l(nruns, 0, 0, benchmarks_l[bs]);
    }
    return 0;
  }

  double exec_time = 0.0, prev_xctime = -1.0;
  for (bs = 0; bs < maskl; bs++) {
    exec_time = 0.0;
    if (benchmask[bs] == '1' && bs < nbench) {
      printf("Execute test for %s ...\n", benchmarks[bs].name);
      exec_time = test1d(nruns, dimx, dimt, benchmarks[bs]);
    } else if (benchmask[bs] == '1') {
      printf("Execute test for %s ...\n", benchmarks2d[bs - nbench].name);
      exec_time = test2d(nruns, dimx, dimy, dimt, benchmarks2d[bs - nbench]);
    }
    if (exec_time < prev_xctime) {
      printf("%s %10.3f better !%s\n", KRED, exec_time, KRESET);
    }
  }
  return 0;
}



/* Get an arguments structure */
struct args_dimt
get2dargs(int dimx, int dimy, int dimt, struct benchspec2d bs)
{
  if (dimx == 0 || dimy == 0 || dimt == 0) {
    dimx = bs.width;
    dimy = bs.height;
    dimt = bs.iters;
  }
  struct args_dimt res;
  res.width = dimx;
  res.height = dimy;
  res.iters = dimt;
  return res;
}



struct args_dimt
get1dargs(int dimx, int dimt, struct benchspec bs)
{
  if (dimx <= 0 || dimt <= 0) {
    dimx = bs.size;
    dimt = bs.iters;
  }
  struct args_dimt res;
  res.width = dimx;
  res.height = 0;
  res.iters = dimt;
  return res;
}



struct args_dimt
get1dargs_l(int dimx, int dimt, struct benchspec1d_l bs)
{
  struct benchspec bsc = {NULL, NULL, NULL, bs.size, bs.iters, 1};
  return get1dargs(dimx, dimt, bsc);
}

void
print_opts()
{
  int i;
  int nopts = sizeof(longopts) / sizeof(struct option) ;
  for(i = 0; i < nopts - 1; i++){
    printf("--%s %s\n", longopts[i].name, opts_msg[i]);
  }
}



double
test1d(int nruns, int dimx, int dimt, struct benchspec benchmark)
{
  int i;
  double t_accu;
  struct args_dimt args = get1dargs(dimx, dimt, benchmark);
  if(verbose_flag > 0) {
    printf("Problem size : \t\t%i \nNumber of iterations : \t%i\n", args.width,
      args.iters);
  }

  double * data_in = aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*data_in) * args.width);
  double * data_out = aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*data_out) * args.width);

  init_data_1d(args.width, data_in);

  t_accu = 0.0;
  for (i = 0; i <= nruns; i ++) {
    printf("%i,", i);
    if (i == 0) {
      benchmark.variant(args, data_in, data_out);
    } else {
      t_accu += benchmark.variant(args, data_in, data_out);
    }
  }
  printf("Done !\n");

  print_test1d_summary(nruns, DISPLAY_VERBOSE, t_accu, benchmark, data_in,
   data_out);

  double *ref_out = aligned_alloc(CACHE_LINE_SIZE,
    sizeof(*ref_out) * args.width);
  djbi1d_sequential(args, data_in, ref_out);

  long diffs;
  if ((diffs = compare_fast(args.width, data_out, ref_out)) > 0) {
    print_check(diffs, args.width, data_out, ref_out);
  }
  free(data_in);
  free(data_out);
  free(ref_out);
  return t_accu;
}



double
test1d_l(int nruns, int dimx, int dimt, struct benchspec1d_l benchmark)
{
  int i;
  double t_accu;
  struct args_dimt args = get1dargs_l(dimx, dimt, benchmark);
  printf("%i\n", args.width);

  long * data_in = aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*data_in) * args.width);
  long * data_out = aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*data_out) * args.width);

  init_data_1d_l(args.width, data_in);

  t_accu = 0.0;
  for (i = 0; i <= nruns; i ++) {
    printf("%i,", i);
    if (i == 0) {
      benchmark.variant(args, data_in, data_out);
    } else {
      t_accu += benchmark.variant(args, data_in, data_out);
    }
  }
  printf("Done !\n");

  print_test1d_l_summary(nruns, DISPLAY_VERBOSE, t_accu, benchmark, data_in,
    data_out);
  /* TODO : benchmark -> specific reference */
  long * ref_out = aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*ref_out) * args.width);
  ljbi1d_sequential(args, data_in, ref_out);
  long diffs;
  if ((diffs = compare_l(data_out, ref_out, args.width))>0) {
    printf("Differences : %li over %i ( %4.2f )\n", diffs, args.width,
      (float) diffs / args.width);
  }
  free(data_in);
  free(data_out);
  return t_accu;
}



double
test2d(int nruns, int dimx, int dimy, int dimt, struct benchspec2d benchmark)
{
  int i;
  double t_accu = 0.0;

  double **data_in, **data_out;

  struct args_dimt args = get2dargs(dimx, dimy, dimt, benchmark);

  data_in = alloc_double_mx(args.width, args.height);
  data_out = alloc_double_mx(args.width, args.height);
  init_data_2d(args.width, args.height, data_in);

  for (i = 0; i <= nruns; i ++) {
    printf("%i,", i);
    if (i == 0) {
      benchmark.variant(args, data_in, data_out);
    } else {
      t_accu += benchmark.variant(args, data_in, data_out);
    }
  }
  printf("Done !\n");

  print_test2d_summary(nruns, DISPLAY_VERBOSE, t_accu, benchmark, data_in,
    data_out);
  free_mx((void **) data_out, dimx);
  free_mx((void **) data_in, dimx);
  return t_accu;
}



void
test_suite_hdiam1d(int num_benchs, int range, int range_iters, char *hdmask,
 struct benchspec * hdiam_benchmarks, FILE * csv_file)
{
  int bm_no, iters_pow, max_threads, n_threads, pow2, run_no;
  double elapsed_time;
  double *data_in, *data_out;

  max_threads = omp_get_max_threads();
  printf("Maximum number of threads : %i\n", max_threads);

  fprintf(csv_file, "%s\n",
    "iterations,size(Kb),algorithm,threads,time");

  for (iters_pow = MIN_ITER_POW; iters_pow < 3 + range_iters; iters_pow ++) {
    printf("%s%i iterations :%s\n", KRED, 1 << iters_pow, KRESET);

    for (bm_no = 0; bm_no < num_benchs; bm_no ++) {
      if (hdmask[bm_no] == '1') {
        for (pow2 = MIN_POW; pow2 < MIN_POW + range; pow2 ++) {

          struct args_dimt args = { (2 << pow2) * KB, 0 , 1 << iters_pow};

          data_in = malloc(CACHE_LINE_SIZE * sizeof(*data_in) * args.width);
          if (data_in == NULL) {
            printf("Allocation of %i failed, aborting...\n", args.width);
            return;
          }
          data_out = malloc(CACHE_LINE_SIZE * sizeof(*data_out) * args.width);
          if (data_out == NULL) {
            free(data_in);
            printf("Allocation of %i failed, aborting...\n", args.width);
            return;
          }

          init_data_1d(args.width, data_in);

          hdiam_benchmarks[bm_no].variant(args, data_in, data_out);

          for(n_threads = 1; n_threads <= max_threads; n_threads++) {

            omp_set_num_threads(n_threads);

            for (run_no = 0; run_no < DEFAULT_NRUNS; run_no ++) {
                elapsed_time = hdiam_benchmarks[bm_no].variant(args, data_in,
                  data_out);
                fprintf(csv_file, "%i,%i,%i,%i,%f\n", args.iters,
                  args.width / KB, bm_no, n_threads, elapsed_time);
            }
          }

          /* Test output */
          free(data_out);
          free(data_in);
        }
      }
    }
  }
}



void
usage(int nbs, int nbs2d, char ** argv, struct benchspec * bs,
  struct benchspec2d * bs2d)
{
  int i;
  printf("Usage : %s <Mask (%i bit wide) > <Test runs / alg.>\n", argv[0],
    nbs + nbs2d);
  printf("Available benchmarks  :\n");
  printf("%s 1D Benchmarks %s\n", KBLU, KRESET);
  for (i = 0; i < nbs; i++) {
    printf("%i - %s\n", i, bs[i].name);
  }
  printf("%s 2D Benchmarks %s\n", KBLU, KRESET);
  for (i = nbs; i < nbs2d + nbs; i++) {
    printf("%i - %s\n", i, bs2d[i - nbs].name);
  }
  printf("%sExample mask to test JACOBI1D_OMP_NAIVE and JACOBI1D_FROM_PLUTO:%s\
          \n\t .\\test 01000001\n", KBLU, KRESET);
  printf("%sChecking correctness with long versions of algorithms : %s\
          \n\t .\\test l\n", KBLU, KRESET);
}
