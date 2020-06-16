#ifndef COVIDSIM_RAND_H_INCLUDED_
#define COVIDSIM_RAND_H_INCLUDED_

#include <inttypes.h>

/* ranf defines */
const int32_t Xm1 = 2147483563;
const int32_t Xm2 = 2147483399;
const int32_t Xa1 = 40014;
const int32_t Xa2 = 40692;
const int32_t Xa1vw = 2082007225;
const int32_t Xa2vw = 784306273;

/* RANDLIB global variables */
extern int **SamplingQueue;
extern int32_t* Xcg1, *Xcg2;
/* RANDLIB functions */
int32_t ignbin(int32_t, double);
int32_t ignpoi(double);
int32_t ignbin_mt(int32_t, double, int);
int32_t ignpoi_mt(double, int);
double ranf(void);
double ranf_mt(int);
void setall(int32_t *, int32_t *);
double sexpo_mt(int);
double sexpo(void);
int32_t mltmod(int32_t, int32_t, int32_t);
double snorm(void);
double snorm_mt(int);
double fsign(double, double);
//added some new beta, gamma generating functions: ggilani 27/11/14
double gen_norm_mt(double, double, int);
double gen_gamma_mt(double, double, int);
//added some new lognormal sampling functions: ggilani 09/02/17
double gen_lognormal(double, double);
void SampleWithoutReplacement(int, int, int);

#endif // COVIDSIM_RAND_H_INCLUDED_
