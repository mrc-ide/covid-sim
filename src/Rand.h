#ifndef COVIDSIM_RAND_H_INCLUDED_
#define COVIDSIM_RAND_H_INCLUDED_

/* ranf defines */
const long Xm1 = 2147483563;
const long Xm2 = 2147483399;
const long Xa1 = 40014;
const long Xa2 = 40692;
const long Xa1vw = 2082007225;
const long Xa2vw = 784306273;

/* RANDLIB global variables */
extern int **SamplingQueue;
extern long* Xcg1, *Xcg2;
/* RANDLIB functions */
long ignbin(long, double);
long ignpoi(double);
long ignbin_mt(long, double, int);
long ignpoi_mt(double, int);
double ranf(void);
double ranf_mt(int);
void setall(long *, long *);
double sexpo_mt(int);
double sexpo(void);
long mltmod(long, long, long);
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
