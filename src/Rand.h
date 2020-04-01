#ifndef SPATIALSIM_RAND_H_INCLUDED_
#define SPATIALSIM_RAND_H_INCLUDED_

/* ranf defines */
#define Xm1 2147483563
#define Xm2 2147483399
#define Xa1 40014
#define Xa2 40692
#define Xa1vw 2082007225
#define Xa2vw 84306273

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
void setall(long, long);
double sexpo_mt(int);
double sexpo(void);
long mltmod(long, long, long);
double snorm(void);
double snorm_mt(int);
double fsign(double, double);
//added some new beta, gamma generating functions: ggilani 27/11/14
double gen_norm(double, double);
double gen_norm_mt(double, double, int);
double gen_gamma(double, double);
double gen_gamma_mt(double, double, int);
//added some new lognormal sampling functions: ggilani 09/02/17
double gen_lognormal(double, double);
double gen_lognormal_mt(double, double, int);
double sgamma(double);
double sgamma_mt(double, int);
void SampleWithoutReplacement(int, int, int);

#endif // SPATIALSIM_RAND_H_INCLUDED_
