#include <cmath>
#include <iostream>
#include "Kernels.h"
#include "Error.h"
#include "Dist.h"
#include "Param.h"

// To speed up calculation of kernel values we provide a couple of lookup
// tables.
//
// nKernel is a P.NKR+1 element table of lookups nKernel[0] is the kernel
// value at a distance of 0, and nKernel[P.NKR] is the kernel value at the
// largest possible distance (diagonal across the bounding box).
//
// nKernelHR is a higher-resolution table of lookups, also of P.NKR+1
// elements.  nKernelHR[n * P.NK_HR] corresponds to nKernel[n] for
// n in [0, P.NKR/P.NK_HR]
//
// Graphically:
//
// Distance 0            ...                              Bound Box diagonal
//          nKernel[0]   ... nKernel[P.NKR / P.NK_HR] ... nKernel[P.NKR]
//          nKernelHR[0] ... nKernelHR[P.NKR]
double *nKernel, *nKernelHR;

void InitKernel(double norm)
{
	double(*Kernel)(double);

	if (P.KernelType == 1)
		Kernel = ExpKernel;
	else if (P.KernelType == 2)
		Kernel = PowerKernel;
	else if (P.KernelType == 3)
		Kernel = GaussianKernel;
	else if (P.KernelType == 4)
		Kernel = StepKernel;
	else if (P.KernelType == 5)
		Kernel = PowerKernelB;
	else if (P.KernelType == 6)
		Kernel = PowerKernelUS;
	else if (P.KernelType == 7)
		Kernel = PowerExpKernel;
	else
		ERR_CRITICAL_FMT("Unknown kernel type %d.\n", P.KernelType);

#pragma omp parallel for schedule(static,500) default(none) \
		shared(P, Kernel, nKernel, nKernelHR, norm)
	for (int i = 0; i <= P.NKR; i++)
	{
		nKernel[i] = (*Kernel)(((double)i) * P.KernelDelta) / norm;
		nKernelHR[i] = (*Kernel)(((double)i) * P.KernelDelta / P.NK_HR) / norm;
	}

#pragma omp parallel for schedule(static,500) default(none) \
		shared(P, CellLookup)
	for (int i = 0; i < P.NCP; i++)
	{
		Cell *l = CellLookup[i];
		l->tot_prob = 0;
		for (int j = 0; j < P.NCP; j++)
		{
			Cell *m = CellLookup[j];
			l->max_trans[j] = (float)numKernel(dist2_cc_min(l, m));
			l->tot_prob += l->max_trans[j] * (float)m->n;
		}
	}
}

//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
//// **** KERNEL DEFINITIONS

double ExpKernel(double r2)
{
	return exp(-sqrt(r2) / P.KernelScale);
}

double PowerKernel(double r2)
{
	double t = -P.KernelShape * log(sqrt(r2) / P.KernelScale + 1);
	return (t < -690) ? 0 : exp(t);
}

double PowerKernelB(double r2)
{
	double t = 0.5 * P.KernelShape * log(r2 / (P.KernelScale * P.KernelScale));
	return (t > 690) ? 0 : (1 / (exp(t) + 1));
}

double PowerKernelUS(double r2)
{
	double t = log(sqrt(r2) / P.KernelScale + 1);
	return (t < -690) ? 0 : (exp(-P.KernelShape * t) + P.KernelP3 * exp(-P.KernelP4 * t)) / (1 + P.KernelP3);
}

double GaussianKernel(double r2)
{
	return exp(-r2 / (P.KernelScale * P.KernelScale));
}

double StepKernel(double r2)
{
	return (r2 > P.KernelScale * P.KernelScale) ? 0 : 1;
}

double PowerExpKernel(double r2)
{
	double d = sqrt(r2);
	double t = -P.KernelShape * log(d / P.KernelScale + 1);
	return (t < -690) ? 0 : exp(t - pow(d / P.KernelP3, P.KernelP4));
}

double numKernel(double r2)
{
	double t = r2 / P.KernelDelta;
	if (t > P.NKR)
	{
		fprintf(stderr, "** %lg  %lg  %lg**\n", r2, P.KernelDelta, t);
		ERR_CRITICAL("r too large in NumKernel\n");
	}

	double s = t * P.NK_HR;
	if (s < P.NKR)
	{
		t = s - floor(s);
		t = (1 - t) * nKernelHR[(int)s] + t * nKernelHR[(int)(s + 1)];
	}
	else
	{
		s = t - floor(t);
		t = (1 - s) * nKernel[(int)t] + s * nKernel[(int)(t + 1)];
	}
	return t;
}
