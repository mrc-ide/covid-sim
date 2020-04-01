#include <cmath>
#include <stdio.h>
#include <cstddef>
#include "Kernels.h"
#include "Error.h"
#include "Dist.h"
#include "Param.h"

double(*Kernel)(double);
double *nKernel, *nKernelHR;
void InitKernel(int DoPlaces, double norm)
{
	int i, j, im;
	std::ptrdiff_t l, m;
	double t2;

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
	t2 = 0;
#pragma omp parallel for private(i) schedule(static,500) //added private i
	for (i = 0; i <= NKR; i++)
	{
		nKernel[i] = (*Kernel)(((double)i) * P.KernelDelta) / norm;
		nKernelHR[i] = (*Kernel)(((double)i) * P.KernelDelta / NK_HR) / norm;
	}

#pragma omp parallel for schedule(static,500) private(i,j,l,m,im)
	for (i = 0; i < P.NCP; i++)
	{
		l = CellLookup[i] - Cells;
		Cells[l].tot_prob = 0;
		for (j = im = 0; j < P.NCP; j++)
		{
			m = CellLookup[j] - Cells;
			Cells[l].tot_prob += (Cells[l].max_trans[j] = (float)numKernel(dist2_cc_min(Cells + l, Cells + m))) * Cells[m].n;
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
	double t;

	t = -P.KernelShape * log(sqrt(r2) / P.KernelScale + 1);

	return (t < -690) ? 0 : exp(t);
}
double PowerKernelB(double r2)
{
	double t;

	t = 0.5 * P.KernelShape * log(r2 / (P.KernelScale * P.KernelScale));

	return (t > 690) ? 0 : (1 / (exp(t) + 1));
}
double PowerKernelUS(double r2)
{
	double t;

	t = log(sqrt(r2) / P.KernelScale + 1);

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
	double d, t;

	d = sqrt(r2);
	t = -P.KernelShape * log(d / P.KernelScale + 1);

	return (t < -690) ? 0 : exp(t - pow(d / P.KernelP3, P.KernelP4));
}
double numKernel(double r2)
{
	double t, s;

	t = r2 / P.KernelDelta;
	if (t > NKR)
	{
		fprintf(stderr, "** %lg  %lg  %lg**\n", r2, P.KernelDelta, t);
		ERR_CRITICAL("r too large in NumKernel\n");
	}
	s = t * NK_HR;
	if (s < NKR)
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
