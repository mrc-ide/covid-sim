#pragma once

#ifndef SPATIALSIM_KERNELS_H_INCLUDED_
#define SPATIALSIM_KERNELS_H_INCLUDED_

extern double *nKernel, *nKernelHR;

void InitKernel(int, double);
double ExpKernel(double);
double PowerKernel(double);
double PowerKernelB(double);
double PowerKernelUS(double);
double PowerExpKernel(double);
double GaussianKernel(double);
double StepKernel(double);
double numKernel(double);

#endif // SPATIALSIM_KERNELS_H_INCLUDED_
