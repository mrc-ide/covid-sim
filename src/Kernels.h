#pragma once

extern double *nKernel, *nKernelHR;

void InitKernel(double);
double ExpKernel(double);
double PowerKernel(double);
double PowerKernelB(double);
double PowerKernelUS(double);
double PowerExpKernel(double);
double GaussianKernel(double);
double StepKernel(double);
double numKernel(double);
