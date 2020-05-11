#ifndef COVIDSIM_KERNELS_H_INCLUDED_
#define COVIDSIM_KERNELS_H_INCLUDED_

enum class Kernel_Type {
    KERNEL_EXPONENTIAL = 1,
    KERNEL_POWER = 2,
    KERNEL_GUSSIAN = 3,
    KERNEL_STEP = 4,
    KERNEL_POWERB = 5, // figure out a better name, power b isn't descriptive
    KERNEL_POWERUS = 6,
    KERNEL_EXPONENTIAL_POWER = 7,
};

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

#endif // COVIDSIM_KERNELS_H_INCLUDED_
