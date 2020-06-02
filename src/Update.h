#ifndef COVIDSIM_UPDATE_H_INCLUDED_
#define COVIDSIM_UPDATE_H_INCLUDED_

void DoDetectedCase(int, double, unsigned short int, int);
void DoCase(int, double, unsigned short int, int);
void DoFalseCase(int, double, unsigned short int, int);
void DoDeath(int, int, int); //added int as argument to record run number: ggilani - 23/10/14
void DoPlaceClose(int, int, unsigned short int, int, int);
void DoPlaceOpen(int, int, unsigned short int, int);
void DoTreatCase(int, unsigned short int, int);
void DoProph(int, unsigned short int, int);
void DoProphNoDelay(int, unsigned short int, int, int);
int DoVacc(int, unsigned short int);
void DoVaccNoDelay(int, unsigned short int);
//SEVERITY ANALYSIS
void DoMild(int, int);
void DoSARI(int, int);
void DoCritical(int, int);
void DoRecoveringFromCritical(int, int);
void DoRecover_FromSeverity(int, int);
void DoDeath_FromCriticalorSARIorILI(int, int);

#endif // COVIDSIM_UPDATE_H_INCLUDED_
