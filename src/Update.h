#pragma once

#ifndef SPATIALSIM_UPDATE_H_INCLUDED_
#define SPATIALSIM_UPDATE_H_INCLUDED_

void DoImmune(int);
void DoInfect(int, double, int, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void DoIncub(int, unsigned short int, int, int); //added int as argument to record run number: ggilani - 23/10/14
void DoDetectedCase(int, double, unsigned short int, int);
void DoCase(int, double, unsigned short int, int);
void DoFalseCase(int, double, unsigned short int, int);
void DoRecover(int, int, int); //added int as argument to record run number: ggilani - 23/10/14. Added thread number to record Severity categories in StateT.
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
void DoDeath_FromCriticalorSARI(int, int);

#endif // SPATIALSIM_UPDATE_H_INCLUDED_
