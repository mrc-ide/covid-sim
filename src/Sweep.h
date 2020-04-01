#pragma once

#ifndef SPATIALSIM_SWEEP_H_INCLUDED_
#define SPATIALSIM_SWEEP_H_INCLUDED_

void TravelReturnSweep(double);
void TravelDepartSweep(double);
void InfectSweep(double, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void IncubRecoverySweep(double, int); //added int as argument to record run number: ggilani - 15/10/14
int TreatSweep(double);
//void HospitalSweep(double); //added hospital sweep function: ggilani - 10/11/14
void DigitalContactTracingSweep(double); // added function to update contact tracing number

#endif // SPATIALSIM_SWEEP_H_INCLUDED_
