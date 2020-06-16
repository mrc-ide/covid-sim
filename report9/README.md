# Main UK runs for Imperial Report 9 Tables 3,4,5 and A1

Results generated as indicated below should exactly match those in the `*T8_NR10*` output files in the `GB_suppress` and `GB_mitigation` folders (see CODECHECK certificate 2020-010
https://doi.org/10.5281/zenodo.3865491). 

Results will be close to but not identical to those in the original report since: (a) the results produced here average over 10 stochastic realisations; (b) the population dataset has changed to Worldpop (open-source) rather than Landscan (closed); (c) the multi-threaded algorithm to create the simulation's household-to-place network has been modified to be single threaded and thus deterministic. Covidsim is now deterministic across platforms (linux, Mac and Windows) and across compilers (gcc, clang, intel and msvc) for a specified number of threads and fixed random number seeds.

## Instructions common to both Suppression and Mitigation runs

Note: `.ps1` files are Microsoft Powershell scripts and `.sh` are bash scripts. 

- Get latest master of https://github.com/mrc-ide/covid-sim and compile with OpenMP enabled and floating point handling set to 'precise' (the cmake file in the repo does this).

- Run `runonce.ps1` or `runonce.sh` in either the `GB_mitigation` or `GB_suppress` folder (no need to run this from both folders). This will create (in the `population` folder) a binary population file `GB_pop2018.bin` from the text population file `GB_pop2018_nhs.txt` and the household-place network file `NetworkGB_8T.bin`. It will also run a simulation for R0=2.4 with no interventions. 

The `GB_pop2018.bin` and `NetworkGB_8T.bin` files produced by this run are required by both the suppression and mitigation runs. Note that other parameter files are not identical between the `GB_suppress` and `GB_mitigation` folders (due to differences in policy duration and triggering). 

##  GB on/off trigger suppression policies in Imperial College NPI report (Report 9 Tables 4 and 5)

- Work in the `GB_suppress` folder

- Run `batch.ps1` (on Windows) or `batch.sh` (on Linux) (script will require editing to substitute your favourite batch job scheduler to run on a compute cluster)

- Run `summariseSup.r` in R (changing cur_path at top of file first)

- Copy first 11 columns of `meanT8_NR10\stats_contain.csv` into the `stats_contain` sheet of a copy of `stats_contain_meanT8_NR10.xlsx`

- Refresh the pivot table in sheet `Tables` of that spreadsheet.

Original results are in `stats_contain_orig.xlsx`. Example results for linux and Windows are in files `stats_contain_meanT8_NR10_*.xlsx`. We also include results for 50 realisations run on 16 threads in `stats_contain_meanT16_NR50`

## Main GB runs for fixed duration mitigation policies in Imperial College  (Report 9 Tables 3 and A1)

- Work in the `GB_mitigation` folder

- Run `batch.ps1` (on Windows) or `batch.sh` (on Linux) (script will require editing to substitute your favourite batch job scheduler to run on a compute cluster)

- Run `summariseMit.r` in R (changing cur_path at top of file first)

- Copy first 11 columns of `meanT8_NR10\stats_mitigation.csv` into the `stats` sheet of a copy of `stats_mitigation_meanT8_NR10.xlsx`

- Refresh the pivot tables in other sheets of the spreadsheet

Original results are in `stats_mitigation_orig.xlsx`. Example results for a range of compilers and platforms are in files `stats_mitigation_meanT8_NR10_*.xlsx`. We also include results for 50 realisations run on 16 threads in `stats_mitigation_meanT16_NR50_msvc.xlsx`.

## Other notes

- The number of model realisations averaged over can by changed by editing the `/NR:10` entry on the command line invocation of the model.
Note that results will only match those provided in the xlsx files if the code is run with `/NR:10`, 8 threads and the current random number seeds.

- The number of threads used by the simulation can be changed by altering the command line option `/c:8`. 
A new network file will need to be generated in that case (via the `runonce` script).
 
- To output every realisation rather than just the average, change the entry for `[Output every realisation]` from 0 to 1 in `preGB_R0=2.0.txt`.

- The last two random number seeds on the command line (which govern the random generated using the epidemic simulation) can be altered without
needing to recreate the network file. The first two random number seeds govern generation of the synthetic population and network file. If they are
changed, a new network file will need to be generated.
