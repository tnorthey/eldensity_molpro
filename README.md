# eldensity_molpro

## Usage

Main code: elDensity3.m  

#### Inputs
  
*mldfile* (string): Path to Molden file output from Molpro.  
*N* (integer): Cubic box symmetric about the origin is size(N,N,N). Default is N=79 (somewhat tested value).  

## Description

Calculates electron density from a Molden file output from Molpro. The electron density is calculated within a cubic box symmetric about the origin with side-length (maximum atom-atom distance + 7.0) Bohr, or (maximum atom-origin distance) + 7.0 Bohr (whichever is larger); or exactly 7.0 Bohr for single atoms. The cubic grid size is (N,N,N).
