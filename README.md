# LoopExtrusion A-GRMC calculations

The supplied code runs the A-GRMC calculation on MATLAB. There are two main programs. The LoopExtrusionBothCondensins.m and MakeStructures.m.

The `LoopExtrusionBothCondensins.m` runs an A-GRMC calculation when supplied with the number of condensins, number of bound condensins, extrusion rate etc. The code
```
delete('LoopListsCond1/*')
delete('LoopListsCond2/*')
```
deletes the old files containing the loop locations. The file `pofL_cdf_inv.mat` contains a numerical the step-size distribution function for loop extrusion. The `SamplePofL.m` helper function samples this distribution function. The parameters controlling the calculation are written at the beginning of the file, as follows:

```
monomerSize = 100; %bps
delT = 1; %seconds
pauseTime = 7; %seconds
residenceTimeRealCond1 = 2*60; %seconds
residenceTimeRealCond2 = 6*60; %seconds
nSites = 1000000; %Each monomer is 100 bp

nCond1Total = 2550;
nCond1Bound = 1210;

nCond2Total = 343;
nCond2Bound = 211;

extrusionRateReal = 1500; %bp/s
```
Additional parameters controlling the simulation print frequency and maximum simulation time respectively are set as follows:
```
freq = 360;
tMaxSim = 60*60/delT;
```
To perform the calculation setup a directory as shown as below.
![SetupScreenshot](https://user-images.githubusercontent.com/13065170/221009541-aa2708ad-814c-479b-bdbe-adff04991508.png)

Then one can run the code `LoopExtrusionBothCondensins.m` to generate the loop locations.

The `MakeStructures.m` code is responsible for accepting loop location output from the previous code and generating structures. The parameters are as follows

```
len = 10000; %1 bead = 1 kbp;
b = 1;
rc = 1;
kBT=1;
sqrtkBT = sqrt(kBT);
relLoopStrength = 1;
k = 3/b^2*kBT;
nSamples = 100;

nLoopsCond1 = 2550;
nLoopsCond2 = 343;
nTimes = 9;
runMax = 1;
fileID = fopen('structureTest.xyz','w');
```
 
These two codes can be used to generate any of the results in our paper. Additional codes are placed in the folder named `Additional Codes`
