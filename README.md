# LoopExtrusionGRMC

The supplied code runs the A-GRMC calculation on MATLAB. There are two main programs. The LoopExtrusionBothCondensins.m and MakeStructures.m.

The LoopExtrusionBothCondensins.m runs an A-GRMC calculation when supplied with the number of condensins, number of bound condensins, extrusion rate etc. The parameters controlling the calculation are written at the beginning of the file, as follows:

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



The MakeStructures.m code is responsible for accepting loop location output from the previous code and generating structures. 
 
These two codes can be used to generate any of the results in our paper.
