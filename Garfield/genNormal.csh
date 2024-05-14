#!/bin/csh

setenv DCSimNtrack 100
setenv DCSimtrackang 0.0

setenv DCSimtrackx 1.0
setenv DCSimOutFile DCsimX1p0A0p0.root
./smalljet

setenv DCSimtrackx 2.0
setenv DCSimOutFile DCsimX2p0A0p0.root
./smalljet

setenv DCSimtrackx 3.0
setenv DCSimOutFile DCsimX3p0A0p0.root
./smalljet 

setenv DCSimtrackx 4.0
setenv DCSimOutFile DCsimX4p0A0p0.root
./smalljet 

setenv DCSimtrackx 5.0
setenv DCSimOutFile DCsimX5p0A0p0.root
./smalljet 

setenv DCSimtrackx 6.0
setenv DCSimOutFile DCsimX6p0A0p0.root
./smalljet 

setenv DCSimtrackx 7.0
setenv DCSimOutFile DCsimX7p0A0p0.root
./smalljet


