#!/bin/bash

./plotLimits.py --combineOutputPrefix fullChain --outputSuffix fullChain --minGluinoMass 1000.0 --maxGluinoMass 1750.0 --plotObservedLimits
./plotLimits.py --combineOutputPrefix fullChain --outputSuffix fullChain_r_0_1 --minGluinoMass 1000.0 --maxGluinoMass 1750.0 --plotObservedLimits --contour_signalStrength 0.1
./plotLimits.py --combineOutputPrefix fullChain --outputSuffix fullChain --minGluinoMass 1000.0 --maxGluinoMass 1750.0
./plotLimits.py --combineOutputPrefix fullChain --outputSuffix fullChain_r_0_1 --minGluinoMass 1000.0 --maxGluinoMass 1750.0 --contour_signalStrength 0.1
