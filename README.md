# Genetic Algorithm Synthesizer

Implements the following paper:

    M. Macret, P. Pasquier, and T. Smyth, “Automatic Calibration of Modified FM Synthesis to Harmonic Sounds using Genetic Algorithms,” Proceedings of the 9th Sound and Music Computing conference (SMC 2012), pp. 387–394, 07 2012.

Uses sounds from the following database:

    “Music instrument samples database,” last visited: April 2023. [Online]. Available: http://theremin.music.uiowa.edu/MIS.html

## Requirements

* MATLAB Signal Processing Toolbox
* MATLAB Audio Toolbox
* MATLAB Global Optimization Toolbox

## Running

* First, sounds specified in `ModFmSynth.m` must be placed in `sounds/`
* Next, run `ModFmSynth.m` to run the genetic algorithm over all Nc and sound files
* Finally save off results to a `.mat` file and read in using `PlotConvergence.m`, `PlotHarmonics.m`, or `PlotSpectra.m`
