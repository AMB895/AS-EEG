# Developmental interactions between oscillatory and spontaneous neural activity support improvements in inhibitory control from adolescence into adulthood

Abstract: *Adolescence is a critical period of development marked by improvements in cognitive control driven by maturation of executive brain systems. Neural mechanisms supporting inhibitory control—the ability to suppress a prepotent response in favor of generating a goal-directed response—mature into adulthood in parallel with refinements in neural frequency-specific neural oscillations and broadband spontaneous activity, supporting more accurate and reliable performance on response inhibition tasks, such as the antisaccade (AS) oculomotor task. Prior work has established that the maturation of neural circuitry (e.g. synaptic pruning and myelination) contributes to age-related decreases in spontaneous activity, supporting higher fidelity transmission of task-relevant activity, such as frequency-specific neural oscillations which are crucial for AS performance. Using longitudinal electroencephalography (EEG) data from 146 participants (72F, 10-34yos, 253 total sessions) during the AS task, we sought to understand the interplay between frontal oscillatory and spontaneous activity during cognitive processing and how it supports developmental improvements in inhibitory control. We found that spontaneous activity developmentally decreased during the preparation of an inhibitory response which was associated with reduced trial-to-trial variability in response latency indicating reductions in spontaneous activity may increase reliability of inhibitory control processing. Furthermore, increased suppression of preparatory high-beta oscillations was associated with decreased latency variability and lower broadband spontaneous activity, with an age interaction indicating this modulation became more pronounced through development. Together, these results illustrate how developmental enhancements in the coordination of frequency-specific oscillations and spontaneous activity during inhibitory response preparation support improvements in inhibitory control, providing a potential mechanism through which inhibitory control becomes more accurate and reliable from adolescence into adulthood.*

## Project Lead
Abigail M. Beatty

## Faculty Lead
Beatriz Luna 

## Project Start Date
September 2024

## Current Project Status
In progress 

## Datasets
LNCD 7T

## Github Repository
https://github.com/AMB895/AS-EEG

## Code Documentation

Script keys:
    
- `C_` : compile scripts. Use to compile data generated into one matrix/ table across all subjects. (e.g. `C_ersp.m` compiles ERSP data into large 4D matrix)
- `G_` : generate scripts. Use to generate data (e.g. `G_ersp.m` generates ERSP data file for each subject)
- `F_` : functions. Functions to run computations (e.g. `F_ersp.m` is a function that computes ERSP)

**Pre-processing**

Pre-preprocessing scripts can be found in `/AS-EEG/preprocessingFunctions/`
1. Raw EEG files were pulled in and trigger values were adjusted to simplier values using [remarkAS.m] for AS trials or [remark.m] for VGS trials
2. Pre-processing pipeline was run with [run_preprocessing_pipelineAS.m] for AS trials or [run_preprocessing_pipelineVGS.m] for VGS trials
    1. Referenced data to mastoids
    2. Bandpass filtered data from 0.5 to 70 Hz
    3. Downsampled data to 150 Hz


**Antisaccade trial scoring**

Antisaccade scoring was done with [score_anti.m]
1. Read in each subjects raw EOG channels [all_anti.m]
 

**Event-related spectral perturbation**



**Spontaneous and evoked power**

**Statistical analyses**



