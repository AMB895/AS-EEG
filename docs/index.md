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
**Outline**


| Step | Outputs |
|-----|-----|
| 1. Pre-processing | Cleaned EEG data ready for analysis |
| 2. Score antisaccade trials | Data table of behavioral performance |
| 3. Event-related spectral perturbation | Data table of average ERSP from significant time-frequency clusters |
| 4. Spontaneous/ evoked power | Data table of mean broadband spontaneous/evoked power |
| 5. Statistical Analyses | R files with model outputs |
| 6. Paper figures | Generate figures for paper and supplementary material |

**Pre-processing**

Pre-preprocessing scripts can be found in `/AS-EEG/preprocessingFunctions/`
1. Raw EEG files were pulled in and trigger values were adjusted to simplier values using [remarkAS.m] for AS trials or [remark.m] for VGS trials
2. Pre-processing pipeline was run with [run_preprocessing_pipelineAS.m] for AS trials or [run_preprocessing_pipelineVGS.m] for VGS trials
    1. Referenced data to mastoids
    2. Bandpass filtered data from 0.5 to 70 Hz
    3. Downsampled data to 150 Hz


**Antisaccade trial scoring**

Antisaccade scoring was done with [score_anti.m](https://github.com/AMB895/AS-EEG/blob/main/eog_scoring/score_anti.m)
1. Read in each subjects raw EOG channels [all_anti.m](https://github.com/AMB895/AS-EEG/blob/main/eog_scoring/all_anti.m)
 
### MATLAB script keys:    
- `C_` : compile data scripts. Use to compile data generated into one matrix/ table across all subjects. (e.g. `C_ersp.m` compiles ERSP data into large 4D matrix)
- `G_` : generate data scripts. Use to generate data (e.g. `G_ersp.m` generates an ERSP data file for each subject)
- `F_` : functions. Functions to run computations (e.g. `F_ersp.m` is a function that computes ERSP)

**Event-related spectral perturbation**

The ERSP for each subject at every channel was computed with [G_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/G_ersp.m) using the functions [F_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/F_ersp.m) and [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m). A 4D matrix (subjects x channels x time x frequency) containing ERSP data for every subject across channels was compiled with [C_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/C_ersp.m)
* [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m) uses Morlet wavelets (## to ## cycles) to compute the time-resolved power spectrum from 3-70 Hz.
* [F_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/F_ersp.m) uses the power spectrums generated from `F_tfpower.m` to compute the ERSP

Significant ERSP time-frequency clusters for frontal and parietal electrodes were found first by running a linear mixed effects model to determine main effects of inverse age and condition (AS/VGS) and their interaction at each ERSP time-frequency point using [G_erspstats.m](https://github.com/AMB895/AS-EEG/blob/main/G_erspstats.m). [G_tfceersp.m](https://github.com/AMB895/AS-EEG/blob/main/G_tfceersp.m) uses the t-values from the linear mixed effects models computed with `G_erspstats.m` as inputs to the threshold-free cluster enhancement algorithm (CITATION) to find significant ERSP clusters from each term in the linear mixed effects models. Finally, [C_clusterersp.m](https://github.com/AMB895/AS-EEG/blob/main/C_clusterersp.m) compiles mean ERSP from significant clusters for each subject into one data table for future analyses. 

**Spontaneous and evoked power**

The spontaneous and evoked power spectrums were computed with [GC_sponevopower.m](https://github.com/AMB895/AS-EEG/blob/main/GC_sponevopower.m) which uses [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m) then averages across time to get broadband spontaneous/ evoked power spectrums. `GC_sponevopower.m` also compiles average broadband power values into data table for future analyses. 

**Statistical analyses**



