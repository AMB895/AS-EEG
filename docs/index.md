# Developmental interactions between oscillatory and spontaneous neural activity support improvements in inhibitory control from adolescence into adulthood

Abstract:

*Adolescence is a critical period of development marked by improvements in cognitive control driven by maturation of executive brain systems. Neural mechanisms supporting inhibitory control—the ability to suppress a prepotent response in favor of generating a goal-directed response—mature into adulthood in parallel with refinements in neural frequency-specific neural oscillations and broadband spontaneous activity, supporting more accurate and reliable performance on response inhibition tasks, such as the antisaccade (AS) oculomotor task. Prior work has established that the maturation of neural circuitry (e.g. synaptic pruning and myelination) contributes to age-related decreases in spontaneous activity, supporting higher fidelity transmission of task-relevant activity, such as frequency-specific neural oscillations which are crucial for AS performance. Using longitudinal electroencephalography (EEG) data from 146 participants (72F, 10-34yos, 253 total sessions) during the AS task, we sought to understand the interplay between frontal oscillatory and spontaneous activity during cognitive processing and how it supports developmental improvements in inhibitory control. We found that spontaneous activity developmentally decreased during the preparation of an inhibitory response which was associated with reduced trial-to-trial variability in response latency indicating reductions in spontaneous activity may increase reliability of inhibitory control processing. Furthermore, increased suppression of preparatory high-beta oscillations was associated with decreased latency variability and lower broadband spontaneous activity, with an age interaction indicating this modulation became more pronounced through development. Together, these results illustrate how developmental enhancements in the coordination of frequency-specific oscillations and spontaneous activity during inhibitory response preparation support improvements in inhibitory control, providing a potential mechanism through which inhibitory control becomes more accurate and reliable from adolescence into adulthood.*

### Project Lead
Abigail M. Beatty

### Faculty Lead
Beatriz Luna 

### Project Start Date
September 2024

### Current Project Status
In progress 

### Datasets
LNCD 7T

### Github Repository
[https://github.com/AMB895/AS-EEG](https://github.com/AMB895/AS-EEG/tree/main)

## CODE DOCUMENTATION
**Outline**

| Step | Outputs |
|------|---------|
| 1. Pre-processing | Cleaned EEG data ready for analysis |
| 2. Score antisaccade trials | Data table of AS behavioral performance |
| 3. Event-related spectral perturbation | Data table of average ERSP from significant time-frequency clusters |
| 4. Spontaneous/ evoked power | Data table of mean broadband spontaneous/evoked power |
| 5. Statistical analyses | R files with model outputs |
| 6. Paper figures | Generate figures for paper and supplementary material |

**Pre-processing**

Specific details on pre-preprocessing scripts found in [EEG_preprocessing](https://github.com/LabNeuroCogDevel/EEG_preprocessing) repository from [LabNeuroCogDevel](https://github.com/LabNeuroCogDevel). Brief overview of pre-processing pipeline:
1. Pre-processing pipeline initiated with [run_preprocessing_V2.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/run_preprocessing_V2.m)
    - Inputs are `task` (e.g. anti, vgs, rest), `maindir` (where things are being saved), and `rawdir` (where raw EEG files are)
Within [run_preprocessing_V2.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/run_preprocessing_V2.m): 
2. Raw EEG files were loaded in and trigger values were adjusted to simplier values using remark functions ([remark_AS_7t.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/functions/remark_functions/remark_AS_7T.m) for AS trials and [remark_vgs.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/functions/remark_functions/remark_vgs.m) for VGS trials)
3. Pre-processing pipeline was run with [preprocessing_pipeline_V2.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/preprocessing_pipeline_V2.m) which involved the following steps:
    1. Referenced data to mastoids [run_filtering.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_filtering.m)
    2. Bandpass filtered data from 0.5 to 70 Hz [run_filtering.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_filtering.m)
    2. Downsampled data to 150 Hz [run_filtering.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_filtering.m) (for current study downsampled to 150 Hz, preprocessing pipeline has since been updated to downsampling to 512 Hz)
    4. Rejected noisy channels [run_channel_rejection.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_channel_rejection.m)
    5. Interpolated rejected channels [run_channel_interpolation](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_channel_interpolation.m)
    6. ICA weights generated
    7. ICA weghts applied to data
    8. Manual ICA rejection of components corresponding to eye blinks, horizontal/ vertical eye movements, ECG activity, channel noise, or movement artifacts [run_manualICArejection.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_channel_interpolation.m)
    9. Notch filtered data at 60 Hz [run_filtering.m](https://github.com/LabNeuroCogDevel/EEG_preprocessing/blob/main/preprocessingTools/steps/run_filtering.m) (current study did this at step 9, updated pre-processing pipeline does this during bandpass filtering step)
    10. Homogenized channel locations, for caps that had 128 channels for comparison to caps with 64 channels
4. Data was segmented into epochs, time-locked to the onset of the preparatory period (for AS task only preparatory period of correct AS trials)


### MATLAB script naming conventions:
- `C_` : compile data scripts. Use to compile data generated into one matrix/ table across all subjects. (e.g. `C_ersp.m` compiles ERSP data into large 4D matrix)
- `G_` : generate data scripts. Use to generate data (e.g. `G_ersp.m` generates an ERSP data file for each subject)
- `F_` : functions. Functions to run computations (e.g. `F_ersp.m` is a function that computes ERSP)


**Event-related spectral perturbation**

The ERSP for each subject at every channel was computed with [G_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/G_ersp.m) using the functions [F_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/F_ersp.m) and [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m). A 4D matrix (subjects x channels x time x frequency) containing ERSP data for every subject across channels was compiled with [C_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/C_ersp.m)
* [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m) uses Morlet wavelets (## to ## cycles) to compute the time-resolved power spectrum from 3-70 Hz.
* [F_ersp.m](https://github.com/AMB895/AS-EEG/blob/main/F_ersp.m) uses the power spectrums generated from [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m) to compute the ERSP

Significant ERSP time-frequency clusters for frontal and parietal electrodes were found first by running a linear mixed effects model to determine main effects of inverse age and condition (AS/VGS) and their interaction at each ERSP time-frequency point using [G_erspstats.m](https://github.com/AMB895/AS-EEG/blob/main/G_erspstats.m). [G_tfceersp.m](https://github.com/AMB895/AS-EEG/blob/main/G_tfceersp.m) uses the t-values from the linear mixed effects models computed with [G_erspstats.m](https://github.com/AMB895/AS-EEG/blob/main/G_erspstats.m) as inputs to the threshold-free cluster enhancement algorithm (CITATION) to find significant ERSP clusters from each term in the linear mixed effects models. Finally, [C_clusterersp.m](https://github.com/AMB895/AS-EEG/blob/main/C_clusterersp.m) compiles mean ERSP from significant clusters for each subject into one data table for future analyses. 

**Spontaneous and evoked power**

The spontaneous and evoked power spectrums were computed with [GC_sponevopower.m](https://github.com/AMB895/AS-EEG/blob/main/GC_sponevopower.m) which uses [F_tfpower.m](https://github.com/AMB895/AS-EEG/blob/main/F_tfpower.m) then averages across time to get broadband spontaneous/ evoked power spectrums. [GC_sponevopower.m](https://github.com/AMB895/AS-EEG/blob/main/GC_sponevopower.m) also compiles average broadband power values into data table for future analyses. 

**Statistics**

All statistical analyses done for the main paper are in [main_sponevo_ersp_analysis.Rmd](https://github.com/AMB895/AS-EEG/blob/main/main_sponevo_ersp_analysis.Rmd) and supplemental statistical analyses are found in [supplement_sponevo_ersp_analysis.Rmd](https://github.com/AMB895/AS-EEG/blob/main/supplement_sponevo_ersp_analysis.Rmd ). 

**Figures**

Main paper figures were generated using [Paper_figs.Rmd](https://github.com/AMB895/AS-EEG/blob/main/Paper_figs.Rmd) and supplemental figures were generated using [Supplemental_figs.Rmd](https://github.com/AMB895/AS-EEG/blob/main/Supplemental_figs.Rmd).

### Project Software 


