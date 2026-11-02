# AS-EEG repository 
Code used to run analyses and generate figures for [insert paper here] 
See [Github page](https://amb895.github.io/AS-EEG/) for detailed project description and code documentation

## Getting started
**Dependencies**
- [EEGLAB 2024.2](https://sccn.ucsd.edu/eeglab/download.php)
- [MATLAB 2021a](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html) or newer
- LIMO tools repository (forked from [LIMO-EEG-Toolbox/limo_tools](https://github.com/LIMO-EEG-Toolbox/limo_tools)

**Modifications**

Some files/folders need to be updated/ made for script to run:
- In `G_ersp.m`:
    - `maindir` = directory to all pre-processed EEG data
    - `merged7t_eeg` = should be a table of AS scores for each subject
    - `epochdirAS` and `epochdirVGS` = directories to EEG epoch data for AS and VGS tasks
    - `subtasksavepath` = folder to save individual data structures for each subject
- In `C_ersp.m`:
    - `maindir` = directory where ERSP data was saved
    - ``taskdir_ersp` = where to save large 4D matrix for each task (AS or VGS)
- In `C_clusterersp.m`:
    - `resultsdir` = where to save table generated
    - `tfcedir` = where to save TFCE outputs
    - `erspdir` = where to find ERSP data
- In `GC_sponevopower.m`
    - `maindirpreproc` = directory to preprocessed EEG data
    - `epochdirAS` and `epochdirVGS` = directories to EEG epoch data for AS and VGS tasks
 
## Authors
Abigail Beatty, University of Pittsburgh

## License
