# HMM Melancholic Depression

This repository contains the scripts to run the HMM in the current sample
of Melancholic Depression and Non-Melancholic Depression

The main file is do_all_hmm_analyses.m

This file does the following:

- load in timeseries data, ROI-extracted.
- run the HMM several times, in order to
    -- choose a suitable K (number of states)
    -- choose a suitable iteration (most consistent one)

- whenever a HMM is run, ALL information of that calculation is stored in a .mat file.

- this allows the user to just select 1 .mat file, and do further analyses with it:
   -- check Viberbi paths (and summary stats of viterbi paths)
   -- check the FO, Dwell Time, etc
   -- you can also check the state transitions.

