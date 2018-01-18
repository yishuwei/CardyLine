[This is the Beta edition. I am still working on user interface and documentation.]

Requirements: MATLAB R2013b or a later version, Signal Processing Toolbox. 

CardyLine is a toolbox enabling one-liner heart rate variability (HRV) analysis 
directly from electrocardiogram (ECG). After you load continuous ECG data as a 
numeric vector in MATLAB (called `ECG_SIG` below), HRV analysis is only one line of 
code with CardyLine: (syntax may change for future editions)
```
   CardyLine(ECG_SIG, Fs, MAINS_FREQ, hrv_analysis_period)
```
where `Fs` is the sampling rate in Hz, `MAINS_FREQ` is the mains (power line) 
frequency in Hz (50 or 60), and `hrv_analysis_period` sets a duration in seconds 
such that HRV is analyzed for consecutive periods/epochs of this duration (leave it 
out or set it to `[]` if you do not want to divide data into epochs); you get as 
output a struct with various fields including the detected heartbeat instants, as 
well as mean inter-beat interval and time- and frequency-domain HRV measures per 
epoch. (See more details with the command `help CardyLine`.)

CardyLine has built-in functionality for noise/artifact rejection. Nonetheless, if 
you already know there are corrupted samples in the ECG data (e.g., large chunk of 
clipped samples), you can indicate them with a logical vector and pass it as an 
additional argument so that CardyLine can incorporate the information.

The two subroutines `extract_heartbeat_sequence` and `extract_hrv_features` can also 
be run independently. The former implements robust heartbeat detection from 
continuous ECG and the latter calculates HRV measures from an inter-beat interval 
time series. See their usage in detail with commands `help extract_heartbeat_sequence` 
and `help extract_hrv_features`. You can, for example, import results from another 
heartbeat detection program or device and analyze HRV with `extract_hrv_features`. 
The main `CardyLine` routine is nothing but a wrapper that invokes the two 
subroutines back-to-back.

Note: CardyLine needs around 6 times as much memory as the input ECG signal takes. 
Memory usage can be reduced by a third or so, if you replace `filtfilt.m` in Signal 
Processing Toolbox with an optimized version (easy to find on the internet). As 
long as the amount of memory is available, CardyLine is able to efficiently process 
long datasets. Per my testing, up to multiple days of continuous recording with a 
high sampling rate could be processed in one go with CardyLine within half an hour 
on an 8-core machine.

This is a working, but temporary edition of CardyLine. You are welcome to check back 
in the future for updates. If you encounter any problem using the toolbox please do 
not hesitate to contact me.
