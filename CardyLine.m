function [FEAT_STRUCT,ECG_SIG] = CardyLine(ECG_SIG,Fs,MAINS_FREQ,hrv_analysis_period,bad_sample)
% CardyLine enables one-liner heart rate variability (HRV) analysis
% directly from electrocardiogram (ECG).
%
%     CardyLine(ECG_SIG, Fs, MAINS_FREQ, hrv_analysis_period, bad_sample)
%
% computes time- and frequency-domain heart rate variability features over
% time from the raw ECG signal.
%
% Input arguments:
%     ECG_SIG is the raw ECG signal (1-lead).
%
%     Fs is the sampling rate in Hz of ECG_SIG.
%
%     MAINS_FREQ is the mains (power line) frequency in Hz (50 or 60).
%         The frequency is notch-filtered before processing.
%
%     hrv_analysis_period (optional) sets a duration in seconds such that
%         HRV is analyzed for consecutive periods/epochs of this duration.
%         For example, to compute HRV features for every 5 minutes, set
%         hrv_analysis_period to 300.
%
%     bad_sample (optional) is a logical vector of the same length as that
%         of ECG_SIG indicating corrupted samples.
%
% The output is a struct with several fields:
%     sample_rate                 - the original sampling rate in Hz
%     heartbeat_sample            - sample indices corresponding to the
%                                   detected heartbeat instants
%     filtered_interbeat_interval - inter-beat intervals in milliseconds,
%                                   with outliers replaced by NaN
%     epoch_Nsample               - number of samples in each epoch
%     epoch_start_sample          - sample indices corresponding to the
%                                   start of each epoch
%     signal_quality              - scores indicating the signal quality
%                                   within each epoch
%     NN                          - mean normal-to-normal interval in
%                                   milliseconds within each epoch
%     SDNN, RMSSD                 - time-domain measures in milliseconds
%                                   of HRV within each epoch
%     LF, HF, LFHFratio           - frequency-domain measures in square
%                                   milliseconds of HRV within each epoch
%
% See also EXTRACT_HEARTBEAT_SEQUENCE, EXTRACT_HRV_FEATURES
%
%
% Author: Yishul Wei. All rights reserved.
%
% CardyLine is intended to be an academic software toolbox. Permission to
% use, copy, modify, and distribute the software and its documentation for
% not-for-profit purposes is granted to any person obtaining a copy of the
% source code, provided that this permission notice appear in all copies.
% For other uses, please contact the author (Y. Wei).
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

if ~(isnumeric(ECG_SIG) && isvector(ECG_SIG))
    error('Please enter a 1-d ECG signal.');
end
Nsamp = numel(ECG_SIG);

if (numel(MAINS_FREQ)~=1) || (MAINS_FREQ < 47) || (MAINS_FREQ > 63)
    error('Please set MAINS_FREQ to either 50 or 60 Hz.');
end

if (nargin < 4) || isempty(hrv_analysis_period)
    hrv_analysis_period = Nsamp;
else
    hrv_analysis_period = min([Nsamp, round(hrv_analysis_period*Fs)]);
end

if nargin < 5
    bad_sample = [];
end

[HEARTBEAT_SEQ, IBI, outlier, ECG_SIG] = ...
    extract_heartbeat_sequence(ECG_SIG, Fs, MAINS_FREQ, bad_sample, 1800);

IBI(outlier) = nan;

%% Analyze HRV per epoch
ebeg = 1:hrv_analysis_period:Nsamp;
ibeg = arrayfun(@(x)min([(numel(HEARTBEAT_SEQ)+1), find(HEARTBEAT_SEQ >= x, 1,'first')]), ebeg);
iend = arrayfun(@(x)max([1, find(HEARTBEAT_SEQ < (x+hrv_analysis_period), 1,'last')]), ebeg);

hrv_feat = struct2table(arrayfun(@(bi,ei)extract_hrv_features( ...
    HEARTBEAT_SEQ(bi:ei)*(1000/Fs), ...
    IBI(bi:ei), ...
    {'NN','SDNN','RMSSD','LFHFratio'}, ...
    hrv_analysis_period*(1000/Fs)), ibeg, iend));

FEAT_STRUCT = struct;
FEAT_STRUCT.sample_rate = Fs;
FEAT_STRUCT.heartbeat_sample = HEARTBEAT_SEQ;
FEAT_STRUCT.filtered_interbeat_interval = IBI;
FEAT_STRUCT.epoch_Nsample = hrv_analysis_period;
FEAT_STRUCT.epoch_start_sample = ebeg;
for f=hrv_feat.Properties.VariableNames; FEAT_STRUCT.(char(f))=hrv_feat.(char(f)).'; end

end
