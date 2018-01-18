function FEAT_STRUCT = extract_hrv_features(TIMESTAMP,IBI,feature_set,signal_duration)
% This program is part of CardyLine, a toolbox enabling one-liner heart
% rate variability (HRV) analysis directly from electrocardiogram (ECG).
%
%     EXTRACT_HRV_FEATURES(TIMESTAMP, IBI, feature_set, signal_duration)
%
% computes time- and frequency-domain heart rate variability features from
% an inter-beat interval time series.
%
% Input arguments:
%     TIMESTAMP is an array indicating the timestamp in milliseconds of
%         each inter-beat interval.
%
%     IBI is an array of inter-beat intervals in milliseconds.
%
%     feature_set is a cell array of strings specifying the features to be
%         calculated. Recognized features are 'NN', 'SDNN', 'RMSSD', 'LF',
%         'HF', and 'LFHFratio'.
%
%     signal_duration (optional) is the duration in milliseconds of the
%         signal from which IBI is sampled. Specify if TIMESTAMP does not
%         cover the entire signal range.
%
% The output is a struct with several fields depending on feature_set:
%     signal_quality    - a score indicating the signal quality
%     NN                - mean normal-to-normal interval in milliseconds
%     SDNN, RMSSD       - time-domain measures of HRV in milliseconds
%     LF, HF, LFHFratio - frequency-domain measures of HRV in square
%                         milliseconds
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

if numel(TIMESTAMP)~=numel(IBI)
    error('Lengths of TIMESTAMP and IBI do not match.');
end
IBI = IBI(isfinite(TIMESTAMP)); TIMESTAMP = TIMESTAMP(isfinite(TIMESTAMP));

if any(diff(TIMESTAMP) <= 0) % enforce sorting
    [TIMESTAMP,ord] = unique(TIMESTAMP); IBI = IBI(ord);
    [TIMESTAMP,ord] = sort(TIMESTAMP,'ascend'); IBI = IBI(ord);
end
IBI(abs(IBI-nanmedian(IBI)) > (7.25*mad(IBI,1))) = nan;
NN = nanmean(IBI); SDNN = nanstd(IBI); RMSSD = sqrt(nanmean(diff(IBI).^2));
time_range = range([TIMESTAMP nan]);

if nargin < 3
    feature_set = {};
end

if (nargin < 4) || (signal_duration < time_range)
    signal_duration = time_range + NN;
end

TIMESTAMP = TIMESTAMP - min(TIMESTAMP) + ((signal_duration - time_range)/2);
TIMESTAMP = TIMESTAMP(isfinite(IBI)); IBI = IBI(isfinite(IBI));

FEAT_STRUCT = struct;
FEAT_STRUCT.signal_quality = ...
    (((NN/SDNN) < 200) && (RMSSD < 280) && ((time_range/signal_duration) > 0.6)) * (sum(IBI)/signal_duration);

if ismember('NN', feature_set)
    FEAT_STRUCT.NN = NN;
end
if ismember('SDNN', feature_set)
    FEAT_STRUCT.SDNN = SDNN;
end
if ismember('RMSSD', feature_set)
    FEAT_STRUCT.RMSSD = RMSSD;
end
if any(ismember({'LF','HF','LFHFratio'}, feature_set))
    % Lomb-Scargle periodogram, resolution 1/200 Hz
    TIMESTAMP = TIMESTAMP(:)/1000; IBI = IBI(:).' - NN;
    spec = zeros(1,80);
    for i=8:80
        wt  = 2*pi*(i/200)*TIMESTAMP;
        s = sin(wt); c = cos(wt);
        
        wtau  = 0.5 * atan2(2*s'*c, (c-s)'*(c+s));
        swtau = sin(wtau); cwtau = cos(wtau);
        ss = s*cwtau - c*swtau; cc = c*cwtau + s*swtau;
        
        spec(i) = ((abs(IBI*cc)^2)/(cc'*cc)) + ((abs(IBI*ss)^2)/(ss'*ss));
    end
    spec = spec * (NN/1000);
    
    if any(ismember({'LF','LFHFratio'}, feature_set))
        FEAT_STRUCT.LF = sum([(spec(8)/2) spec(9:29) (spec(30)/2)]) /200;
    end
    if any(ismember({'HF','LFHFratio'}, feature_set))
        FEAT_STRUCT.HF = sum([(spec(30)/2) spec(31:79) (spec(80)/2)]) /200;
    end
    if ismember('LFHFratio', feature_set)
        FEAT_STRUCT.LFHFratio = FEAT_STRUCT.LF / FEAT_STRUCT.HF;
    end
end
end
