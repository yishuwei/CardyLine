function [HEARTBEAT_SEQ,IBI,outlier,ECG_SIG] = extract_heartbeat_sequence(ECG_SIG,Fs,MAINS_FREQ,input_bad_sample,proc_period)
% This program is part of CardyLine, a toolbox enabling one-liner heart
% rate variability (HRV) analysis directly from electrocardiogram (ECG).
%
%     [HEARTBEAT_SEQ,IBI,outlier] = EXTRACT_HEARTBEAT_SEQUENCE(ECG_SIG, Fs, MAINS_FREQ, input_bad_sample, proc_period)
%
% detects and aligns heartbeat in the raw ECG signal.
%
% Input arguments:
%     ECG_SIG is the raw ECG signal (1-lead).
%
%     Fs is the sampling rate in Hz of ECG_SIG.
%
%     MAINS_FREQ is the mains (power line) frequency in Hz (50 or 60).
%         The frequency is notch-filtered before processing.
%
%     input_bad_sample (optional) is a logical vector of the same length as
%         that of ECG_SIG indicating corrupted samples.
%
%     proc_period (optional) sets a duration in seconds such that ECG_SIG
%         is processed in consecutive periods of this duration.
%
% Output:
%     HEARTBEAT_SEQ is an array of sample indices corresponding to the
%         detected heartbeat instants.
%
%     IBI is an array of inter-beat intervals in milliseconds calculated
%         from HEARTBEAT_SEQ.
%
%     outlier is a logical array indicating whether each entry in IBI is an
%         outlier.
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
ECG_SIG = ECG_SIG(:).';
Nsamp = numel(ECG_SIG);

if (numel(MAINS_FREQ)~=1) || (MAINS_FREQ < 47) || (MAINS_FREQ > 63)
    error('Please set MAINS_FREQ to either 50 or 60 Hz.');
end

if (nargin < 4) || (numel(input_bad_sample)~=Nsamp)
    input_bad_sample = false(1,Nsamp);
else
    input_bad_sample = logical(input_bad_sample(:)).';
end

if nargin < 5
    proc_period = Nsamp;
else
    proc_period = min([Nsamp, round(proc_period*Fs)]);
end

%% Filter signal
fprintf('Initiating...');

[b1,a1] = butter(3, (MAINS_FREQ+[-2 2])*2/Fs,'stop');
[b2,a2] = butter(3, [5 45]*2/Fs);
b = conv(b1,b2); a = conv(a1,a2);
ECG_SIG = filtfilt(b,a, ECG_SIG);
fprintf('..');

[b,a] = butter(3, [5 15]*2/Fs);
ECG_SIGF = filtfilt(b,a, ECG_SIG);
fprintf('..');

%% Heartbeat detection using Pan-Tompkins algorithm
fprintf('\nRunning heartbeat detection');

sbeg = [(1:proc_period:Nsamp) (Nsamp+1)];

seg_beat = cell(1,numel(sbeg)); seg_maybe_beat = cell(1,numel(sbeg));
for seg=1:(numel(sbeg)-1)
    bad_sample = deviant_indicator(ECG_SIGF(sbeg(seg):(sbeg(seg+1)-1)), Fs, 0.50, input_bad_sample(sbeg(seg):(sbeg(seg+1)-1)));
    
    [~,wbeg,wend] = get_bouts(~bad_sample); wbeg = wbeg+sbeg(seg)-1; wend = wend+sbeg(seg)-1;
    beat = arrayfun(@(bi,ei)[ (PanTompkins(ECG_SIG(bi:ei),Fs,nan,ECG_SIGF(bi:ei))+bi-1) ], wbeg, wend, 'UniformOutput',false);
    seg_beat{seg} = cell2mat(beat);
    
    [xspan,xbeg,xend] = get_bouts(bad_sample); xbeg = xbeg(xspan>(2.1*Fs))+sbeg(seg)-1; xend = xend(xspan>(2.1*Fs))+sbeg(seg)-1;
    beat = arrayfun(@(bi,ei)[ (PanTompkins(ECG_SIG(bi:ei),Fs,nan,ECG_SIGF(bi:ei))+bi-1) ], xbeg, xend, 'UniformOutput',false);
    beat( arrayfun(@(bi,ei)any(input_bad_sample(bi:ei)), xbeg, xend) ) = {nan(1,0)};
    seg_maybe_beat{seg} = cell2mat(beat);
    
    fprintf('.');
end

%% Heartbeat marker alignment
fprintf('\nRunning heartbeat alignment');

for seg=1:(numel(sbeg)-1)
    % 5 sec padding each end
    wbeg = max(1, sbeg(seg)-round(5*Fs));
    wend = min(Nsamp, sbeg(seg+1)+round(5*Fs)-1);
    
    if seg==1
        beat = seg_beat{1};
    else
        ipad = find(seg_beat{seg-1} >= wbeg, 1,'first');
        beat = [seg_beat{seg-1}(ipad:end) seg_beat{seg}];
    end
        ipad = find(seg_beat{seg+1} <= wend, 1,'last');
        beat = [beat seg_beat{seg+1}(1:ipad)];
    
    beat = align_heartbeat_markers(ECG_SIG(wbeg:wend), Fs, beat-wbeg+1, seg_maybe_beat{seg}-wbeg+1) + wbeg - 1;
    beat(input_bad_sample(beat)) = [];
    
    ibeg = find(beat >= sbeg(seg), 1,'first');
    iend = find(beat < sbeg(seg+1), 1,'last');
    seg_beat{seg} = beat(ibeg:iend);
    
    fprintf('.');
end

[HEARTBEAT_SEQ, IBI, outlier] = calibrate_heartbeat_sequence(cell2mat(seg_beat), Fs);

fprintf('\n');

end

function bad_sample = deviant_indicator(ECG_SIGF,Fs,smooth_window,input_bad_sample)

ECG_SIGI = conv(ECG_SIGF, [-1 -2 0 2 1]/8);
ECG_SIGI = abs(hilbert(ECG_SIGI)).^2;

smooth_window = uint32(smooth_window*Fs);
if (idivide(smooth_window,2)*2==smooth_window); smooth_window = smooth_window + 1; end
    
ECG_SIGI = conv(ECG_SIGI, ones(1,smooth_window)/double(smooth_window));

trail = 2 + idivide(smooth_window - 1, 2,'floor');
ECG_SIGI = ECG_SIGI((trail+1):(end-trail));

% find large peaks
[peak_amp,peak_idx] = findpeaks(ECG_SIGI, 'MinPeakDistance',floor(smooth_window));%'MinPeakHeight',median(cumsum(ECG_SIGI)./(1:Nsamp)));
[peak_idx,ord] = sort(peak_idx,'ascend'); peak_amp = peak_amp(ord);

badpeak = input_bad_sample(peak_idx);
badpeak = badpeak | (peak_amp > (median(peak_amp(~badpeak))+5*mad(peak_amp(~badpeak),1)));

bad_sample = input_bad_sample;
for p=sort(find(badpeak),'ascend')
    if p==1
        ibeg = 1;
    elseif badpeak(p-1)
        ibeg = peak_idx(p-1);
    else
        ibeg = find(ECG_SIGI(peak_idx(p-1):peak_idx(p)) > peak_amp(p-1), 1,'first') + peak_idx(p-1) - 1;
    end
    if p==numel(peak_idx)
        iend = numel(bad_sample);
    elseif badpeak(p+1)
        iend = peak_idx(p+1);
    else
        iend = find(ECG_SIGI(peak_idx(p):peak_idx(p+1)) > peak_amp(p+1), 1,'last') + peak_idx(p) - 1;
    end
    bad_sample(ibeg:iend) = true;
end

[goodspan,ibeg,iend] = get_bouts(~bad_sample);
for b=1:numel(ibeg); if goodspan(b)<(2*Fs); bad_sample(ibeg(b):iend(b)) = true; end;end
end

function [boutspan,ibeg,iend] = get_bouts(OnOffSeq)

beg_end = diff([0 OnOffSeq 0]);
ibeg = sort(find(beg_end > 0),'ascend');
iend = sort(find(beg_end < 0),'ascend');
boutspan = iend - ibeg;
iend = iend - 1;

end

function [beat,PEAK_SIGN] = PanTompkins(ECG_SIG,Fs,PEAK_SIGN,ECG_SIGF)
% Pan & Tompkins (1985) thresholding algorithm, with some adjustments

if isnan(PEAK_SIGN)
    % use first 10 sec to figure out peak sign
    [~,PEAK_SIGN] = PanTompkins(ECG_SIG(1:min(10*Fs,end)),Fs,0,ECG_SIGF(1:min(10*Fs,end)));
end

smooth_window = round(0.150*Fs);

ECG_SIGI = conv(ECG_SIGF, [-1 -2 0 2 1]/8);
ECG_SIGI = ECG_SIGI.^2;
ECG_SIGI = conv(ECG_SIGI, ones(1,smooth_window)/smooth_window);

[peak_amp,peak_idx] = findpeaks(ECG_SIGI, 'MinPeakDistance',floor(0.20*Fs));
[peak_idx,ord] = sort(peak_idx,'ascend'); peak_amp = peak_amp(ord);

%% Initialize thresholds
SPKI = max(ECG_SIGI(1:min(2*Fs,end)));
NPKI = 0.5*mean(ECG_SIGI(1:min(2*Fs,end)));
THRESHOLD_I1 = 0.75 * NPKI + 0.25 * SPKI;

SPKF = max(abs(ECG_SIGF(1:min(2*Fs,end))));
NPKF = 0.5*mean(abs(ECG_SIGF(1:min(2*Fs,end))));
THRESHOLD_F1 = 0.75 * NPKF + 0.25 * SPKF;

%% Adaptive thresholding
if PEAK_SIGN > 0.5
    locate_peak = @max;
elseif PEAK_SIGN < -0.5
    locate_peak = @min;
else
    locate_peak = @(x)max(abs(x));
end

beatI =[]; beatF =[];
stable_ibi = -1; mean_ibi = -1; stable = true;
p = 1;
while p <= numel(peak_idx)
    
    PEAKI = peak_amp(p); peakI_idx = peak_idx(p);
    
    found_R_wave = false;
    search_back = false;
    
    %% Locate peak in BP signal
    if peakI_idx > smooth_window
        [PEAKF,peakF_idx] = locate_peak(ECG_SIGF((peakI_idx-smooth_window):min(peakI_idx,end)));
        peakF_idx = peakI_idx - smooth_window + peakF_idx - 1;
    else
        [PEAKF,peakF_idx] = locate_peak(ECG_SIGF(1:peakI_idx));
    end
    
    %% Search back
    ref_ibi = -1;
    if stable_ibi > 0
        ref_ibi = stable_ibi;
    elseif mean_ibi > 0
        ref_ibi = mean_ibi;
    end
    
    if (ref_ibi > 0) && ((peakF_idx-beatF(end)) > (1.66*ref_ibi))
        [PEAKI_sb,peakI_idx_sb] = max(ECG_SIGI((beatI(end)+round(0.360*Fs)):(peakI_idx-round(0.200*Fs))));
        peakI_idx_sb = beatI(end) + round(0.360*Fs) + peakI_idx_sb - 1;
        
        if PEAKI_sb > (THRESHOLD_I1/2)
            % search with a more lenient range than in the normal case
            search_upper_bound = max(beatF(end)+round(ref_ibi), peakI_idx_sb+round(0.150*Fs));
            [PEAKF_sb,peakF_idx_sb] = locate_peak(ECG_SIGF((peakI_idx_sb-smooth_window):min(search_upper_bound,peakF_idx)));
            peakF_idx_sb = peakI_idx_sb - smooth_window + peakF_idx_sb - 1;
            
            if abs(PEAKF_sb) > (THRESHOLD_F1/2)
                beatI = [beatI peakI_idx_sb];
                beatF = [beatF peakF_idx_sb];
                
                SPKI = 0.25 * PEAKI_sb + 0.75 * SPKI;
                SPKF = 0.25 * abs(PEAKF_sb) + 0.75 * SPKF;
                search_back = true;
                found_R_wave = true;
            end
        end
    end
    
    %% Determine T or R wave
    if (~found_R_wave) && (PEAKI > THRESHOLD_I1)
        
        found_T_wave = false;
        
        if ~isempty(beatF) && ((peakF_idx-beatF(end)) < (0.360*Fs))
            if abs(ECG_SIGI(peakI_idx)-ECG_SIGI(peakI_idx-round(0.075*Fs))) < (0.5*abs(ECG_SIGI(beatI(end))-ECG_SIGI(beatI(end)-round(0.075*Fs))))
                found_T_wave = true;
            end
        end
        
        if (~found_T_wave) && (abs(PEAKF) > THRESHOLD_F1)
            beatI = [beatI peakI_idx];
            beatF = [beatF peakF_idx];
            
            SPKI = 0.125 * PEAKI + 0.875 * SPKI;
            SPKF = 0.125 * abs(PEAKF) + 0.875 * SPKF;
            found_R_wave = true;
        end
    end
    
    if ~found_R_wave
        NPKI = 0.125 * PEAKI + 0.875 * NPKI;
        NPKF = 0.125 * abs(PEAKF) + 0.875 * NPKF;
    else
        %% NEW: Delete false positive when unstable
        if (stable_ibi > 0) && ~stable
            last_ibi = beatF(end) - beatF(end-1);
            previous_ibi = beatF(end-1) - beatF(end-2);
            
            crit = abs(last_ibi + previous_ibi - mean_ibi);
            
            % merge two epochs if the result deviates less from mean
            if (abs(last_ibi-mean_ibi) > crit) && (abs(previous_ibi-mean_ibi) > crit)
                beatI(end-1) = [];
                beatF(end-1) = [];
            end
        end
        
        %% Update mean interbeat interval
        if numel(beatF) > 2
            last_ibi = beatF(end) - beatF(end-1);
            stable = ~((last_ibi < (0.92*mean_ibi)) || (last_ibi > (1.16*mean_ibi)));
            mean_ibi = mean(diff(beatF(max(end-8,1):end)));
            
            if stable
                stable_ibi = mean_ibi;
            end
        elseif numel(beatF)==2
            mean_ibi = diff(beatF);
        end
    end
    
    %% Adapt thresholds
    THRESHOLD_I1 = 0.75 * NPKI + 0.25 * SPKI;
    THRESHOLD_F1 = 0.75 * NPKF + 0.25 * SPKF;
    if ~stable
        THRESHOLD_I1 = 0.5 * THRESHOLD_I1;
        THRESHOLD_F1 = 0.5 * THRESHOLD_F1;
    end
    
    if ~search_back
        p = p + 1;
    end
end

%% Locate peaks in raw signal
beatF = unique(beatF);
Nbeat = numel(beatF);
beat = zeros(1,Nbeat);
for p=1:Nbeat
    peakF_idx = beatF(p);

    if peakF_idx > round(0.075*Fs)
        [~,peak_idx] = locate_peak(ECG_SIG((peakF_idx-round(0.075*Fs)):min(peakF_idx+round(0.075*Fs),end)));
        peak_idx = peakF_idx - round(0.075*Fs) + peak_idx - 1;
    else
        [~,peak_idx] = locate_peak(ECG_SIG(1:(peakF_idx+round(0.075*Fs))));
    end
    
    beat(p) = peak_idx;
end
beat = sort(unique(beat),'ascend');

if abs(PEAK_SIGN) < 0.5
    PEAK_SIGN = sign(median(ECG_SIGF(beatF)));
end
end

function [beat,ECG_template] = align_heartbeat_markers(ECG_SIG,Fs,beat,maybe_beat)

Nsamp = numel(ECG_SIG);

%% Correct for false positives
[maybe_beat, ibi, oibi] = calibrate_heartbeat_sequence([beat maybe_beat],Fs);
ref_ibi = nanmedian(ibi(~oibi)*(Fs/1000)); if ~isfinite(ref_ibi); ref_ibi = 60*Fs/72; end
crit_ibi = max(0.3*Fs, ref_ibi-7.5*mad(ibi(~oibi),1));

ibi = ibi(isfinite(maybe_beat)); oibi = oibi(isfinite(maybe_beat));
dibi = diff(ibi); odibi = dibi < (prctile(dibi,25)-1.5*iqr(dibi));

maybe_beat = maybe_beat(isfinite(maybe_beat));
maybe_beat( oibi & (([false odibi] & [odibi false]) | (ibi < crit_ibi)) ) = [];

%% Make ECG template
beat = beat(ismember(beat, maybe_beat));

pre_R = min([round(0.2*Fs), (max(beat)-1)]); post_R = min([round(min(0.6*Fs, ref_ibi-0.1*Fs)), (Nsamp-min(beat))]);
beat = beat((beat > pre_R) & ((beat+post_R) <= Nsamp));

Ntemplate = pre_R + post_R + 1;
ECG_template = nanmean(ECG_SIG(repmat(beat(:),1,Ntemplate) + repmat((-pre_R):post_R,numel(beat),1)), 1);
ECG_template = ECG_template - mean(ECG_template); ECG_template = ECG_template / norm(ECG_template);

%% Fit ECG template to raw signal
beat = maybe_beat((maybe_beat > ref_ibi) & ((maybe_beat+post_R) <= Nsamp));

if any(isnan(ECG_template))
    beat = [];
else
    %% Efficient computation of sliding-window Pearson correlation
    mvm = 0; mvv = zeros(1,Nsamp+post_R+1);
    for i=1:Ntemplate
        m = mvm; mvm = mvm + (ECG_SIG(i)/Ntemplate);
        mvv(i+1) = mvv(i) + ECG_SIG(i)*(ECG_SIG(i)-mvm-m);
    end
    for i=(Ntemplate+1):Nsamp
        m = mvm; d = ECG_SIG(i)-ECG_SIG(i-Ntemplate); mvm = mvm + (d/Ntemplate);
        mvv(i+1) = mvv(i) + d*(ECG_SIG(i)-mvm+ECG_SIG(i-Ntemplate)-m);
    end
    for i=(Nsamp+1):(Nsamp+post_R)
        m = mvm; mvm = mvm - (ECG_SIG(i-Ntemplate)/Ntemplate);
        mvv(i+1) = mvv(i) - ECG_SIG(i-Ntemplate)*(ECG_SIG(i-Ntemplate)-mvm-m);
    end
    C = fliplr(xcorr(fliplr(ECG_SIG), fliplr(ECG_template)));
    C = C((post_R+1):(Nsamp+post_R)) ./ sqrt(mvv((post_R+2):end));
    
    C = abs(C);
    
    ref_ibi = floor(ref_ibi); CP = 1;
    for p=1:numel(beat)
        wbeg = max(CP + pre_R, beat(p)-ref_ibi);
        
        if p < numel(beat)
            wend = min(beat(p)+ref_ibi, beat(p+1)-pre_R);
        else
            wend = min(beat(p)+ref_ibi, Nsamp-post_R);
        end
        
        [~,CP] = max([nan C(wbeg:wend)]); CP = CP + wbeg - 1;
        
        beat(p) = CP;
    end
    
    %% Correct for false negatives
    maybe_beat = [0 beat (Nsamp+1)]; oibi = diff(maybe_beat) > (1.66*ref_ibi);
    
    ins = arrayfun(@(bi)insert_markers(C,maybe_beat(bi),maybe_beat(bi+1),pre_R,1.66*ref_ibi), find(oibi), 'UniformOutput',false);
    
    beat = round(sort(horzcat(beat,ins{:}),'ascend'));
    
    beat(C(beat) <= 0.5) = [];
end
end

function ins = insert_markers(C,ins_beg,ins_end,min_interval,max_interval)

[~,CP] = max(C((ins_beg+min_interval):(ins_end-min_interval)));
CP = CP + ins_beg + min_interval - 1;

ins = CP;

if (CP - ins_beg) > max_interval
    ins = [insert_markers(C,ins_beg,CP,min_interval,max_interval) ins];
end
if (ins_end - CP) > max_interval
    ins = [ins insert_markers(C,CP,ins_end,min_interval,max_interval)];
end
end

function [sorted_beat,orig_ibi,outlier] = calibrate_heartbeat_sequence(beat,Fs)

if any(diff(beat(isfinite(beat))) <= 0) % enforce sorting
    sorted_beat = sort(unique(beat),'ascend');
else
    sorted_beat = beat(find(isfinite(beat), 1,'first'):end);
end

ibi = diff(sorted_beat)*(1000/Fs); beat = sorted_beat(2:end);
ibi = ibi(isfinite(beat)); beat = beat(isfinite(beat));
orig_ibi = ibi; dorig_ibi = diff(ibi);

ibi((ibi < 300) | (ibi > 1700)) = nan; Noutlier = -1;

while (numel(ibi) > 1) && (sum(isnan(ibi)) > Noutlier)
    Noutlier = sum(isnan(ibi)); dibi = diff(ibi);
    
    outlier = (dorig_ibi > (prctile(dibi,75)+1.5*iqr(dibi))) | (dorig_ibi < (prctile(dibi,25)-1.5*iqr(dibi))) | isnan(dorig_ibi);
    outlier = [true outlier] & [outlier true];
    
    for odd=find(outlier)
        vcn = odd+[-2 -1 1 2]; vcn = vcn((vcn > 0) & (vcn <= numel(ibi)));
        
        vcn = nanmedian(ibi(vcn));
        if ((ibi(odd) >= (0.8*vcn)) && (ibi(odd) <= (1.2*vcn))); outlier(odd)=false; end
    end
    
    ibi(outlier) = nan;
    
    outlier = abs(dibi-nanmedian(dibi)) > (9*mad(dibi,1));
    outlier = [false outlier] | [outlier false];
    outlier = outlier | ( [true isnan(dibi)] & [isnan(dibi) true] & (abs(ibi-nanmedian(ibi)) > (7.5*mad(ibi,1))) );
    ibi(outlier) = nan;
end

outlier = isnan(ibi);

if ~isempty(sorted_beat)
    sorted_beat = [sorted_beat(1) beat];
    orig_ibi = [nan orig_ibi];
    outlier = [true outlier];
end
end
