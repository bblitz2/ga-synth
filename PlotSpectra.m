clear variables;
close all;

modfm = @(t,fm,fc,I) exp(I.*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t); % Normalize? -- should be handled by w
classicfm = @(t,fm,fc,I) cos(2*pi*fc.*t + I.*sin(2*pi*fm.*t));
modfm_synth = @(t,fm,k,I,w) sum(w.*modfm(t,fm,k.*fm,I));
classicfm_synth = @(t,fm,k,I,w) sum(w.*classicfm(t,fm,k.*fm,I));

fs = 44100;

%% Process target sounds

% Hyperparameter definitions
params = [];
params.fs = 44100;  % Sample rate
params.Nharm = 10;  % Number of harmonics in analysis
params.Lw = 0.010*params.fs;  % Analysis window size in samples
params.zpf = 4;  % Zero padding factor


params1 = params;
x1 = audioread(fullfile("sounds", "trumpet", "Trumpet.novib.ff.C6.stereo.aif"));
if size(x1,2) == 2
    x1 = (x1(:,1) + x1(:,2)) / 2;  % Convert from stereo to mono
end
params1.f0 = median(pitch(x1, fs, "Range", [50 2000], "Method", "SRH"));
params1.fm = params1.f0;  % Modulation frequency
[T_trumpet, spec_trumpet] = harmonic_analysis(x1,params1.fs,params1.f0,params1.Lw,params1.zpf,params1.Nharm);

params2 = params;
x2 = audioread(fullfile("sounds", "viola", "Viola.arco.ff.sulG.C4.stereo.aif"));
if size(x2,2) == 2
    x2 = (x2(:,1) + x2(:,2)) / 2;  % Convert from stereo to mono
end
params2.f0 = median(pitch(x2, fs, "Range", [50 2000], "Method", "SRH"));
params2.fm = params2.f0;  % Modulation frequency
[T_viola, spec_viola] = harmonic_analysis(x2,params2.fs,params2.f0,params2.Lw,params2.zpf,params2.Nharm);

%% Synthesize best candidates
load("results_4_27.mat");
rtab = struct2table(results);
viola_chrom = rtab(rtab.fidx == 2 & rtab.fm_method_idx == 2 & rtab.Nc == 4, :).best_chrom{1};
trumpet_chrom = rtab(rtab.fidx == 1 & rtab.fm_method_idx == 2 & rtab.Nc == 4, :).best_chrom{1};

t1 = 0:1/fs:length(x1)/fs;
[~, Wavg] = evaluate(x1, T_trumpet, modfm, trumpet_chrom, params1);
x3 = synthesize(length(x1), Wavg, modfm_synth, trumpet_chrom, params1);

t2 = 0:1/fs:length(x2)/fs;
[~, Wavg] = evaluate(x2, T_viola, modfm, viola_chrom, params2);
x4 = synthesize(length(x2), Wavg, modfm_synth, viola_chrom, params2);

%% Plot spectrograms

nfft = params.Lw*(1+params.zpf);

figure;

subplot(2,2,1);
spectrogram(x1,params.Lw,floor(params.Lw/2),nfft,fs);
colormap("jet");
xlim([0 9]);
title("Trumpet target spectrum");

subplot(2,2,2);
spectrogram(x2,params.Lw,floor(params.Lw/2),nfft,fs);
colormap("jet");
ylim([0 3.5]);
xlim([0 6.5]);
title("Viola target spectrum");

subplot(2,2,3);
spectrogram(x3,params.Lw,floor(params.Lw/2),nfft,fs);
colormap("jet");
xlim([0 9]);
title("Trumpet best candidate spectrum");

subplot(2,2,4);
spectrogram(x4,params.Lw,floor(params.Lw/2),nfft,fs);
colormap("jet");
ylim([0 3.5]);
xlim([0 6.5]);
title("Viola best candidate spectrum");

%% Write audio
audiowrite("sounds/trumpet_target.mp4",x1,fs);
audiowrite("sounds/viola_target.mp4",x2,fs);
audiowrite("sounds/trumpet_best_candidate.mp4",x3,fs);
audiowrite("sounds/viola_best_candidate.mp4",x4./max(x4),fs);
