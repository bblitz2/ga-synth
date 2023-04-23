clear variables;
close all;

modfm = @(t,fm,fc,I) exp(I.*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t); % Normalize? -- should be handled by w
classicfm = @(t,fm,fc,I) cos(2*pi*fc + I.*sin(2*pi*fm.*t));
synth = @(t,fm,k,I,w) sum(w.*modfm(t,fm,k.*fm,I));

% Load the sound file
% x = audioread("Oboe.ff.A4.stereo.aiff");
x = audioread("Trumpet.novib.ff.C6.stereo.aiff");

x = (x(:,1) + x(:,2)) / 2;  % Convert from stereo to mono
% TODO: better way of extracting f0?
[pxx,fwelch] = pwelch(x,44100*0.01,220,10*44100,44100);
[~,maxloc] = max(pxx);
figure; plot(pxx(1:20000));

% Hyperparameter definitions
params = [];
% Trumpet: 1046.5
% Oboe: 439.5
params.f0 = 1046.5;  % Determined using pwelch to nearest 0.1Hz
params.fm = params.f0;  % Modulation frequency
params.Nc = 4;  % Number of carrier frequencies
params.fs = 44100;  % Sample rate
params.Nharm = 10;  % Number of harmonics in analysis
params.Lw = 0.010*params.fs;  % Analysis window size in samples
params.zpf = 4;  % Zero padding factor

% Harmonic analysis of target spectrum
T = harmonic_analysis(x,params.fs,params.f0,params.Lw,params.zpf,params.Nharm);

% Genetic algorithm configuration
options = optimoptions("ga");
options.CrossoverFcn = "crossoversinglepoint";
options.CrossoverFraction = 0.8;
options.EliteCount = 2;
options.FunctionTolerance = 10e-10;
options.MaxGenerations = 300;
options.MaxStallGenerations = 50;
options.MutationFcn = "mutationgaussian";
options.PlotFcn = "gaplotbestf"; %"gaplotbestindiv"; %"gaplotbestf";
options.PopulationSize = 100;
options.SelectionFcn = "selectiontournament";
options.StallTest = "geometricWeighted";

% Run genetic algorithm to find best chromosome
best_chrom = ga( ...
    @(chrom) evaluate(x, T, modfm, chrom, params), ...
    2*params.Nc, ...
    [], ...
    [], ...
    [], ...
    [], ...
    zeros(1,2*params.Nc), ...
    [10 20 10 20 10 20 10 20 10 20 10 20], ...
    [], ...
    1:2:2*params.Nc, ...
    options ...
);

[err, Wavg, That] = evaluate(x, T, modfm, best_chrom, params);
y = synthesize(length(x), Wavg, synth, best_chrom, params);
sound(y,params.fs,16);

% Plot reconstruction
figure; imagesc(squeeze((T - That).^2)); colorbar;
Ty = harmonic_analysis(y,params.fs,params.f0,params.Lw,params.zpf,params.Nharm);
figure; imagesc(squeeze(T)); colorbar;
figure; imagesc(squeeze(Ty)); colorbar;
figure; plot(y);



