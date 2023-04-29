clear variables;
close all;

modfm = @(t,fm,fc,I) exp(I.*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t); % Normalize? -- should be handled by w
classicfm = @(t,fm,fc,I) cos(2*pi*fc.*t + I.*sin(2*pi*fm.*t));
synth = @(t,fm,k,I,w) sum(w.*modfm(t,fm,k.*fm,I));

fs = 44100;

% Genetic algorithm configuration
options = optimoptions("ga");
options.CrossoverFcn = "crossoversinglepoint";
options.CrossoverFraction = 0.8;
options.EliteCount = 2;
options.FunctionTolerance = 1e-10;
options.MaxGenerations = 300;
options.MaxStallGenerations = 50;
options.MutationFcn = "mutationgaussian";
options.PlotFcn = "gaplotbestf"; %"gaplotbestindiv"; %"gaplotbestf";
options.PopulationSize = 100;
options.SelectionFcn = "selectiontournament";
options.StallTest = "geometricWeighted";

sound_files = dir("sounds/*.aif");

filenames = {...
    "Trumpet.novib.ff.C6.stereo.aiff", ...
    "Viola.arco.ff.sulG.C4.stereo.aif", ...
    "Oboe.ff.A4.stereo.aiff", ...
    "Horn.ff.C4.stereo.aif", ...
    "AltoSax.NoVib.ff.C4.stereo.aif", ...
    "Bass.arco.ff.sulD.C3.stereo.aif", ...
    "Bassoon.ff.C3.stereo.aif", ...
    "Cello.arco.ff.sulG.C4.stereo.aif", ...
    "EbClarinet.ff.C4.stereo.aif", ...
    "Flute.nonvib.ff.C5.stereo.aif", ...
    "TenorTrombone.ff.C4.stereo.aif", ...
    "Tuba.ff.C3.stereo.aif", ...
    "Violin.arco.ff.sulG.C5.stereo.aif", ...
};

results = [];

for fidx = 10
    x = audioread(fullfile("sounds", filenames{fidx}));
    if size(x,2) == 2
        x = (x(:,1) + x(:,2)) / 2;  % Convert from stereo to mono
    end
    f0 = median(pitch(x, fs, "Range", [50 2000], "Method", "SRH"));
    
    for Nc = [2 4 6]
        % Hyperparameter definitions
        params = [];
        params.f0 = f0;  % Determined using pwelch to nearest 0.1Hz
        params.fm = params.f0;  % Modulation frequency
        params.Nc = Nc;  % Number of carrier frequencies
        params.fs = 44100;  % Sample rate
        params.Nharm = 10;  % Number of harmonics in analysis
        params.Lw = 0.010*params.fs;  % Analysis window size in samples
        params.zpf = 4;  % Zero padding factor
        params.f0 = f0;
        T = harmonic_analysis(x,params.fs,params.f0,params.Lw,params.zpf,params.Nharm);

        fm_methods = {classicfm, modfm};
        for fm_method_idx = 1:2
            tic
            % Run genetic algorithm to find best chromosome
            [best_chrom,fval,exitflag,output] = ga( ...
                @(chrom) evaluate(x, T, fm_methods{fm_method_idx}, chrom, params), ...
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
            toc
            result = [];
            result.fidx = fidx;
            result.Nc = Nc;
            result.fm_method_idx = fm_method_idx;
            result.best_chrom = best_chrom;
            result.fval = fval;
            result.output = output;
            results = [results result];
        end
    end
end



