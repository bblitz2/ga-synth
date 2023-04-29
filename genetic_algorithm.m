% Run the genetic algorithm until convergence
function [x,fval,exitflag,output] = genetic_algorithm(x, params)

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
[x,fval,exitflag,output] = ga( ...
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
