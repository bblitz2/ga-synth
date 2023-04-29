clear variables;
close all;

modfm = @(t,fm,fc,I) exp(I.*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t); % Normalize? -- should be handled by w
classicfm = @(t,fm,fc,I) cos(2*pi*fc.*t + I.*sin(2*pi*fm.*t));
modfm_synth = @(t,fm,k,I,w) sum(w.*modfm(t,fm,k.*fm,I));
classicfm_synth = @(t,fm,k,I,w) sum(w.*classicfm(t,fm,k.*fm,I));

fs = 44100;

% Hyperparameter definitions
params = [];
params.fs = 44100;  % Sample rate
params.Nharm = 10;  % Number of harmonics in analysis
params.Lw = 0.010*params.fs;  % Analysis window size in samples
params.zpf = 4;  % Zero padding factor

%% Plot results
load("results_4_27.mat");
rtab = struct2table(results);
viola_chroms = rtab(rtab.fidx == 2 & rtab.fm_method_idx == 2, :).best_chrom;
trumpet_chroms = rtab(rtab.fidx == 1 & rtab.fm_method_idx == 2, :).best_chrom;
chrom_lists = [viola_chroms trumpet_chroms];
instruments = ["Viola" "Trumpet"];
files = [
    fullfile("sounds", "Viola.arco.ff.sulG.C4.stereo.aif")
    fullfile("sounds", "Trumpet.novib.ff.C6.stereo.aiff")
];
Nc = [2 4 6];
harmonics = [1 4 7];
carrier_colors = ['#176e2e' '#ff0000' '#000000'];
instr_xlim_max = [700 1000];

figure;
colororder(["#0000ff" "#176e2e" "#ff0000" "#000000"])
for instr_idx = 1:2
    x = audioread(files(instr_idx));
    x = (x(:,1) + x(:,2)) / 2;
    params.f0 = median(pitch(x, fs, "Range", [50 2000], "Method", "SRH"));
    params.fm = params.f0;
    T = harmonic_analysis(x,params.fs,params.f0,params.Lw,params.zpf,params.Nharm);
    for harm_idx = 1:3
        subplot(2,3,(instr_idx-1)*3+harm_idx);
        y = abs(squeeze(T(harmonics(harm_idx)+1,:,:)));
        plot(y);
        xlim([0 instr_xlim_max(instr_idx)]);
        hold on;
        xlabel("Analysis Window");
        ylabel("Magnitude");
        title(['Harmonic' num2str(harmonics(harm_idx))]);
        for nc_idx = 1:3
            chrom = chrom_lists{nc_idx,instr_idx};
            params.Nc = Nc(nc_idx);
            [err, Wavg, T_hat] = evaluate(x, T, modfm, chrom, params);
            y = abs(squeeze(T_hat(harmonics(harm_idx)+1,:,:)));
            plot(y);
            xlim([0 instr_xlim_max(instr_idx)]);
            hold on;
        end
    end
end
legend(["target", "2 carriers", "4 carriers", "6 carriers"]);


