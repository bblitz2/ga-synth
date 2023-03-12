clear variables;
close all;

modfm = @(t,fm,fc,I) exp(I.*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t);
synth = @(t,fm,k,I,w) sum(w.*modfm(t,fm,k.*fm,I));

% Load the sound file
x = audioread("Trumpet.novib.ff.C6.stereo.aiff");
x = (x(:,1) + x(:,2)) / 2;  % Convert from stereo to mono
% TODO: better way of extracting f0?
% [pxx,fwelch] = pwelch(sound,fs*window_s,fs*overlap_s,10*fs,fs);
% [~,maxloc] = max(pxx);
% f0 = fwelch(maxloc);

% Hyperparameter definitions
f0 = 1046.4;  % Determined using pwelch to nearest 0.1Hz
fm = f0;  % Modulation frequency
Nc = 4;  % Number of carrier frequencies
fs = 44000;  % Sample rate
Nharm = 10;  % Number of harmonics in analysis
Lw = 0.010*fs;  % Analysis window size in samples
zpf = 4;  % Zero padding factor
fpass = 22;  % Passband frequency for envelope filter
Lx = length(x);
t = (0:Lx-1)*(1/fs);

% Chromosome: Nc pairs of (ki, I)
chromosome = [
    0   1
    1   4
    2   2
    4   1
];

% Parameter definitions and example synthesis

k = chromosome(:,1);
I = chromosome(:,2);
y = synth(t,fm,k,I,ones(Nc,Lx));

% Solve for W
T = harmonic_analysis(x,fs,f0,Lw,zpf,Nharm);
A = harmonic_analysis(y,fs,k*f0,Lw,zpf,Nharm);
Atrans = permute(A,[2 1 3]);
W = pagemldivide(pagemtimes(Atrans,A),pagemtimes(Atrans,T));
W = fillmissing(W,'nearest',3);

% Average overlaps in W
Wavg = W;
Wavg(:,:,1:2:end) = (W(:,:,1:2:end)+W(:,:,2:2:end))/2;
Wavg(:,:,2:2:end-1) = (W(:,:,2:2:end-1)+W(:,:,3:2:end))/2;

% Determine least squares error (fitness)
E = T - pagemtimes(A,Wavg);
err = sum(pagemtimes(permute(E,[2 1 3]),E));

% Sound reconstruction
w = lowpass(repelem(squeeze(Wavg),1,Lw/2)',fpass,fs)';
y = synth(t,fm,k,I,w(:,1:length(t)));
sound(y,fs,16)

% Harmonic analysis
function A = harmonic_analysis(x,fs,fcs,Lw,zpf,Nharm)
    Lx = length(x);
    nfft = Lw*zpf;
    nwindows = 2*ceil(Lx/Lw);
    windowsl = zeros(Lw,nwindows/2);
    windowsr = windowsl;
    windowsl(1:Lx) = x;
    windowsr(1:Lx-(Lw/2)) = x(Lw/2+1:end);
    windows = zeros(nfft,nwindows);
    windows(1:Lw,1:2:nwindows) = windowsl;
    windows(1:Lw,2:2:nwindows) = windowsr;
    spec = fft(windows);
    specfreqs = (0:nfft-1)*fs/nfft;
    A = zeros(Nharm,length(fcs),nwindows);
    % TODO: speed up the interpolation
    for i = 1:nwindows
        bins = spec(:,i);
        for j = 1:length(fcs)
            for h = 1:Nharm
%                 A(h,j,i) = abs(bins(round(fcs(j)*h/nfft)))
                A(h,j,i) = abs(interp1(specfreqs,bins,fcs(j)*h,'linear'));
            end
        end
    end
end
