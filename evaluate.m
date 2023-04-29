% Evaluate fitness and determine envelopes
function [err, Wavg, That] = evaluate(x, T, synthmodel, chromosome, params)
    f0 = params.f0;  % Determined using pwelch to nearest 0.1Hz
    fm = params.fm;  % Modulation frequency
    fs = params.fs;  % Sample rate
    Nharm = params.Nharm;  % Number of harmonics in analysis
    Lw = params.Lw;  % Analysis window size in samples
    zpf = params.zpf;  % Zero padding factor
    Lx = length(x);
    t = (0:Lx-1)*(1/fs);
    chromosome = reshape(chromosome,2,[])';
    k = chromosome(:,1);
    I = chromosome(:,2);
    Nc = size(chromosome,1);  % Number of carrier frequencies
    
    % Solve for W using least squares
    A = zeros(Nharm,Nc,size(T,3));
    for c = 1:Nc
        xc = synthmodel(t,fm,k(c)*f0,I(c));
        A(:,c,:) = harmonic_analysis(xc,fs,f0,Lw,zpf,Nharm);
    end
    Atrans = permute(A,[2 1 3]);
    W = pagemldivide(pagemtimes(Atrans,A),pagemtimes(Atrans,T));
    W = fillmissing(W,'nearest',3);
    
    % Average overlaps in W
    Wavg = W;
    Wavg(:,:,1:2:end-1) = (W(:,:,1:2:end-1)+W(:,:,2:2:end))/2;
    Wavg(:,:,2:2:end-1) = (W(:,:,2:2:end-1)+W(:,:,3:2:end))/2;

    % Determine least squares error (fitness)
    That = pagemtimes(A,Wavg);
    E = T - That;
    err = sum(pagemtimes(permute(E,[2 1 3]),E));
end
