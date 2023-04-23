% Sound reconstruction
function y = synthesize(Lx, Wavg, synth, chromosome, params)
    fpass = 22;  % Passband frequency for envelope filter
    t = (0:Lx-1)*(1/params.fs);
    chromosome = reshape(chromosome,2,[])';
    k = chromosome(:,1);
    I = chromosome(:,2);
    Wfilt = lowpass(squeeze(Wavg)',fpass,params.fs/params.Lw);
    w = resample(Wfilt,Lx,size(Wfilt,1))';
    y = synth(t,params.fm,k,I,w(:,1:length(t)));
end

