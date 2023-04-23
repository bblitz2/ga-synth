% Harmonic analysis
function A = harmonic_analysis(x,fs,fc,Lw,zpf,Nharm)
    nfft = Lw*(1+zpf);

    % Do STFT for all windows
    spec = abs(spectrogram(x,Lw,floor(Lw/2),nfft));

    % Find frequency locations between bins for each harmonic
    h = (1:Nharm);
    startfreqs = fc.*h;
    startbin = floor(startfreqs.*nfft./fs)+1;
    binoffs = (fc.*h - (startbin-1).*fs./nfft).*nfft./fs;
    stopbin = startbin+1;

    % Interpolate harmonic weight at each frequency location
    Ac = spec(startbin,:)+(spec(stopbin,:)-spec(startbin,:)).*binoffs';
    A(:,1,:) = Ac;
end
