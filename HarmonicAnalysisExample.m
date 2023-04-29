clear variables;
close all;

modfm = @(t,fc,fm,I) exp(I*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t);
classicfm = @(t,fc,fm,I) cos(2*pi*fc.*t + I*sin(2*pi*fm.*t));
trapezoid = @(x,a,b,c,d) max(min(min((x-a)./(b-a),1),(d-x)./(d-c)),0);

fc = 440;
fm = 440;
I = 2;

nharm = 10;
fs = 44100;
dur_sec = 2.0;

t = 0:1/fs:dur_sec-1/fs;
N = length(t);


I = 5;
figure;

c1 = classicfm(t,fc,fm,I);
m1 = modfm(t,fc,fm,I).*exp(-I);
[cpxx, ~] = periodogram(c1,hann(N),441*5,fs);
[mpxx, f] = periodogram(m1,hann(N),441*5,fs);

p1 = plot(f, 10*log10(cpxx), "o:");
p1.LineWidth = 1;
p1.MarkerSize = 5;
hold on;
p2 = plot(f, 10*log10(mpxx), ".:");
p2.LineWidth = 2;
p2.MarkerSize = 20;
xlabel("Frequency (kHz)");
ylabel("Power/frequency (dB/Hz)");
title("Power Spectral Density, f_m=440Hz, k=1, I=5");
legend(["ClassicFM" "ModFM"])
xlim([0 10000]);


%%
% I = 10;
% figure;
% periodogram(m2,hann(N),441*5,fs);
% title("Power Spectral Density of m_2");
% xlim([0 10]);

I = 2;
w1 = trapezoid(t,0,0.1,0.3,1.7);
m1 = modfm(t,1*fm,fm,I).*exp(-I);
figure; spectrogram(m1,441,220,441*5,fs);

I = 8;
w2 = trapezoid(t,0.3,1.7,1.9,2.0);
m2 = modfm(t,5*fm,fm,I).*exp(-I);
figure; spectrogram(m2,441,220,441*5,fs);

%%
figure;
spectrogram(sin(2*pi*fm.*t),441,220,441*5,fs);
xlim([0 8]);
title("Spectrogram of Modulator Tone: sin(2*pi*f_m.*t), f_m=440Hz");
audiowrite("demo1.mp4",sin(2*pi*fm.*t),fs);


figure;
spectrogram(m1,441,220,441*5,fs);
xlim([0 8]);
title(["Spectrogram of x_1(t) = e^{I_1 cos(2\pi f_m t)} cos(2\pi k_1 f_m t),",
"f_m=440Hz, k_1=1, I_1=2"]);
audiowrite("demo2.mp4",m1,fs);

figure;
spectrogram(m2,441,220,441*5,fs);
xlim([0 8]);
title(["Spectrogram of x_2(t) = e^{I_2 cos(2\pi f_m t)} cos(2\pi k_2 f_m t),",
"f_m=440Hz, k_2=5, I_2=8"]);
audiowrite("demo3.mp4",m2,fs);

figure; plot(t, [w1' w2']);
legend(["w_1(t)", "w_2(t)"]);
xlabel("Time (s)");
ylim([0, 1.5]);
ylabel("Envelope");

m1m2 = m1.*w1 + m2.*w2;
figure; spectrogram(m1m2,441,220,441*5,fs);
title("Spectrogram of Synthesized Sound: s(t) = w_1(t)x_1(t)+w_2(t)x_2(t)");
xlim([0 8]);
sound(m1m2,fs)
audiowrite("demo4.mp4",m1m2,fs);


% I = 5;
% figure;
% subplot(2,1,1);
% stem(f/fc,normspec(modfm(t,fc,fm,I)));
% subplot(2,1,2);
% stem(f/fc,normspec(classicfm(t,fc,fm,I)));
% 
% 
% I = 10;
% figure;
% subplot(2,1,1);
% stem(f/fc,normspec(modfm(t,fc,fm,I)));
% subplot(2,1,2);
% stem(f/fc,normspec(classicfm(t,fc,fm,I)));

% fs = 44000;
% duration = 2;
% N = fs*duration;
% t = (0:N-1)*(1/fs);
% f = fs*(0:(N/2))/N;

% example = modfm(t,440,440,20)/exp(20)/3 ...
%     + modfm(t,880,440,3)/exp(3)/3 ...
%     + modfm(t,1320,440,4)/exp(4)/3;
% figure; plot(example);
% figure; stem(f/fm,normspec(example));
% sound(example,fs);


function spec = normspec(x)
    spec = fft(x);
    spec = spec(1:length(x)/2+1);
    spec = spec./sum(spec);
end