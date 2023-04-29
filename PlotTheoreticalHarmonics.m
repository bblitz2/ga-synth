clear variables;
close all;

modfm = @(t,fc,fm,I) exp(I*cos(2*pi*fm.*t)).*cos(2*pi*fc.*t);
classicfm = @(t,fc,fm,I) cos(2*pi*fc.*t + I*sin(2*pi*fm.*t));

fc = 440;
fm = 440;
I = 5;

nharm = 10;
N = (nharm-1)*2;
fs = N*fc;
t = (0:N-1)*(1/fs);
f = fs*(0:(N/2))/N;

I = 5;
figure;
subplot(2,1,1);
stem(f/fc,normspec(modfm(t,fc,fm,I)));
subplot(2,1,2);
stem(f/fc,normspec(classicfm(t,fc,fm,I)));


I = 10;
figure;
subplot(2,1,1);
stem(f/fc,normspec(modfm(t,fc,fm,I)));
subplot(2,1,2);
stem(f/fc,normspec(classicfm(t,fc,fm,I)));


function spec = normspec(x)
    spec = fft(x);
    spec = spec(1:length(x)/2+1);
    spec = spec./sum(spec);
end