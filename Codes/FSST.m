%% FSST
clear;

n=0:0.0001:1-1/10/(3)^(1/2);
fs=4096; %Hz
f=1./((1-n).^2);

xn=sin(2*pi*f);
figure(1);
plot(xn);

fsst(xn);