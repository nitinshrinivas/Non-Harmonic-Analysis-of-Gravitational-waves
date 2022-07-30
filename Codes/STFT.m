%% STFT
clear;

n=0:0.0001:1-1/10/(3)^(1/2);
fs=4096; %Hz
f=1./((1-n).^2);
% xn=chirp(n,1,1-1/10/(3)^(1/2),300,'quadratic');

xn=sin(2*pi*f);
figure(1);
plot(xn);

stft(xn); 
 [s,f1,t1]=stft(xn);          % f1-> contains frequencies at which stft calculated
                                % t1-> mid points of all the windows
% [val,f2]=max(abs(s(:,(1:length(t1)))));
 figure(2);
% plot(t1/10000,abs((f2-64)*32));
% hold on;
plot(t1/10000,f(t1));
% ylim([0,300]);
% hold off;

% rms_error_stft=((sum(abs((f2-64)*32)-f(t1)).^2)/length(t1))^(1/2);