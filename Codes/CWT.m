clear;

n=0:0.0001:1-1/10/(3)^(1/2);
fs=4096; %Hz
f=1./((1-n).^2);

xn=sin(2*pi*f);
figure(1);
plot(xn);
cwt(xn,fs);
% figure(20);
% [wt,f1]=cwt(xn,fs);
% 
% for z=1:length(n)
%     [max_val,max_ind(z)]=max(abs(wt(:,z)));
%     fz(z)=f1(max_ind(z));
% end
% plot(abs(fz));
% hold on;
% plot(f,'r');
% hold off;
% 
% rms_error_CWT=((sum((abs(fz)-f).^2))/length(f))^(1/2);
% 
