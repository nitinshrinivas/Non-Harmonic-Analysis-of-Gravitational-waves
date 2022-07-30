clear;
dt=.0001;
% n=0:dt:1-1/10/(3)^(1/2);
n=0:dt:1-1/20;
fs=4096; %Hz
f=1./((1-n).^2);

xn=sin(2*pi*f);
figure(1);
plot(xn);

% F(A, f , ϕ) = 1/N*(sum_{n=0}^{n=N-1}{x(n*dt) − A*cos(2*pi*f*n*dt + phi)}^2)

% dF/df= 2*[x(n*dt)-A*cos(2*pi*f*n*dt+phi)]*A*sin(2*pi*f*n*dt + phi)*2*pi*n*dt

% d2F/df2=
% 2*A*(2*pi*n*dt)^2*[{x(n*dt)-A*cos(2*pi*n*f*dt+phi)}*cos(2*pi*f*n*dt+ phi)+
% A*sin(2*pi*f*n*dt + phi)^2

% d2F/df dphi=
% 2*A*(2*pi*n*dt)*[{x(n*dt)-A*cos(2*pi*f*n*dt+phi)}*cos(2*pi*f*n*dt+ phi)+
% A*sin(2*pi*f*n*dt + phi)^2

% dF/dphi= 2*[x(n*dt)-A*cos(2*pi*f*n*dt+phi)]*A*sin(2*pi*f*n*dt + phi)

% d2F/dphi2=
% 2*A*[{x(n*dt)-A*cos(2*pi*n*f*dt+phi)}*cos(2*pi*f*n*dt+ phi)+
% A*sin(2*pi*f*n*dt + phi)^2

%dF/dA=-2*(xn-A*cos(2*pi*f*n*dt+phi))*cos(2*pi*f*n*dt+phi);

% x=[f phi];

% dF/df--->  d11
% dF/dphi--> d21
% dF/dA---> d31
% d2F/df2--> d12
% d2F/dphi2->d22
% d2F/df dphi-> d3

nm=1/3;
% f0=300;
f0=400;
for k=0:15
    fk(k+1)=(1/(nm+1))^k*f0;

end
tk(1)=0;
for k=2:17
    tk(k)=1/fk(17-k+1);
end

no_of_windows=16;

for k=1:16
win_length(k)=floor(tk(k+1)/dt);
end

f_ans=zeros(1,no_of_windows);           %stores optimum values of f, phi and A for given windows
% phi_ans=zeros(1,no_of_windows);
% A=zeros(1,no_of_windows);

win_no=1;

% d11=zeros(1,N); %d011;
% d21=zeros(1,N); %d021;
% d12=zeros(1,N); %d012;
% d22=zeros(1,N); %d022;
% d3 =zeros(1,N); %d03;
% d31=zeros(1,N); %d031;
% 
correction=[1,1];
const=0;                %stores numbers of samples covered
while(win_no-1<no_of_windows)
    win_len=win_length(win_no);
    if(const+win_len>length(xn))
        break;
    end
    xk=fft(xn(const+1:const+win_len),win_len);
    [max_val, max_ind]=max(abs(xk));
    x=[(max_ind)/win_len/dt, 0];
    A=max_val;
    
    d11=zeros(1,win_len); %d011;
    d21=zeros(1,win_len); %d021;
    d12=zeros(1,win_len); %d012;
    d22=zeros(1,win_len); %d022;
    d3 =zeros(1,win_len); %d03;
    d31=zeros(1,win_len); %d031;
    
%     for q=1:100
while(abs(sum(correction))>1e-10)
   for l=1:win_len
    d11(l)=2*(xn(const+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*A*sin(2*pi*x(1)*l*dt+x(2))*2*pi*l*dt;
    
    d21(l)=2*(xn(const+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*A*sin(2*pi*x(1)*l*dt+x(2));
  
    d12(l)=2*(2*pi*l*dt)^2*A*((xn(const+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*cos(2*pi*l*dt*x(1)+x(2))+A*sin(2*pi*x(1)*l*dt+x(2))^2);

    d22(l)=2*A*((xn(const+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*cos(2*pi*l*dt*x(1)+x(2))+A*sin(2*pi*x(1)*l*dt+x(2))^2);

    d3(l)=2*A*(2*pi*l*dt)*((xn(const+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*cos(2*pi*l*dt*x(1)+x(2))+A*sin(2*pi*x(1)*l*dt+x(2))^2);
   
    d31(l)=-2*(xn(const+l)-A*cos(2*pi*x(1)*l*dt+x(2)))*cos(2*pi*x(1)*l*dt+x(2));
   end
    d011=1/win_len*sum(d11);
    d021=1/win_len*sum(d21);
    d012=1/win_len*sum(d12);
    d022=1/win_len*sum(d22);
    d03=1/win_len*sum(d3);
    d031=1/win_len*sum(d31);
    
    J=d021 * d022-d03^2;
    
    neu=0.01;            %weighting coeff
    
    j1=d011*d022-d03*d012;
    j2=d021*d012-d03*d011;
    
    correction=[j1,j2];
    
    x=x-neu/J*correction;
    A=A-neu*d031;
    
end
    f_ans(win_no)=x(1);
%     phi_ans(win_no)=x(2);
    win_no=win_no+1;
    const=const+win_len;
end

figure(2);
plot(abs(f_ans));

% This plot of multiwindow NHA might not be very accurate but we did get
% good plots for NHA


