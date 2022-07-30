clear;
dt=.0001;
n=0:dt:1-1/10/(3)^(1/2);
fs=4096;                                 %Hz
f=1./((1-n).^2);

xn=sin(2*pi*f);
% xn=chirp(n,1,1-1/10/(3)^(1/2),300,'quadratic');
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

N=512;               % window length
no_of_windows=floor(length(xn)/N);

% f_ans=zeros(1,no_of_windows*N);           %stores optimum values of f for given windows
f_ans=zeros(1,no_of_windows);           


win_no=1;                                   % initialising the win_no

d11=zeros(1,N); %d011;
d21=zeros(1,N); %d021;
d12=zeros(1,N); %d012;
d22=zeros(1,N); %d022;
d3 =zeros(1,N); %d03;
d31=zeros(1,N); %d031;
correction=[1,1];                               % initialising the correction term

while(win_no-1<no_of_windows)

    
xk=fft(xn((win_no-1)*N+1:win_no*N),5*N);       % finding the fft for the given window
    [max_val ,max_ind]=max(abs(xk));           % finding the freq with max component
    x=[(max_ind-1)/N/dt,0];

    A=abs(max_val);
%     for q=1:100
while(sum(correction)>1e-10)
   for l=1:N
    d11(l)=2*(xn((win_no-1)*N+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*A*sin(2*pi*x(1)*l*dt+x(2))*2*pi*l*dt;
    
    d21(l)=2*(xn((win_no-1)*N+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*A*sin(2*pi*x(1)*l*dt+x(2));
  
    d12(l)=2*(2*pi*l*dt)^2*A*((xn((win_no-1)*N+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*cos(2*pi*l*dt*x(1)+x(2))+A*sin(2*pi*x(1)*l*dt+x(2))^2);

    d22(l)=2*A*((xn((win_no-1)*N+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*cos(2*pi*l*dt*x(1)+x(2))+A*sin(2*pi*x(1)*l*dt+x(2))^2);

    d3(l)=2*A*(2*pi*l*dt)*((xn((win_no-1)*N+l)-A*cos(2*pi*l*dt*x(1)+x(2)))*cos(2*pi*l*dt*x(1)+x(2))+A*sin(2*pi*x(1)*l*dt+x(2))^2);
   
    d31(l)=-2*(xn((win_no-1)*N+l)-A*cos(2*pi*x(1)*l*dt+x(2)))*cos(2*pi*x(1)*l*dt+x(2));
   end
    d011=1/N*sum(d11);
    d021=1/N*sum(d21);
    d012=1/N*sum(d12);
    d022=1/N*sum(d22);
    d03 =1/N*sum(d3);
    d031=1/N*sum(d31);
    
    J=d021 * d022-d03^2;
    
    neu=0.2;            %weighting coeff
    
    j1=d011*d022-d03*d012;
    j2=d021*d012-d03*d011;
    
    correction=[j1,j2];
    
    x=x-neu/J*correction;
    A=A-neu*d031;
    
end
%     for k=1:N
%     f_ans((win_no-1)*N+k)=k*x(1);
%     phi((win_no-1)*N+k)=k*x(2);
%     end
    f_ans(win_no)=x(1);
    phi(win_no)=x(2);
    
    win_no=win_no+1;
    
end
figure(3);
plot(f_ans);





