%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Developed by Ehsan Mohebianfar in 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptionPrice=FDEuropeanput6() % FDM for European put option
format long e
% clear all
clc
%%%%%%%%%%%%%%%%%%%%%% Inputing Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=100; % number of centers
c=1:0.01:8; % the value of constant shape parameter
M=100; % Number of time steps
r=0.05;
sigma=0.2;
K=10;   
T=0.5;
Smax=30;
h=Smax/(N-1);
S=0:h:Smax;
deltato=T/M;
to=0:deltato:T;% to = T-t     to(n)=n*deltato
%perturb=0.5*((sigma^2)/r);
%retion=c/h;
Exact= evalin('base','exact');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=zeros(1,N);
for i=2:N-1
    a(i)=.5*deltato*r*S(i);
end
b=zeros(1,N);
for i=2:N-1
    b(i)=.25*deltato*(sigma^2)*(S(i)^2);
end
U=zeros(2,M+1); % including boundary values
%%%%%%%% 
for n=1:M+1
    U(1,n)=K*exp(-r*to(n));
    U(2,n)=0;
end
%%%%%%%%
V=zeros(N-2,M+1); % option values in all time steps. Note: V(:,M+1)= desired values 
%%%%%%%% 
for i=1:N-2
    V(i,1)=max(K-S(i+1),0);
end
%%%%%%%%%%%%%%%%%%%%%%
RMSE=zeros(1,701);
MaxE=zeros(1,701);
%%%%%%%%%%%%%%%%%%%%%%
for p=1:701
A1=zeros(3,N-2); % Constant Matrix
A2=zeros(3,N-2); % Constant Matrix
w=(1/(2*h));
z=(1/(h^2));
for j=1:N-2
    A1(1,j)=b(j+1)*z-a(j+1)*w;
    A1(2,j)=-2*b(j+1)*z-(1+0.5*deltato*r);
    A1(3,j)=b(j+1)*z+a(j+1)*w;
end
for j=1:N-2
    A2(1,j)=A1(1,j);
    A2(2,j)=-2*b(j+1)*z+(1-0.5*deltato*r);
    A2(3,j)=A1(3,j);
end
A=zeros(N-2,N-2); % Constant Matrix
B=zeros(N-2,N-2); % Constant Matrix
for j=1:N-2
    A(j,j)=A1(2,j);
    B(j,j)=A2(2,j);
end
for j=1:N-3
    A(j+1,j)=A1(1,j+1);
    B(j+1,j)=A2(1,j+1);
end
for j=1:N-3
    A(j,j+1)=A1(3,j);
    B(j,j+1)=A2(3,j);
end
%%%%%%%%
RHS=zeros(N-2,M);
for n=1:M      %%%% Main Part
    RHS(:,n)=-B*V(:,n);
    RHS(1,n)=RHS(1,n)-A2(1,1)*U(1,n)-A1(1,1)*U(1,n+1);
    RHS(N-2,n)=RHS(N-2,n)-A2(3,N-2)*U(2,n)-A1(3,N-2)*U(2,n+1);
    %%%%%%%%% system solving %%%%%%%%%%%%
    V(:,n+1)=A\RHS(:,n);
end
App=V(:,M+1);
sum=0;
for i=1:N-2
    sum=sum+(Exact(i)-App(i))^2;
end
RMSE(p)=(1/(N-2))*sqrt(sum);
MaxE(p)=max(abs(Exact-App));
end
OptionPrice= MaxE