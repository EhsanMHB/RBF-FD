%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Developed by Ehsan Mohebianfar in 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptionPrice=RBFFDEuropeanput()
format long e
clear all
clc
r=0.05;
sigma=0.2;
X=10;
Smin=1;
Smax=30;
T=0.5;
M=6;
N=121;
epsilon=1.5;
S=Smin:(Smax-Smin)/(N-1):Smax;
deltat=T/(M-1);
t=T:-deltat:0;
L=zeros(N,N);
Alpha=zeros(N,M);
A=zeros(N,M);
V=zeros(N,M);
rondeV=zeros(N,M);
rond2V=zeros(N,M);
rondeL=zeros(N,N);
rond2L=zeros(N,N);
phi=zeros(1,N);
for i=1:N
    for j=1:N
        L(i,j)=sqrt(1+((epsilon^2)*(S(i)-S(j))^2));
    end
end
for i=1:N
    for j=1:N
        rondeL(i,j)=((epsilon^2)*(S(i)-S(j)))/sqrt(1+((epsilon^2)*(S(i)-S(j))^2));
    end
end
for i=1:N
    for j=1:N
        rond2L(i,j)=(epsilon^2)/((1+((epsilon^2)*(S(i)-S(j))^2))^(3/2));
    end
end
for f=1:N
for i=1:N
    if X>S(f) 
        A(i,1)=max(X-S(i),0)+r*S(i)*deltat+r*deltat*max(X-S(i),0);
    else
        A(i,1)=max(X-S(i),0)+r*deltat*max(X-S(i),0);
    end
end
Alpha(:,2)=L\A(:,1);
for k=2:M-1
    A(1,k)=(1+r*deltat)*X*exp(-r*(T-t(k)));
    for i=2:N
        V(i,k)=L(i,:)*Alpha(:,k);
        rondeV(i,k)=rondeL(i,:)*Alpha(:,k);
        rond2V(i,k)=rond2L(i,:)*Alpha(:,k);
        A(i,k)=V(i,k)-(.5*(sigma^2)*(S(i)^2)*deltat*rond2V(i,k))-r*S(i)*deltat*rondeV(i,k)+r*deltat*V(i,k);
        Alpha(:,k+1)=L\A(:,k);
    end
end
sum=0;
for i=1:N
    phi(i)=sqrt(1+((epsilon^2)*(S(f)-S(i))^2));
end 
solution=phi*Alpha(:,M);
d1=(log(S(f)/10)+(.05+.5*.2^2)*.5)/(.2*sqrt(.5));
d2=d1-.2*sqrt(.5);
N1 = 0.5*(1+erf(-d1/sqrt(2)));
N2 = 0.5*(1+erf(-d2/sqrt(2)));
ExactValue=10.*exp(-.05*.5).*N2-S(f).*N1;
sum=sum+abs(ExactValue-solution);
end
RalativeError=(1/N)*sum;
OptionPrice=RalativeError


