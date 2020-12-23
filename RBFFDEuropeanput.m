%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Developed by Ehsan Mohebianfar in 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptionPrice=RBFFDEuropeanput() % put option
format long e
% clear all
clc
%%%%%%%%%%%%%%%%%%%%%% Inputing Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.05;
sigma=0.2;
X=10;   
Smin=1;
Smax=30;
T=0.5;
%%%%%%%%%%%%%%%%%%%%%% Defination Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=2;    % Index of time steps
M=6;
N=121;
h=(Smax-Smin)/(N-1);  
S=Smin:h:Smax;
deltato=T/(M-1);
to=0:deltato:T;% to = T-t
to(M+1)=T+deltato;
c=1; % the value of constant shape parameter
V=zeros(N,M+1);
U=zeros(N-2,M+1);
%%%%%%%%%%%%%%%%%%% Assign elements of 1'th colume of V using initial condition %%%%%%%%%%%%%%%%%%%%%%%%
for i=2:N-1
    V(i,1)=max(X-S(i),0);
    U(i-1,1)=V(i,1);
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Assign elements of 1'th and N'th rows of V using boundray conditions %%%%%%%%%%%%%%%
for j=1:M+1
    V(1,j)=X*exp(-r*to(j));
    V(N,j)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(N-2,N-2);  % coefficient matrix, A will be produced one time in the algorithm
w=(1/(2*h))*(1+((h^2)/(2*(c^2)))); % producing first differential coefficients
z=(1/(h^2))*(1+((h^2)/(c^2)));     % producing two differential coefficients
%%%%%%%%%%%%%%%%%%% Assign elements of A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(1,1)= .5*(sigma^2)*(S(2)^2)*(-2*z)-r-(2/deltato);
A(1,2)= .5*(sigma^2)*(S(2)^2)*z+r*S(2)*w;
for k=2:N-3  % k'th row of the coefficient matrix A
    A(k,k-1)=.5*(sigma^2)*(S(k+1)^2)*z+r*S(k+1)*(-w);
    A(k,k)=.5*(sigma^2)*(S(k+1)^2)*(-2*z)-r-(2/deltato);
    A(k,k+1)=.5*(sigma^2)*(S(k+1)^2)*z+r*S(k+1)*w;   
end
A(N-2,N-3)=.5*(sigma^2)*(S(N-1)^2)*z+r*S(N-1)*(-w);
A(N-2,N-2)=.5*(sigma^2)*(S(N-1)^2)*(-2*z)-r-(2/deltato);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RHS=zeros(N-2,M);
%%%%%%%%%%%%%%%%%%% Assign elements of RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:M
    RHS(1,k)=-(.5*(sigma^2)*(S(2)^2)*z+r*S(2)*(-w))*(V(1,k)+V(1,k+1))-... % k+1>M if k=M
            (.5*(sigma^2)*(S(2)^2)*(-2*z)-r+(2/deltato))*V(2,k)-...
            (.5*(sigma^2)*(S(2)^2)*z+r*S(2)*w)*V(3,k);
    for f=2:N-3
        RHS(f,k)=-(.5*(sigma^2)*(S(f+1)^2)*z+r*S(f+1)*(-w))*V(f,k)-...
                (.5*(sigma^2)*(S(f+1)^2)*(-2*z)-r+(2/deltato))*V(f+1,k)-...
                (.5*(sigma^2)*(S(f+1)^2)*z+r*S(f+1)*w)*V(f+2,k);
    end
    RHS(N-2,k)=-(.5*(sigma^2)*(S(N-1)^2)*z+r*S(N-1)*(-w))*(V(N-2,k))-...
             (.5*(sigma^2)*(S(N-1)^2)*(-2*z)-r+(2/deltato))*V(N-1,k)-...
             (.5*(sigma^2)*(S(N-1)^2)*z+r*S(N-1)*w)*(V(N,k)+V(N,k+1));    % k+1>M if k=M
    U(:,m)=A\U(:,m-1);
    for i=2:N-1
        V(i,m)=U(i-1,m);
    end
    m=m+1;
end
OptionPrice=A



