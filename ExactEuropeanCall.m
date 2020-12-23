%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Developed by Ehsan Mohebianfar in 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptionPriceExactEuropeanCall=ExactEuropeanCall()
format long e
clc
% fprintf('Input the risk-free interest rate= \n');
%     r= input('');
% fprintf('Input the volatility of stock price= \n');
%     ro= input('');
% fprintf('Input the Strike Price of Stock Price= \n');
%     SP= input('');    
%fprintf('Input the Stock Price= \n');
%    S= input('');    
% fprintf('Input the Maturity Time= \n');
%     T= input('');
N=1024;
r=.05;
ro=.2;
SP=10;
T=.5;
Smax=30;
h=Smax/(N-1);
S=h:h:Smax-h;
d1=(log(S/SP)+(r+.5*ro^2)*T)/(ro*sqrt(T));
d2=d1-ro*sqrt(T);
N1 = 0.5*(1+erf(d1/sqrt(2)));
N2 = 0.5*(1+erf(d2/sqrt(2)));
ExactValue=S.*N1-SP.*exp(-r*T).*N2;
OptionPriceExactEuropeanCall= ExactValue