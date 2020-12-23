%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Developed by Ehsan Mohebianfar in 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptionPriceExactEuropeanPut=ExactEuropeanPut() % exact valumes for put option
format long e
tic
%clc
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
N=120;
SP=10;
Smax=30;
h=Smax/(N-1);
S=h:h:Smax-h;
r=0.05; % risk-free interest rate
ro=0.2; % volatility of stock price
 % Strike Price of Stock Price
T=.5; % maturity
    d1=(log(S/SP)+(r+.5*ro^2)*T)/(ro*sqrt(T));
    d2=d1-ro*sqrt(T);
    N1 = 0.5*(1+erf(-d1/sqrt(2)));
    N2 = 0.5*(1+erf(-d2/sqrt(2)));
    ExactValue=SP.*exp(-r*T).*N2-S.*N1;
    OptionPriceExactEuropeanPut= ExactValue
toc