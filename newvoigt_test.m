clear all; close all; clc;

%Define constants
m = 1.99442e-26;            %mass CIII ion
rC = 229.687;           %CIII central wavelenght, nm
c = 3e8;                    %speed of light, m/s
W = 2.87e-2;                %electron impact parameter, angstroms
centroidGuess = 229.7;
wave = 226.628:.005:232.76;
dx = 0.005;


rNf = [];
rTf = [];
rCf = [];
rAf = [];

options = optimset('MaxFunEvals',1e6,'MaxIter',1e5);
weights = 'off';

%Test temp and density
N = 1.512E16;
T = 800;

V = fitVoigtConv(W,m,c,wave, dx, N,T,rC);

dat = V./max(V) + 0.1*rand(1,length(V));

[result, fval, exitflag, output]= fminsearch (@(P) fitVoigtfittingfunc(m ,c , W, wave, dat, ...
rC, weights,P) ,[14e16 100 centroidGuess 1] , options ); %

rNf = result(1);            %density
rTf = result(2);           %background temp
rCf = result(4);            %centroid
rAf = result(3);            %amplitude

V_out = fitVoigtConv(W,m,c,wave,dx,rNf,rTf,rC);

figure(1)
plot(wave,dat)
hold on
plot(wave,V_out./max(V_out))
hold off
