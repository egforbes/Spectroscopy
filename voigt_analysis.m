% clear all; close all; clc;

%TITLE: Voigt Fit for ZaP-HD Data
%AUTHOR: Eleanor Forbes, University of Washington
%LAST EDIT: 24 JULY 2017
%OBJECTIVE: Fit a Voigt profile to ICCD data using a Gaussian sum to
%account for background (cold) radiation in the experiment.
%INPUTS: Shotnumber of interest.
%FUNCTIONS:

%Define constants
m = 1.99442e-26;            %mass CIII ion
center = 229.687;           %CIII central wavelenght, nm
c = 3e8;                    %speed of light, m/s
W = 2.8e-2;                %electron impact parameter, angstroms
centroidGuess = 229.7;
vdata = newdata;
xfine = min(x):0.005:max(x);

%Voigt Profile Data
rNf = [];                   %density
rTcf = [];                  %background temperature
rThf = [];                  %core plasma temperature
rCf = [];                   %centroid
rAf = [];                   %amplitude


%Use fminsearch to find parameters P(1), P(2), P(3), P(4), P(5) (density,
%background temp, core temp, centroid, amplitude) to fit data.
options = optimset('MaxFunEvals',1e6,'MaxIter',1e5);
[result fval exitflag output] = fminsearch(@(P) fitVoigtG2funct(m,c,W,x, vdata,center,P),[2e17 20 2e3 centroidGuess 1],options);

%%
rNf = result(1);            %density
rTcf = result(2);           %background temp
rThf = result(3);           %core temp
rCf = result(4);            %centroid
rAf = result(5);            %amplitude

DX = mean(diff(xfine));

voigtFit = rAf*fitvoigtG2Convolve(W,m,c,xfine,DX,rNf,rThf,rTcf,rCf);

figure()
plot(xfine,voigtFit)
% plot(x,voigtFit);
hold on
plot(x,newdata);
legend('Fit','Data')
%%
fig1 = figure();
