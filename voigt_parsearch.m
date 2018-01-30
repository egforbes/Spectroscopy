clear all; close all; clc;

%Search the plasma parameter space to determine the error of the new voigt
%fitting function

%Set parameters
m = 1.99e-26;  %mass CIII kg
c = 3e8;       %speed of light m/s
wave = 226.628:.005:232.76;  %wavelength range
W = 2.8E-2;        %impact parameter
dx = 0.005;
rC = 229.687; %center wavelength


%Parameter search
N = 1E16:1E17:1E18;
T = 10:100:2000;



%Initial guess for other function
par0 =[229.7; 1; 0.3; 0.003];
dat(:,1) = wave;
%%
for jdens = 1:length(N)
    for jtemp = 1:length(T)
    V = fitVoigtConv(W,m,c,wave,dx,N(jdens),T(jtemp),rC);
    dat(:,2) = V./max(V);
    [parmin,resnorm,res,exitflag] = fit2voigt(dat,par0);
    fit=voigt(dat(:,1),parmin);
    results(:,jtemp,jdens) = parmin;
%     figure(1)
%     subplot(2,1,1)
%     plot(dat(:,1),dat(:,2),'b',dat(:,1),fit,'r');
%     legend('data','fit')
%     subplot(2,1,2);
%     plot(dat(:,1),res,'k'); 
%     legend('residual','location','north');
%     pause(0.1)
    end
end

%%
%Compute the error in temperature and density
fGw(:,:) = results(3,:,:)*2;    %computed gaussian width
fLw(:,:) = results(4,:,:)*2;    %computed lorentzian width

fT = (fGw./(rC*10^-9)).^2.*(m*c^2)./(8*log(2));
fN = fLw.*10./(W*2e-16);
% 
% figure(1)
% pcolor(fT), shading interp, colorbar
% figure(2)
% pcolor(fN), shading interp, colorbar


for jder = 1:length(N)
    for jter = 1:length(T)
    lsqer(jter,jder) = abs((T(jter) - fT(jter,jder))*(N(jder) - fN(jter,jder)))/(T(jter)*N(jder));
    end
end

contour(N,T,lsqer,500), colorbar
% xlim([3E17 9E17])
% caxis([0.1 0.3])
set(gca,'Fontsize',14)
xlabel('Density (cm^-3)')
xlabel('Density (cm^-^3)')
ylabel('Temperature (eV)')
title('Least Squares Error for Noise in Voigt Fit')

%%
%Now do the same thing but with noise

for jdens = 1:length(N)
    for jtemp = 1:length(T)
    V = fitVoigtConv(W,m,c,wave,dx,N(jdens),T(jtemp),rC);
    dat(:,2) = V./max(V)+0.1.*rand(1,length(V));
    [parmin,resnorm,res,exitflag] = fit2voigt(dat,par0);
    fit=voigt(dat(:,1),parmin);
    results(:,jtemp,jdens) = parmin;
%     figure(1)
%     subplot(2,1,1)
%     plot(dat(:,1),dat(:,2),'b',dat(:,1),fit,'r');
%     legend('data','fit')
%     subplot(2,1,2);
%     plot(dat(:,1),res,'k'); 
%     legend('residual','location','north');
%     pause(0.1)
    end
end
