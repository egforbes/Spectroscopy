function voigtFit = fitVoigtConv(W,m,c,rx,Dx,rN,rT,rC)

Lw = (2e-16)*W*rN/10;                     %Lorentzian Width
Gw = sqrt(8*log(2)*rT/(m*c^2))*rC*10^(-9);   %Gaussian Width
G = exp((-4*log(2)*(rx-rC*ones(size(rx))).^2)/(Gw^2));
L = Lw^2./((4*(rx-rC*ones(size(rx))).^2+Lw^2));
voigtFit=convn(G,L,'same')*Dx;
