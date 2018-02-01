function f = fitVoigtfittingfunc (m ,c ,W ,x , vdata , centerNom , weights , P)

boxX =5;
boxY =0.1e16 ;

% wavelength range
%xx = -0.3:0.001:0.3;
xx = min(x):0.005: max(x);
%xx=xx ’+ centerNom .* ones ( length (xx) ,1);
% lambdaC = center * ones ( length (xx) ,1);

% evaluate lorentzian and gaussian as function
%of P (2) = Temperature and P (1) = density
%P (3) = centroid , P (4) = amplitude
fLw =(2e-16)^2*W^2*P(1)^2/10^2;
fGw =(sqrt(8*P(2)*log(2)/(m*c^2))*centerNom*(10^(-9))).^2;
fL = fLw./((4*(xx-P(3)).^2+fLw));
fG = exp((-4*log(2)*(xx-P(3)).^2)/fGw);

% convolve L and G to get Voigt , that has same size as L ( and G)
fdx = mean(diff(xx));
voigt = P(4)*convn(fG,fL,'same')*fdx;
Voigt = interp1(xx,voigt,x);


if strcmp(weights,'off')
    % Least squares fit , no weighing
    f =sum(abs(Voigt - vdata).^2)/length(x);
    if P(2) < 1
        f =1e7;
    end
    % if P (1) > 3.0 e16 || P (1) < 3.3 e16 || P (2) >11 || P (2) < 10
    % f=1 e7;
    % end
else % with weighing
    for j =1:length(x)
        if j >49 && j <70 % || j > length (x) -50) %10 - 31 , or 49 -70
            a(j)= abs(Voigt(j) - vdata(j))^2;
        else
            a(j) =10/(Voigt(j)+10)*abs(Voigt(j) - vdata(j))^2;
        end
    end
    f =sum(a)/length(x);
    if P(2)< 1
        f =Inf;
    end
end