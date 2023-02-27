function [Z,Zerr]=calc_Z(rhoa,rhoaerr,phi,phierr,f)
%
%[Z,Zvar]=calc_Z(rhoa,rhoaerr,phi,phierr,f)
%
%    [Z,Zvar]=calc_Z(rhoa,rhoaerr,phi,phierr,f)
%    calculates the impedance tensor and variance 
%    for frequencies f based on the apparent
%    resistivity rhoa and the phase phi.
%    The apparent resistivity and phase have to be
%    2x2xN arrays where N corresponds to the number
%    of frequencies.
%   JF_20220324: Unit in the outputed Z and Zvar are [mV/km]/nT

% Specifying constants
mu = 4*pi*10^-7; 

% Calculating the impedance tensor
Z = zeros(2,2,size(rhoa,3));
Zerr = zeros(2,2,size(rhoaerr,3));
for j = 1:size(Z,3)
    Z(:,:,j) = sign(phi(:,:,j)).*(2*pi*f(j)*mu*rhoa(:,:,j)./...
        (1+tand(phi(:,:,j)).^2)).^(1/2)+...
        1i*sign(phi(:,:,j)).*(2*pi*f(j)*mu*rhoa(:,:,j)./...
        (1+cotd(phi(:,:,j)).^2)).^(1/2);
%     Zvar(:,:,j) = sqrt(2*pi*f(j)*mu*(rhoa(:,:,j)+rhoaerr(:,:,j)))-abs(Z(:,:,j));
    % convert the unit of impedance to [mV/km]/[nT]
    Z(:,:,j) = Z(:,:,j)*1/(4*pi)*10^4;
%     Zvar(:,:,j) = Zvar(:,:,j)*1/(4*pi)*10^4;
    % JF in 2018-03-16, Calculating the impedance tensor variance(error)
%     Zerr(:,:,j) = sqrt((5/4)*f(j)*rhoaerr(:,:,j).^2./rhoa(:,:,j) + 5*f(j)*rhoa(:,:,j).*(phierr(:,:,j)*pi/180).^2); %[mV/km]/[nT]  field unit
%     Zerr(:,:,j) = sqrt(rhoaerr(:,:,j).^2./(0.4*(1/f(j)).*rhoa(:,:,j)));  %[mV/km]/[nT] field unit
    Zerr(:,:,j) = sqrt((abs(Z(:,:,j)+conj(Z(:,:,j))).^4).*(phierr(:,:,j)*pi/180).^2./(8*abs(Z(:,:,j)).^2.*cosd(phi(:,:,j)).^4));  %[mV/km]/[nT] field unit
end
