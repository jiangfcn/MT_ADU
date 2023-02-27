function [rhoa,rhoaerr,phi,phierr]=calc_MT(Z,Zvar,f)
%
%[rhoa,rhoaerr,phi,phierr]=calc_MT(Z,Zvar,f)
%
%    [rhoa,rhoaerr,phi,phierr]=calc_MT(Z,Zvar,f) calculates the apparent
%    resifreqtivity and phase including corresponding errors for
%    frequencies f based on the impedance tensor Z. The impedance tensor
%    has to be a 2x2xN array where N corresponds to the number of
%    frequencies.

% NOTE: Please check the unit of the impedance when calling this function
% which should be SI unit and NOT be field measurement unit. Field: E/B
% [mV/km]/nT) SI: E/H [V/m]/[A/m], corresponding \Omegam

% Specifying constants
mu = 4*pi*10^-7;

% Calculating apparent resistivities and phases
nf=size(Z,3);
rhoa = zeros(2,2,nf);
phi = zeros(2,2,nf);
for ifreq = 1:nf
    rhoa(:,:,ifreq) = 1/(2*pi*f(ifreq)*mu)*abs(Z(:,:,ifreq)).^2; % using the SI unit to calculate the rhoa and phase.
    phi(:,:,ifreq) = atan2(imag(Z(:,:,ifreq)),real(Z(:,:,ifreq)))*360/(2*pi); % four quadrant.
    
    % Calculating apparent resifreqtivity and phase errors (using Ersan's
    % rules)
%     erZZ(:,:,ifreq) =
%     Zvar(:,:,ifreq)./sqrt(Z(:,:,ifreq).*conj(Z(:,:,ifreq)));
%     rhoaerr(:,:,ifreq) =
%     rhoa(:,:,ifreq).*((2*erZZ(:,:,ifreq))+(erZZ(:,:,ifreq)).^2);
%     phierr(:,:,ifreq) = rhoaerr(:,:,ifreq)./rhoa(:,:,ifreq) *100*0.29;
    
    %JF, using Chen's rules
    Z(:,:,ifreq) = Z(:,:,ifreq)*(1/(4*pi*10^-4)); % !!!!!converting the SI unit to field unit when calculate the error using Chen's rules.
    Zvar(:,:,ifreq) = Zvar(:,:,ifreq)*(1/(4*pi*10^-4));
    rhoavar(:,:,ifreq) = 0.4./f(ifreq).*rhoa(:,:,ifreq).*Zvar(:,:,ifreq);
    phivar(:,:,ifreq) = 8*(cosd(phi(:,:,ifreq)).^4).*Zvar(:,:,ifreq).*abs(Z(:,:,ifreq)).^2./(abs(Z(:,:,ifreq)+conj(Z(:,:,ifreq))).^4);
    rhoaerr(:,:,ifreq) = sqrt(rhoavar(:,:,ifreq));
    phierr(:,:,ifreq) = sqrt(phivar(:,:,ifreq)).*180/pi;
end

% Otherwise make errors zero
%rhoaerr = zeros(2,2,size(Z,3)); phierr = zeros(2,2,size(Z,3));

% Need to include error propagation