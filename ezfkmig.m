function [migRF,z] = ezfkmig(RF,dt,pitch,ee)
% RF: Processed B-Scan
% dt: time-step(seconds)
% pitch: spatial step (meters)
% ee : relative permittivity 
% migRF: reconstructed image after migration
% z : depth axis to be used for plotting (imagesc(x,z,miRF))
    
fs=1/dt;
c=3*10^8/sqrt(ee);
[nt0,nx0] = size(RF);

% Zero-padding
nt = 2^(nextpow2(nt0)+1);
nx = 2*nx0;

% Exploding Reflector Model velocity
ERMv = c/sqrt(4);

% FFT
fftRF = fftshift(fft2(RF,nt,nx));

% Linear interpolation
f = (-nt/2:nt/2-1)*fs/nt;
kx = (-nx/2:nx/2-1)/pitch/nx;
[kx,f] = meshgrid(kx,f);
fkz = ERMv*sign(f).*sqrt(kx.^2+f.^2/ERMv^2);
fftRF = interp2(kx,f,fftRF,kx,fkz,'linear',0);

% Jacobian (optional)
kz = (-nt/2:nt/2-1)'/ERMv/fs/nt;
fftRF = bsxfun(@times,fftRF,kz)./(fkz+eps);

% IFFT & Migrated RF
migRF = ifft2(ifftshift(fftRF),'symmetric');
migRF = migRF(1:nt0,1:nx0);
%
z=(1:nt0)*dt*c/2;
end







