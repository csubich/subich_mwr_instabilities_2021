function [modes, phases] = label_eigvecs(eigvec,Nx,Lx)
dk = 2*pi/Lx; % Smallest increment of wavenumber
kvec = dk*([0:floor(Nx/2) ceil((-Nx+1)/2):-1])'; % Wavenumber vector

if (nargout > 1 && size(eigvec,1) ~= 2*Nx)
    error(sprintf('eigvec matrix must have %d rows (has %d)',2*Nx,size(eigvec,1)))
end

f_phi = fft(eigvec(1:Nx,:),[],1);

[~,maxidx] = max(abs(f_phi),[],1);
modes = kvec(maxidx);

if (nargout > 1)
    phi_amp = f_phi(sub2ind([Nx,2*Nx],maxidx,1:2*Nx));
    f_u = fft(eigvec(Nx+(1:Nx),:),[],1);
    u_amp = f_u(sub2ind([Nx,2*Nx],maxidx,1:2*Nx));
    phases = u_amp./phi_amp;
end