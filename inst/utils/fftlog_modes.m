function m = fftlog_modes(N)
% FFTLOG_MODES  Fourier mode ordering matching MATLAB fft convention.

    if mod(N,2) == 0
        m = [0:(N/2), (-(N/2)+1):-1].';
    else
        m = [0:((N-1)/2), (-(N-1)/2):-1].';
    end
end