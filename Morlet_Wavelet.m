% This function calculates the wavelet transform
% inputs:
% Sig : time domain of signal for transformation it should be N trials *
% time points
% 
% freqs2use: center frequencies for wavelent transform
% Fs : sampling frequency 
% output:
% Analytic_Sig : the complex value analytic signal   

function Analytic_Sig = Morlet_Wavelet(Sig, freqs2use, num_cycles, Fs)

padding_pnt = 1000; %% zero Padding

Sig2 = [Sig, zeros(size(Sig,1),padding_pnt)];
Sig2 = transpose(Sig2);
Sig2 = Sig2(:);
Sig2 = [zeros(padding_pnt,1); Sig2];

pnts=size(Sig,2)+padding_pnt;
Ntrial=size(Sig,1);

% the decompistion performed on frequency domain for simplicity because of the fourier
% property that the convlution in time is equal to multiplication at
% frequency domain

% design the wavelet function
Time          = -1:1/Fs:1;
half_wavelet  = (length(Time)-1)/2;

% determine the number of cycle for creating the wavelet function
% num_cycles    = 7*ones(1,length(freqs2use));
% num_cycles    = logspace(log10(3),log10(10),length(freqs2use));

N_wavelet     = length(Time);
N_data        = pnts*Ntrial + padding_pnt;
N_convolution = N_wavelet+N_data-1;

% data FFTs
data_fft1 = fft(Sig2,N_convolution);
data_fft1 = transpose(data_fft1);

%% Core: calculate the decomposition for each frquecy specified in freqs2use
Analytic_Sig = zeros(Ntrial,length(freqs2use),size(Sig,2));

for fi=1:length(freqs2use)

    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*Time) .* exp(-Time.^2./(2*(s^2))) ,N_convolution);

    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,N_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);

    % remove the padding
    convolution_result_fft = convolution_result_fft(padding_pnt+1:end);
    A_sig_h = reshape(convolution_result_fft,[pnts,Ntrial]);
    A_sig_h = A_sig_h(1:end-padding_pnt,:);

    Analytic_Sig(:,fi,:)=(transpose(A_sig_h));


end % end frequency loop


end