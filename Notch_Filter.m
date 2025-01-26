% This function gets the data and removes the power line noise
% Input : Sig [trial * sample] , fn( frequency of power line noise), 
% Fs (sampling frequency in Hz)
% the output is the clearn Signal ("Sig")

function sig_ = Notch_Filter(Sig,fn,Fs)

HD = designfilt('bandstopiir','FilterOrder',20, ...
    'HalfPowerFrequency1',fn-2,'HalfPowerFrequency2',fn+2, ...
    'SampleRate',Fs);
       
% vectorize the signal:
sig_ = transpose(Sig);
sig_ = sig_(:);

%- zero_padiing  : cause of edge effect problem during filtering
sig_ = [zeros(500,1); sig_ ; zeros(500,1)];

sig_ = filtfilt(HD,sig_);
sig_=sig_(501:end-500); % remove padding

%- reshaping the filtered signal back to trial * sample
sig_ = reshape(sig_,[size(Sig,2),size(Sig,1)]);
sig_ = transpose(sig_);

end