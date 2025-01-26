% The main code for calculating the pac (MI) 
% last edit date: '23-Nov-2024'
% Ashkan Farrokhi, Mohammad Aliramezani
% IPM
%% description: this code calculate the PAC on a sample LFP data 
%- Input: "Data" is a structure including "LFP", "Data.ChannelName", ...
%- and "Data.ChannelNumber" related to the main data.
%- LFP: a matrix [trial * sample] for a channel
%- the main results on whole recording data published in : ...
%- "Aliramezani, M., et al., Delta-alpha/beta coupling as a signature of ...
%- visual working memory in the prefrontal cortex. iScience, 2024. 27(8)."
%- this version prepared for the "star protocols" journal
%% load sample data:
clear;
close all hidden
clc
load LFP_Sample.mat  % load sample LFP data
Sig = Data.LFP ;
fprintf('Sample Data Loaded   Channel Name: %s \n',Data.ChannelName)
fprintf('----------------------------------------------------- \n')

%% parameter:
Fs= 500; % sampling frequency

%- set the frequency sets use for PAC calculation:
fl = 1:1:15;          % for phase component
fh = 1:1:80;         % for Amplitude component

freqs2use=1:1:100;  % freqeuncies that signal decomposed into it using wavelet
% - note: all frequencies that selected in fl and fh should be exist in "freqs2use"

% determine the number of cycle for creating the wavelet function
% num_cycles    = 7*ones(1,length(freqs2use));
num_cycles    = logspace(log10(3),log10(10),length(freqs2use));


%- set N_shuffle as the size of the surrogate distribution for PAC calculation
N_shuffle = 100;

%% remove Artifcat 
sigma=3; % use for calculating the threshold 
artif = Remove_Artifact(Sig, sigma); % calculate the polluted trials
Sig =Sig(~artif,:); % remove polluted trials

fprintf('Artifact removed     Number of Removed Trials: %d \n',sum(artif))
fprintf('----------------------------------------------------- \n')
%% Notch filter
fn =60;   % frequency of power line noise (Hz)
Sig = Notch_Filter(Sig,fn,Fs);
fprintf('Power Line Noise filtered \n')
fprintf('----------------------------------------------------- \n')
%% Normalizaton
% using zscore to normalize the data of electrode
Sig = zscore (Sig,[],'all');
fprintf('Signal Normalized \n')
fprintf('----------------------------------------------------- \n')
%% Time-frequency decompositon using Complex Wavelet
Analytic_Sig = Morlet_Wavelet(Sig, freqs2use, num_cycles, Fs) ;
fprintf('Wavelet Transform completed \n')
fprintf('----------------------------------------------------- \n')
%% Select the time period for analysis
%- the time period seletcer based on the time interval we want to calculat
%- the PAC in it. for more details reffer to the our paper published at
%- "iScience", 2024

%- interwals description:
%Inter-trial(0-500) Fixation(501-1000) S1(1001-1250) D1(1251-2000) S2(2001-2250) D2(2251-3000) target(3000-3500)

% T1 = 1001; T2 = 1250;  %%S1
% period='S1';

% T1 = 1251; T2 = 2000;  %%D1
% period='D1';

% T1 = 2001; T2 = 2250;  %%S2
% period='S2';
 
% T1 = 2251; T2 = 3000;  %%D2
% period='D2-';

% T1 = 1; T2 = 500;  %%Inter-trial
% period='Inter-trial';


%- set the interval:
T1 = 501; T2 = 1000;  %%Fixation
period='Fixation';


%- Trimed the data based on selected time interval
Analytic_Sig = Analytic_Sig(:,:,T1:T2);

fprintf('Time interval selected: %s  ,Sample from %d to %d \n',period,T1,T2)
fprintf('----------------------------------------------------- \n')  
%% Calculation of phase amplitude coupling based on Tort et al 2010
%- measure of PAC: mudulation index (MI)

fprintf('PAC calculation process: \n')
 

for f_1 = 1:length(fh)      %for amplitude component

    %find the index of current frequency fh(f_1) in decomposed signal
    f_amp = fh(f_1);
    f_amp = find(freqs2use==f_amp);

    % get the amplitude component
    Amp = squeeze(Analytic_Sig(:,f_amp,:));
    Amp = abs(Amp);

    for f_2= 1:length(fl) %for phase component
        

        Mgs = sprintf(['current analysis: \n' ,...
            'Amplitude fre. = %dHz(%d-%d) \nPhase fre. = %dHz(%d-%d) '] ,...
            fh(f_1),f_1,length(fh),fl(f_2),f_2,length(fl) );
        fprintf(Mgs)

        %find the index of current frequency fl(f_2) in decomposed signal
        f_phs = fl(f_2);
        f_phs = find(freqs2use==f_phs);

        % get the phase component
        Phs = squeeze(Analytic_Sig(:,f_phs,:));
        Phs = angle(Phs);


        %- calculate the original MI for current [f_amp,f_phs]
        mi = PACcalculator_MI(Phs,Amp);
        PAC.MI(f_1,f_2) = mi;

        %- calculate the shuffled MI for current [f_amp,f_phs]
        %- at each itteration the trials shuffled for AMP signal at
        %- AMP_shuffle. then mi_shuffled calculated using AMP_shuffle and
        %- Phs components at specified frequency
        mi_shuff = zeros(N_shuffle,1);
        for n = 1:N_shuffle
           
            Cond=1;
            while Cond==1
                ind =randi (size(Amp,1),1);
                if ind <3 || ind>size(Amp,1)-3 
                    % do nothing
                else
                    Cond=0;
                end
            end % while
            
            Amp_shuffled = [Amp(ind+1:size(Amp,1),:); Amp(1:ind,:)];
            mi_shuff (n) =  PACcalculator_MI(Phs,Amp_shuffled);

        end


        % put resutls in PAC as Output
        PAC.MI_Shuffle_mean(f_1,f_2) = mean(mi_shuff);
        PAC.MI_Normal(f_1,f_2) =  mi - mean(mi_shuff) ;
    
        fprintf(repmat('\b',1,length(Mgs)))
    end
end

fprintf('PAC calculation completed\n')
fprintf('----------------------------------------------------- \n') 

%% plot the resutls:

figure(1)
h=pcolor(fl,fh,PAC.MI);
h.EdgeColor = 'none';
colormap(jet); 
shading interp
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')    
set(gca, 'fontsize', 40, 'fontweight', 'bold');
title('Non-normalized PAC');
caxis([0,0.005])
xticks([1 5:5:15])
xticklabels([1 5:5:15])
yticks([1 5:5:100])
yticklabels([1 5:5:100])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
h=pcolor(fl,fh,PAC.MI_Shuffle_mean);
h.EdgeColor = 'none';
colormap(jet); 
shading interp
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')    
set(gca, 'fontsize', 40, 'fontweight', 'bold');
title('PAC Shuffled');
caxis([0,0.001])
xticks([1 5:5:15])
xticklabels([1 5:5:15])
yticks([1 5:5:100])
yticklabels([1 5:5:100])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
h=pcolor(fl,fh,PAC.MI_Normal);
h.EdgeColor = 'none';
colormap(jet); 
shading interp
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')    
set(gca, 'fontsize', 40, 'fontweight', 'bold');
title('Normalized PAC');
caxis([0,0.005])
xticks([1 5:5:15])
xticklabels([1 5:5:15])
yticks([1 5:5:100])
yticklabels([1 5:5:100])
%%