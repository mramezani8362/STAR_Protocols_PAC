% This function computes phase amplitude coupling   
% Tort, A. B., Komorowski, R., Eichenbaum, H., & Kopell, N. (2010). 
% Measuring phase-amplitude coupling between neuronal oscillations of 
% different frequencies.
% Journal of neurophysiology, 104(2), 1195-1210.
%
% inputs:
% phase : phase values
%
% power: power values
%  
% output:
% MI : the value of mutual information btween pahse and power

function MI = PACcalculator_MI(phase,power)

n_hist_bins = 18;
phase_edges=linspace(-pi,pi,n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);
for i=1:n_hist_bins
    amp_by_phases(i) = ...
        mean(power(phase>phase_edges(i) & phase<phase_edges(i+1)));
end
P_r=amp_by_phases/sum(amp_by_phases);
MI=(1*sum(P_r.*log(P_r))+log(length(P_r)))/log(length(P_r));


end 