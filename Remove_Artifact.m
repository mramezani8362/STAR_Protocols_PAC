% This function gets the data and identify those trial in Sig which contain
% artifact based thresholding on amplitude fluctuations
% Input : Sig [trial * sample] , Sigma
% Sigma will use for calculating the threshold (th_it)
% the output is a vector ("artif") with 0 and 1 as elements, represents ...
% if the related trial include artifact (1) or clean (0)

function artif = Remove_Artifact(Sig, sigma)

val_max=max(Sig,[],2);
val_min=min(Sig,[],2);
val_=val_max-val_min;
th_it=median(val_)+sigma*std(val_);

val_=max(Sig,[],2)-min(Sig,[],2);
artif=zeros(size(Sig,1),1);
artif(abs(val_)>th_it)=1;
artif(abs(val_)==0)=1;

end