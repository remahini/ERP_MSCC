

function [MK_ERP]=MKMS_ERP(ERP_Subj,Subj,k,G)

for g=1:G
    for s=1:Subj % subjects
        
        x=squeeze(ERP_Subj(:,:,s,g)); % samples x channels x subjects x groups
 
        [A_opt,c_idx,Res] = modkmeans(x',k); % chan x samples
        
        
        MK_ERP(:,s,g)=c_idx';
        
    end
end
end