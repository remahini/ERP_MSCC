
function [AHC_ERP]=AAHC_ERP(ERP_Subj,Subj,k,G)

for g=1:G
    for s=1:Subj % subjects
        
        x=squeeze(ERP_Subj(:,:,s,g)); % samples x channels x subjects x groups


        [A_all,c_idx] = raahc(x',k); % AAHC chan x samples 

        
        AHC_ERP(:,s,g)=c_idx{1,1}';
        
    end
end
end