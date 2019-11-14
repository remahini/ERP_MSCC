% ---------K-means Clustering for ERP data different methods--------------

function [K_ERP]=k_means_ERP(ERP_Subj,Subj,k,G)

for g=1:G
    for s=1:Subj % subjects
        
        x=squeeze(ERP_Subj(:,:,s,g)); % samples x channels x subjects x groups
        
        %         [c_idx,cen]=kmeans(x,K); % normal kmeans toolbox
        %         [c_idx,cen]=kmeans(x,K,'distance','cityblock');
        %         [c_idx,cen]=kmeans(x,K,'distance','cosine');
        [c_idx,cen]=kmeans(x,k,'distance','correlation');
        %         [c_idx,cen]=kmeans(x,K,'distance','Hamming');
        
        K_ERP(:,s,g)=c_idx;
        
    end
end
end




