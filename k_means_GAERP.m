% ---------K-means Clustering for ERP data different methods--------------

function [K_ERP]=k_means_GAERP(x,k)


% x is grand average data e.g. 4800 (observation )x30 (feature)


% % for g=1:G
% %     for s=1:Subj % subjects

% %         x=squeeze(x(:,:)); % samples x channels x subjects x groups

%         [c_idx,cen]=kmeans(x,K); % normal kmeans toolbox
%         [c_idx,cen]=kmeans(x,K,'distance','cityblock');
%         [c_idx,cen]=kmeans(x,K,'distance','cosine');
[c_idx,cen]=kmeans(x,k,'distance','correlation');
%         [c_idx,cen]=kmeans(x,K,'distance','Hamming');

K_ERP=c_idx;
% %
% %     end
% % end
end




