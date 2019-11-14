%  --------------------Hirerichal Clustering ERP data ------------------
% The ERP data has well structure for clustering time sammples and
% channnels 500x58 we applied Matlab clustering function on this dataset

function [H_ERP]=Hierarchical_GAERP(x,k)

% x is grand average data e.g. 4800 (observation )x30 (feature)

% % for g=1:G
% %     for s=1:Subj % subjects


%         x=squeeze(ERP_Subj(:,:,s,g));

%  clustering with diferent similarity functions---------------------------

%         cl=clusterdata(x,'linkage','complete','distance','euclidean','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','seuclidean','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','cityblock','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','minkowski','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','chebychev','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','mahalanobis','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','cosine','maxclust',k);
cl=clusterdata(x,'linkage','complete','distance','correlation','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','spearman','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','hamming','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','jaccard','maxclust',k);
H_ERP=cl;

% %     end
% % end
end

%----------------- The end of Hierarchical clustering ---------------------