% Written by Reza Mahini 2019

function index=Indexing_M2(clust_idx,k,nSam)

for i=1:k
    cl_inf(i).data= find(clust_idx==i);
end
% index(1,:)=[1 c_idx(1)]; % may need this
index=[];
index(1,:)=[1 clust_idx(1)];
j1=2;

for j=1:nSam-1
    if abs(clust_idx(j)-clust_idx(j+1))>=1
        index(j1,:)=[j+1 clust_idx(j+1)]; % Start end of each clusters
        j1=j1+1;
    end
end

index(j1,:)=[j+1 clust_idx(j)];
end
