%  using of SOM (self-orgnizating map)neural network aim clustering
%  one stimuli as demo. the motivation of this demo we have used matlab
%  SOM toolbox for clustering ERP data (for one stimuli 58 channel and 500
%  time points)

function [S_ERP]=SOM_GAERP(x,k)

% x is grand average data e.g. 4800 (observation )x30 (feature)

x=x';
m=1; n=k;

% % for g=1:G
% %     for s=1:Subj % subjects

% %         x=squeeze(ERP_Subj(:,:,s,g))';

%----------------------Clustering with SOM ----------------------

net = selforgmap([m n]); % 2-dimension layer of 6 neurons arranged in an 2x3
% hexagonal grid for this example.

[net,tr] = train(net,x); % Now the network is ready to be optimized with *train*.
% % nntraintool

y = net(x);
cl_idx = (vec2ind(y)); % cluster indexes by SOM
S_ERP=cl_idx';

% --------------------------- membership plot --------------------------
figure;
imagesc(cl_idx');
title(['The Membership for all Stimulus Grand Averaged Data']);
set(gca,'fontsize',13)

close all;
% %     end % for stim
% % end
nntraintool('close');
end



