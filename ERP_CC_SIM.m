% This toolbox provides an opportunity to process spatio-temporal ERP data
% utilizing a cluster analysis methodology.
 
% A novel consensus clustering and time-window selection mechanism has been
% developed for measuring best possible time-window for measuring ERPs of interest.
 
% Additionally, a repeated ANOVA measurement included to enhance the statistical analysis
% based on mean latency amplitude for the ERP of interest.
 
% We have provided a demo version for the toolbox in a simulated ERP data.
% The first version By Reza Mahini June 2017, later updated on Aug 2019.
 
% Please add MST1.0 to matlab path to use microstae analysis functions
% provided by Andreas Trier Poulsen and cite the addressed researches in below as well:

% Download from:
% https://github.com/DTUComputeCognitiveSystems/Microstate-EEGlab-toolbox/tree/master/MST1.0
 

%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
%  Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
%  Andreas Trier Poulsen, atpo@dtu.dk
%  Technical University of Denmark, Cognitive systems - February 2017
 
% Here is how to cite the current toolbox:
 
% Two manuscripts about this toolbox are under review please kindly cite
% this toolbox as:
 
% 1-R.Mahini et al., "Determination of Time-Window of Event-Related Potential
% Using Multiple-Set Consensus Clustering" , Under review in the Journal of Neuroscience Methods 
% 2- R.Mahini et al. ,"Optimal Number of Clusters by Measuring Similarity among
% Topographies for Spatio-temporal ERP Analysis", under review in Brain
% Topography journal.
% 3- Reza Mahini. (2019, July 22). Opt_NC_ERP (Version v1.0.0). Zenodo. http://doi.org/10.5281/zenodo.3345259
 
 
% We are sorry to state that, some functions are encrypted temporary,
% We will update the current version as we could publish the related works.
% We hope you kindly understand us and patiently wait for our final version update.
 
 
% Inputs:
% Information about the ERP data such as: number of groups, Number of Stimuli,
% number of subjects, number of components, number of time samples, start time (ms)
% , End time (ms), and more important experimental measurement time-window.
 
% Outputs:
% Time-windows for conditions/groups, Topographical maps and waveform plots,
% Statistical power analysis results.
 
% Updated in Aug 2019



clc;
clear;
close all;

delete K_lg_eig_M2.mat % To make sure the dataset is replaced with new in DSC

% -------------------------Input ERP (all)----------------------------

load chanlocs;
load Data_5;

chanLoc=chanlocs;
size(Data_5)

tic

% ------------------------------------------------------------------------
%  **************************** Initializing ******************************

global G St Sa Subj startEph endEph nSam k Comp ERP_Subj inDaGA_M1 ...
    twStart twEnd stimSet compSet chanLoc IndAnlz plotonoff InSim_Thr...
    minSamThr

G=1;
St=2;
Sa=300;
Subj=5;
startEph=-100;
endEph=600;
Chan=65;
nSam=St*Sa;
Comp=2;

stbtst=0;  % Stability and method selection test keyes
SelMeth=0; % Method selection activetion

stimSet={'Cond1','Cond2'}; % the conditions name
compSet={'N2','P3'}; % components

k=5; % the optimal number of clusters obtained from to Opt_NC toolbox
Maxcount=1; % maximum number of running code

% **** Inner-similarity threshold, it lets to include map with ...
% higher than this threshold to join candidate map list (Sensetivity parameter)
minSamThr=20; % min number of samples for TW
InSim_Thr=0.80;  % the innersimilarity threshold

% ------------------------------------------------------------------------
% ****************************** Data preparing  **************************

inData=DimPrep(Data_5,Chan,Sa,St,Subj,G);
% Channel x Sample x Stim x Subject x Group

[ERP_Subj,inDaGA_M1]=Data_Preparing(inData,Subj,St,Sa,G); % subjects data and grand average data modeled data

x=inDaGA_M1; % the grand average data ""samples x channels x group""


% ------------------------------------------------------------------------
% **************** Initializing for consensus clustering ****************

clmethod={'K-means','Hierachial','FCM','SOM ','DFS','MKMS', 'AAHC'};

% "1"= kmeans,"3"= FCM, "4"=SOM, "5"=DFS, "6"=MKMS, and "7"=AAHC

M_list=[1 2 6]; % the selected clustering methods
Stb_list=[1 6]; % the clustering methods for stabilization
rep=[3 3];

twStart=[201 265]; % (ms) The experimental interesting time course
twEnd=[266 357]; % (ms)

% [G,k,start_ms,end_ms,twStart,twEnd]=Init_ERP_CC();

CSPA_f_result=[];
CSPA_Cl_pow_res=[];


% Stability test (if it's necessary) --------------------------------------

if stbtst==1
    Maxrep=20;
    Stab=MthStb(inDaGA_M1,k,Maxrep,Stb_list); % access by --> Stab(m).OptRep
end

tic

stb=input('Do we need stablilization? Enter (1=yes, 0=No)?');
IndAnlz=input('Do you need individual subject analysis results (No =0, Yes=1) ? ');
plotonoff=input(' Do you need the component plot for each subject (No =0, Yes=1) ?');


% *********************   Clustering with multiple-method    *******************


MethLabs_M2=Labeling_AllCls(ERP_Subj,Subj,k,G,Sa,Maxcount,rep,stb,M_list,Stb_list);


% Selecting methods for furthemore procedure (if it's necessary)------------

if SelMeth==1
    Base=find(M_list==5);
    Nmeth=length(M_list);
    [Simlist]=SimLblBtw_MKMS(MethLabs_M2,Subj,Base,Nmeth);
    disp(Simlist);
end


for count=1:Maxcount
    
    [count]
    
    % -------------------------------------------------------------------------
    % ***********************  Consensus Clustering L1, L2 ********************
    
    inLabel=MethLabs_M2(count).data;
    
    [CSPA_H, compSubj_L1]=CC_L1_SIM(inLabel,Subj,k,G); % Consensus clustering Level 1
    
    if ~isempty(compSubj_L1)
        compSubj(count)=compSubj_L1;
    end
    close all;
    
    [CSPA_Label_H_L2]=CC_L2_SIM(x,CSPA_H,k,G); % Consensus clustering Level 2
    
    CC_idx(:,:,count)=CSPA_Label_H_L2;
    
    %%
    
    for com=1:Comp
        
        [v,w]=time_conv_ts(Sa,twStart(com),twEnd(com),startEph,endEph);
        
        % Channels' address for processing
        
        if com==1
            selChan={'P2','P6','PO4'};
            ch_loc=[49   56    62];
        else
            selChan={'CP2','CPz','Cz'};
            ch_loc=[42    58    65];
        end
        
        for g=1:G
            
            x1=squeeze(x(:,:,g));
            
            Clu_idx=CSPA_Label_H_L2(:,g);
            
            index=Indexing_M2(Clu_idx,k,nSam);
            
            
            PlotAmp_M1_SIM(x1,index,Sa,startEph,endEph,St,ch_loc,selChan,com,compSet,stimSet);
            
            
            [CSPA_f_result,comp_pow,innerCorr,winnID,InnSimWcl]=...
                Comp_detect_ERP_CC_Upd(Clu_idx,x1,chanlocs,k,Sa,St,v,w,com,stimSet,compSet,InSim_Thr,minSamThr);
            
            [selected_TW,TWs_ms,selTWs_ms,sel_innerCorr,InnSim,selPower_amp]=...
                Sel_TW_Upd(CSPA_f_result,innerCorr,v,w,St,g,Sa,startEph,endEph,winnID,InnSimWcl,comp_pow); % TWs selection algorithm
            
            
            for st=1:St
                
                compCorr=[];
                compCorr=sel_innerCorr(st).data;
                n=size(compCorr,1);
                meanRow=sum(sum(compCorr,1))-n;
                innSim(st,g)=meanRow/(n^2-n);
                
                % %                 SC=selected_TW(st,1,g);
                meantop_amp(st,:)=selPower_amp(st).data;
                
            end
            
            % Plot the component1 --------------------------------------------
            %             if plot_On
            
            for st=1:St
                
                % % %                     SC=selected_TW(st,1,g);
                % % %                     WI=winnID(st,SC,g);
                WI=selected_TW(st,1);
                
                figure('Renderer', 'painters', 'Position', [10 10 750 350])
                
                %         set(gcf,'outerposition',get(0,'screensize'));
                subplot(1,2,1);
                % % %                     topoplot(squeeze(power_amp(st,SC,:,g)),chanLoc)
                topoplot(squeeze(meantop_amp(st,:)),chanLoc)
                
                title(['Topography map, ClustNo.', int2str(WI),', ', stimSet{st},]);
                set(gca,'fontsize',12);
                colorbar;
                caxis([-2 2]);
                
                subplot(1,2,2)
                imagesc(sel_innerCorr(st).data);
                title(['Samples Correlation']);
                xlabel('Sample #');
                ylabel('Sample #');
                set(gca,'fontsize',12);
                colorbar;
                caxis([-1 1]);
                
            end
            
            %             end
            
        end
        
        selTWs(com).data=selTWs_ms;
        disp(compSet{com})
        disp(selTWs(com).data)
        
        disp(['Inner-similarity in ',compSet{com}, 'C1     C2'])
        disp(InnSim)
        
        compGroup_CC(count).comp(com).innSimm(k).data=InnSim; %  access innSim(st,g), st =stimulus, g= group
        % % %         compGroup_CC(count).comp(com).Corr(k).data=sel_innerCorr; % innerCorr
        compGroup_CC(count).comp(com).sel_TW(k).data=selected_TW;
        compGroup_CC(count).comp(com).idx(k).data=Clu_idx;
        compGroup_CC(count).comp(com).sel_TW_ms(k).data=selTWs(com).data; % ms
        compGroup_CC(count).comp(com).meantop(k).data=meantop_amp; % access meantop_amp(st,:,g)
        
        %______________ Update
        
        for g=1:G % this is necessary
            selected_TW(:,:,g)=compGroup_CC(count).comp(com).sel_TW(k).data;
        end
        
        % Obtained TW based analysis ---------------------------------------------
        
        [SPSS_tab_avg]=ERP_statTable2_100s(selected_TW,inData,ch_loc,Subj,St,G);
        
        [ranova_tbl]=ranova_ERP_100(SPSS_tab_avg);
        
        stim_pvalue(count,com,:)=ranova_tbl.pValue(3);
        chan_pvalue(count,com,:)=ranova_tbl.pValue(5);
        intStCh_pvalue(count,com,:)=ranova_tbl.pValue(7);
        
        ranova_tot.count(count).comp(com).data=ranova_tbl;
        
    end
    
    % % %            pause
end

disp('Statistical analysis results: ');
StaTab=table(stim_pvalue(1),stim_pvalue(2),chan_pvalue(1),chan_pvalue(2),intStCh_pvalue(1),...
    intStCh_pvalue(2),'VariableNames',{'Stim_N2','Stim_P3','Chan_N2','Chan_P3','intStChN2','intStChP3'})

toc

% Saving most important results

% save compGroup_CC_KMSHCMKMS  compGroup_CC;



%% FUNCTIONS in this procedure

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


%% Statistical power analysis for within factors (Stim x Chann) -----------

function [ranova_tbl]=ranova_ERP_100(SPSS_tab_avg)

Rdata=SPSS_tab_avg(:,2:7);

Group={'G1';'G1';'G1';'G1';'G1'};

varNames={'Group','St1_ch1','St1_ch2','St1_ch3','St2_ch1','St2_ch2','St2_ch3'};

tbl = table(Group,Rdata(:,1),Rdata(:,2),Rdata(:,3),Rdata(:,4),...
    Rdata(:,5),Rdata(:,6),'VariableNames',varNames);

factNames = {'Stim','Chan'};

within_R = table({'St1';'St1';'St1';'St2';'St2';'St2'},{'Ch1';'Ch2';'Ch3';'Ch1';'Ch2';'Ch3'},'VariableNames',factNames);

rm = fitrm(tbl,'St1_ch1-St2_ch3~1','WithinDesign',within_R);

[ranova_tbl] = ranova(rm, 'WithinModel','Stim*Chan');

end




