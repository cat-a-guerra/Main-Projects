%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% BRAIN RESPONSE FUNCTIONS AND NEUROVASCULAR COUPLING IN TYPE 2 DIABETES:
% INSIGHTS FROM FMRI
% 
%                       Catarina Guerra | 2015240209
%                               December 2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all
clc

load('ROI_regions.mat');
load('HRF_parameters.mat');
load('headers.mat');
load('covariates_HRFdata.mat');
load('diamarker_summary.mat');


%% CONSTANTS:

n_subjects = 141;
n_T2DM = 64;
n_CNT = 77;
n_rois = 22;
n_parameters = 11;


%% STATISTICS:

% In this section, we estimate the statistics of the HRF parameters of
% interest from every subject in each group and condition and in every ROI. 
% First, we apply a Shapiro-Wilk test in order to assess if the data has a 
% normal distribution per group and condition. Then, we assess the 
% significant difference for the same HRF parameters between groups. To do 
% so, we apply a Two-sample T-test or a Wilcoxon ranksum test, depending on 
% whether the data regarding the same HRF parameter in both groups of the 
% same condition has a normal distribution or not, respectively.
% Afterwards, we apply a one-way ANOVA test in order to see if the age is a
% factor which has some sort of influence on the HRF parameters in each 
% group.


% HRF parameters data from each group and condition in every ROI
CNT_Thr = zeros(n_CNT,n_rois*length(HRF_parameters));
T2DM_Thr = zeros(n_T2DM,n_rois*length(HRF_parameters));
CNT_Sub = zeros(n_CNT,n_rois*length(HRF_parameters));
T2DM_Sub = zeros(n_T2DM,n_rois*length(HRF_parameters));

% Shapiro-Wilk test results of the HRF parameters from each group and 
% condition in every ROI 
sw_CNT_Thr = zeros(3,n_rois*length(HRF_parameters));
sw_T2DM_Thr = zeros(3,n_rois*length(HRF_parameters));
sw_CNT_Sub = zeros(3,n_rois*length(HRF_parameters));
sw_T2DM_Sub = zeros(3,n_rois*length(HRF_parameters));


% Column numbers for the HRF parameters with normal distribution from each 
% group and condition in every ROI
avg_std_indexes_CNT_Thr = [];
avg_std_indexes_T2DM_Thr = [];
avg_std_indexes_CNT_Sub = [];
avg_std_indexes_T2DM_Sub = [];

% Column numbers for the HRF parameters with non-normal distribution from 
% each group and condition in every ROI
median_iqr_indexes_CNT_Thr = [];
median_iqr_indexes_T2DM_Thr = [];
median_iqr_indexes_CNT_Sub = [];
median_iqr_indexes_T2DM_Sub = [];


% Average + standard-deviation of the normal HRF parameters // Median + 
% interquartile range of the non-normal HRF parameters from each group and
% condition in every ROI
avg_parameters_CNT_Thr = [];
avg_parameters_T2DM_Thr = [];
std_parameters_CNT_Thr = [];
std_parameters_T2DM_Thr = [];

median_parameters_CNT_Thr = [];
median_parameters_T2DM_Thr = [];
iqr_parameters_CNT_Thr = [];
iqr_parameters_T2DM_Thr = [];

avg_parameters_CNT_Sub = [];
avg_parameters_T2DM_Sub = [];
std_parameters_CNT_Sub = [];
std_parameters_T2DM_Sub = [];

median_parameters_CNT_Sub = [];
median_parameters_T2DM_Sub = [];
iqr_parameters_CNT_Sub = [];
iqr_parameters_T2DM_Sub = [];


% Column numbers for the HRF parameters with normal // non-normal 
% distribution from different groups but same condition in every ROI
ttest_indexes_thr = [];
ttest_indexes_sub = [];

wilcoxon_indexes_thr = [];
wilcoxon_indexes_sub = [];


% T-test and Wilcoxon test results of the HRF parameters from each 
% condition in every ROI
ttest_thr = [];
ttest_sub = [];

wilcoxon_thr = [];
wilcoxon_sub = [];


for w=2:width(covariates_HRFdata)
    
    
    % --------------------- THRESHOLD CONDITION -------------------------
        
    
    CNT_Thr(1:n_CNT,w-1) = table2array(covariates_HRFdata(1:2:154,w));                                    
    
    [h,p,stats] = swtest(table2array(covariates_HRFdata(1:2:154,w)), 0.05, true);                                             % estimates h and p of the Shapiro-Wilk test for each HRF parameter of the CNT subjects in every ROI
    sw_CNT_Thr(:,w-1) = [h;p;stats];
    
    if h==0                                                                                                             % in case of a normal distribution --> group average of the HRF parameters
        avg_std_indexes_CNT_Thr = [avg_std_indexes_CNT_Thr (w-1)];
        avg_parameters_CNT_Thr = [avg_parameters_CNT_Thr mean(table2array(covariates_HRFdata(1:2:154,w)),1)];
        std_parameters_CNT_Thr = [std_parameters_CNT_Thr std(table2array(covariates_HRFdata(1:2:154,w)),0,1)];
    else                                                                                                                % in case of a non-normal distribution --> group median of the HRF parameters
        median_iqr_indexes_CNT_Thr = [median_iqr_indexes_CNT_Thr (w-1)];
        median_parameters_CNT_Thr = [median_parameters_CNT_Thr median(table2array(covariates_HRFdata(1:2:154,w)),1)];
        iqr_parameters_CNT_Thr = [iqr_parameters_CNT_Thr iqr(table2array(covariates_HRFdata(1:2:154,w)),1)];
    end
            

    T2DM_Thr(1:n_T2DM,w-1) = table2array(covariates_HRFdata(155:2:end,w));
    
    [h,p,stats] = swtest(table2array(covariates_HRFdata(155:2:end,w)), 0.05, true);                                           % estimates h and p of the Shapiro-Wilk test for each HRF parameter of the T2DM subjects in every ROI                               
    sw_T2DM_Thr(:,w-1) = [h;p;stats];
    
    if h==0
        avg_std_indexes_T2DM_Thr = [avg_std_indexes_T2DM_Thr (w-1)];
        avg_parameters_T2DM_Thr = [avg_parameters_T2DM_Thr mean(table2array(covariates_HRFdata(155:2:end,w)),1)];
        std_parameters_T2DM_Thr = [std_parameters_T2DM_Thr std(table2array(covariates_HRFdata(155:2:end,w)),0,1)];
    else
        median_iqr_indexes_T2DM_Thr = [median_iqr_indexes_T2DM_Thr (w-1)];
        median_parameters_T2DM_Thr = [median_parameters_T2DM_Thr median(table2array(covariates_HRFdata(155:2:end,w)),1)];
        iqr_parameters_T2DM_Thr = [iqr_parameters_T2DM_Thr iqr(table2array(covariates_HRFdata(155:2:end,w)),1)];
    end
    
    
    % Testing significant differences for the same HRF parameters between 
    % groups
    if sw_CNT_Thr(1,w-1)==0 && sw_T2DM_Thr(1,w-1)==0                                                                    % in case of a normal distribution in the same HRF parameter from different groups --> Two sample T-test
        ttest_indexes_thr = [ttest_indexes_thr (w-1)];
        [h,p,~,stats] = ttest2(CNT_Thr(:,w-1),T2DM_Thr(:,w-1));
        ttest_thr = [ttest_thr [h;p;stats.tstat]];
    else                                                                                                                % in case of a non-normal distribution in the same HRF parameter from different groups --> Wilcoxon ranksum test
        wilcoxon_indexes_thr = [wilcoxon_indexes_thr (w-1)];
        [p,h,stats] = ranksum(CNT_Thr(:,w-1),T2DM_Thr(:,w-1));
        wilcoxon_thr = [wilcoxon_thr [h;p;stats.zval]];
    end
        
    
    % --------------------- SUBMAXIMUM CONDITION ------------------------
    
    
    CNT_Sub(1:n_CNT,w-1) = table2array(covariates_HRFdata(2:2:154,w));                                    
    
    [h,p,stats] = swtest(table2array(covariates_HRFdata(2:2:154,w)), 0.05, true);                                             
    sw_CNT_Sub(:,w-1) = [h;p;stats];
    
    if h==0                                                                                                             
        avg_std_indexes_CNT_Sub = [avg_std_indexes_CNT_Sub (w-1)];
        avg_parameters_CNT_Sub = [avg_parameters_CNT_Sub mean(table2array(covariates_HRFdata(2:2:154,w)),1)];
        std_parameters_CNT_Sub = [std_parameters_CNT_Sub std(table2array(covariates_HRFdata(2:2:154,w)),0,1)];
    else                                                                                                                
        median_iqr_indexes_CNT_Sub = [median_iqr_indexes_CNT_Sub (w-1)];
        median_parameters_CNT_Sub = [median_parameters_CNT_Sub median(table2array(covariates_HRFdata(2:2:154,w)),1)];
        iqr_parameters_CNT_Sub = [iqr_parameters_CNT_Sub iqr(table2array(covariates_HRFdata(2:2:154,w)),1)];
    end  
    
    
    T2DM_Sub(1:n_T2DM,w-1) = table2array(covariates_HRFdata(156:2:end,w));

    [h,p,stats] = swtest(table2array(covariates_HRFdata(156:2:end,w)), 0.05, true);                                           
    sw_T2DM_Sub(:,w-1) = [h;p;stats];
    
    if h==0
        avg_std_indexes_T2DM_Sub = [avg_std_indexes_T2DM_Sub (w-1)];
        avg_parameters_T2DM_Sub = [avg_parameters_T2DM_Sub mean(table2array(covariates_HRFdata(156:2:end,w)),1)];
        std_parameters_T2DM_Sub = [std_parameters_T2DM_Sub std(table2array(covariates_HRFdata(156:2:end,w)),0,1)];
    else
        median_iqr_indexes_T2DM_Sub = [median_iqr_indexes_T2DM_Sub (w-1)];
        median_parameters_T2DM_Sub = [median_parameters_T2DM_Sub median(table2array(covariates_HRFdata(156:2:end,w)),1)];
        iqr_parameters_T2DM_Sub = [iqr_parameters_T2DM_Sub iqr(table2array(covariates_HRFdata(156:2:end,w)),1)];
    end
    
    
    % Testing significant differences for the same HRF parameters between 
    % groups
    if sw_CNT_Sub(1,w-1)==0 && sw_T2DM_Sub(1,w-1)==0                        
        ttest_indexes_sub = [ttest_indexes_sub (w-1)];
        [h,p,~,stats] = ttest2(CNT_Sub(:,w-1),T2DM_Sub(:,w-1));
        ttest_sub = [ttest_sub [h;p;stats.tstat]];
    else                                                                    
        wilcoxon_indexes_sub = [wilcoxon_indexes_sub (w-1)];
        [p,h,stats] = ranksum(CNT_Sub(:,w-1),T2DM_Sub(:,w-1));
        wilcoxon_sub = [wilcoxon_sub [h;p;stats.zval]];
    end
end



% **************************** One-way ANOVA *****************************


anova_thr = zeros(2,n_rois*n_parameters);
anova_sub = zeros(2,n_rois*n_parameters);

anova_headers = cell(1,n_rois*n_parameters);
anova_group = cell(n_subjects,1);


% Deletes original database information regarding subject no.66, who not 
% included in the image analysis due to problems in co-registration
diamarker_summary(66,:) = [];

% Assigns each subject's age at the time of the scan according to its group
% and then assembles all ages in one array
ages_CNT = table2array(diamarker_summary(n_T2DM+1:end,3));
ages_T2DM = table2array(diamarker_summary(1:n_T2DM,3));
total_ages = [ages_CNT; ages_T2DM];


% Assigns the group of each subject
for s=1:n_subjects
    
    if s<78
        anova_group{s}='CNT';
    else
        anova_group{s}='T2DM';
    end
end


for r=1:n_rois
    for par=1:n_parameters

        % Automatically sets the table headers
        anova_headers{par+n_parameters*(r-1)} = strcat(strrep(ROI_regions{r},' ','_'), strrep(HRF_parameters{par},' ','_'));                  % Replaces spaces by '_'
        
        [~,stats]=aoctool(table2array(covariates_HRFdata(1:2:end,1+par+n_parameters*(r-1))),total_ages,anova_group,0.05,strip(HRF_parameters{1}),'Age','Group','off');
        anova_thr(:,par+n_parameters*(r-1))=[stats{2,6}; stats{2,5}];
        
        [~,stats]=aoctool(table2array(covariates_HRFdata(2:2:end,1+par+n_parameters*(r-1))),total_ages,anova_group,0.05,strip(HRF_parameters{1}),'Age','Group','off');
        anova_sub(:,par+n_parameters*(r-1))=[stats{2,6}; stats{2,5}];
    end
    
end



%% ADJUSTED P-VALUES

% In this section, we estimate the adjusted p-values according to the 
% Benjamini & Hochberg (1995) approach for the Two-sample T-test, Wilcoxon
% ranksum test and One-way ANOVA test, using a False Discovery Rate (FDR) 
% of 0.10. 

% The results are stored in the following order:
%       1st round of adjustments: 
%                   1. Two-sample T-test Thr
%                   2. Two-sample T-test Sub
%                   3. Wilcoxon ranksum test Thr
%                   4. Wilcoxon ranksum test Sub

%       2nd round of adjustments:
%                   1. One-way ANOVA Thr
%                   2. One-way ANOVA Sub


if (isempty(ttest_thr)==0 || isempty(wilcoxon_thr)==0) && ((isempty(ttest_sub)==0 || isempty(wilcoxon_sub)==0))                                                       % if one of the statistical tests takes place in both conditions
    if (isempty(ttest_thr)==0 && isempty(wilcoxon_thr)==0) && (isempty(ttest_sub)==0 && isempty(wilcoxon_sub)==0)                                                     % if both statistical tests take place in both conditions
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_thr(2,:) ttest_sub(2,:) wilcoxon_thr(2,:) wilcoxon_sub(2,:)]',0.10,'pdep','yes');                            % 'yes' in order to see the report
    elseif (isempty(ttest_thr)==0 && isempty(wilcoxon_thr)==0) && (isempty(ttest_sub)==0 && isempty(wilcoxon_sub)==1)												  % if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Sub condition
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_thr(2,:) ttest_sub(2,:) wilcoxon_thr(2,:)]',0.10,'pdep','yes');                                                
    elseif (isempty(ttest_thr)==0 && isempty(wilcoxon_thr)==0) && (isempty(ttest_sub)==1 && isempty(wilcoxon_sub)==0)												  % if both statistical tests take place in the Thr condition, but only Wilcoxon ranksum test takes place in the Sub condition
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_thr(2,:) wilcoxon_thr(2,:) wilcoxon_sub(2,:)]',0.10,'pdep','yes');                                             
    elseif (isempty(ttest_thr)==0 && isempty(wilcoxon_thr)==1) && (isempty(ttest_sub)==0 && isempty(wilcoxon_sub)==0)												  % if only Two-sample T-test takes place in the Thr condition, but both statistical tests take place in the Sub condition
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_thr(2,:) ttest_sub(2,:) wilcoxon_sub(2,:)]',0.10,'pdep','yes');                                                                                                               
    elseif (isempty(ttest_thr)==0 && isempty(wilcoxon_thr)==1) && (isempty(ttest_sub)==0 && isempty(wilcoxon_sub)==1)										          % if only Two-sample T-test takes place in the both conditions
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_thr(2,:) ttest_sub(2,:)]',0.10,'pdep','yes');
    elseif (isempty(ttest_thr)==0 && isempty(wilcoxon_thr)==1) && (isempty(ttest_sub)==1 && isempty(wilcoxon_sub)==0)												  % if only Two-sample T-test takes place in the Thr condition, but only Wilcoxon ranksum test takes place in the Sub condition
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_thr(2,:) wilcoxon_sub(2,:)]',0.10,'pdep','yes');   
    elseif (isempty(ttest_thr)==1 && isempty(wilcoxon_thr)==0) && (isempty(ttest_sub)==0 && isempty(wilcoxon_sub)==0)												  % if only Wilcoxon ranksum test takes place in the Thr condition, but both statistical tests take place in the Sub condition
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_sub(2,:) wilcoxon_thr(2,:) wilcoxon_sub(2,:)]',0.10,'pdep','yes');                                                                                                        
    elseif (isempty(ttest_thr)==1 && isempty(wilcoxon_thr)==0) && (isempty(ttest_sub)==0 && isempty(wilcoxon_sub)==1)												  % if only Wilcoxon ranksum test takes place in the Thr condition, but only Two-sample T-test takes place in the Sub condition
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([ttest_sub(2,:) wilcoxon_thr(2,:)]',0.10,'pdep','yes');
    elseif (isempty(ttest_thr)==1 && isempty(wilcoxon_thr)==0) && (isempty(ttest_sub)==1 && isempty(wilcoxon_sub)==0)												  % if only Wilcoxon ranksum test takes place in the both conditions
        [bh_h,crit_p,adj_ci,adj_p]=fdr_bh([wilcoxon_thr(2,:) wilcoxon_sub(2,:)]',0.10,'pdep','yes');         
    end
end


[bh_h_anova,crit_p_anova,adj_ci_anova,adj_p_anova]=fdr_bh([anova_thr(1,:) anova_sub(1,:)]',0.10,'pdep','yes');

 
%% FINAL STORAGING:

% In this section, we create tables regarding the HRF parameters of
% interest, as well as the Shapiro-Wilk test, Two-sample T-test and 
% Wilcoxon ranksum test results for them and the FDR correction for the
% last two mentioned tests.


% Column with the statistical test's parameter names to identify its 
% corresponding results    
sw_test_values = array2table({'Hypothesis test result', 'p-value','W-value'}','VariableNames', {'Value'});  
ttest_test_values = array2table({'Hypothesis test result', 'p-value','t-value'}','VariableNames', {'Value'});  
wilcoxon_test_values = array2table({'Hypothesis test result', 'p-value','z-value'}','VariableNames', {'Value'});  


% ************* HRF parameters data and Shapiro-Wilk tables **************
    

% ««««««««««««««««««««« HRF parameters data tables »»»»»»»»»»»»»»»»»»»»»»» 
    

subjects_CNT_Thr = cell(n_CNT,1);
subjects_CNT_Sub = cell(n_CNT,1);
subjects_T2DM_Thr = cell(n_T2DM,1);
subjects_T2DM_Sub = cell(n_T2DM,1);


for c = 1:n_CNT
    
    % Identifies each subject in each group according to its condition
    subjects_CNT_Thr{c} = strcat('Subject', {' '}, num2str(c), {' '}, 'CNT', {' '}, 'Thr');
    subjects_CNT_Sub{c} = strcat('Subject', {' '}, num2str(c), {' '}, 'CNT', {' '}, 'Sub');
    
end

for t = 1:n_T2DM
    
    % Identifies each subject in each group according to its condition
    subjects_T2DM_Thr{t} = strcat('Subject', {' '}, num2str(t), {' '}, 'T2DM', {' '}, 'Thr');
    subjects_T2DM_Sub{t} = strcat('Subject', {' '}, num2str(t), {' '}, 'T2DM', {' '}, 'Sub');
    
end


% Forms tables with the HRF parameters' data from every subject in each 
% group and condition and in every ROI as well as its headers 
CNT_Thr_table = array2table(CNT_Thr,'VariableNames', headers);              
CNT_Sub_table = array2table(CNT_Sub,'VariableNames', headers);          
T2DM_Thr_table = array2table(T2DM_Thr,'VariableNames', headers);         
T2DM_Sub_table = array2table(T2DM_Sub,'VariableNames', headers);  

subjects_CNT_Thr_table = array2table(subjects_CNT_Thr,'VariableNames', {'Subject'});
subjects_CNT_Sub_table = array2table(subjects_CNT_Sub,'VariableNames', {'Subject'});
subjects_T2DM_Thr_table = array2table(subjects_T2DM_Thr,'VariableNames', {'Subject'});
subjects_T2DM_Sub_table = array2table(subjects_T2DM_Sub,'VariableNames', {'Subject'});


% Merges the HRF parameters' table with the corresponding subjects' number
% for each group and condition
CNT_Thr_table = [subjects_CNT_Thr_table CNT_Thr_table];              
CNT_Sub_table = [subjects_CNT_Sub_table CNT_Sub_table];          
T2DM_Thr_table = [subjects_T2DM_Thr_table T2DM_Thr_table];          
T2DM_Sub_table = [subjects_T2DM_Sub_table T2DM_Sub_table]; 

    
    
% «««««««««««««««««««««««« Shapiro-Wilk tables »»»»»»»»»»»»»»»»»»»»»»»»»»» 


% Forms tables with the Shapiro-Wilk test results for the HRF parameters 
% data from every subject in each group and condition and in every ROI 
% as well as its headers 
sw_CNT_Thr_table = array2table(sw_CNT_Thr,'VariableNames', headers);               
sw_CNT_Sub_table = array2table(sw_CNT_Sub,'VariableNames', headers);       
sw_T2DM_Thr_table = array2table(sw_T2DM_Thr,'VariableNames', headers);      
sw_T2DM_Sub_table = array2table(sw_T2DM_Sub,'VariableNames', headers);       
                
% Merges tables
sw_CNT_Thr_table = [sw_test_values sw_CNT_Thr_table];
sw_CNT_Sub_table = [sw_test_values sw_CNT_Sub_table];
sw_T2DM_Thr_table = [sw_test_values sw_T2DM_Thr_table];
sw_T2DM_Sub_table = [sw_test_values sw_T2DM_Sub_table];
    
    

% ********* Two sample T-test and Wilcoxon ranksum test tables ***********
    

% ««««««««««««««««««««««««« Two sample T-test »»»»»»»»»»»»»»»»»»»»»»»»»»» 
    

if isempty(ttest_thr)==0                                                                    % when a HRF parameter from different groups in the Thr condition has normal distribution  
        
    ttest_thr_headers = cell(1,size(ttest_thr,2));                                          % headers for the table

    for l = 1:length(ttest_indexes_thr)
        ttest_thr_headers(l) = headers(ttest_indexes_thr(l));                               % gets the headers of the normal Thr parameters
    end

    % Forms tables with the Two sample T-test results for the HRF
    % parameters with normal distribution from different groups in the Thr 
    % condition as well as its headers 
    ttest_thr_table = array2table(ttest_thr,'VariableNames', ttest_thr_headers);            

    % Merges the previous table with the statistical test's parameter names 
    ttest_thr_table = [ttest_test_values ttest_thr_table];
end
    
    
if isempty(ttest_sub)==0                                                                    % when a HRF parameter from different groups in the Sub condition has normal distribution       

    ttest_sub_headers = cell(1,size(ttest_sub,2));                                          

    for l = 1:length(ttest_indexes_sub)
        ttest_sub_headers(l) = headers(ttest_indexes_sub(l));                               % gets the headers of the normal Sub parameters
    end
    
    % Forms tables with the Two sample T-test results for the HRF
    % parameters with normal distribution from different groups in the Sub 
    % condition as well as its headers
    ttest_sub_table = array2table(ttest_sub,'VariableNames', ttest_sub_headers);            

    % Merges the previous table with the statistical test's parameter names 
    ttest_sub_table = [ttest_test_values ttest_sub_table];
end
    
    
% ««««««««««««««««««««««« Wilcoxon ranksum test »»»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(wilcoxon_thr)==0                                                                     % when a HRF parameter from different groups in the Thr condition has non-normal distribution 

    wilcoxon_thr_headers = cell(1,size(wilcoxon_thr,2));                                       

    for l = 1:length(wilcoxon_indexes_thr)
        wilcoxon_thr_headers(l) = headers(wilcoxon_indexes_thr(l));                             % gets the headers of non-normal Thr HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the HRF
    % parameters with non-normal distribution from different groups in the 
    % Thr condition as well as its headers
    wilcoxon_thr_table = array2table(wilcoxon_thr,'VariableNames', wilcoxon_thr_headers);       

    % Merges the previous table with the statistical test's parameter names 
    wilcoxon_thr_table = [wilcoxon_test_values wilcoxon_thr_table];
end


if isempty(wilcoxon_sub)==0                                                                     % when a HRF parameter from different groups in the Sub condition has non-normal distribution  

    wilcoxon_sub_headers = cell(1,size(wilcoxon_sub,2));                                        

    for l = 1:length(wilcoxon_indexes_sub)
        wilcoxon_sub_headers(l) = headers(wilcoxon_indexes_sub(l));                             % gets the headers of non-normal Sub HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the HRF
    % parameters with non-normal distribution from different groups in the 
    % Sub condition as well as its headers
    wilcoxon_sub_table = array2table(wilcoxon_sub,'VariableNames', wilcoxon_sub_headers);       

    % Merges the previous table with the statistical test's parameter names 
    wilcoxon_sub_table = [wilcoxon_test_values wilcoxon_sub_table];
end
    
    

% ************************** HRF Parameters ******************************
    

% «««««««««««««««««««« Average and standard deviation »»»»»»»»»»»»»»»»»»»»


% -------------------------------- CNT -----------------------------------


if isempty(avg_parameters_CNT_Thr)==0                                                                   % when a CNT Thr HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI

    avg_headers_CNT_Thr = cell(1,size(avg_parameters_CNT_Thr,2));                                       % headers for the table

    for l = 1:length(avg_std_indexes_CNT_Thr)
        avg_headers_CNT_Thr(l) = headers(avg_std_indexes_CNT_Thr(l));                                   % gets the headers of normal CNT Thr HRF parameters
    end

    % Forms tables with the average and standard deviation for the CNT Thr
    % HRF parameters with normal distribution as well as its headers
    avg_parameters_CNT_Thr_table = array2table(avg_parameters_CNT_Thr,'VariableNames', avg_headers_CNT_Thr);       
    std_parameters_CNT_Thr_table = array2table(std_parameters_CNT_Thr,'VariableNames', avg_headers_CNT_Thr);              
end


if isempty(avg_parameters_CNT_Sub)==0                                                                   % when a CNT Sub HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI      

    avg_headers_CNT_Sub = cell(1,size(avg_parameters_CNT_Sub,2));                                         

    for l = 1:length(avg_std_indexes_CNT_Sub)
        avg_headers_CNT_Sub(l) = headers(avg_std_indexes_CNT_Sub(l));                                   % gets the headers of normal CNT Sub HRF parameters
    end

    % Forms tables with the average and standard deviation for the CNT Sub
    % HRF parameters with normal distribution as well as its headers
    avg_parameters_CNT_Sub_table = array2table(avg_parameters_CNT_Sub,'VariableNames', avg_headers_CNT_Sub);       
    std_parameters_CNT_Sub_table = array2table(std_parameters_CNT_Sub,'VariableNames', avg_headers_CNT_Sub);         
end
    
    
% -------------------------------- T2DM ----------------------------------


if isempty(avg_parameters_T2DM_Thr)==0                                                                  % when a T2DM Thr HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI       

    avg_headers_T2DM_Thr = cell(1,size(avg_parameters_T2DM_Thr,2));                                       

    for l = 1:length(avg_std_indexes_T2DM_Thr)
        avg_headers_T2DM_Thr(l) = headers(avg_std_indexes_T2DM_Thr(l));                                 % gets the headers of normal T2DM Thr HRF parameters
    end
    
    % Forms tables with the average and standard deviation for the T2DM Thr
    % HRF parameters with normal distribution as well as its headers
    avg_parameters_T2DM_Thr_table = array2table(avg_parameters_T2DM_Thr,'VariableNames', avg_headers_T2DM_Thr);    
    std_parameters_T2DM_Thr_table = array2table(std_parameters_T2DM_Thr,'VariableNames', avg_headers_T2DM_Thr);      
end


if isempty(avg_parameters_T2DM_Sub)==0                                                                  % when a T2DM Sub HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI      

    avg_headers_T2DM_Sub = cell(1,size(avg_parameters_T2DM_Sub,2));                                       

    for l = 1:length(avg_std_indexes_T2DM_Sub)
        avg_headers_T2DM_Sub(l) = headers(avg_std_indexes_T2DM_Sub(l));                                 % gets the headers of normal T2DM Sub HRF parameters
    end

    % Forms tables with the average and standard deviation for the T2DM Sub
    % HRF parameters with normal distribution as well as its headers
    avg_parameters_T2DM_Sub_table = array2table(avg_parameters_T2DM_Sub,'VariableNames', avg_headers_T2DM_Sub);    
    std_parameters_T2DM_Sub_table = array2table(std_parameters_T2DM_Sub,'VariableNames', avg_headers_T2DM_Sub);      
end



% ««««««««««««««««««« Median and inter-quartile range »»»»»»»»»»»»»»»»»»»»


% -------------------------------- CNT -----------------------------------


if isempty(median_parameters_CNT_Thr)==0                                                                            % when a CNT Thr HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI    

    median_headers_CNT_Thr = cell(1,size(median_parameters_CNT_Thr,2));                                             % headers for the table

    for l = 1:length(median_iqr_indexes_CNT_Thr)
        median_headers_CNT_Thr(l) = headers(median_iqr_indexes_CNT_Thr(l));                                         % gets the headers of non-normal CNT Thr HRF parameters
    end
    
    % Forms tables with the median and interquartile range for the CNT Thr
    % HRF parameters with non-normal distribution as well as its headers
    median_parameters_CNT_Thr_table = array2table(median_parameters_CNT_Thr,'VariableNames', median_headers_CNT_Thr);          
    iqr_parameters_CNT_Thr_table = array2table(iqr_parameters_CNT_Thr,'VariableNames', median_headers_CNT_Thr);                            
end

if isempty(median_parameters_CNT_Sub)==0                                                                            % when a CNT Sub HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI              

    median_headers_CNT_Sub = cell(1,size(median_parameters_CNT_Sub,2));                                            

    for l = 1:length(median_iqr_indexes_CNT_Sub)
        median_headers_CNT_Sub(l) = headers(median_iqr_indexes_CNT_Sub(l));                                         % gets the headers of non-normal CNT Sub HRF parameters
    end
    
    % Forms tables with the median and interquartile range for the CNT Sub
    % HRF parameters with non-normal distribution as well as its headers
    median_parameters_CNT_Sub_table = array2table(median_parameters_CNT_Sub,'VariableNames', median_headers_CNT_Sub);          
    iqr_parameters_CNT_Sub_table = array2table(iqr_parameters_CNT_Sub,'VariableNames', median_headers_CNT_Sub);               
end


% -------------------------------- T2DM ----------------------------------


if isempty(median_parameters_T2DM_Thr)==0                                                                           % when a T2DM Thr HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI        

    median_headers_T2DM_Thr = cell(1,size(median_parameters_T2DM_Thr,2));                                          

    for l = 1:length(median_iqr_indexes_T2DM_Thr)
        median_headers_T2DM_Thr(l) = headers(median_iqr_indexes_T2DM_Thr(l));                                       % gets the headers of non-normal T2DM Thr HRF parameters
    end

    % Forms tables with the median and interquartile range for the T2DM Thr
    % HRF parameters with non-normal distribution as well as its headers
    median_parameters_T2DM_Thr_table = array2table(median_parameters_T2DM_Thr,'VariableNames', median_headers_T2DM_Thr);       
    iqr_parameters_T2DM_Thr_table = array2table(iqr_parameters_T2DM_Thr,'VariableNames', median_headers_T2DM_Thr);       
end

if isempty(median_parameters_T2DM_Sub)==0                                                                           % when a T2DM Sub HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI        

    median_headers_T2DM_Sub = cell(1,size(median_parameters_T2DM_Sub,2));                                          

    for l = 1:length(median_iqr_indexes_T2DM_Sub)
        median_headers_T2DM_Sub(l) = headers(median_iqr_indexes_T2DM_Sub(l));                                       % gets the headers of non-normal T2DM Sub HRF parameters
    end

    % Forms tables with the median and interquartile range for the T2DM Sub
    % HRF parameters with non-normal distribution as well as its headers
    median_parameters_T2DM_Sub_table = array2table(median_parameters_T2DM_Sub,'VariableNames', median_headers_T2DM_Sub);       
    iqr_parameters_T2DM_Sub_table = array2table(iqr_parameters_T2DM_Sub,'VariableNames', median_headers_T2DM_Sub);       
end


% ************************ One-way ANOVA tables *************************


anova_values = array2table({'p-value', 'F-value'}','VariableNames', {'Value'});  


% Forms tables
anova_thr_table = array2table(anova_thr, 'VariableNames', anova_headers);
anova_sub_table = array2table(anova_sub, 'VariableNames', anova_headers);

% Merges tables
anova_thr_table = [anova_values anova_thr_table];
anova_sub_table = [anova_values anova_sub_table];



% *********************** FDR Correction tables *************************
   

% Column with the statistical test's parameter names to identify its 
% corresponding results    
bh_test_values = array2table({'Hypothesis test result', 'Adjusted p-value'}','VariableNames', {'Value'});  


% ««««««««««««««««««««««««« Two sample T-test »»»»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(ttest_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the HRF parameters in the Threshold 
    % condition as well as its headers
    bh_ttest_thr_table = array2table([bh_h(1:size(ttest_thr,2),1) adj_p(1:size(ttest_thr,2),1)]','VariableNames', ttest_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    bh_ttest_thr_table = [bh_test_values bh_ttest_thr_table];
end


if isempty(ttest_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the HRF parameters in the Submaximum 
    % condition as well as its headers
    bh_ttest_sub_table = array2table([bh_h(1+size(ttest_thr,2):size(ttest_thr,2)+size(ttest_sub,2),1) adj_p(1+size(ttest_thr,2):size(ttest_thr,2)+size(ttest_sub,2),1)]','VariableNames', ttest_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    bh_ttest_sub_table = [bh_test_values bh_ttest_sub_table];
end


% «««««««««««««««««««««««« Wilcoxon ranksum test »»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(wilcoxon_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the HRF parameters in 
    % the Threshold condition as well as its headers
    bh_wilcoxon_thr_table = array2table([bh_h(1+size(ttest_thr,2)+size(ttest_sub,2):size(ttest_thr,2)+size(ttest_sub,2)+size(wilcoxon_thr,2),1) adj_p(1+size(ttest_thr,2)+size(ttest_sub,2):size(ttest_thr,2)+size(ttest_sub,2)+size(wilcoxon_thr,2),1)]','VariableNames', wilcoxon_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    bh_wilcoxon_thr_table = [bh_test_values bh_wilcoxon_thr_table];
end


if isempty(wilcoxon_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the HRF parameters in 
    % the Submaximum condition as well as its headers
    bh_wilcoxon_sub_table = array2table([bh_h(1+size(ttest_thr,2)+size(ttest_sub,2)+size(wilcoxon_thr,2):end,1) adj_p(1+size(ttest_thr,2)+size(ttest_sub,2)+size(wilcoxon_thr,2):end,1)]','VariableNames', wilcoxon_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    bh_wilcoxon_sub_table = [bh_test_values bh_wilcoxon_sub_table];
end



% «««««««««««««««««««««««««««« One-way ANOVA »»»»»»»»»»»»»»»»»»»»»»»»»»»» 

    
% Forms tables with the Benjamini & Hochberg approach to control the FDR of 
% the One-way ANOVA for the HRF parameters in both conditions as well as 
% its headers
bh_anova_thr_table = array2table([bh_h_anova(1:size(anova_thr,2),1) adj_p_anova(1:size(anova_thr,2),1)]','VariableNames', anova_headers);
bh_anova_sub_table = array2table([bh_h_anova(1+size(anova_thr,2):size(anova_thr,2)+size(anova_sub,2),1) adj_p_anova(1+size(anova_thr,2):size(anova_thr,2)+size(anova_sub,2),1)]','VariableNames', anova_headers);

% Merges the previous table with the statistical test's parameter names
bh_anova_thr_table = [bh_test_values bh_anova_thr_table];
bh_anova_sub_table = [bh_test_values bh_anova_sub_table];
