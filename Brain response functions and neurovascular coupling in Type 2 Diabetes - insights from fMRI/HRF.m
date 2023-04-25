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

load('ROI_betas.mat');
load('ROI_SEs.mat');
load('ROI_regions.mat');


%% CONSTANTS:

n_points = 8;                             % data points
time = (0:n_points-1)* 2.5;               % volumes per subject * TR (2.5 s)
n_subjects = 141;
n_T2DM = 64; 
n_CNT = 77;
n_rois = 22;
n_parameters = 11;


%% RAW DATA EXTRACTION AND STORAGING:

% In this section, we retrieve and then store the raw data (betas and 
% standard errors) from both sides across subjects of each group and 
% condition


betas_T2DM_Thr = cell(n_T2DM*2*n_points,n_rois+1);
betas_T2DM_Sub = cell(n_T2DM*2*n_points,n_rois+1);

SE_T2DM_Thr = cell(n_T2DM*2*n_points,n_rois+1);
SE_T2DM_Sub = cell(n_T2DM*2*n_points,n_rois+1);

betas_CNT_Thr = cell(n_CNT*2*n_points,n_rois+1);
betas_CNT_Sub = cell(n_CNT*2*n_points,n_rois+1);

SE_CNT_Thr = cell(n_CNT*2*n_points,n_rois+1);
SE_CNT_Sub = cell(n_CNT*2*n_points,n_rois+1);


% Data from the same condition and subject has a 16-row step, hence the 
% n_points*2 (=16) expression in the storaging arrays.
% Data from different subjects has a 32-row step, hence the n_points*4
% (=32) expression.

% Overall storaging order:
%       1. Betas / SE from the left side (Thr / Sub)
%       2. Betas / SE from the right side (Thr / Sub)


% -------------------------------- CNT -----------------------------------

%   n_points*4*n_T2DM+n_points*2+1 = 2065 --> row from which left 
%   hemisphere activation data of CNT subjects in the Threshold condition
%   begins
%   n_points*4*n_T2DM+1 = 2049 --> row from which right hemisphere 
%   activation data of CNT subjects in the Threshold condition begins
%   n_points*4*n_T2DM+n_points*3+1 = 2073 --> row from which left 
%   hemisphere activation data of CNT subjects in the Submaximum condition 
%   begins
%   n_points*4*n_T2DM+n_points+1 = 2057 --> row from which right hemisphere 
%   activation data of CNT subjects in the Submaximum condition begins

for c = 1:n_CNT
    betas_CNT_Thr((1+n_points*2*(c-1):n_points*2+n_points*2*(c-1)),:) = [ROIbetas((n_points*4*n_T2DM+n_points*2+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points*3)+((c-1)*n_points*4),:); ROIbetas((n_points*4*n_T2DM+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points)+((c-1)*n_points*4),:)];
    betas_CNT_Sub((1+n_points*2*(c-1):n_points*2+n_points*2*(c-1)),:) = [ROIbetas((n_points*4*n_T2DM+n_points*3+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points*4)+((c-1)*n_points*4),:); ROIbetas((n_points*4*n_T2DM+n_points+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points*2)+((c-1)*n_points*4),:)];
    SE_CNT_Thr((1+n_points*2*(c-1):n_points*2+n_points*2*(c-1)),:) = [ROISEs((n_points*4*n_T2DM+n_points*2+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points*3)+((c-1)*n_points*4),:); ROISEs((n_points*4*n_T2DM+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points)+((c-1)*n_points*4),:)];
    SE_CNT_Sub((1+n_points*2*(c-1):n_points*2+n_points*2*(c-1)),:) = [ROISEs((n_points*4*n_T2DM+n_points*3+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points*4)+((c-1)*n_points*4),:); ROISEs((n_points*4*n_T2DM+n_points+1)+((c-1)*n_points*4): (n_points*4*n_T2DM+n_points*2)+((c-1)*n_points*4),:)];    
end


% ------------------------------- T2DM -----------------------------------

%   n_points*2+1 = 17 --> row from which left hemisphere activation data 
%   of T2DM subjects in the Threshold condition begins
%   1 --> row from which right hemisphere activation data of T2DM subjects 
%   in the Threshold condition begins
%   n_points*3+1 = 25 --> row from which left hemisphere activation data 
%   of T2DM subjects in the Submaximum condition begins
%   n_points+1 = 9 --> row from which right hemisphere activation data of
%   T2DM subjects in the Submaximum condition begins

for t = 1:n_T2DM
    betas_T2DM_Thr((1+n_points*2*(t-1):n_points*2+n_points*2*(t-1)),:) = [ROIbetas((n_points*2+1)+((t-1)*n_points*4): (n_points*3)+((t-1)*n_points*4),:); ROIbetas(1+((t-1)*n_points*4): n_points+((t-1)*n_points*4),:)];
    betas_T2DM_Sub((1+n_points*2*(t-1):n_points*2+n_points*2*(t-1)),:) = [ROIbetas((n_points*3+1)+((t-1)*n_points*4): (n_points*4)+((t-1)*n_points*4),:); ROIbetas((n_points+1)+((t-1)*n_points*4): (n_points*2)+((t-1)*n_points*4),:)];
    SE_T2DM_Thr((1+n_points*2*(t-1):n_points*2+n_points*2*(t-1)),:) = [ROISEs((n_points*2+1)+((t-1)*n_points*4): (n_points*3)+((t-1)*n_points*4),:); ROISEs(1+((t-1)*n_points*4): n_points+((t-1)*n_points*4),:)];
    SE_T2DM_Sub((1+n_points*2*(t-1):n_points*2+n_points*2*(t-1)),:) = [ROISEs((n_points*3+1)+((t-1)*n_points*4): (n_points*4)+((t-1)*n_points*4),:); ROISEs((n_points+1)+((t-1)*n_points*4): (n_points*2)+((t-1)*n_points*4),:)];
end


% -----------------------------------------------------------------------


CNT_info = cell(length(betas_CNT_Thr),1);
T2DM_info = cell(length(betas_T2DM_Thr),1);

headers_ROIs = cell(1,length(ROI_regions));


% Identifies each datapoint from each CNT subject in each condition
for c = 1:n_CNT
    for side = 1:2
        for dp = 1:n_points
            if side == 1
                CNT_info{dp+(n_points*(side-1))+(n_points*2*(c-1))} = {sprintf('%s %d %s %d %s', 'Subject', c, 'datapoint', dp, 'Left')};
            else
                CNT_info{dp+(n_points*(side-1))+(n_points*2*(c-1))} = {sprintf('%s %d %s %d %s', 'Subject', c, 'datapoint', dp, 'Right')};
            end
        end
    end
end

% Identifies each datapoint from each T2DM subject in each condition
for t = 1:n_T2DM
    for side = 1:2
        for dp = 1:n_points
            if side == 1
                T2DM_info{dp+(n_points*(side-1))+(n_points*2*(t-1))} = {sprintf('%s %d %s %d %s', 'Subject', t, 'datapoint', dp, 'Left')};
            else
                T2DM_info{dp+(n_points*(side-1))+(n_points*2*(t-1))} = {sprintf('%s %d %s %d %s', 'Subject', t, 'datapoint', dp, 'Right')};
            end
        end
    end
end


% Automatically sets the table headers (ROI names)
for r = 1:length(ROI_regions)
    headers_ROIs{r} = char(strrep(strcat(ROI_regions(r)),' ','_'));         % concatenates strings, replaces spaces by _ and converts to char
end


% Forms tables
raw_betas_CNT_Thr = array2table(betas_CNT_Thr(:,2:end),'VariableNames',headers_ROIs);
raw_betas_CNT_Sub = array2table(betas_CNT_Sub(:,2:end),'VariableNames',headers_ROIs);
raw_SE_CNT_Thr = array2table(SE_CNT_Thr(:,2:end),'VariableNames',headers_ROIs);
raw_SE_CNT_Sub = array2table(SE_CNT_Sub(:,2:end),'VariableNames',headers_ROIs);

raw_betas_T2DM_Thr = array2table(betas_T2DM_Thr(:,2:end),'VariableNames',headers_ROIs);
raw_betas_T2DM_Sub = array2table(betas_T2DM_Sub(:,2:end),'VariableNames',headers_ROIs);
raw_SE_T2DM_Thr = array2table(SE_T2DM_Thr(:,2:end),'VariableNames',headers_ROIs);
raw_SE_T2DM_Sub = array2table(SE_T2DM_Sub(:,2:end),'VariableNames',headers_ROIs);

headers_CNT = array2table(CNT_info,'VariableNames', {'Subject_datapoint_side'});           
headers_T2DM = array2table(T2DM_info,'VariableNames', {'Subject_datapoint_side'});           


% Merges tables
table_raw_betas_CNT_Thr = [headers_CNT raw_betas_CNT_Thr];
table_raw_betas_CNT_Sub = [headers_CNT raw_betas_CNT_Sub];
table_raw_SE_CNT_Thr = [headers_CNT raw_SE_CNT_Thr];
table_raw_SE_CNT_Sub = [headers_CNT raw_SE_CNT_Sub];

table_raw_betas_T2DM_Thr = [headers_T2DM raw_betas_T2DM_Thr];
table_raw_betas_T2DM_Sub = [headers_T2DM raw_betas_T2DM_Sub];
table_raw_SE_T2DM_Thr = [headers_T2DM raw_SE_T2DM_Thr];
table_raw_SE_T2DM_Sub = [headers_T2DM raw_SE_T2DM_Sub];


%% AVERAGE AND MEDIAN HRF CURVES:

% In this section, we estimate the average and median HRF curves in each 
% condition and group per Regions of Interest (ROIs).


beta_thr_both_sides = zeros(n_points,2);                         % 2 since each column represents a stimulation side - left and right
beta_sub_both_sides = zeros(n_points,2);

beta_per_subject_CNT_Thr = zeros(n_points,n_CNT);
beta_per_subject_CNT_Sub = zeros(n_points,n_CNT);
beta_per_subject_T2DM_Thr = zeros(n_points,n_T2DM);
beta_per_subject_T2DM_Sub = zeros(n_points,n_T2DM);

total_avg_HRF_CNT_Thr = zeros(n_points,n_rois);
total_avg_HRF_CNT_Sub = zeros(n_points,n_rois);
total_avg_HRF_T2DM_Thr = zeros(n_points,n_rois);
total_avg_HRF_T2DM_Sub = zeros(n_points,n_rois);

total_std_HRF_CNT_Thr = zeros(n_points,n_rois);
total_std_HRF_CNT_Sub = zeros(n_points,n_rois);
total_std_HRF_T2DM_Thr = zeros(n_points,n_rois);
total_std_HRF_T2DM_Sub = zeros(n_points,n_rois);

total_median_HRF_CNT_Thr = zeros(n_points,n_rois);
total_median_HRF_CNT_Sub = zeros(n_points,n_rois);
total_median_HRF_T2DM_Thr = zeros(n_points,n_rois);
total_median_HRF_T2DM_Sub = zeros(n_points,n_rois);

total_iqr_HRF_CNT_Thr = zeros(n_points,n_rois);
total_iqr_HRF_CNT_Sub = zeros(n_points,n_rois);
total_iqr_HRF_T2DM_Thr = zeros(n_points,n_rois);
total_iqr_HRF_T2DM_Sub = zeros(n_points,n_rois);


for r=2:n_rois+1
    for i=1:n_subjects
        for j=1:2
            beta_thr_both_sides(:,j) = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));    % each column is a side in each subject 
            beta_sub_both_sides(:,j) = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r)); 
        end                                            

        if i<65                                                                     % until subject 65 --> T2DM 
            beta_per_subject_T2DM_Thr(:,i) = mean(beta_thr_both_sides,2);           % mean of the beta values from each side per datapoint of T2DM subjects in the Thr condition in a ROI
            beta_per_subject_T2DM_Sub(:,i) = mean(beta_sub_both_sides,2);
        else                                                                        % after subject 65 --> CNT
            beta_per_subject_CNT_Thr(:,i-n_T2DM) = mean(beta_thr_both_sides,2);     % mean of the beta values from each side per datapoint of CNT subjects in the Thr condition in a ROI 
            beta_per_subject_CNT_Sub(:,i-n_T2DM) = mean(beta_sub_both_sides,2);
        end
    
        % Deletes to not duplicate or interfere with data
        beta_thr_both_sides = zeros(n_points,2);
        beta_sub_both_sides = zeros(n_points,2);
    end  
    
    
    % --------------------- THRESHOLD CONDITION -------------------------
    
    % Average and median HRF per datapoint in each subject group of the Thr 
    % condition in a ROI and its corresponding standard deviation and 
    % interquartile range 
    total_avg_HRF_CNT_Thr(:,r-1) = mean(beta_per_subject_CNT_Thr,2);            % stores the average beta values per datapoint of CNT subjects in the Thr condition across ROIs
    total_avg_HRF_T2DM_Thr(:,r-1) = mean(beta_per_subject_T2DM_Thr,2);          % stores the average beta values per datapoint of T2DM subjects in the Thr condition across ROIs
    total_std_HRF_CNT_Thr(:,r-1) = std(beta_per_subject_CNT_Thr,0,2);           % stores the standard deviation of the beta values per datapoint of CNT subjects in the Thr condition across ROIs
    total_std_HRF_T2DM_Thr(:,r-1) = std(beta_per_subject_T2DM_Thr,0,2);         % stores the standard deviation of the beta values per datapoint of T2DM subjects in the Thr condition across ROIs
    
    total_median_HRF_CNT_Thr(:,r-1) = median(beta_per_subject_CNT_Thr,2);       % stores the median beta values per datapoint of CNT subjects in the Thr condition across ROIs
    total_median_HRF_T2DM_Thr(:,r-1) = median(beta_per_subject_T2DM_Thr,2);     % stores the median beta values per datapoint of T2DM subjects in the Thr condition across ROIs
    total_iqr_HRF_CNT_Thr(:,r-1) = iqr(beta_per_subject_CNT_Thr,2);             % stores the interquartile range of the beta values per datapoint of CNT subjects in the Thr condition across ROIs
    total_iqr_HRF_T2DM_Thr(:,r-1) = iqr(beta_per_subject_T2DM_Thr,2);           % stores the interquartile range of the beta values per datapoint of T2DM subjects in the Thr condition across ROIs
        
        
    % ---------------------- SUBMAXIMUM CONDITION ------------------------
   
    % Average and median HRF per datapoint in each subject group of the Sub 
    % condition in a ROI and its corresponding standard deviation and 
    % interquartile range             
    total_avg_HRF_CNT_Sub(:,r-1) = mean(beta_per_subject_CNT_Sub,2);         
    total_avg_HRF_T2DM_Sub(:,r-1) = mean(beta_per_subject_T2DM_Sub,2);      
    total_std_HRF_CNT_Sub(:,r-1) = std(beta_per_subject_CNT_Sub,0,2);            
    total_std_HRF_T2DM_Sub(:,r-1) = std(beta_per_subject_T2DM_Sub,0,2);           
    
    total_median_HRF_CNT_Sub(:,r-1) = median(beta_per_subject_CNT_Sub,2);         
    total_median_HRF_T2DM_Sub(:,r-1) = median(beta_per_subject_T2DM_Sub,2);      
    total_iqr_HRF_CNT_Sub(:,r-1) = iqr(beta_per_subject_CNT_Sub,2);                  
    total_iqr_HRF_T2DM_Sub(:,r-1) = iqr(beta_per_subject_T2DM_Sub,2);               
       
        
    % Deletes to not duplicate or interfere with data
    beta_per_subject_CNT_Thr = zeros(n_points,n_CNT);
    beta_per_subject_CNT_Sub = zeros(n_points,n_CNT);
    beta_per_subject_T2DM_Thr = zeros(n_points,n_T2DM);
    beta_per_subject_T2DM_Sub = zeros(n_points,n_T2DM);
end


%% PLOTTING THE AVERAGE AND MEDIAN HRF CURVES: 

% In this section, we plot the average and median HRF curve per group and
% condition in each ROI, separating positive from negative signal change
% regions (regions where the HRF curve increases after an activation, 
% and regions where the HRF curve decreases after an activation,
% respectively).


% Plots the average HRF according to each set of ROIs
figure('Name', sprintf('Average HRF - Positive signal change ROIs'))
for a=1:10                                                                  % 10 positive signal change regions of interest
    subplot(2,5,a);
    shadedErrorBar(time, total_avg_HRF_CNT_Thr(:,a), (total_std_HRF_CNT_Thr(:,a))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_avg_HRF_T2DM_Thr(:,a), (total_std_HRF_T2DM_Thr(:,a))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_avg_HRF_CNT_Sub(:,a), (total_std_HRF_CNT_Sub(:,a))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_avg_HRF_T2DM_Sub(:,a), (total_std_HRF_T2DM_Sub(:,a))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
    title(ROI_regions{a});
    xlabel('Time (s)')
    ylabel('Beta values')    
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Average HRF - Negative signal change ROIs'))
for d=1:12                                                                  % 12 negative signal change regions of interest
    subplot(3,4,d);
    shadedErrorBar(time, total_avg_HRF_CNT_Thr(:,10+d), (total_std_HRF_CNT_Thr(:,10+d))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_avg_HRF_T2DM_Thr(:,10+d), (total_std_HRF_T2DM_Thr(:,10+d))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)    
    shadedErrorBar(time, total_avg_HRF_CNT_Sub(:,10+d), (total_std_HRF_CNT_Sub(:,10+d))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_avg_HRF_T2DM_Sub(:,10+d), (total_std_HRF_T2DM_Sub(:,10+d))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)    
    title(ROI_regions{10+d});
    xlabel('Time (s)')
    ylabel('Beta values')  
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


% Plots the median HRF according to each set of ROIs
figure('Name', sprintf('Median HRF - Positive signal change ROIs'))
for a=1:10                                                       
    subplot(2,5,a);
    shadedErrorBar(time, total_median_HRF_CNT_Thr(:,a), (total_iqr_HRF_CNT_Thr(:,a))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_median_HRF_T2DM_Thr(:,a), (total_iqr_HRF_T2DM_Thr(:,a))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)    
    shadedErrorBar(time, total_median_HRF_CNT_Sub(:,a), (total_iqr_HRF_CNT_Sub(:,a))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_median_HRF_T2DM_Sub(:,a), (total_iqr_HRF_T2DM_Sub(:,a))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)    
    title(ROI_regions{a});
    xlabel('Time (s)')
    ylabel('Beta values')    
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Median HRF - Negative signal change ROIs'))
for d=1:12                                                       
    subplot(3,4,d);
    shadedErrorBar(time, total_median_HRF_CNT_Thr(:,10+d), (total_iqr_HRF_CNT_Thr(:,10+d))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_median_HRF_T2DM_Thr(:,10+d), (total_iqr_HRF_T2DM_Thr(:,10+d))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_median_HRF_CNT_Sub(:,10+d), (total_iqr_HRF_CNT_Sub(:,10+d))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_median_HRF_T2DM_Sub(:,10+d), (total_iqr_HRF_T2DM_Sub(:,10+d))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
    title(ROI_regions{10+d});
    xlabel('Time (s)')
    ylabel('Beta values')  
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


%% GRAND AVERAGE / MEDIAN ANALYSIS - PLOTS:

% In this section, we plot the overall average / median of the average /
% median HRF curves per group and condition in each ROI for each set of 
% ROIs 


% Plots the average HRF according to each set of ROIs
figure('Name', sprintf('Grand Average HRF - Positive signal change ROIs'))
shadedErrorBar(time, mean(total_avg_HRF_CNT_Thr(:,1:10),2), (std(total_avg_HRF_CNT_Thr(:,1:10),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_HRF_T2DM_Thr(:,1:10),2), (std(total_avg_HRF_T2DM_Thr(:,1:10),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_HRF_CNT_Sub(:,1:10),2), (std(total_avg_HRF_CNT_Sub(:,1:10),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, mean(total_avg_HRF_T2DM_Sub(:,1:10),2), (std(total_avg_HRF_T2DM_Sub(:,1:10),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights') 
title('Grand Average HRF - Positive signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Grand Average HRF - Negative signal change ROIs'))
shadedErrorBar(time, mean(total_avg_HRF_CNT_Thr(:,11:end),2), (std(total_avg_HRF_CNT_Thr(:,11:end),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_HRF_T2DM_Thr(:,11:end),2), (std(total_avg_HRF_T2DM_Thr(:,11:end),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_HRF_CNT_Sub(:,11:end),2), (std(total_avg_HRF_CNT_Sub(:,11:end),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, mean(total_avg_HRF_T2DM_Sub(:,11:end),2), (std(total_avg_HRF_T2DM_Sub(:,11:end),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights')  
title('Grand Average HRF - Negative signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


% Plots the median HRF according to each set of ROIs
figure('Name', sprintf('Grand Median HRF - Positive signal change ROIs'))
shadedErrorBar(time, median(total_median_HRF_CNT_Thr(:,1:10),2), (std(total_median_HRF_CNT_Thr(:,1:10),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_HRF_T2DM_Thr(:,1:10),2), (std(total_median_HRF_T2DM_Thr(:,1:10),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_HRF_CNT_Sub(:,1:10),2), (std(total_median_HRF_CNT_Sub(:,1:10),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, median(total_median_HRF_T2DM_Sub(:,1:10),2), (std(total_median_HRF_T2DM_Sub(:,1:10),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights')  
title('Grand Median HRF - Positive signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Grand Median HRF - Negative signal change ROIs'))
shadedErrorBar(time, median(total_median_HRF_CNT_Thr(:,11:end),2), (std(total_median_HRF_CNT_Thr(:,11:end),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_HRF_T2DM_Thr(:,11:end),2), (std(total_median_HRF_T2DM_Thr(:,11:end),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_HRF_CNT_Sub(:,11:end),2), (std(total_median_HRF_CNT_Sub(:,11:end),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, median(total_median_HRF_T2DM_Sub(:,11:end),2), (std(total_median_HRF_T2DM_Sub(:,11:end),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights')  
title('Grand Median HRF - Negative signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


%% COEFFICIENT OF VARIATION:

% In this section, we estimate the coefficient of variation for the HRF 
% peak amplitude and peak latency in each ROI


beta_thr_both_sides = zeros(n_points,2);                         % 2 since each column represents a stimulation side - left and right
beta_sub_both_sides = zeros(n_points,2);

peak_per_subject_CNT_Thr = zeros(n_CNT,n_rois);
peak_per_subject_CNT_Sub = zeros(n_CNT,n_rois);
peak_per_subject_T2DM_Thr = zeros(n_T2DM,n_rois);
peak_per_subject_T2DM_Sub = zeros(n_T2DM,n_rois);

peak_latency_per_subject_CNT_Thr = zeros(n_CNT,n_rois);
peak_latency_per_subject_CNT_Sub = zeros(n_CNT,n_rois);
peak_latency_per_subject_T2DM_Thr = zeros(n_T2DM,n_rois);
peak_latency_per_subject_T2DM_Sub = zeros(n_T2DM,n_rois);
    
total_CV_peak_CNT_Thr = zeros(1,n_rois);
total_CV_peak_CNT_Sub = zeros(1,n_rois);
total_CV_peak_T2DM_Thr = zeros(1,n_rois);
total_CV_peak_T2DM_Sub = zeros(1,n_rois);

total_CV_peak_latency_CNT_Thr = zeros(1,n_rois);
total_CV_peak_latency_CNT_Sub = zeros(1,n_rois);
total_CV_peak_latency_T2DM_Thr = zeros(1,n_rois);
total_CV_peak_latency_T2DM_Sub = zeros(1,n_rois);


for r=2:(n_rois+1)
    for i=1:n_subjects
        for j=1:2
            beta_thr_both_sides(:,j) = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));    % each column is a side in each subject 
            beta_sub_both_sides(:,j) = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));
        end

        mean_thr = mean(beta_thr_both_sides,2);                                         % mean of the beta values from each side per datapoint of a subject in the Thr condition in a ROI
        mean_sub = mean(beta_sub_both_sides,2);                                               

        if i<65                                                                         % until subject 65 --> T2DM
            [peak_T2DM_Thr, peak_latency_T2DM_Thr] = max(mean_thr(2:6,1));              % gets the peak amplitude and peak latency within the sectioned time range (2.5 to 12.5 s) in T2DM subjects of the Thr condition
            [peak_T2DM_Sub, peak_latency_T2DM_Sub] = max(mean_sub(2:6,1));                           
            peak_per_subject_T2DM_Thr(i,r-1) = peak_T2DM_Thr;                           % stores the HRF peak amplitude value of T2DM subjects in the Thr condition in each ROI 
            peak_per_subject_T2DM_Sub(i,r-1) = peak_T2DM_Sub;
            peak_latency_per_subject_T2DM_Thr(i,r-1) = peak_latency_T2DM_Thr*2.5;       % stores the HRF peak latency value of T2DM subjects in the Thr condition in each ROI
            peak_latency_per_subject_T2DM_Sub(i,r-1) = peak_latency_T2DM_Sub*2.5;
        else                                                                            % after subject 65 --> CNT
            [peak_CNT_Thr, peak_latency_CNT_Thr] = max(mean_thr(2:6,1));                % gets the peak amplitude and peak latency within the sectioned time range (2.5 to 12.5 s) in CNT individuals of the Thr condition
            [peak_CNT_Sub, peak_latency_CNT_Sub] = max(mean_sub(2:6,1));
            peak_per_subject_CNT_Thr(i-n_T2DM,r-1) = peak_CNT_Thr;                      % stores the HRF peak amplitude value of CNT subjects in the Thr condition in each ROI
            peak_per_subject_CNT_Sub(i-n_T2DM,r-1) = peak_CNT_Sub; 
            peak_latency_per_subject_CNT_Thr(i-n_T2DM,r-1) = peak_latency_CNT_Thr*2.5;  % stores the HRF peak latency value of CNT subjects in the Thr condition in each ROI
            peak_latency_per_subject_CNT_Sub(i-n_T2DM,r-1) = peak_latency_CNT_Sub*2.5;
        end    
        
        % Delete to not duplicate or interfere with data
        beta_thr_both_sides = zeros(n_points,2);
        beta_sub_both_sides = zeros(n_points,2);
    end    
end


for r=1:22
  
    % Estimates the Coefficient of Variation (standard deviation / average) 
    % of the HRF peak amplitude and peak latency in each group and 
    % condition per ROI
    total_CV_peak_CNT_Thr(1,r) = std(peak_per_subject_CNT_Thr(:,r),0,1)/mean(peak_per_subject_CNT_Thr(:,r),1);                  
    total_CV_peak_CNT_Sub(1,r) = std(peak_per_subject_CNT_Sub(:,r),0,1)/mean(peak_per_subject_CNT_Sub(:,r),1);
    total_CV_peak_T2DM_Thr(1,r) = std(peak_per_subject_T2DM_Thr(:,r),0,1)/mean(peak_per_subject_T2DM_Thr(:,r),1);
    total_CV_peak_T2DM_Sub(1,r) = std(peak_per_subject_T2DM_Sub(:,r),0,1)/mean(peak_per_subject_T2DM_Sub(:,r),1);

    total_CV_peak_latency_CNT_Thr(1,r) = std(peak_latency_per_subject_CNT_Thr(:,r),0,1)/mean(peak_latency_per_subject_CNT_Thr(:,r),1);
    total_CV_peak_latency_CNT_Sub(1,r) = std(peak_latency_per_subject_CNT_Sub(:,r),0,1)/mean(peak_latency_per_subject_CNT_Sub(:,r),1);
    total_CV_peak_latency_T2DM_Thr(1,r) = std(peak_latency_per_subject_T2DM_Thr(:,r),0,1)/mean(peak_latency_per_subject_T2DM_Thr(:,r),1);
    total_CV_peak_latency_T2DM_Sub(1,r) = std(peak_latency_per_subject_T2DM_Sub(:,r),0,1)/mean(peak_latency_per_subject_T2DM_Sub(:,r),1);     

end


% Plots the Coefficient of Variation graphs with data concerning each
% condition and group per ROI

% ««««««««««««««««««««««««« Peak Amplitude »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

figure
scatter(1:22,total_CV_peak_CNT_Thr,30,'MarkerEdgeColor',[0 .5 .5],...
                                      'MarkerFaceColor',[0 .7 .7],...
                                      'LineWidth',1)
hold on
scatter(1:22,total_CV_peak_T2DM_Thr,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                       'MarkerFaceColor',[0.87 0.19 0.19],...
                                       'LineWidth',1)
hold on
scatter(1:22,total_CV_peak_CNT_Sub,30,'MarkerEdgeColor',[0 .5 .5],...
                                      'LineWidth',0.75)
hold on
scatter(1:22,total_CV_peak_T2DM_Sub,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                       'LineWidth',0.75)
                            
set(gca, 'XTick', 1:22, 'XTickLabel', ROI_regions);                             % labels the x axis as the ROI names
xtickangle(45)                                                                  % rotates the labels by 45º
title('Coefficient of Variation - Peak Amplitude');
xlabel('ROI')
ylabel('Coefficient of Variation')  
legend ('CNT Thr', 'T2DM Thr', 'CNT Sub', 'T2DM Sub')


% «««««««««««««««««««««««««« Peak Latency »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

figure
scatter(1:22,total_CV_peak_latency_CNT_Thr,30,'MarkerEdgeColor',[0 .5 .5],...
                                              'MarkerFaceColor',[0 .7 .7],...
                                              'LineWidth',1)
hold on
scatter(1:22,total_CV_peak_latency_T2DM_Thr,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                               'MarkerFaceColor',[0.87 0.19 0.19],...
                                               'LineWidth',1)
hold on                            
scatter(1:22,total_CV_peak_latency_CNT_Sub,30,'MarkerEdgeColor',[0 .5 .5],...
                                              'LineWidth',0.75)
hold on
scatter(1:22,total_CV_peak_latency_T2DM_Sub,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                               'LineWidth',0.75)
                            
set(gca, 'XTick', 1:22, 'XTickLabel', ROI_regions);
xtickangle(45)
title('Coefficient of Variation - Peak Latency');
xlabel('ROI')
ylabel('Coefficient of Variation')  
legend ('CNT Thr', 'T2DM Thr', 'CNT Sub', 'T2DM Sub')