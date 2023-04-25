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


%% HRF PARAMETER TABLES PER CONDITION AND SET OF ROIS

% In this section, we assign the HRF parameters of each subject to its
% corresponding condition (Thr or Sub) and set of ROIs (positive and 
% negative signal change ROIs)


% Gets the set of subjects for each condition
subjects_Thr = table2array(covariates_HRFdata(1:2:end,1));
subjects_Sub = table2array(covariates_HRFdata(2:2:end,1));


% Assigns table headers and HRF parameters according to the subject's 
% condition and positive or negative signal change ROIs
psc_headers = headers(:,1:10*n_parameters);
psc_Thr_parameters = table2array(covariates_HRFdata(1:2:end,2:10*n_parameters+1));
psc_Sub_parameters = table2array(covariates_HRFdata(2:2:end,2:10*n_parameters+1));

nsc_headers = headers(:,10*n_parameters+1:end);
nsc_Thr_parameters = table2array(covariates_HRFdata(1:2:end,10*n_parameters+2:n_rois*n_parameters+1));
nsc_Sub_parameters = table2array(covariates_HRFdata(2:2:end,10*n_parameters+2:n_rois*n_parameters+1));


% Forms tables
psc_Thr_parameters_data = array2table(psc_Thr_parameters,'VariableNames', psc_headers);                   % table with the HRF parameters' data of each subject in the Thr condition in the positive signal change ROIs as well as its headers      
psc_Thr_subjects = array2table(subjects_Thr,'VariableNames', {'Subject'});                

psc_Sub_parameters_data = array2table(psc_Sub_parameters,'VariableNames', psc_headers);                   % table with the HRF parameters' data of each subject in the Sub condition in the positive signal change ROIs as well as its headers     
psc_Sub_subjects = array2table(subjects_Sub,'VariableNames', {'Subject'});   

nsc_Thr_parameters_data = array2table(nsc_Thr_parameters,'VariableNames', nsc_headers);                   % table with the HRF parameters' data of each subject in the Thr condition in the negative signal change ROIs as well as its headers      
nsc_Thr_subjects = array2table(subjects_Thr,'VariableNames', {'Subject'});                

nsc_Sub_parameters_data = array2table(nsc_Sub_parameters,'VariableNames', nsc_headers);                   % table with the HRF parameters' data of each subject in the Sub condition in the negative signal change ROIs as well as its headers      
nsc_Sub_subjects = array2table(subjects_Sub,'VariableNames', {'Subject'}); 


% Merges tables
psc_Thr_parameters_data = [psc_Thr_subjects psc_Thr_parameters_data];
psc_Sub_parameters_data = [psc_Sub_subjects psc_Sub_parameters_data];

nsc_Thr_parameters_data = [nsc_Thr_subjects nsc_Thr_parameters_data];
nsc_Sub_parameters_data = [nsc_Sub_subjects nsc_Sub_parameters_data];


%% HRF PARAMETER TABLES PER SET OF ROIS

% In this section, we estimate the HRF parameters per average condition
% across subjects (average of Thr and Sub condition data) and assign them 
% to each set of ROIs (positive and negative signal change ROIs)


psc_parameters_CNT = zeros(n_CNT, n_parameters*10);
psc_parameters_T2DM = zeros(n_T2DM, n_parameters*10);

nsc_parameters_CNT = zeros(n_CNT, n_parameters*12);
nsc_parameters_T2DM = zeros(n_T2DM, n_parameters*12);

subjects = cell(n_subjects,1);


% Estimates the HRF parameters per average condition across CNT subjects in
% both sets of ROIs
for c=1:n_CNT
    psc_parameters_CNT(c,:) = mean([psc_Thr_parameters(c,:); psc_Sub_parameters(c,:)],1);
    nsc_parameters_CNT(c,:) = mean([nsc_Thr_parameters(c,:); nsc_Sub_parameters(c,:)],1);
end

% Estimates the HRF parameters per average condition across T2DM subjects 
% in both sets of ROIs
for t=1:n_T2DM
    psc_parameters_T2DM(t,:) = mean([psc_Thr_parameters(n_CNT+t,:); psc_Sub_parameters(n_CNT+t,:)],1);
    nsc_parameters_T2DM(t,:) = mean([nsc_Thr_parameters(n_CNT+t,:); nsc_Sub_parameters(n_CNT+t,:)],1);
end


% Identifies each subject from each group *without* conditions
for l=1:n_subjects
    subjects(l,:) = regexprep(subjects_Thr{l}, ' Thr', '');                         % removes the ' Thr' string of each element of that cell array
end


% Forms tables
psc_parameters_data = array2table([psc_parameters_CNT; psc_parameters_T2DM],'VariableNames', psc_headers);                   % table with the HRF parameters of the average condition across subjects in the positive signal change ROIs as well as its headers     
psc_subjects = array2table(subjects,'VariableNames', {'Subject'});

nsc_parameters_data = array2table([nsc_parameters_CNT; nsc_parameters_T2DM],'VariableNames', nsc_headers);                   % table with the HRF parameters of the average condition across subjects in the negative signal change ROIs as well as its headers      
nsc_subjects = array2table(subjects,'VariableNames', {'Subject'});
  

% Merges tables
psc_parameters_data = [psc_subjects psc_parameters_data];
nsc_parameters_data = [nsc_subjects nsc_parameters_data];


%% AVERAGE HRF PARAMETER TABLES PER SET OF ROIS

% In this section, we estimate the average HRF parameters per average 
% condition and sets of ROIs across subjects


% Gets the HRF parameters per average condition across subjects of each 
% group in each set of ROIs
psc_parameters = [psc_parameters_CNT; psc_parameters_T2DM];
nsc_parameters = [nsc_parameters_CNT; nsc_parameters_T2DM];


avg_psc_parameters = zeros(n_subjects,n_parameters);
avg_nsc_parameters = zeros(n_subjects,n_parameters);

avg_headers = cell(1,n_parameters);


for par=1:n_parameters
    
    % Estimates the average HRF parameters per average condition and sets 
    % of ROIs across subjects
    avg_psc_parameters(:,par) = mean(psc_parameters(:,par:n_parameters:end),2);
    avg_nsc_parameters(:,par) = mean(nsc_parameters(:,par:n_parameters:end),2);
    
    % Setting automatically the table headers (name of the HRF parameters)
    avg_headers{par} = char(strrep(strip(HRF_parameters(par)),' ','_'));                          % removes the first space of this cell array's elements, replaces the remaining spaces by _ and converts to char
end


% Forms tables
avg_psc_parameters_data = array2table(avg_psc_parameters,'VariableNames', avg_headers);                   % table with the average HRF parameters per average condition and positive signal change ROIs across subjects as well as its headers       
avg_psc_subjects = array2table(subjects,'VariableNames', {'Subject'});

avg_nsc_parameters_data = array2table(avg_nsc_parameters,'VariableNames', avg_headers);                   % table with the average HRF parameters per average condition and negative signal change ROIs across subjects as well as its headers       
avg_nsc_subjects = array2table(subjects,'VariableNames', {'Subject'});
  

% Merges tables
avg_psc_parameters_data = [avg_psc_subjects avg_psc_parameters_data];
avg_nsc_parameters_data = [avg_nsc_subjects avg_nsc_parameters_data];


%% AVERAGE HRF PARAMETER TABLES PER SET OF ROIS AND GROUP

% In this section, we estimate the average and standard deviation of the 
% HRF parameters per average condition, sets of ROIs and group (T2DM or
% CNT)


% Identifies each group of subjects
group = {'CNT', 'T2DM'}';


% Gets the average of the HRF parameters per average condition, sets of
% ROIs and group in each set of ROIs and group
avg_psc_parameters_CNT = mean(avg_psc_parameters(1:n_CNT,:),1);
avg_psc_parameters_T2DM = mean(avg_psc_parameters(n_CNT+1:end,:),1);

avg_nsc_parameters_CNT = mean(avg_nsc_parameters(1:n_CNT,:),1);
avg_nsc_parameters_T2DM = mean(avg_psc_parameters(n_CNT+1:end,:),1);


% Gets the standard deviation of the HRF parameters per average condition, 
% sets of ROIs and group in each set of ROIs and group
std_psc_parameters_CNT = std(avg_psc_parameters(1:n_CNT,:),0,1);
std_psc_parameters_T2DM = std(avg_psc_parameters(n_CNT+1:end,:),0,1);

std_nsc_parameters_CNT = std(avg_nsc_parameters(1:n_CNT,:),0,1);
std_nsc_parameters_T2DM = std(avg_psc_parameters(n_CNT+1:end,:),0,1);


% Forms tables
group_avg_psc_parameters_data = array2table([avg_psc_parameters_CNT; avg_psc_parameters_T2DM],'VariableNames', avg_headers);                   % table with the average of the HRF parameters per average condition, positive signal change ROIs and group as well as its headers       
group_avg_psc = array2table(group,'VariableNames', {'Group'});

group_avg_nsc_parameters_data = array2table([avg_nsc_parameters_CNT; avg_nsc_parameters_T2DM],'VariableNames', avg_headers);                   % table with the average of the HRF parameters per average condition, negative signal change ROIs and group as well as its headers       
group_avg_nsc = array2table(group,'VariableNames', {'Group'});
  
group_std_psc_parameters_data = array2table([std_psc_parameters_CNT; std_psc_parameters_T2DM],'VariableNames', avg_headers);                   % table with the standard deviation of the HRF parameters per average condition, positive signal change ROIs and group as well as its headers  
group_std_psc = array2table(group,'VariableNames', {'Group'});

group_std_nsc_parameters_data = array2table([std_nsc_parameters_CNT; std_nsc_parameters_T2DM],'VariableNames', avg_headers);                   % table with the standard deviation of the HRF parameters per average condition, negative signal change ROIs and group as well as its headers
group_std_nsc = array2table(group,'VariableNames', {'Group'});


% Merges tables
group_avg_psc_parameters_data = [group_avg_psc group_avg_psc_parameters_data];
group_avg_nsc_parameters_data = [group_avg_nsc group_avg_nsc_parameters_data];

group_std_psc_parameters_data = [group_std_psc group_std_psc_parameters_data];
group_std_nsc_parameters_data = [group_std_nsc group_std_nsc_parameters_data];


%% INFORMATIVE TABLES PER CONDITION

% In this section, we create three different tables with data regarding
% each subject from each group in three different conditions: Thr, Sub and
% the average of both

% The tables have the following disposition:
% Group || Age || HRF parameters in all ROIs || average HRF parameters per set of ROIs       
%  (1)  || (1) ||           (22x11)          ||                (11x2)
%
% 1. The HRF parameter columns show the positive signal change ROIs first
% and then the negative ones.
% 2. Data from CNT subjects is shown first, and then data regarding T2DM 
% subjects.


% Deletes original database information regarding subject no.66, who not 
% included in the image analysis due to problems in co-registration
diamarker_summary(66,:) = [];

% Assigns each subject's age at the time of the scan according to its group
ages_CNT = table2array(diamarker_summary(n_T2DM+1:end,3));
ages_T2DM = table2array(diamarker_summary(1:n_T2DM,3));


avg_psc_parameters_Thr = zeros(n_subjects,n_parameters);
avg_psc_parameters_Sub = zeros(n_subjects,n_parameters);
avg_nsc_parameters_Thr = zeros(n_subjects,n_parameters);
avg_nsc_parameters_Sub = zeros(n_subjects,n_parameters);

avg_psc_headers = cell(1,n_parameters);
avg_nsc_headers = cell(1,n_parameters);


for par=1:n_parameters
    
    % Average of the HRF parameters per sets of ROIs in the Thr and Sub 
    % conditions across subjects
    avg_psc_parameters_Thr(:,par) = mean(psc_Thr_parameters(:,par:n_parameters:end),2);    
    avg_psc_parameters_Sub(:,par) = mean(psc_Sub_parameters(:,par:n_parameters:end),2);
    
    avg_nsc_parameters_Thr(:,par) = mean(nsc_Thr_parameters(:,par:n_parameters:end),2);
    avg_nsc_parameters_Sub(:,par) = mean(nsc_Sub_parameters(:,par:n_parameters:end),2);
    
    
    % Sets automatically the table headers for the average HRF parameters  
    % per set of ROIs across subjects    
    avg_psc_headers{par} = strcat('avg_psc_', avg_headers{par});
    avg_nsc_headers{par} = strcat('avg_nsc_', avg_headers{par});
end


% Forms tables
subject_data = array2table(subjects,'VariableNames', {'Subject'});                                                  % table with the subjects' number per group - the first column of the three tables, which allows to identify data - as well as its headers
ages_data = array2table([ages_CNT; ages_T2DM],'VariableNames', {'Age'});                                            % table with each subjects' age as well as its headers

avg_psc_parameters_Thr_data = array2table(avg_psc_parameters_Thr,'VariableNames', avg_psc_headers);                 % table with each subject's average of the HRF parameters in the positive signal change ROIs in the Thr condition as well as its headers
avg_psc_parameters_Sub_data = array2table(avg_psc_parameters_Sub,'VariableNames', avg_psc_headers);                 % table with each subject's average of the HRF parameters in the positive signal change ROIs in the Sub condition as well as its headers
informative_avg_psc_parameters_data = array2table(avg_psc_parameters,'VariableNames', avg_psc_headers);                   % table with each subject's average of the HRF parameters in the positive signal change ROIs in the average condition as well as its headers       

avg_nsc_parameters_Thr_data = array2table(avg_nsc_parameters_Thr,'VariableNames', avg_nsc_headers);                 % table with each subject's average of the HRF parameters in the negative signal change ROIs in the Thr condition as well as its headers
avg_nsc_parameters_Sub_data = array2table(avg_nsc_parameters_Sub,'VariableNames', avg_nsc_headers);                 % table with each subject's average of the HRF parameters in the negative signal change ROIs in the Sub condition as well as its headers
informative_avg_nsc_parameters_data = array2table(avg_nsc_parameters,'VariableNames', avg_nsc_headers);                   % table with each subject's average of the HRF parameters in the positive signal change ROIs in the average condition as well as its headers       


% Merges tables
informative_table_avg = [subject_data ages_data psc_parameters_data(:,2:end) nsc_parameters_data(:,2:end) informative_avg_psc_parameters_data(:,:) informative_avg_nsc_parameters_data(:,:)];
informative_table_Thr = [subject_data ages_data psc_Thr_parameters_data(:,2:end) nsc_Thr_parameters_data(:,2:end) avg_psc_parameters_Thr_data(:,:) avg_nsc_parameters_Thr_data(:,:)];
informative_table_Sub = [subject_data ages_data psc_Sub_parameters_data(:,2:end) nsc_Sub_parameters_data(:,2:end) avg_psc_parameters_Sub_data(:,:) avg_nsc_parameters_Sub_data(:,:)];
