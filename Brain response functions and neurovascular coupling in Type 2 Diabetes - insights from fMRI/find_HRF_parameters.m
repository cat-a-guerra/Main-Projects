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


%% CONSTANTS:

n_subjects = 141;
n_T2DM = 64;
n_CNT = 77;
n_rois = 22;
n_parameters = 11;
n_points = 8;
time = (0:n_points-1)* 2.5;                              % volumes per subject * TR (2.5 s)

ROI_regions = {'L IPL BA40', 'L Insula BA13', 'L Precuneus BA7', 'R IFG BA9', 'R MFG BA8', 'R MFG BA46', 'R MT BA19', ... 
               'R SFG BA6', 'R SPL BA7', 'R V2 BA18', 'L AC BA32', 'L CG BA31', 'L PC BA30', 'L PrcL BA5', 'L PrhG BA36', ...
               'R CG BA24', 'R Insula BA13', 'R PC BA23', 'R PrecG BA4', 'R PrecG BA43', 'R PrimSens BA1', 'R STG BA39'};


%% INITIALIZATION:

% In this section, we get every mean value of beta per subject in each
% group and condition for each ROI


thr_both_sides = zeros(n_points,2);                         % 2 since each column represents a stimulation side - left and right
sub_both_sides = zeros(n_points,2);

mean_per_subject_CNT_Thr = zeros(n_points,n_CNT*n_rois);
mean_per_subject_CNT_Sub = zeros(n_points,n_CNT*n_rois);
mean_per_subject_T2DM_Thr = zeros(n_points,n_T2DM*n_rois);
mean_per_subject_T2DM_Sub = zeros(n_points,n_T2DM*n_rois);


for r=2:n_rois+1
    for i=1:n_subjects
        for j=1:2
            thr_both_sides(:,j) = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));    % each column is a side in each subject 
            sub_both_sides(:,j) = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));
        end

        if i<65                                                                                 % until subject 65 --> T2DM 
            mean_per_subject_T2DM_Thr(:,i+n_T2DM*(r-2)) = mean(thr_both_sides,2);               % mean of the beta values from each side per datapoint of T2DM subjects in the Thr condition in a ROI
            mean_per_subject_T2DM_Sub(:,i+n_T2DM*(r-2)) = mean(sub_both_sides,2);
        else                                                                                    % after subject 65 --> CNT
            mean_per_subject_CNT_Thr(:,i-n_T2DM+n_CNT*(r-2)) = mean(thr_both_sides,2);          % mean of the beta values from each side per datapoint of CNT subjects in the Thr condition in a ROI 
            mean_per_subject_CNT_Sub(:,i-n_T2DM+n_CNT*(r-2)) = mean(sub_both_sides,2);
        end
    
        % Deletes to not duplicate or interfere with data
        thr_both_sides = zeros(n_points,2);
        sub_both_sides = zeros(n_points,2);
    end  
end


%% PARAMETER ESTIMATION:

% In this section, we estimate and store the HRF parameters of interest 
% per subject in each group and condition for each ROI


total_parameters_thr_CNT = zeros(n_CNT,n_parameters*n_rois);
total_parameters_sub_CNT = zeros(n_CNT,n_parameters*n_rois);
total_parameters_thr_T2DM = zeros(n_T2DM,n_parameters*n_rois);
total_parameters_sub_T2DM = zeros(n_T2DM,n_parameters*n_rois);


% ********************************* CNT **********************************
        
for c = 1:n_CNT
    for r = 1:n_rois                                                                                                                    % 22 ROI
        
        % Stores the HRF curves per CNT subjects in each condition for each 
        % ROI 
        data.CNT.subject(c).ROI(r).Threshold.HRFCurve = mean_per_subject_CNT_Thr(1:n_points,r+n_rois*(c-1));
        data.CNT.subject(c).ROI(r).Submaximum.HRFCurve = mean_per_subject_CNT_Sub(1:n_points,r+n_rois*(c-1));

        
        % -------------------- THRESHOLD CONDITION ----------------------

        [peak,volume] = max(data.CNT.subject(c).ROI(r).Threshold.HRFCurve(2:6));                                                        % peak: maximum within the sectioned time range (2.5 to 12.5 s - 2nd to 6th index since indexing starts at 1)
        latency = volume * 2.5;                                                                                                         % peak latency: time to peak (volume where the maximum takes place * TR)
        slope_to_peak = minus(peak,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1))/latency;                                           % relative slope to peak: (peak-1st point)/peak latency
        total_area = trapz((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points)- ... 
                            min(data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points)));                                            % area under the curve (trapezoidal method)
        pos_area = positive_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points));                         % area of the positive sections of the HRF curve
        first_pos_area = positive_area((0:2)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:3));                                   % area of the first positive section of the HRF curve (from 0 to 5 s)
        second_pos_area = positive_area((2:4)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(3:5));                                  % area of the second section of the HRF curve (from 5 to 10 s)
        third_pos_area = positive_area((4:7)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(5:n_points));                            % area of the third section of the HRF curve (from 10 to 17.5 s)
        neg_area = negative_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points));                         % area of the negative sections of the HRF curve 
        initial_dip_area = negative_area((0:volume)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:(volume+1)));                   % area of the initial dip of the HRF curve
        undershoot_area = negative_area((volume:7)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve((volume+1):n_points));             % area of the undershoot of the HRF curve
                
        % Gathers and stores every HRF parameter per CNT subjects in the 
        % Thr condition for each ROI
        data.CNT.subject(c).ROI(r).Threshold.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];
        total_parameters_thr_CNT(c,1+n_parameters*(r-1):n_parameters+n_parameters*(r-1)) = data.CNT.subject(c).ROI(r).Threshold.HRFParameters;
       
        
        % -------------------- SUBMAXIMUM CONDITION ---------------------

        [peak,volume] = max(data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(2:6));                             
        latency = volume * 2.5;                                                                             
        slope_to_peak = minus(peak,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1))/latency;                                  
        total_area = trapz((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points)- ... 
                            min(data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points)));    
        pos_area = positive_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points));      
        first_pos_area = positive_area((0:2)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:3));                   
        second_pos_area = positive_area((2:4)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(3:5));                  
        third_pos_area = positive_area((4:7)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(5:n_points));            
        neg_area = negative_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points));      
        initial_dip_area = negative_area((0:volume)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:(volume+1)));                   
        undershoot_area = negative_area((volume:7)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve((volume+1):n_points));           
        
        % Gathers and stores every HRF parameter per CNT subjects in the 
        % Sub condition for each ROI
        data.CNT.subject(c).ROI(r).Submaximum.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];       
        total_parameters_sub_CNT(c,1+n_parameters*(r-1):n_parameters+n_parameters*(r-1)) = data.CNT.subject(c).ROI(r).Submaximum.HRFParameters;
    end
    
    % Average HRF parameters per condition in each ROI across subjects
    avg_parameters_thr_CNT = mean(total_parameters_thr_CNT,1);
    avg_parameters_sub_CNT = mean(total_parameters_sub_CNT,1);   
end


% ******************************** T2DM **********************************


for t=1:n_T2DM
    for r = 1:n_rois
        
        % Stores the HRF curves per T2DM subjects in each condition for 
        % each ROI 
        data.T2DM.subject(t).ROI(r).Threshold.HRFCurve = mean_per_subject_T2DM_Thr(1:n_points,r+n_rois*(t-1));
        data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve = mean_per_subject_T2DM_Sub(1:n_points,r+n_rois*(t-1));
        
        % -------------------- THRESHOLD CONDITION ----------------------

        [peak,volume] = max(data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(2:6));                              
        latency = volume * 2.5;                                                                             
        slope_to_peak = minus(peak,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1))/latency;
        total_area = trapz((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points)- ... 
                            min(data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points)));  
        pos_area = positive_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points));
        first_pos_area = positive_area((0:2)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:3));                   
        second_pos_area = positive_area((2:4)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(3:5));                  
        third_pos_area = positive_area((4:7)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(5:n_points));            
        neg_area = negative_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points));         
        initial_dip_area = negative_area((0:volume)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:(volume+1)));                  
        undershoot_area = negative_area((volume:7)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve((volume+1):n_points));           
        
        % Gathers and stores every HRF parameter per T2DM subjects in the 
        % Thr condition for each ROI
        data.T2DM.subject(t).ROI(r).Threshold.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];        
        total_parameters_thr_T2DM(t,1+n_parameters*(r-1):n_parameters+n_parameters*(r-1)) = data.T2DM.subject(t).ROI(r).Threshold.HRFParameters;

        
        % -------------------- SUBMAXIMUM CONDITION ---------------------

        [peak,volume] = max(data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(2:6));                             
        latency = volume * 2.5;                                                                             
        slope_to_peak = minus(peak,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1))/latency;                                  
        total_area = trapz((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points)- ... 
                            min(data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points)));                         
        pos_area = positive_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points));      
        first_pos_area = positive_area((0:2)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:3));                   
        second_pos_area = positive_area((2:4)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(3:5));                  
        third_pos_area = positive_area((4:7)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(5:n_points)); 
        neg_area = negative_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points));      
        initial_dip_area = negative_area((0:volume)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:(volume+1)));                   
        undershoot_area = negative_area((volume:7)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve((volume+1):n_points));       
        
        % Gathers and stores every HRF parameter per T2DM subjects in the 
        % Sub condition for each ROI 
        data.T2DM.subject(t).ROI(r).Submaximum.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];
        total_parameters_sub_T2DM(t,1+n_parameters*(r-1):n_parameters+n_parameters*(r-1)) = data.CNT.subject(t).ROI(r).Submaximum.HRFParameters;

    end
    
    % Average HRF parameters per condition in each ROI across subjects
    avg_parameters_thr_T2DM = mean(total_parameters_thr_T2DM,1);
    avg_parameters_sub_T2DM = mean(total_parameters_sub_T2DM,1);
end


%% GRAND AVERAGE / MEDIAN ANALYSIS - HRF PARAMETERS:

% In this section, we do an overall average / median of the HRF parameters
% of interest of each subject across ROIs for each condition and group


psc_avg_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
psc_avg_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
psc_avg_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
psc_avg_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);

psc_std_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
psc_std_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
psc_std_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
psc_std_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);

psc_median_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
psc_median_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
psc_median_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
psc_median_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);

psc_iqr_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
psc_iqr_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
psc_iqr_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
psc_iqr_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);


nsc_avg_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
nsc_avg_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
nsc_avg_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
nsc_avg_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);

nsc_std_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
nsc_std_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
nsc_std_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
nsc_std_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);

nsc_median_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
nsc_median_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
nsc_median_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
nsc_median_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);

nsc_iqr_parameters_CNT_Thr = zeros(n_CNT,n_parameters);                                     
nsc_iqr_parameters_CNT_Sub = zeros(n_CNT,n_parameters);                                     
nsc_iqr_parameters_T2DM_Thr = zeros(n_T2DM,n_parameters);                                     
nsc_iqr_parameters_T2DM_Sub = zeros(n_T2DM,n_parameters);


% Overall average and standard deviation / median and interquartile range
% of the HRF parameters of interest
        
for par = 1:n_parameters
    

    % ***************** Positive signal change ROIs *********************
    
    psc_avg_parameters_CNT_Thr(:,par) = mean(total_parameters_thr_CNT(:,par:n_parameters:n_parameters*10),2);                                     
    psc_avg_parameters_CNT_Sub(:,par) = mean(total_parameters_sub_CNT(:,par:n_parameters:n_parameters*10),2);                                     
    psc_avg_parameters_T2DM_Thr(:,par) = mean(total_parameters_thr_T2DM(:,par:n_parameters:n_parameters*10),2);                                     
    psc_avg_parameters_T2DM_Sub(:,par) = mean(total_parameters_sub_T2DM(:,par:n_parameters:n_parameters*10),2);  

    psc_std_parameters_CNT_Thr(:,par) = std(total_parameters_thr_CNT(:,par:n_parameters:n_parameters*10),0,2);                                     
    psc_std_parameters_CNT_Sub(:,par) = std(total_parameters_sub_CNT(:,par:n_parameters:n_parameters*10),0,2);                                     
    psc_std_parameters_T2DM_Thr(:,par) = std(total_parameters_thr_T2DM(:,par:n_parameters:n_parameters*10),0,2);                                     
    psc_std_parameters_T2DM_Sub(:,par) = std(total_parameters_sub_T2DM(:,par:n_parameters:n_parameters*10),0,2); 
    
    psc_median_parameters_CNT_Thr(:,par) = median(total_parameters_thr_CNT(:,par:n_parameters:n_parameters*10),2);                                     
    psc_median_parameters_CNT_Sub(:,par) = median(total_parameters_sub_CNT(:,par:n_parameters:n_parameters*10),2);                                     
    psc_median_parameters_T2DM_Thr(:,par) = median(total_parameters_thr_T2DM(:,par:n_parameters:n_parameters*10),2);                                     
    psc_median_parameters_T2DM_Sub(:,par) = median(total_parameters_sub_T2DM(:,par:n_parameters:n_parameters*10),2);
    
    psc_iqr_parameters_CNT_Thr(:,par) = iqr(total_parameters_thr_CNT(:,par:n_parameters:n_parameters*10),2);                                     
    psc_iqr_parameters_CNT_Sub(:,par) = iqr(total_parameters_sub_CNT(:,par:n_parameters:n_parameters*10),2);                                     
    psc_iqr_parameters_T2DM_Thr(:,par) = iqr(total_parameters_thr_T2DM(:,par:n_parameters:n_parameters*10),2);                                     
    psc_iqr_parameters_T2DM_Sub(:,par) = iqr(total_parameters_sub_T2DM(:,par:n_parameters:n_parameters*10),2);
    

    % ***************** Negative signal change ROIs *********************

    nsc_avg_parameters_CNT_Thr(:,par) = mean(total_parameters_thr_CNT(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_avg_parameters_CNT_Sub(:,par) = mean(total_parameters_sub_CNT(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_avg_parameters_T2DM_Thr(:,par) = mean(total_parameters_thr_T2DM(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_avg_parameters_T2DM_Sub(:,par) = mean(total_parameters_sub_T2DM(:,n_parameters*10+par:n_parameters:end),2); 

    nsc_std_parameters_CNT_Thr(:,par) = std(total_parameters_thr_CNT(:,n_parameters*10+par:n_parameters:end),0,2);                                     
    nsc_std_parameters_CNT_Sub(:,par) = std(total_parameters_sub_CNT(:,n_parameters*10+par:n_parameters:end),0,2);                                     
    nsc_std_parameters_T2DM_Thr(:,par) = std(total_parameters_thr_T2DM(:,n_parameters*10+par:n_parameters:end),0,2);                                     
    nsc_std_parameters_T2DM_Sub(:,par) = std(total_parameters_sub_T2DM(:,n_parameters*10+par:n_parameters:end),0,2);

    nsc_median_parameters_CNT_Thr(:,par) = median(total_parameters_thr_CNT(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_median_parameters_CNT_Sub(:,par) = median(total_parameters_sub_CNT(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_median_parameters_T2DM_Thr(:,par) = median(total_parameters_thr_T2DM(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_median_parameters_T2DM_Sub(:,par) = median(total_parameters_sub_T2DM(:,n_parameters*10+par:n_parameters:end),2);
    
    nsc_iqr_parameters_CNT_Thr(:,par) = iqr(total_parameters_thr_CNT(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_iqr_parameters_CNT_Sub(:,par) = iqr(total_parameters_sub_CNT(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_iqr_parameters_T2DM_Thr(:,par) = iqr(total_parameters_thr_T2DM(:,n_parameters*10+par:n_parameters:end),2);                                     
    nsc_iqr_parameters_T2DM_Sub(:,par) = iqr(total_parameters_sub_T2DM(:,n_parameters*10+par:n_parameters:end),2);
    
    
end


%% FINAL STORAGING:

% In this section, we create tables regarding the HRF parameters of 
% interest and its information per group and condition in each subject

% Storaging order:        ROI1 p.1 || ROI1 p.2 || ... || ROI22 p.11
%                                       1. CNT Thr
%                                       2. CNT Sub
%                                       3. T2DM Thr
%                                       4. T2DM Sub


% ************************** Individual data ****************************

% «««««««««««««««««««««««««««« Table data »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

group = cell(n_subjects*2,1);                                                               

CNT_total_parameters = zeros(n_CNT*2,n_parameters*n_rois);                                      % gets the HRF parameters of all CNT subjects in both conditions and in every ROI
T2DM_total_parameters = zeros(n_T2DM*2,n_parameters*n_rois);


for c = 1:n_CNT
    
    % Gathers every parameter of interest of every CNT subject in each condition and in every ROI
    for r = 1:n_rois
        CNT_total_parameters((c*2)-1:c*2,1+n_parameters*(r-1):n_parameters+n_parameters*(r-1)) = [data.CNT.subject(c).ROI(r).Threshold.HRFParameters; data.CNT.subject(c).ROI(r).Submaximum.HRFParameters];
    end
    
    % Identifies each subject in each group according to its condition
    for l = 1:2
        if l==1
            group{l+2*(c-1)} = strcat('Subject', {' '}, num2str(c), {' '}, 'CNT', {' '}, 'Thr');
        else
            group{l+2*(c-1)} = strcat('Subject', {' '}, num2str(c), {' '}, 'CNT', {' '}, 'Sub');
        end
    end
end


for t = 1:n_T2DM
    
    % Gathers every parameter of interest of every T2DM subject in each condition and in every ROI
    for r = 1:n_rois
        T2DM_total_parameters((t*2)-1:t*2,1+n_parameters*(r-1):n_parameters+n_parameters*(r-1)) = [data.T2DM.subject(t).ROI(r).Threshold.HRFParameters; data.T2DM.subject(t).ROI(r).Submaximum.HRFParameters];
    end
    
    % Identifies each subject in each group according to its condition
    for l = 1:2
        if l==1
            group{(n_CNT*2)+(l+2*(t-1))} = strcat('Subject', {' '}, num2str(t), {' '}, 'T2DM', {' '}, 'Thr');
        else
            group{(n_CNT*2)+(l+2*(t-1))} = strcat('Subject', {' '}, num2str(t), {' '}, 'T2DM', {' '}, 'Sub');
        end
    end
end

% Gathers every parameter of interest of all subjects in each condition and in every ROI
HRF_data = [CNT_total_parameters; T2DM_total_parameters];


% ««««««««««««««««««««««««««« Table setting »»»»»»»»»»»»»»»»»»»»»»»»»»»»»

headers = cell(1,n_rois*n_parameters);                                                                                                                                                                     % cell array for the table headings
HRF_parameters = {' peak amplitude', ' peak latency', ' relative slope to peak', ' AUC', ' APCS', ' fPCSA', ' sPCSA', ' tPCSA', ' ANCS', ' initial dip area', ' undershoot area'};         % parameters of interest

for r = 1:n_rois

    % Automatically sets the table headers
    for par = 1:length(HRF_parameters)
        headers{par+n_parameters*(r-1)} = char(strrep(strcat(ROI_regions(r),HRF_parameters(par)),' ','_'));   % concatenates strings, replaces spaces by _ and converts to char
    end
end

% Forms tables
covariates_HRFdata = array2table(HRF_data,'VariableNames',headers);                                 % table with the HRF data from every subject in each condition and in every ROI as well as its headers 
covariates_group = array2table(group,'VariableNames',{'Subject'});                                  % table with subject number, group and condition

% Merges tables
covariates_HRFdata = [covariates_group covariates_HRFdata];

% Saves data
save 'ROI_regions.mat' ROI_regions;
save 'headers.mat' headers;
save 'covariates_HRFdata.mat' covariates_HRFdata;
save 'HRF_parameters.mat' HRF_parameters;


% *************************** Average data *****************************

% «««««««««««««««««««««««««««« Table data »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

avg_HRF_data = [avg_parameters_thr_CNT; avg_parameters_sub_CNT; avg_parameters_thr_T2DM; avg_parameters_sub_T2DM];

% ««««««««««««««««««««««««««« Table setting »»»»»»»»»»»»»»»»»»»»»»»»»»»»»

avg_headers = headers;
avg_group = {'CNT','CNT','T2DM','T2DM'}';

% Forms tables
avg_covariates_HRFdata = array2table(avg_HRF_data,'VariableNames', avg_headers);                   % table with the average HRF data in each condition and in every ROI as well as its headers       
avg_covariates_group = array2table(avg_group,'VariableNames', {'Group'});                

% Merges tables
avg_covariates_HRFdata = [avg_covariates_group avg_covariates_HRFdata];

% Saves data
save 'avg_covariates_HRFdata.mat' avg_covariates_HRFdata;


% *********************** Grand Average/Median *************************

grand_headers = cell(1,n_parameters);

subjects_CNT = cell(n_CNT,1);
subjects_T2DM = cell(n_T2DM,1); 


% Setting automatically the table headers (name of the HRF parameters)
for par = 1:n_parameters
    grand_headers{par} = char(strrep(strip(HRF_parameters(par)),' ','_'));                          % removes the first space of HRF parameters, replaces the remaining spaces by _ and converts to char
end


% Identifies each subject from each group
for c = 1:n_CNT
    subjects_CNT{c} = strcat('Subject', {' '}, num2str(c));   
end

for t = 1:n_T2DM
    subjects_T2DM{t} = strcat('Subject', {' '}, num2str(t));   
end


% Sets the first column of the table, which corresponds to subjects from 
% each group
subject_CNT_column = array2table(subjects_CNT,'VariableNames',{'Subject'});
subject_T2DM_column = array2table(subjects_T2DM,'VariableNames',{'Subject'});


% ««««««««««««««««««««« Positive Section Curve ROIs »»»»»»»»»»»»»»»»»»»»»»»

% Forms tables
psc_avg_parameters_CNT_Thr_data = array2table(psc_avg_parameters_CNT_Thr,'VariableNames',grand_headers);                 
psc_avg_parameters_CNT_Sub_data = array2table(psc_avg_parameters_CNT_Sub,'VariableNames',grand_headers);
psc_avg_parameters_T2DM_Thr_data = array2table(psc_avg_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
psc_avg_parameters_T2DM_Sub_data = array2table(psc_avg_parameters_T2DM_Sub,'VariableNames',grand_headers);                 

psc_std_parameters_CNT_Thr_data = array2table(psc_std_parameters_CNT_Thr,'VariableNames',grand_headers);  
psc_std_parameters_CNT_Sub_data = array2table(psc_std_parameters_CNT_Sub,'VariableNames',grand_headers);                 
psc_std_parameters_T2DM_Thr_data = array2table(psc_std_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
psc_std_parameters_T2DM_Sub_data = array2table(psc_std_parameters_T2DM_Sub,'VariableNames',grand_headers); 

psc_median_parameters_CNT_Thr_data = array2table(psc_median_parameters_CNT_Thr,'VariableNames',grand_headers);                 
psc_median_parameters_CNT_Sub_data = array2table(psc_median_parameters_CNT_Sub,'VariableNames',grand_headers);
psc_median_parameters_T2DM_Thr_data = array2table(psc_median_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
psc_median_parameters_T2DM_Sub_data = array2table(psc_median_parameters_T2DM_Sub,'VariableNames',grand_headers);                 

psc_iqr_parameters_CNT_Thr_data = array2table(psc_iqr_parameters_CNT_Thr,'VariableNames',grand_headers);  
psc_iqr_parameters_CNT_Sub_data = array2table(psc_iqr_parameters_CNT_Sub,'VariableNames',grand_headers);                 
psc_iqr_parameters_T2DM_Thr_data = array2table(psc_iqr_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
psc_iqr_parameters_T2DM_Sub_data = array2table(psc_iqr_parameters_T2DM_Sub,'VariableNames',grand_headers);


% Merges tables
psc_avg_parameters_CNT_Thr_data = [subject_CNT_column psc_avg_parameters_CNT_Thr_data];
psc_avg_parameters_CNT_Sub_data = [subject_CNT_column psc_avg_parameters_CNT_Sub_data];
psc_avg_parameters_T2DM_Thr_data = [subject_T2DM_column psc_avg_parameters_T2DM_Thr_data];                 
psc_avg_parameters_T2DM_Sub_data = [subject_T2DM_column psc_avg_parameters_T2DM_Sub_data];   

psc_std_parameters_CNT_Thr_data = [subject_CNT_column psc_std_parameters_CNT_Thr_data];  
psc_std_parameters_CNT_Sub_data = [subject_CNT_column psc_std_parameters_CNT_Sub_data];                 
psc_std_parameters_T2DM_Thr_data = [subject_T2DM_column psc_std_parameters_T2DM_Thr_data];                 
psc_std_parameters_T2DM_Sub_data = [subject_T2DM_column psc_std_parameters_T2DM_Sub_data]; 

psc_median_parameters_CNT_Thr_data = [subject_CNT_column psc_median_parameters_CNT_Thr_data];                 
psc_median_parameters_CNT_Sub_data = [subject_CNT_column psc_median_parameters_CNT_Sub_data];
psc_median_parameters_T2DM_Thr_data = [subject_T2DM_column psc_median_parameters_T2DM_Thr_data];                 
psc_median_parameters_T2DM_Sub_data = [subject_T2DM_column psc_median_parameters_T2DM_Sub_data];                 

psc_iqr_parameters_CNT_Thr_data = [subject_CNT_column psc_iqr_parameters_CNT_Thr_data];  
psc_iqr_parameters_CNT_Sub_data = [subject_CNT_column psc_iqr_parameters_CNT_Sub_data];                 
psc_iqr_parameters_T2DM_Thr_data = [subject_T2DM_column psc_iqr_parameters_T2DM_Thr_data];                 
psc_iqr_parameters_T2DM_Sub_data = [subject_T2DM_column psc_iqr_parameters_T2DM_Sub_data];


% ««««««««««««««««««««« Negative Section Curve ROIs »»»»»»»»»»»»»»»»»»»»»»»

% Forms tables
nsc_avg_parameters_CNT_Thr_data = array2table(nsc_avg_parameters_CNT_Thr,'VariableNames',grand_headers);                 
nsc_avg_parameters_CNT_Sub_data = array2table(nsc_avg_parameters_CNT_Sub,'VariableNames',grand_headers);
nsc_avg_parameters_T2DM_Thr_data = array2table(nsc_avg_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
nsc_avg_parameters_T2DM_Sub_data = array2table(nsc_avg_parameters_T2DM_Sub,'VariableNames',grand_headers);                 

nsc_std_parameters_CNT_Thr_data = array2table(nsc_std_parameters_CNT_Thr,'VariableNames',grand_headers);  
nsc_std_parameters_CNT_Sub_data = array2table(nsc_std_parameters_CNT_Sub,'VariableNames',grand_headers);                 
nsc_std_parameters_T2DM_Thr_data = array2table(nsc_std_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
nsc_std_parameters_T2DM_Sub_data = array2table(nsc_std_parameters_T2DM_Sub,'VariableNames',grand_headers); 

nsc_median_parameters_CNT_Thr_data = array2table(nsc_median_parameters_CNT_Thr,'VariableNames',grand_headers);                 
nsc_median_parameters_CNT_Sub_data = array2table(nsc_median_parameters_CNT_Sub,'VariableNames',grand_headers);
nsc_median_parameters_T2DM_Thr_data = array2table(nsc_median_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
nsc_median_parameters_T2DM_Sub_data = array2table(nsc_median_parameters_T2DM_Sub,'VariableNames',grand_headers);                 

nsc_iqr_parameters_CNT_Thr_data = array2table(nsc_iqr_parameters_CNT_Thr,'VariableNames',grand_headers);  
nsc_iqr_parameters_CNT_Sub_data = array2table(nsc_iqr_parameters_CNT_Sub,'VariableNames',grand_headers);                 
nsc_iqr_parameters_T2DM_Thr_data = array2table(nsc_iqr_parameters_T2DM_Thr,'VariableNames',grand_headers);                 
nsc_iqr_parameters_T2DM_Sub_data = array2table(nsc_iqr_parameters_T2DM_Sub,'VariableNames',grand_headers); 


% Merges tables
nsc_avg_parameters_CNT_Thr_data = [subject_CNT_column nsc_avg_parameters_CNT_Thr_data];                 
nsc_avg_parameters_CNT_Sub_data = [subject_CNT_column nsc_avg_parameters_CNT_Sub_data];
nsc_avg_parameters_T2DM_Thr_data = [subject_T2DM_column nsc_avg_parameters_T2DM_Thr_data];                 
nsc_avg_parameters_T2DM_Sub_data = [subject_T2DM_column nsc_avg_parameters_T2DM_Sub_data];                 

nsc_std_parameters_CNT_Thr_data = [subject_CNT_column nsc_std_parameters_CNT_Thr_data];  
nsc_std_parameters_CNT_Sub_data = [subject_CNT_column nsc_std_parameters_CNT_Sub_data];                 
nsc_std_parameters_T2DM_Thr_data = [subject_T2DM_column nsc_std_parameters_T2DM_Thr_data];                 
nsc_std_parameters_T2DM_Sub_data = [subject_T2DM_column nsc_std_parameters_T2DM_Sub_data]; 

nsc_median_parameters_CNT_Thr_data = [subject_CNT_column nsc_median_parameters_CNT_Thr_data];                 
nsc_median_parameters_CNT_Sub_data = [subject_CNT_column nsc_median_parameters_CNT_Sub_data];
nsc_median_parameters_T2DM_Thr_data = [subject_T2DM_column nsc_median_parameters_T2DM_Thr_data];                 
nsc_median_parameters_T2DM_Sub_data = [subject_T2DM_column nsc_median_parameters_T2DM_Sub_data];                 

nsc_iqr_parameters_CNT_Thr_data = [subject_CNT_column nsc_iqr_parameters_CNT_Thr_data];  
nsc_iqr_parameters_CNT_Sub_data = [subject_CNT_column nsc_iqr_parameters_CNT_Sub_data];                 
nsc_iqr_parameters_T2DM_Thr_data = [subject_T2DM_column nsc_iqr_parameters_T2DM_Thr_data];                 
nsc_iqr_parameters_T2DM_Sub_data = [subject_T2DM_column nsc_iqr_parameters_T2DM_Sub_data];