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
load('ROI_regions.mat');
load('HRF_parameters.mat');


%% CONSTANTS:

n_points = 8;                             % number of data points
n_subjects = 141;
n_T2DM = 64; 
n_CNT = 77;
n_rois = 22;
n_parameters = 11;                  

           
%% HRF PARAMETERS FROM THE AVERAGE AND MEDIAN HRF CURVES:

% In this section, we estimate the HRF parameters of interest from the 
% average and median HRF curves in each condition (Thr, Sub and average of 
% both) and group per ROI.


beta_thr_both_sides = zeros(n_points,2);                        % 2 since each column represents a stimulation side - left and right
beta_sub_both_sides = zeros(n_points,2);

beta_per_subject_CNT_Thr = zeros(n_points,n_CNT);
beta_per_subject_CNT_Sub = zeros(n_points,n_CNT);
beta_per_subject_CNT_Avg = zeros(n_points,n_CNT);
beta_per_subject_T2DM_Thr = zeros(n_points,n_T2DM);
beta_per_subject_T2DM_Sub = zeros(n_points,n_T2DM);
beta_per_subject_T2DM_Avg = zeros(n_points,n_T2DM);

total_average_HRF_parameters_CNT_Thr = zeros(n_parameters,n_rois);
total_average_HRF_parameters_CNT_Sub = zeros(n_parameters,n_rois);
total_average_HRF_parameters_CNT_Avg = zeros(n_parameters,n_rois);
total_average_HRF_parameters_T2DM_Thr = zeros(n_parameters,n_rois);
total_average_HRF_parameters_T2DM_Sub = zeros(n_parameters,n_rois);
total_average_HRF_parameters_T2DM_Avg = zeros(n_parameters,n_rois);

total_median_HRF_parameters_CNT_Thr = zeros(n_parameters,n_rois);
total_median_HRF_parameters_CNT_Sub = zeros(n_parameters,n_rois);
total_median_HRF_parameters_CNT_Avg = zeros(n_parameters,n_rois);
total_median_HRF_parameters_T2DM_Thr = zeros(n_parameters,n_rois);
total_median_HRF_parameters_T2DM_Sub = zeros(n_parameters,n_rois);
total_median_HRF_parameters_T2DM_Avg = zeros(n_parameters,n_rois);


for r=2:n_rois+1
    for i=1:n_subjects
        for j=1:2
            beta_thr_both_sides(:,j) = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));                     % each column is a side in each subject 
            beta_sub_both_sides(:,j) = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));
        end

        if i<65                                                                     % until subject 65 --> T2DM 
            beta_per_subject_T2DM_Thr(:,i) = mean(beta_thr_both_sides,2);           % mean of the beta values from each side per datapoint of T2DM subjects in the Thr condition in a ROI 
            beta_per_subject_T2DM_Sub(:,i) = mean(beta_sub_both_sides,2);
            beta_per_subject_T2DM_Avg(:,i) = mean([mean([beta_thr_both_sides(:,1) beta_sub_both_sides(:,1)],2) mean([beta_thr_both_sides(:,2) beta_sub_both_sides(:,2)],2)],2);         % average of the average beta values in each side --> avg(avg(Left Thr + Left Sub) + avg(Right Thr + Right Sub))
        else                                                                        % after subject 65 --> CNT
            beta_per_subject_CNT_Thr(:,i-n_T2DM) = mean(beta_thr_both_sides,2);     % mean of the beta values from each side per datapoint of CNT subjects in the Thr condition in a ROI 
            beta_per_subject_CNT_Sub(:,i-n_T2DM) = mean(beta_sub_both_sides,2);
            beta_per_subject_CNT_Avg(:,i-n_T2DM) = mean([mean([beta_sub_both_sides(:,1) beta_thr_both_sides(:,1)],2) mean([beta_sub_both_sides(:,2) beta_thr_both_sides(:,2)],2)],2);   
        end
    
        % Deletes to not duplicate or interfere with data
        beta_thr_both_sides = zeros(n_points,2);
        beta_sub_both_sides = zeros(n_points,2);
    end  
    
    
    % --------------------- THRESHOLD CONDITION -------------------------
    
    % Average and median HRF per datapoint in each subject group in the Thr 
    % condition in a ROI 
    avg_HRF_CNT_Thr = mean(beta_per_subject_CNT_Thr,2);                        
    avg_HRF_T2DM_Thr = mean(beta_per_subject_T2DM_Thr,2);                                        
    
    median_HRF_CNT_Thr = median(beta_per_subject_CNT_Thr,2);                    
    median_HRF_T2DM_Thr = median(beta_per_subject_T2DM_Thr,2);                  
        
    
        % «««««««««««««««««««««««« PARAMETERS »»»»»»»»»»»»»»»»»»»»»»»»»»»
        
        % HRF parameters of interest per group in the Thr condition for 
        % each ROI based on the average / median HRF 
        
        
        % **************************** CNT ******************************

        [peak,volume] = max(avg_HRF_CNT_Thr(2:6,1));                                                        % peak amplitude: maximum within the sectioned time range (2.5 to 12.5 s - 2nd to 6th index since indexing starts at 1)
        peak_latency = volume * 2.5;                                                                        % peak latency: time to peak (volume where the maximum takes place * TR)
        slope_to_peak = minus(peak,avg_HRF_CNT_Thr(1))/peak_latency;                                        % slope to peak: (peak-1st point)/latency
        total_area = trapz((0:n_points-1)*2.5,avg_HRF_CNT_Thr - min(avg_HRF_CNT_Thr));                      % area under the curve (trapezoidal method, from 0 to 17.5 s)
        pos_area = positive_area((0:n_points-1)*2.5,avg_HRF_CNT_Thr);                                       % area of the positive sections of the HRF curve
        first_pos_area = positive_area((0:2)*2.5,avg_HRF_CNT_Thr(1:3,1));                                   % area of the first positive section of the HRF curve (from 0 to 5 s)
        second_pos_area = positive_area((2:4)*2.5,avg_HRF_CNT_Thr(3:5,1));                                  % area of the second section of the HRF curve (from 5 to 10 s)
        third_pos_area = positive_area((4:n_points-1)*2.5,avg_HRF_CNT_Thr(5:n_points,1));                   % area of the third section of the HRF curve (from 10 to 17.5 s)
        neg_area = negative_area((0:n_points-1)*2.5,avg_HRF_CNT_Thr);                                       % area of the negative sections of the HRF curve
        initial_dip_area = negative_area((0:volume)*2.5,avg_HRF_CNT_Thr(1:(volume+1),1));                   % area of the initial dip of the HRF curve
        undershoot_area = negative_area((volume:n_points-1)*2.5,avg_HRF_CNT_Thr((volume+1):n_points,1));    % area of the undershoot of the HRF curve
        total_average_HRF_parameters_CNT_Thr(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        
        [peak,volume] = max(median_HRF_CNT_Thr(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_HRF_CNT_Thr(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_HRF_CNT_Thr - min(median_HRF_CNT_Thr));
        pos_area = positive_area((0:n_points-1)*2.5,median_HRF_CNT_Thr);
        first_pos_area = positive_area((0:2)*2.5,median_HRF_CNT_Thr(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_HRF_CNT_Thr(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,median_HRF_CNT_Thr(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_HRF_CNT_Thr);
        initial_dip_area = negative_area((0:volume)*2.5,median_HRF_CNT_Thr(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,median_HRF_CNT_Thr((volume+1):n_points,1));
        total_median_HRF_parameters_CNT_Thr(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];

        
        % *************************** T2DM ******************************
        
        [peak,volume] = max(avg_HRF_T2DM_Thr(2:6,1));                                       
        peak_latency = volume * 2.5;                                                
        slope_to_peak = minus(peak,avg_HRF_T2DM_Thr(1))/peak_latency;                       
        total_area = trapz((0:n_points-1)*2.5,avg_HRF_T2DM_Thr - min(avg_HRF_T2DM_Thr));        
        pos_area = positive_area((0:n_points-1)*2.5,avg_HRF_T2DM_Thr);                      
        first_pos_area = positive_area((0:2)*2.5,avg_HRF_T2DM_Thr(1:3,1));                  
        second_pos_area = positive_area((2:4)*2.5,avg_HRF_T2DM_Thr(3:5,1));                 
        third_pos_area = positive_area((4:n_points-1)*2.5,avg_HRF_T2DM_Thr(5:n_points,1));                  
        neg_area = negative_area((0:n_points-1)*2.5,avg_HRF_T2DM_Thr);                      
        initial_dip_area = negative_area((0:volume)*2.5,avg_HRF_T2DM_Thr(1:(volume+1),1));  
        undershoot_area = negative_area((volume:n_points-1)*2.5,avg_HRF_T2DM_Thr((volume+1):n_points,1));   
        total_average_HRF_parameters_T2DM_Thr(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];

        [peak,volume] = max(median_HRF_T2DM_Thr(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_HRF_T2DM_Thr(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_HRF_T2DM_Thr - min(median_HRF_T2DM_Thr));
        pos_area = positive_area((0:n_points-1)*2.5,median_HRF_T2DM_Thr);
        first_pos_area = positive_area((0:2)*2.5,median_HRF_T2DM_Thr(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_HRF_T2DM_Thr(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,median_HRF_T2DM_Thr(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_HRF_T2DM_Thr);
        initial_dip_area = negative_area((0:volume)*2.5,median_HRF_T2DM_Thr(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,median_HRF_T2DM_Thr((volume+1):n_points,1));
        total_median_HRF_parameters_T2DM_Thr(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];

        
    % ---------------------- SUBMAXIMUM CONDITION ------------------------
   
    % Average and median HRF per datapoint in each subject group in the Sub 
    % condition in a ROI
    avg_HRF_CNT_Sub = mean(beta_per_subject_CNT_Sub,2);                      
    avg_HRF_T2DM_Sub = mean(beta_per_subject_T2DM_Sub,2);                                      

    median_HRF_CNT_Sub = median(beta_per_subject_CNT_Sub,2);                  
    median_HRF_T2DM_Sub = median(beta_per_subject_T2DM_Sub,2);                                                 
       
    
        % «««««««««««««««««««««««« PARAMETERS »»»»»»»»»»»»»»»»»»»»»»»»»»»

        % HRF parameters of interest per group in the Sub condition for 
        % each ROI based on the average / median HRF  
        
        
        % **************************** CNT ******************************
        
        [peak,volume] = max(avg_HRF_CNT_Sub(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,avg_HRF_CNT_Sub(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,avg_HRF_CNT_Sub - min(avg_HRF_CNT_Sub));
        pos_area = positive_area((0:n_points-1)*2.5,avg_HRF_CNT_Sub);
        first_pos_area = positive_area((0:2)*2.5,avg_HRF_CNT_Sub(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,avg_HRF_CNT_Sub(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,avg_HRF_CNT_Sub(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,avg_HRF_CNT_Sub);
        initial_dip_area = negative_area((0:volume)*2.5,avg_HRF_CNT_Sub(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,avg_HRF_CNT_Sub((volume+1):n_points,1));
        total_average_HRF_parameters_CNT_Sub(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];

        [peak,volume] = max(median_HRF_CNT_Sub(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_HRF_CNT_Sub(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_HRF_CNT_Sub - min(median_HRF_CNT_Sub));
        pos_area = positive_area((0:n_points-1)*2.5,median_HRF_CNT_Sub);
        first_pos_area = positive_area((0:2)*2.5,median_HRF_CNT_Sub(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_HRF_CNT_Sub(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,median_HRF_CNT_Sub(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_HRF_CNT_Sub);
        initial_dip_area = negative_area((0:volume)*2.5,median_HRF_CNT_Sub(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,median_HRF_CNT_Sub((volume+1):n_points,1));
        total_median_HRF_parameters_CNT_Sub(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        
        
        % *************************** T2DM ******************************

        [peak,volume] = max(avg_HRF_T2DM_Sub(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,avg_HRF_T2DM_Sub(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,avg_HRF_T2DM_Sub - min(avg_HRF_T2DM_Sub));
        pos_area = positive_area((0:n_points-1)*2.5,avg_HRF_T2DM_Sub);
        first_pos_area = positive_area((0:2)*2.5,avg_HRF_T2DM_Sub(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,avg_HRF_T2DM_Sub(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,avg_HRF_T2DM_Sub(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,avg_HRF_T2DM_Sub);
        initial_dip_area = negative_area((0:volume)*2.5,avg_HRF_T2DM_Sub(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,avg_HRF_T2DM_Sub((volume+1):n_points,1));
        total_average_HRF_parameters_T2DM_Sub(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        
        [peak,volume] = max(median_HRF_T2DM_Sub(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_HRF_T2DM_Sub(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_HRF_T2DM_Sub - min(median_HRF_T2DM_Sub));
        pos_area = positive_area((0:n_points-1)*2.5,median_HRF_T2DM_Sub);
        first_pos_area = positive_area((0:2)*2.5,median_HRF_T2DM_Sub(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_HRF_T2DM_Sub(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,median_HRF_T2DM_Sub(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_HRF_T2DM_Sub);
        initial_dip_area = negative_area((0:volume)*2.5,median_HRF_T2DM_Sub(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,median_HRF_T2DM_Sub((volume+1):n_points,1));
        total_median_HRF_parameters_T2DM_Sub(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];      

           
    % ----------------------- AVERAGE CONDITION -------------------------
    
    % Average and median HRF per datapoint in each subject group in the Avg 
    % condition in a ROI 
    avg_HRF_CNT_Avg = mean(beta_per_subject_CNT_Avg,2);                        
    avg_HRF_T2DM_Avg = mean(beta_per_subject_T2DM_Avg,2);                                        
    
    median_HRF_CNT_Avg = median(beta_per_subject_CNT_Avg,2);                    
    median_HRF_T2DM_Avg = median(beta_per_subject_T2DM_Avg,2);                  
    
    
        % «««««««««««««««««««««««« PARAMETERS »»»»»»»»»»»»»»»»»»»»»»»»»»»

        % HRF parameters of interest per group in the Avg condition for 
        % each ROI based on the average / median HRF 


        % **************************** CNT ******************************

        [peak,volume] = max(avg_HRF_CNT_Avg(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,avg_HRF_CNT_Avg(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,avg_HRF_CNT_Avg - min(avg_HRF_CNT_Avg));
        pos_area = positive_area((0:n_points-1)*2.5,avg_HRF_CNT_Avg);
        first_pos_area = positive_area((0:2)*2.5,avg_HRF_CNT_Avg(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,avg_HRF_CNT_Avg(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,avg_HRF_CNT_Avg(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,avg_HRF_CNT_Avg);
        initial_dip_area = negative_area((0:volume)*2.5,avg_HRF_CNT_Avg(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,avg_HRF_CNT_Avg((volume+1):n_points,1));
        total_average_HRF_parameters_CNT_Avg(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];

        [peak,volume] = max(median_HRF_CNT_Avg(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_HRF_CNT_Avg(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_HRF_CNT_Avg - min(median_HRF_CNT_Avg));
        pos_area = positive_area((0:n_points-1)*2.5,median_HRF_CNT_Avg);
        first_pos_area = positive_area((0:2)*2.5,median_HRF_CNT_Avg(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_HRF_CNT_Avg(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,median_HRF_CNT_Avg(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_HRF_CNT_Avg);
        initial_dip_area = negative_area((0:volume)*2.5,median_HRF_CNT_Avg(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,median_HRF_CNT_Avg((volume+1):n_points,1));
        total_median_HRF_parameters_CNT_Avg(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];


        % *************************** T2DM ******************************

        [peak,volume] = max(avg_HRF_T2DM_Avg(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,avg_HRF_T2DM_Avg(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,avg_HRF_T2DM_Avg - min(avg_HRF_T2DM_Avg));
        pos_area = positive_area((0:n_points-1)*2.5,avg_HRF_T2DM_Avg);
        first_pos_area = positive_area((0:2)*2.5,avg_HRF_T2DM_Avg(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,avg_HRF_T2DM_Avg(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,avg_HRF_T2DM_Avg(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,avg_HRF_T2DM_Avg);
        initial_dip_area = negative_area((0:volume)*2.5,avg_HRF_T2DM_Avg(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,avg_HRF_T2DM_Avg((volume+1):n_points,1));
        total_average_HRF_parameters_T2DM_Avg(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];

        [peak,volume] = max(median_HRF_T2DM_Avg(2:6,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_HRF_T2DM_Avg(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_HRF_T2DM_Avg - min(median_HRF_T2DM_Avg));
        pos_area = positive_area((0:n_points-1)*2.5,median_HRF_T2DM_Avg);
        first_pos_area = positive_area((0:2)*2.5,median_HRF_T2DM_Avg(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_HRF_T2DM_Avg(3:5,1));
        third_pos_area = positive_area((4:n_points-1)*2.5,median_HRF_T2DM_Avg(5:n_points,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_HRF_T2DM_Avg);
        initial_dip_area = negative_area((0:volume)*2.5,median_HRF_T2DM_Avg(1:(volume+1),1));
        undershoot_area = negative_area((volume:n_points-1)*2.5,median_HRF_T2DM_Avg((volume+1):n_points,1));
        total_median_HRF_parameters_T2DM_Avg(:,r-1) = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];      


    % Deletes to not duplicate or interfere with data
    beta_per_subject_CNT_Thr = zeros(n_points,n_CNT);
    beta_per_subject_CNT_Sub = zeros(n_points,n_CNT);
    beta_per_subject_CNT_Avg = zeros(n_points,n_CNT);
    beta_per_subject_T2DM_Thr = zeros(n_points,n_T2DM);
    beta_per_subject_T2DM_Sub = zeros(n_points,n_T2DM);
    beta_per_subject_T2DM_Avg = zeros(n_points,n_T2DM);
end


%% STATISTICS:

% In this section, we estimate the statistics of the HRF parameters of
% interest from the average and median HRF curves in each condition (Thr, 
% Sub, average) and group per set of ROIs (positive and negative signal 
% change ROIs).
% First, we apply a Shapiro-Wilk test in order to assess if the data has a 
% normal distribution per group in the same condition. Then, we assess if 
% there are significant differences for the same HRF parameters between 
% groups. To do so, we apply a Two-sample T-test or a Wilcoxon ranksum 
% test, depending on whether the data regarding the same HRF parameter in 
% both groups of the same condition has a normal distribution or not, 
% respectively.


sw_psc_avg_HRF_CNT_Thr = zeros(3,n_parameters);
sw_psc_avg_HRF_CNT_Sub = zeros(3,n_parameters);
sw_psc_avg_HRF_CNT_Avg = zeros(3,n_parameters);
sw_psc_avg_HRF_T2DM_Thr = zeros(3,n_parameters);
sw_psc_avg_HRF_T2DM_Sub = zeros(3,n_parameters);
sw_psc_avg_HRF_T2DM_Avg = zeros(3,n_parameters);

sw_nsc_avg_HRF_CNT_Thr = zeros(3,n_parameters);
sw_nsc_avg_HRF_CNT_Sub = zeros(3,n_parameters);
sw_nsc_avg_HRF_CNT_Avg = zeros(3,n_parameters);
sw_nsc_avg_HRF_T2DM_Thr = zeros(3,n_parameters);
sw_nsc_avg_HRF_T2DM_Sub = zeros(3,n_parameters);
sw_nsc_avg_HRF_T2DM_Avg = zeros(3,n_parameters);

sw_psc_median_HRF_CNT_Thr = zeros(3,n_parameters);
sw_psc_median_HRF_CNT_Sub = zeros(3,n_parameters);
sw_psc_median_HRF_CNT_Avg = zeros(3,n_parameters);
sw_psc_median_HRF_T2DM_Thr = zeros(3,n_parameters);
sw_psc_median_HRF_T2DM_Sub = zeros(3,n_parameters);
sw_psc_median_HRF_T2DM_Avg = zeros(3,n_parameters);

sw_nsc_median_HRF_CNT_Thr = zeros(3,n_parameters);
sw_nsc_median_HRF_CNT_Sub = zeros(3,n_parameters);
sw_nsc_median_HRF_CNT_Avg = zeros(3,n_parameters);
sw_nsc_median_HRF_T2DM_Thr = zeros(3,n_parameters);
sw_nsc_median_HRF_T2DM_Sub = zeros(3,n_parameters);
sw_nsc_median_HRF_T2DM_Avg = zeros(3,n_parameters);


psc_avg_std_indexes_avg_HRF_CNT_Thr = [];
psc_avg_std_indexes_avg_HRF_CNT_Sub = []; 
psc_avg_std_indexes_avg_HRF_CNT_Avg = []; 
psc_avg_std_indexes_avg_HRF_T2DM_Thr = [];
psc_avg_std_indexes_avg_HRF_T2DM_Sub = [];
psc_avg_std_indexes_avg_HRF_T2DM_Avg = [];

nsc_avg_std_indexes_avg_HRF_CNT_Thr = [];
nsc_avg_std_indexes_avg_HRF_CNT_Sub = []; 
nsc_avg_std_indexes_avg_HRF_CNT_Avg = []; 
nsc_avg_std_indexes_avg_HRF_T2DM_Thr = [];
nsc_avg_std_indexes_avg_HRF_T2DM_Sub = [];
nsc_avg_std_indexes_avg_HRF_T2DM_Avg = [];

psc_avg_std_indexes_median_HRF_CNT_Thr = [];
psc_avg_std_indexes_median_HRF_CNT_Sub = []; 
psc_avg_std_indexes_median_HRF_CNT_Avg = []; 
psc_avg_std_indexes_median_HRF_T2DM_Thr = [];
psc_avg_std_indexes_median_HRF_T2DM_Sub = [];
psc_avg_std_indexes_median_HRF_T2DM_Avg = [];

nsc_avg_std_indexes_median_HRF_CNT_Thr = [];
nsc_avg_std_indexes_median_HRF_CNT_Sub = []; 
nsc_avg_std_indexes_median_HRF_CNT_Avg = []; 
nsc_avg_std_indexes_median_HRF_T2DM_Thr = [];
nsc_avg_std_indexes_median_HRF_T2DM_Sub = [];
nsc_avg_std_indexes_median_HRF_T2DM_Avg = [];


psc_median_iqr_indexes_avg_HRF_CNT_Thr = [];
psc_median_iqr_indexes_avg_HRF_CNT_Sub = [];
psc_median_iqr_indexes_avg_HRF_CNT_Avg = [];
psc_median_iqr_indexes_avg_HRF_T2DM_Thr = [];
psc_median_iqr_indexes_avg_HRF_T2DM_Sub = [];
psc_median_iqr_indexes_avg_HRF_T2DM_Avg = [];

nsc_median_iqr_indexes_avg_HRF_CNT_Thr = [];
nsc_median_iqr_indexes_avg_HRF_CNT_Sub = [];
nsc_median_iqr_indexes_avg_HRF_CNT_Avg = [];
nsc_median_iqr_indexes_avg_HRF_T2DM_Thr = [];
nsc_median_iqr_indexes_avg_HRF_T2DM_Sub = [];
nsc_median_iqr_indexes_avg_HRF_T2DM_Avg = [];

psc_median_iqr_indexes_median_HRF_CNT_Thr = [];
psc_median_iqr_indexes_median_HRF_CNT_Sub = [];
psc_median_iqr_indexes_median_HRF_CNT_Avg = [];
psc_median_iqr_indexes_median_HRF_T2DM_Thr = [];
psc_median_iqr_indexes_median_HRF_T2DM_Sub = [];
psc_median_iqr_indexes_median_HRF_T2DM_Avg = [];

nsc_median_iqr_indexes_median_HRF_CNT_Thr = [];
nsc_median_iqr_indexes_median_HRF_CNT_Sub = [];
nsc_median_iqr_indexes_median_HRF_CNT_Avg = [];
nsc_median_iqr_indexes_median_HRF_T2DM_Thr = [];
nsc_median_iqr_indexes_median_HRF_T2DM_Sub = [];
nsc_median_iqr_indexes_median_HRF_T2DM_Avg = [];


psc_avg_parameters_avg_HRF_CNT_Thr = [];
psc_avg_parameters_avg_HRF_CNT_Sub = [];
psc_avg_parameters_avg_HRF_CNT_Avg = [];
psc_avg_parameters_avg_HRF_T2DM_Thr = [];
psc_avg_parameters_avg_HRF_T2DM_Sub = [];
psc_avg_parameters_avg_HRF_T2DM_Avg = [];

nsc_avg_parameters_median_HRF_CNT_Thr = []; 
nsc_avg_parameters_median_HRF_CNT_Sub = []; 
nsc_avg_parameters_median_HRF_CNT_Avg = []; 
nsc_avg_parameters_median_HRF_T2DM_Thr = [];
nsc_avg_parameters_median_HRF_T2DM_Sub = [];
nsc_avg_parameters_median_HRF_T2DM_Avg = [];

psc_avg_parameters_median_HRF_CNT_Thr = [];
psc_avg_parameters_median_HRF_CNT_Sub = [];
psc_avg_parameters_median_HRF_CNT_Avg = [];
psc_avg_parameters_median_HRF_T2DM_Thr = [];
psc_avg_parameters_median_HRF_T2DM_Sub = [];
psc_avg_parameters_median_HRF_T2DM_Avg = [];

nsc_avg_parameters_avg_HRF_CNT_Thr = []; 
nsc_avg_parameters_avg_HRF_CNT_Sub = []; 
nsc_avg_parameters_avg_HRF_CNT_Avg = []; 
nsc_avg_parameters_avg_HRF_T2DM_Thr = [];
nsc_avg_parameters_avg_HRF_T2DM_Sub = [];
nsc_avg_parameters_avg_HRF_T2DM_Avg = [];


psc_std_parameters_avg_HRF_CNT_Thr = [];
psc_std_parameters_avg_HRF_CNT_Sub = [];
psc_std_parameters_avg_HRF_CNT_Avg = [];
psc_std_parameters_avg_HRF_T2DM_Thr = [];
psc_std_parameters_avg_HRF_T2DM_Sub = [];
psc_std_parameters_avg_HRF_T2DM_Avg = [];

nsc_std_parameters_avg_HRF_CNT_Thr = [];
nsc_std_parameters_avg_HRF_CNT_Sub = [];
nsc_std_parameters_avg_HRF_CNT_Avg = [];
nsc_std_parameters_avg_HRF_T2DM_Thr = [];
nsc_std_parameters_avg_HRF_T2DM_Sub = [];
nsc_std_parameters_avg_HRF_T2DM_Avg = [];

psc_std_parameters_median_HRF_CNT_Thr = [];
psc_std_parameters_median_HRF_CNT_Sub = [];
psc_std_parameters_median_HRF_CNT_Avg = [];
psc_std_parameters_median_HRF_T2DM_Thr = [];
psc_std_parameters_median_HRF_T2DM_Sub = [];
psc_std_parameters_median_HRF_T2DM_Avg = [];

nsc_std_parameters_median_HRF_CNT_Thr = [];
nsc_std_parameters_median_HRF_CNT_Sub = [];
nsc_std_parameters_median_HRF_CNT_Avg = [];
nsc_std_parameters_median_HRF_T2DM_Thr = [];
nsc_std_parameters_median_HRF_T2DM_Sub = [];
nsc_std_parameters_median_HRF_T2DM_Avg = [];


psc_median_parameters_avg_HRF_CNT_Thr = [];
psc_median_parameters_avg_HRF_CNT_Sub = [];
psc_median_parameters_avg_HRF_CNT_Avg = [];
psc_median_parameters_avg_HRF_T2DM_Thr = [];
psc_median_parameters_avg_HRF_T2DM_Sub = [];
psc_median_parameters_avg_HRF_T2DM_Avg = [];

nsc_median_parameters_avg_HRF_CNT_Thr = [];
nsc_median_parameters_avg_HRF_CNT_Sub = [];
nsc_median_parameters_avg_HRF_CNT_Avg = [];
nsc_median_parameters_avg_HRF_T2DM_Thr = [];
nsc_median_parameters_avg_HRF_T2DM_Sub = [];
nsc_median_parameters_avg_HRF_T2DM_Avg = [];

psc_median_parameters_median_HRF_CNT_Thr = [];
psc_median_parameters_median_HRF_CNT_Sub = [];
psc_median_parameters_median_HRF_CNT_Avg = [];
psc_median_parameters_median_HRF_T2DM_Thr = [];
psc_median_parameters_median_HRF_T2DM_Sub = [];
psc_median_parameters_median_HRF_T2DM_Avg = [];

nsc_median_parameters_median_HRF_CNT_Thr = [];
nsc_median_parameters_median_HRF_CNT_Sub = [];
nsc_median_parameters_median_HRF_CNT_Avg = [];
nsc_median_parameters_median_HRF_T2DM_Thr = [];
nsc_median_parameters_median_HRF_T2DM_Sub = [];
nsc_median_parameters_median_HRF_T2DM_Avg = [];


psc_iqr_parameters_avg_HRF_CNT_Thr = [];
psc_iqr_parameters_avg_HRF_CNT_Sub = [];
psc_iqr_parameters_avg_HRF_CNT_Avg = [];
psc_iqr_parameters_avg_HRF_T2DM_Thr = [];
psc_iqr_parameters_avg_HRF_T2DM_Sub = [];
psc_iqr_parameters_avg_HRF_T2DM_Avg = [];

nsc_iqr_parameters_avg_HRF_CNT_Thr = [];
nsc_iqr_parameters_avg_HRF_CNT_Sub = [];
nsc_iqr_parameters_avg_HRF_CNT_Avg = [];
nsc_iqr_parameters_avg_HRF_T2DM_Thr = [];
nsc_iqr_parameters_avg_HRF_T2DM_Sub = [];
nsc_iqr_parameters_avg_HRF_T2DM_Avg = [];

psc_iqr_parameters_median_HRF_CNT_Thr = [];
psc_iqr_parameters_median_HRF_CNT_Sub = [];
psc_iqr_parameters_median_HRF_CNT_Avg = [];
psc_iqr_parameters_median_HRF_T2DM_Thr = [];
psc_iqr_parameters_median_HRF_T2DM_Sub = [];
psc_iqr_parameters_median_HRF_T2DM_Avg = [];

nsc_iqr_parameters_median_HRF_CNT_Thr = [];
nsc_iqr_parameters_median_HRF_CNT_Sub = [];
nsc_iqr_parameters_median_HRF_CNT_Avg = [];
nsc_iqr_parameters_median_HRF_T2DM_Thr = [];
nsc_iqr_parameters_median_HRF_T2DM_Sub = [];
nsc_iqr_parameters_median_HRF_T2DM_Avg = [];


psc_ttest_indexes_avg_HRF_thr = [];
psc_ttest_indexes_avg_HRF_sub = [];
psc_ttest_indexes_avg_HRF_avg = [];
nsc_ttest_indexes_avg_HRF_thr = [];
nsc_ttest_indexes_avg_HRF_sub = [];
nsc_ttest_indexes_avg_HRF_avg = [];

psc_wilcoxon_indexes_avg_HRF_thr = [];
psc_wilcoxon_indexes_avg_HRF_sub = [];
psc_wilcoxon_indexes_avg_HRF_avg = [];
nsc_wilcoxon_indexes_avg_HRF_thr = [];
nsc_wilcoxon_indexes_avg_HRF_sub = [];
nsc_wilcoxon_indexes_avg_HRF_avg = [];

psc_ttest_indexes_median_HRF_thr = [];
psc_ttest_indexes_median_HRF_sub = [];
psc_ttest_indexes_median_HRF_avg = [];
nsc_ttest_indexes_median_HRF_thr = [];
nsc_ttest_indexes_median_HRF_sub = [];
nsc_ttest_indexes_median_HRF_avg = [];

psc_wilcoxon_indexes_median_HRF_thr = [];
psc_wilcoxon_indexes_median_HRF_sub = [];
psc_wilcoxon_indexes_median_HRF_avg = [];
nsc_wilcoxon_indexes_median_HRF_thr = [];
nsc_wilcoxon_indexes_median_HRF_sub = [];
nsc_wilcoxon_indexes_median_HRF_avg = [];


psc_ttest_avg_HRF_thr = [];
psc_ttest_avg_HRF_sub = [];
psc_ttest_avg_HRF_avg = [];
nsc_ttest_avg_HRF_thr = [];
nsc_ttest_avg_HRF_sub = [];
nsc_ttest_avg_HRF_avg = [];

psc_ttest_median_HRF_thr = [];
psc_ttest_median_HRF_sub = [];
psc_ttest_median_HRF_avg = [];
nsc_ttest_median_HRF_thr = [];
nsc_ttest_median_HRF_sub = [];
nsc_ttest_median_HRF_avg = [];

psc_wilcoxon_avg_HRF_thr = [];
psc_wilcoxon_avg_HRF_sub = [];
psc_wilcoxon_avg_HRF_avg = [];
nsc_wilcoxon_avg_HRF_thr = [];
nsc_wilcoxon_avg_HRF_sub = [];
nsc_wilcoxon_avg_HRF_avg = [];

psc_wilcoxon_median_HRF_thr = [];
psc_wilcoxon_median_HRF_sub = [];
psc_wilcoxon_median_HRF_avg = [];
nsc_wilcoxon_median_HRF_thr = [];
nsc_wilcoxon_median_HRF_sub = [];
nsc_wilcoxon_median_HRF_avg = [];



for par = 1:n_parameters

    
    % ---------------------- THRESHOLD CONDITION ------------------------

    
    
    % «««««««««««««««««««««««««« Average HRF »»»»»»»»»»»»»»»»»»»»»»»»»»»»
    
    
    % ****************** Positive Signal Change ROIs ********************


    [h,p,stats]=swtest(total_average_HRF_parameters_CNT_Thr(par,1:10),0.05,true);               % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Thr condition for each ROI of the positive signal change ROIs
    sw_psc_avg_HRF_CNT_Thr(:,par) = [h;p;stats];
    
    
    if h==0                                                                                     % in case of a normal distribution --> group average of the HRF parameters
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the  
        % positive signal change ROIs
        psc_avg_std_indexes_avg_HRF_CNT_Thr = [psc_avg_std_indexes_avg_HRF_CNT_Thr par];        
        
        % Average + standard-deviation of the average HRF's parameters with 
        % normal distribution of CNT subjects in the Thr condition in the
        % positive signal change ROIs
        psc_avg_parameters_avg_HRF_CNT_Thr = [psc_avg_parameters_avg_HRF_CNT_Thr mean(total_average_HRF_parameters_CNT_Thr(par,1:10),2)];
        psc_std_parameters_avg_HRF_CNT_Thr = [psc_std_parameters_avg_HRF_CNT_Thr std(total_average_HRF_parameters_CNT_Thr(par,1:10),0,2)];
    
    else                                                                                        % in case of a non-normal distribution --> group median of the HRF parameters

        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of CNT subjects in the Thr condition in 
        % the positive signal change ROIs
        psc_median_iqr_indexes_avg_HRF_CNT_Thr = [psc_median_iqr_indexes_avg_HRF_CNT_Thr par];
                 
        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribuition of CNT subjects in the Thr condition in
        % the positive signal change ROIs
        psc_median_parameters_avg_HRF_CNT_Thr = [psc_median_parameters_avg_HRF_CNT_Thr median(total_average_HRF_parameters_CNT_Thr(par,1:10),2)];
        psc_iqr_parameters_avg_HRF_CNT_Thr = [psc_iqr_parameters_avg_HRF_CNT_Thr iqr(total_average_HRF_parameters_CNT_Thr(par,1:10),2)];

    end
    
    
    [h,p,stats] = swtest(total_average_HRF_parameters_T2DM_Thr(par,1:10), 0.05, true);          % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of T2DM subjects in the Thr condition for each ROI of the positive signal change ROIs       
    sw_psc_avg_HRF_T2DM_Thr(:,par) = [h;p;stats];
    
    
    if h==0                                                                                     
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the  
        % positive signal change ROIs
        psc_avg_std_indexes_avg_HRF_T2DM_Thr = [psc_avg_std_indexes_avg_HRF_T2DM_Thr par];
        
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the
        % positive signal change ROIs
        psc_avg_parameters_avg_HRF_T2DM_Thr = [psc_avg_parameters_avg_HRF_T2DM_Thr mean(total_average_HRF_parameters_T2DM_Thr(par,1:10),2)];
        psc_std_parameters_avg_HRF_T2DM_Thr = [psc_std_parameters_avg_HRF_T2DM_Thr std(total_average_HRF_parameters_T2DM_Thr(par,1:10),0,2)];
 
    else                                                                                                                                                                          
         
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of T2DM subjects in the Thr condition in 
        % the positive signal change ROIs
        psc_median_iqr_indexes_avg_HRF_T2DM_Thr = [psc_median_iqr_indexes_avg_HRF_T2DM_Thr par];

        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Thr condition in
        % the positive signal change ROIs
        psc_median_parameters_avg_HRF_T2DM_Thr = [psc_median_parameters_avg_HRF_T2DM_Thr median(total_average_HRF_parameters_T2DM_Thr(par,1:10),2)];
        psc_iqr_parameters_avg_HRF_T2DM_Thr = [psc_iqr_parameters_avg_HRF_T2DM_Thr iqr(total_average_HRF_parameters_T2DM_Thr(par,1:10),2)];
    
    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_psc_avg_HRF_CNT_Thr(1,par)==0 && sw_psc_avg_HRF_T2DM_Thr(1,par)==0                    % in case of a normal distribution in the same HRF parameter from different groups --> Two sample T-test
        psc_ttest_indexes_avg_HRF_thr = [psc_ttest_indexes_avg_HRF_thr par];                        % indexes for the corresponding HRF parameters from the average HRF with normal distribution in the Thr condition in the positive signal change ROIs
        [h,p,~,stats] = ttest2(total_average_HRF_parameters_CNT_Thr(par,1:10)',total_average_HRF_parameters_T2DM_Thr(par,1:10)');
        psc_ttest_avg_HRF_thr = [psc_ttest_avg_HRF_thr [h;p;stats.tstat]];                          % T-test results of the HRF parameters from the average HRF in the Thr condition in the positive signal change ROIs
    else                                                                                        % in case of a non-normal distribution in the same HRF parameter from different groups --> Wilcoxon ranksum test
        psc_wilcoxon_indexes_avg_HRF_thr = [psc_wilcoxon_indexes_avg_HRF_thr par];                  % indexes for the corresponding HRF parameters from the average HRF with non-normal distribution in the Thr condition in the positive signal change ROIs
        [p,h,stats] = ranksum(total_average_HRF_parameters_CNT_Thr(par,1:10)',total_average_HRF_parameters_T2DM_Thr(par,1:10)');
        psc_wilcoxon_avg_HRF_thr = [psc_wilcoxon_avg_HRF_thr [h;p;stats.zval]];                     % Wilcoxon test results of the HRF parameters from the average HRF in the Thr condition in the positive signal change ROIs
    end
 
    
    
    % ****************** Negative Signal Change ROIs ********************

    
    [h,p,stats]=swtest(total_average_HRF_parameters_CNT_Thr(par,11:end),0.05,true);             % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Thr condition for each ROI of the negative signal change ROIs             
    sw_nsc_avg_HRF_CNT_Thr(:,par) = [h;p;stats];

    
    if h==0                                                                                     
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the
        % negative signal change ROIs
        nsc_avg_std_indexes_avg_HRF_CNT_Thr = [nsc_avg_std_indexes_avg_HRF_CNT_Thr par];
                
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the
        % negative signal change ROIs
        nsc_avg_parameters_avg_HRF_CNT_Thr = [nsc_avg_parameters_avg_HRF_CNT_Thr mean(total_average_HRF_parameters_CNT_Thr(par,11:end),2)];
        nsc_std_parameters_avg_HRF_CNT_Thr = [nsc_std_parameters_avg_HRF_CNT_Thr std(total_average_HRF_parameters_CNT_Thr(par,11:end),0,2)];
        
    else                                                                                        
        
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of CNT subjects in the Thr condition in
        % the negative signal change ROIs        
        nsc_median_iqr_indexes_avg_HRF_CNT_Thr = [nsc_median_iqr_indexes_avg_HRF_CNT_Thr par];
                
        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of CNT subjects in the Thr condition in
        % the negative signal change ROIs        
        nsc_median_parameters_avg_HRF_CNT_Thr = [nsc_median_parameters_avg_HRF_CNT_Thr median(total_average_HRF_parameters_CNT_Thr(par,11:end),2)];
        nsc_iqr_parameters_avg_HRF_CNT_Thr = [nsc_iqr_parameters_avg_HRF_CNT_Thr iqr(total_average_HRF_parameters_CNT_Thr(par,11:end),2)];   

    end
    
    
    [h,p,stats]=swtest(total_average_HRF_parameters_T2DM_Thr(par,11:end),0.05,true);            % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Thr condition for each ROI of the negative signal change ROIs             
    sw_nsc_avg_HRF_T2DM_Thr(:,par) = [h;p;stats];


    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the
        % negative signal change ROIs
        nsc_avg_std_indexes_avg_HRF_T2DM_Thr = [nsc_avg_std_indexes_avg_HRF_T2DM_Thr par];
        
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the
        % negative signal change ROIs        
        nsc_avg_parameters_avg_HRF_T2DM_Thr = [nsc_avg_parameters_avg_HRF_T2DM_Thr mean(total_average_HRF_parameters_T2DM_Thr(par,11:end),2)];
        nsc_std_parameters_avg_HRF_T2DM_Thr = [nsc_std_parameters_avg_HRF_T2DM_Thr std(total_average_HRF_parameters_T2DM_Thr(par,11:end),0,2)];
        
    else                                                                                            
        
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of T2DM subjects in the Thr condition in
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_avg_HRF_T2DM_Thr = [nsc_median_iqr_indexes_avg_HRF_T2DM_Thr par];       
        
        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Thr condition in
        % the negative signal change ROIs  
        nsc_median_parameters_avg_HRF_T2DM_Thr = [nsc_median_parameters_avg_HRF_T2DM_Thr median(total_average_HRF_parameters_T2DM_Thr(par,11:end),2)];
        nsc_iqr_parameters_avg_HRF_T2DM_Thr = [nsc_iqr_parameters_avg_HRF_T2DM_Thr iqr(total_average_HRF_parameters_T2DM_Thr(par,11:end),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_nsc_avg_HRF_CNT_Thr(1,par)==0 && sw_nsc_avg_HRF_T2DM_Thr(1,par)==0                        
        nsc_ttest_indexes_avg_HRF_thr = [nsc_ttest_indexes_avg_HRF_thr par];                            % indexes for the corresponding HRF parameters from the average HRF with normal distribution in the Thr condition in the negative signal change ROIs
        [h,p,~,stats] = ttest2(total_average_HRF_parameters_CNT_Thr(par,11:end)',total_average_HRF_parameters_T2DM_Thr(par,11:end)');
        nsc_ttest_avg_HRF_thr = [nsc_ttest_avg_HRF_thr [h;p;stats.tstat]];                              % T-test results of the HRF parameters from the average HRF in the Thr condition in the negative signal change ROIs
    else                                                                                            
        nsc_wilcoxon_indexes_avg_HRF_thr = [nsc_wilcoxon_indexes_avg_HRF_thr par];                      % indexes for the corresponding HRF parameters from the average HRF with non-normal distribution in the Thr condition in the negative signal change ROIs
        [p,h,stats] = ranksum(total_average_HRF_parameters_CNT_Thr(par,11:end)',total_average_HRF_parameters_T2DM_Thr(par,11:end)');
        nsc_wilcoxon_avg_HRF_thr = [nsc_wilcoxon_avg_HRF_thr [h;p;stats.zval]];                         % Wilcoxon test results of the HRF parameters from the average HRF in the Thr condition in the negative signal change ROIs
    end

        
   % ««««««««««««««««««««««««««« Median HRF »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    
        
   % ******************** Positive Signal Change ROIs ********************

       
    [h,p,stats]=swtest(total_median_HRF_parameters_CNT_Thr(par,1:10),0.05,true);                % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of CNT subjects in the Thr condition for each ROI of the positive signal change ROIs
    sw_psc_median_HRF_CNT_Thr(:,par) = [h;p;stats];
   

    if h==0                                                                                     
        
        % Indexes for the corresponding median HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_median_HRF_CNT_Thr = [psc_avg_std_indexes_median_HRF_CNT_Thr par];

        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the
        % positive signal change ROIs
        psc_avg_parameters_median_HRF_CNT_Thr = [psc_avg_parameters_median_HRF_CNT_Thr mean(total_median_HRF_parameters_CNT_Thr(par,1:10),2)];
        psc_std_parameters_median_HRF_CNT_Thr = [psc_std_parameters_median_HRF_CNT_Thr std(total_median_HRF_parameters_CNT_Thr(par,1:10),0,2)];
        
    else                                                                                        
        
        % Indexes for the corresponding median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Thr condition in
        % the positive signal change ROIs        
        psc_median_iqr_indexes_median_HRF_CNT_Thr = [psc_median_iqr_indexes_median_HRF_CNT_Thr par];

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Thr condition in
        % the positive signal change ROIs        
        psc_median_parameters_median_HRF_CNT_Thr = [psc_median_parameters_median_HRF_CNT_Thr median(total_median_HRF_parameters_CNT_Thr(par,1:10),2)];
        psc_iqr_parameters_median_HRF_CNT_Thr = [psc_iqr_parameters_median_HRF_CNT_Thr iqr(total_median_HRF_parameters_CNT_Thr(par,1:10),2)];

    end
    
    
    [h,p,stats] = swtest(total_median_HRF_parameters_T2DM_Thr(par,1:10), 0.05, true);           % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of T2DM subjects in the Thr condition for each ROI of the positive signal change ROIs                     
    sw_psc_median_HRF_T2DM_Thr(:,par) = [h;p;stats];

    
    if h==0 

        % Indexes for the corresponding median HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_median_HRF_T2DM_Thr = [psc_avg_std_indexes_median_HRF_T2DM_Thr par];
        
        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the
        % positive signal change ROIs        
        psc_avg_parameters_median_HRF_T2DM_Thr = [psc_avg_parameters_median_HRF_T2DM_Thr mean(total_median_HRF_parameters_T2DM_Thr(par,1:10),2)];
        psc_std_parameters_median_HRF_T2DM_Thr = [psc_std_parameters_median_HRF_T2DM_Thr std(total_median_HRF_parameters_T2DM_Thr(par,1:10),0,2)];
  
    else
        
        % Indexes for the corresponding median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Thr condition in
        % the positive signal change ROIs 
        psc_median_iqr_indexes_median_HRF_T2DM_Thr = [psc_median_iqr_indexes_median_HRF_T2DM_Thr par];       
        
        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Thr condition in
        % the positive signal change ROIs     
        psc_median_parameters_median_HRF_T2DM_Thr = [psc_median_parameters_median_HRF_T2DM_Thr median(total_median_HRF_parameters_T2DM_Thr(par,1:10),2)];
        psc_iqr_parameters_median_HRF_T2DM_Thr = [psc_iqr_parameters_median_HRF_T2DM_Thr iqr(total_median_HRF_parameters_T2DM_Thr(par,1:10),2)];

    end


    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_psc_median_HRF_CNT_Thr(1,par)==0 && sw_psc_median_HRF_T2DM_Thr(1,par)==0                    % in case of a normal distribution in the same HRF parameter from different groups --> Two sample T-test
        psc_ttest_indexes_median_HRF_thr = [psc_ttest_indexes_median_HRF_thr par];                        % indexes for the corresponding HRF parameters from the median HRF with normal distribution in the Thr condition in the positive signal change ROIs
        [h,p,~,stats] = ttest2(total_median_HRF_parameters_CNT_Thr(par,1:10)',total_median_HRF_parameters_T2DM_Thr(par,1:10)');
        psc_ttest_median_HRF_thr = [psc_ttest_median_HRF_thr [h;p;stats.tstat]];                          % T-test results of the HRF parameters from the median HRF in the Thr condition in the positive signal change ROIs
    else                                                                                        % in case of a non-normal distribution in the same HRF parameter from different groups --> Wilcoxon ranksum test
        psc_wilcoxon_indexes_median_HRF_thr = [psc_wilcoxon_indexes_median_HRF_thr par];                  % indexes for the corresponding HRF parameters from the median HRF with non-normal distribution in the Thr condition in the positive signal change ROIs
        [p,h,stats] = ranksum(total_median_HRF_parameters_CNT_Thr(par,1:10)',total_median_HRF_parameters_T2DM_Thr(par,1:10)');
        psc_wilcoxon_median_HRF_thr = [psc_wilcoxon_median_HRF_thr [h;p;stats.zval]];                     % Wilcoxon test results of the HRF parameters from the median HRF in the Thr condition in the positive signal change ROIs
    end   
    
    
    
    % ****************** Negative Signal Change ROIs ********************

    
    [h,p,stats]=swtest(total_median_HRF_parameters_CNT_Thr(par,11:end),0.05,true);              % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of CNT subjects in the Thr condition for each ROI of the negative signal change ROIs                           
    sw_nsc_median_HRF_CNT_Thr(:,par) = [h;p;stats];
    
   
    if h==0                                                                                     
        
        % Indexes for the corresponding median HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the
        % negative signal change ROIs
        nsc_avg_std_indexes_median_HRF_CNT_Thr = [nsc_avg_std_indexes_median_HRF_CNT_Thr par];

        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of CNT subjects in the Thr condition in the
        % negative signal change ROIs
        nsc_avg_parameters_median_HRF_CNT_Thr = [nsc_avg_parameters_median_HRF_CNT_Thr mean(total_median_HRF_parameters_CNT_Thr(par,11:end),2)];
        nsc_std_parameters_median_HRF_CNT_Thr = [nsc_std_parameters_median_HRF_CNT_Thr std(total_median_HRF_parameters_CNT_Thr(par,11:end),0,2)];

    else                                                                                        
        
        % Indexes for the corresponding HRF parameters with non-normal 
        % distribution from the average and median HRF of CNT subjects in 
        % the Thr condition in negative signal change ROIs        
        nsc_median_iqr_indexes_median_HRF_CNT_Thr = [nsc_median_iqr_indexes_median_HRF_CNT_Thr par];

        % Median + interquartile range of the HRF parameters with 
        % non-normal distribution from the average and median HRF of CNT 
        % subjects in the Thr condition in negative signal change ROIs        
        nsc_median_parameters_median_HRF_CNT_Thr = [nsc_median_parameters_median_HRF_CNT_Thr median(total_median_HRF_parameters_CNT_Thr(par,11:end),2)];
        nsc_iqr_parameters_median_HRF_CNT_Thr = [nsc_iqr_parameters_median_HRF_CNT_Thr iqr(total_median_HRF_parameters_CNT_Thr(par,11:end),2)];

    end


    [h,p,stats] = swtest(total_median_HRF_parameters_T2DM_Thr(par,11:end), 0.05, true);             % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of T2DM subjects in the Thr condition for each ROI of the negative signal change ROIs                              
    sw_nsc_median_HRF_T2DM_Thr(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of T2DM subjects in the Thr condition in the
        % negative signal change ROIs
        nsc_avg_std_indexes_median_HRF_T2DM_Thr = [nsc_avg_std_indexes_median_HRF_T2DM_Thr par];

        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of T2DM subjects in the Thr condition in the
        % negative signal change ROIs        
        nsc_avg_parameters_median_HRF_T2DM_Thr = [nsc_avg_parameters_median_HRF_T2DM_Thr mean(total_median_HRF_parameters_T2DM_Thr(par,11:end),2)];
        nsc_std_parameters_median_HRF_T2DM_Thr = [nsc_std_parameters_median_HRF_T2DM_Thr std(total_median_HRF_parameters_T2DM_Thr(par,11:end),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding median HRF's parameters with  
        % non-normal distribution of T2DM subjects in the Thr condition in
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_median_HRF_T2DM_Thr = [nsc_median_iqr_indexes_median_HRF_T2DM_Thr par];       

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Thr condition in 
        % the negative signal change ROIs  
        nsc_median_parameters_median_HRF_T2DM_Thr = [nsc_median_parameters_median_HRF_T2DM_Thr median(total_median_HRF_parameters_T2DM_Thr(par,11:end),2)];
        nsc_iqr_parameters_median_HRF_T2DM_Thr = [nsc_iqr_parameters_median_HRF_T2DM_Thr iqr(total_median_HRF_parameters_T2DM_Thr(par,11:end),2)];

    end
    
        
    if sw_nsc_median_HRF_CNT_Thr(1,par)==0 && sw_nsc_median_HRF_T2DM_Thr(1,par)==0                  
        nsc_ttest_indexes_median_HRF_thr = [nsc_ttest_indexes_median_HRF_thr par];                      % indexes for the corresponding HRF parameters from the median HRF with normal distribution in the Thr condition in the negative signal change ROIs
        [h,p,~,stats] = ttest2(total_median_HRF_parameters_CNT_Thr(par,11:end)',total_median_HRF_parameters_T2DM_Thr(par,11:end)');
        nsc_ttest_median_HRF_thr = [nsc_ttest_median_HRF_thr [h;p;stats.tstat]];                        % T-test results of the HRF parameters from the median HRF in the Thr condition in the negative signal change ROIs
    else                                                                                            
        nsc_wilcoxon_indexes_median_HRF_thr = [nsc_wilcoxon_indexes_median_HRF_thr par];                % indexes for the corresponding HRF parameters from the median HRF with non-normal distribution in the Thr condition in the negative signal change ROIs
        [p,h,stats] = ranksum(total_median_HRF_parameters_CNT_Thr(par,11:end)',total_median_HRF_parameters_T2DM_Thr(par,11:end)');
        nsc_wilcoxon_median_HRF_thr = [nsc_wilcoxon_median_HRF_thr [h;p;stats.zval]];                   % Wilcoxon test results of the HRF parameters from the median HRF in the Thr condition in the negative signal change ROIs
    end   
    
    
    % --------------------- SUBMAXIMUM CONDITION ------------------------

    
    % «««««««««««««««««««««««««« Average HRF »»»»»»»»»»»»»»»»»»»»»»»»»»»»
    
    
    % ****************** Positive Signal Change ROIs ********************

    
    [h,p,stats]=swtest(total_average_HRF_parameters_CNT_Sub(par,1:10),0.05,true);                   % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Sub condition for each ROI of the positive signal change ROIs    
    sw_psc_avg_HRF_CNT_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of CNT subjects in the Sub condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_avg_HRF_CNT_Sub = [psc_avg_std_indexes_avg_HRF_CNT_Sub par];        
        
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of CNT subjects in the Sub condition in the
        % positive signal change ROIs  
        psc_avg_parameters_avg_HRF_CNT_Sub = [psc_avg_parameters_avg_HRF_CNT_Sub mean(total_average_HRF_parameters_CNT_Sub(par,1:10),2)];      
        psc_std_parameters_avg_HRF_CNT_Sub = [psc_std_parameters_avg_HRF_CNT_Sub std(total_average_HRF_parameters_CNT_Sub(par,1:10),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of CNT subjects in the Sub condition in 
        % the positive signal change ROIs          
        psc_median_iqr_indexes_avg_HRF_CNT_Sub = [psc_median_iqr_indexes_avg_HRF_CNT_Sub par];

        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of CNT subjects in the Sub condition in 
        % the positive signal change ROIs 
        psc_median_parameters_avg_HRF_CNT_Sub = [psc_median_parameters_avg_HRF_CNT_Sub median(total_average_HRF_parameters_CNT_Sub(par,1:10),2)];        
        psc_iqr_parameters_avg_HRF_CNT_Sub = [psc_iqr_parameters_avg_HRF_CNT_Sub iqr(total_average_HRF_parameters_CNT_Sub(par,1:10),2)];

    end
    
    
    [h,p,stats] = swtest(total_average_HRF_parameters_T2DM_Sub(par,1:10), 0.05, true);              % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of T2DM subjects in the Sub condition for each ROI of the positive signal change ROIs                                                                         
    sw_psc_avg_HRF_T2DM_Sub(:,par) = [h;p;stats];
      
    
    if h==0                                                                                         

        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of T2DM subjects in the Sub condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_avg_HRF_T2DM_Sub = [psc_avg_std_indexes_avg_HRF_T2DM_Sub par];
        
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of T2DM subjects in the Sub condition in the
        % positive signal change ROIs  
        psc_avg_parameters_avg_HRF_T2DM_Sub = [psc_avg_parameters_avg_HRF_T2DM_Sub mean(total_average_HRF_parameters_T2DM_Sub(par,1:10),2)];
        psc_std_parameters_avg_HRF_T2DM_Sub = [psc_std_parameters_avg_HRF_T2DM_Sub std(total_average_HRF_parameters_T2DM_Sub(par,1:10),0,2)];

    else                                                                                            

        % Indexes for the corresponding average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Sub condition in 
        % the positive signal change ROIs  
        psc_median_iqr_indexes_avg_HRF_T2DM_Sub = [psc_median_iqr_indexes_avg_HRF_T2DM_Sub par];

        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Sub condition in
        % the positive signal change ROIs        
        psc_median_parameters_avg_HRF_T2DM_Sub = [psc_median_parameters_avg_HRF_T2DM_Sub median(total_average_HRF_parameters_T2DM_Sub(par,1:10),2)];
        psc_iqr_parameters_avg_HRF_T2DM_Sub = [psc_iqr_parameters_avg_HRF_T2DM_Sub iqr(total_average_HRF_parameters_T2DM_Sub(par,1:10),2)];
        
    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_psc_avg_HRF_CNT_Sub(1,par)==0 && sw_psc_avg_HRF_T2DM_Sub(1,par)==0                        
        psc_ttest_indexes_avg_HRF_sub = [psc_ttest_indexes_avg_HRF_sub par];                            % indexes for the corresponding HRF parameters from the average HRF with normal distribution in the Sub condition in the positive signal change ROIs
        [h,p,~,stats] = ttest2(total_average_HRF_parameters_CNT_Sub(par,1:10)',total_average_HRF_parameters_T2DM_Sub(par,1:10)');
        psc_ttest_avg_HRF_sub = [psc_ttest_avg_HRF_sub [h;p;stats.tstat]];                              % T-test results of the HRF parameters from the average HRF in the Sub condition in the positive signal change ROIs
    else                                                                                            
        psc_wilcoxon_indexes_avg_HRF_sub = [psc_wilcoxon_indexes_avg_HRF_sub par];                      % indexes for the corresponding HRF parameters from the average HRF with non-normal distribution in the Sub condition in the positive signal change ROIs
        [p,h,stats] = ranksum(total_average_HRF_parameters_CNT_Sub(par,1:10)',total_average_HRF_parameters_T2DM_Sub(par,1:10)');
        psc_wilcoxon_avg_HRF_sub = [psc_wilcoxon_avg_HRF_sub [h;p;stats.zval]];                         % Wilcoxon test results of the HRF parameters from the average HRF in the Sub condition in the positive signal change ROIs
    end
    
    
    
    % ****************** Negative Signal Change ROIs ********************

   
    [h,p,stats]=swtest(total_average_HRF_parameters_CNT_Sub(par,11:end),0.05,true);                 % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Sub condition for each ROI of the negative signal change ROIs
    sw_nsc_avg_HRF_CNT_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of CNT subjects in the Sub condition in the
        % negative signal change ROIs
        nsc_avg_std_indexes_avg_HRF_CNT_Sub = [nsc_avg_std_indexes_avg_HRF_CNT_Sub par];
        
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of CNT subjects in the Sub condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_avg_HRF_CNT_Sub = [nsc_avg_parameters_avg_HRF_CNT_Sub mean(total_average_HRF_parameters_CNT_Sub(par,11:end),2)];                
        nsc_std_parameters_avg_HRF_CNT_Sub = [nsc_std_parameters_avg_HRF_CNT_Sub std(total_average_HRF_parameters_CNT_Sub(par,11:end),0,2)];

    else                                                                                             
        
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of CNT subjects in the Sub condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_avg_HRF_CNT_Sub = [nsc_median_iqr_indexes_avg_HRF_CNT_Sub par];
        
        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of CNT subjects in the Sub condition in 
        % the negative signal change ROIs    
        nsc_median_parameters_avg_HRF_CNT_Sub = [nsc_median_parameters_avg_HRF_CNT_Sub median(total_average_HRF_parameters_CNT_Sub(par,11:end),2)];
        nsc_iqr_parameters_avg_HRF_CNT_Sub = [nsc_iqr_parameters_avg_HRF_CNT_Sub iqr(total_average_HRF_parameters_CNT_Sub(par,11:end),2)];
  
    end
    
    
    [h,p,stats] = swtest(total_average_HRF_parameters_T2DM_Sub(par,11:end), 0.05, true);            % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of T2DM subjects in the Sub condition for each ROI of the negative signal change ROIs                              
    sw_nsc_avg_HRF_T2DM_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of T2DM subjects in the Sub condition in the
        % negative signal change ROIs        
        nsc_avg_std_indexes_avg_HRF_T2DM_Sub = [nsc_avg_std_indexes_avg_HRF_T2DM_Sub par];

        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of T2DM subjects in the Sub condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_avg_HRF_T2DM_Sub = [nsc_avg_parameters_avg_HRF_T2DM_Sub mean(total_average_HRF_parameters_T2DM_Sub(par,11:end),2)];
        nsc_std_parameters_avg_HRF_T2DM_Sub = [nsc_std_parameters_avg_HRF_T2DM_Sub std(total_average_HRF_parameters_T2DM_Sub(par,11:end),0,2)];
    
    else                                                                                            
        
        % Indexes for the corresponding average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Sub condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_avg_HRF_T2DM_Sub = [nsc_median_iqr_indexes_avg_HRF_T2DM_Sub par];       

        % Median + interquartile range of the average HRF's parameters with
        % non-normal distribution of T2DM subjects in the Sub condition in 
        % the negative signal change ROIs  
        nsc_median_parameters_avg_HRF_T2DM_Sub = [nsc_median_parameters_avg_HRF_T2DM_Sub median(total_average_HRF_parameters_T2DM_Sub(par,11:end),2)];
        nsc_iqr_parameters_avg_HRF_T2DM_Sub = [nsc_iqr_parameters_avg_HRF_T2DM_Sub iqr(total_average_HRF_parameters_T2DM_Sub(par,11:end),2)];
        
    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_nsc_avg_HRF_CNT_Sub(1,par)==0 && sw_nsc_avg_HRF_T2DM_Sub(1,par)==0                        
        nsc_ttest_indexes_avg_HRF_sub = [nsc_ttest_indexes_avg_HRF_sub par];                            % indexes for the corresponding HRF parameters from the average HRF with normal distribution in the Sub condition in the negative signal change ROIs
        [h,p,~,stats] = ttest2(total_average_HRF_parameters_CNT_Sub(par,11:end)',total_average_HRF_parameters_T2DM_Sub(par,11:end)');
        nsc_ttest_avg_HRF_sub = [nsc_ttest_avg_HRF_sub [h;p;stats.tstat]];                              % T-test results of the HRF parameters from the average HRF in the Sub condition in the negative signal change ROIs
    else                                                                                            
        nsc_wilcoxon_indexes_avg_HRF_sub = [nsc_wilcoxon_indexes_avg_HRF_sub par];                      % indexes for the corresponding HRF parameters from the average HRF with non-normal distribution in the Sub condition in the negative signal change ROIs
        [p,h,stats] = ranksum(total_average_HRF_parameters_CNT_Sub(par,11:end)',total_average_HRF_parameters_T2DM_Sub(par,11:end)');
        nsc_wilcoxon_avg_HRF_sub = [nsc_wilcoxon_avg_HRF_sub [h;p;stats.zval]];                         % Wilcoxon test results of the HRF parameters from the average HRF in the Sub condition in the negative signal change ROIs
    end
    

        
    % «««««««««««««««««««««««««« Median HRF »»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    
        
    % ****************** Positive Signal Change ROIs ********************

    
    [h,p,stats]=swtest(total_median_HRF_parameters_CNT_Sub(par,1:10),0.05,true);                    % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of CNT subjects in the Sub condition for each ROI of the positive signal change ROIs
    sw_psc_median_HRF_CNT_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of CNT subjects in the Sub condition in the positive
        % signal change ROIs
        psc_avg_std_indexes_median_HRF_CNT_Sub = [psc_avg_std_indexes_median_HRF_CNT_Sub par];
        
        
        % Average + standard-deviation of the median's HRF parameters with  
        % normal distribution of CNT subjects in the Sub condition in the
        % positive signal change ROIs  
        psc_avg_parameters_median_HRF_CNT_Sub = [psc_avg_parameters_median_HRF_CNT_Sub mean(total_median_HRF_parameters_CNT_Sub(par,1:10),2)];
        psc_std_parameters_median_HRF_CNT_Sub = [psc_std_parameters_median_HRF_CNT_Sub std(total_median_HRF_parameters_CNT_Sub(par,1:10),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding median HRF's parameters with  
        % non-normal distribution of CNT subjects in the Sub condition in
        % the positive signal change ROIs          
        psc_median_iqr_indexes_median_HRF_CNT_Sub = [psc_median_iqr_indexes_median_HRF_CNT_Sub par];

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Sub condition in
        % the positive signal change ROIs 
        psc_median_parameters_median_HRF_CNT_Sub = [psc_median_parameters_median_HRF_CNT_Sub median(total_median_HRF_parameters_CNT_Sub(par,1:10),2)];
        psc_iqr_parameters_median_HRF_CNT_Sub = [psc_iqr_parameters_median_HRF_CNT_Sub iqr(total_median_HRF_parameters_CNT_Sub(par,1:10),2)];

    end
    
      
    [h,p,stats] = swtest(total_median_HRF_parameters_T2DM_Sub(par,1:10), 0.05, true);               % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of T2DM subjects in the Sub condition for each ROI of the positive signal change ROIs                                                                         
    sw_psc_median_HRF_T2DM_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         

        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of T2DM subjects in the Sub condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_median_HRF_T2DM_Sub = [psc_avg_std_indexes_median_HRF_T2DM_Sub par];
        
        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of T2DM subjects in the Sub condition in the
        % positive signal change ROIs  
        psc_avg_parameters_median_HRF_T2DM_Sub = [psc_avg_parameters_median_HRF_T2DM_Sub mean(total_median_HRF_parameters_T2DM_Sub(par,1:10),2)];
        psc_std_parameters_median_HRF_T2DM_Sub = [psc_std_parameters_median_HRF_T2DM_Sub std(total_median_HRF_parameters_T2DM_Sub(par,1:10),0,2)];

    else                                                                                            

        % Indexes for the corresponding median HRF's parameters with  
        % non-normal distribution of T2DM subjects in the Sub condition in 
        % the positive signal change ROIs  
        psc_median_iqr_indexes_median_HRF_T2DM_Sub = [psc_median_iqr_indexes_median_HRF_T2DM_Sub par];

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Sub condition in
        % the positive signal change ROIs        
        psc_median_parameters_median_HRF_T2DM_Sub = [psc_median_parameters_median_HRF_T2DM_Sub median(total_median_HRF_parameters_T2DM_Sub(par,1:10),2)];
        psc_iqr_parameters_median_HRF_T2DM_Sub = [psc_iqr_parameters_median_HRF_T2DM_Sub iqr(total_median_HRF_parameters_T2DM_Sub(par,1:10),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups  
    if sw_psc_median_HRF_CNT_Sub(1,par)==0 && sw_psc_median_HRF_T2DM_Sub(1,par)==0                  
        psc_ttest_indexes_median_HRF_sub = [psc_ttest_indexes_median_HRF_sub par];                      % indexes for the corresponding HRF parameters from the median HRF with normal distribution in the Sub condition in the positive signal change ROIs
        [h,p,~,stats] = ttest2(total_median_HRF_parameters_CNT_Sub(par,1:10)',total_median_HRF_parameters_T2DM_Sub(par,1:10)');
        psc_ttest_median_HRF_sub = [psc_ttest_median_HRF_sub [h;p;stats.tstat]];                        % T-test results of the HRF parameters from the median HRF in the Sub condition in the positive signal change ROIs
    else                                                                                            
        psc_wilcoxon_indexes_median_HRF_sub = [psc_wilcoxon_indexes_median_HRF_sub par];                % indexes for the corresponding HRF parameters from the median HRF with non-normal distribution in the Sub condition in the positive signal change ROIs
        [p,h,stats] = ranksum(total_median_HRF_parameters_CNT_Sub(par,1:10)',total_median_HRF_parameters_T2DM_Sub(par,1:10)');
        psc_wilcoxon_median_HRF_sub = [psc_wilcoxon_median_HRF_sub [h;p;stats.zval]];                   % Wilcoxon test results of the HRF parameters from the median HRF in the Sub condition in the positive signal change ROIs
    end

       
    
    % ****************** Negative Signal Change ROIs ********************

    
    [h,p,stats]=swtest(total_median_HRF_parameters_CNT_Sub(par,11:end),0.05,true);                  % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of CNT subjects in the Sub condition for each ROI of the negative signal change ROIs
    sw_nsc_median_HRF_CNT_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of CNT subjects in the Sub condition in the negative
        % signal change ROIs
        nsc_avg_std_indexes_median_HRF_CNT_Sub = [nsc_avg_std_indexes_median_HRF_CNT_Sub par];

        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of CNT subjects in the Sub condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_median_HRF_CNT_Sub = [nsc_avg_parameters_median_HRF_CNT_Sub mean(total_median_HRF_parameters_CNT_Sub(par,11:end),2)];
        nsc_std_parameters_median_HRF_CNT_Sub = [nsc_std_parameters_median_HRF_CNT_Sub std(total_median_HRF_parameters_CNT_Sub(par,11:end),0,2)];

    else                                                                                             
        
        % Indexes for the corresponding median HRF's parameters with  
        % non-normal distribution of CNT subjects in the Sub condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_median_HRF_CNT_Sub = [nsc_median_iqr_indexes_median_HRF_CNT_Sub par];

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Sub condition in 
        % the negative signal change ROIs    
        nsc_median_parameters_median_HRF_CNT_Sub = [nsc_median_parameters_median_HRF_CNT_Sub median(total_median_HRF_parameters_CNT_Sub(par,11:end),2)];
        nsc_iqr_parameters_median_HRF_CNT_Sub = [nsc_iqr_parameters_median_HRF_CNT_Sub iqr(total_median_HRF_parameters_CNT_Sub(par,11:end),2)];
  
    end
    
      
    [h,p,stats] = swtest(total_median_HRF_parameters_T2DM_Sub(par,11:end), 0.05, true);             % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of T2DM subjects in the Sub condition for each ROI of the negative signal change ROIs                               
    sw_nsc_median_HRF_T2DM_Sub(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of T2DM subjects in the Sub condition in the
        % negative signal change ROIs        
        nsc_avg_std_indexes_median_HRF_T2DM_Sub = [nsc_avg_std_indexes_median_HRF_T2DM_Sub par];

        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of T2DM subjects in the Sub condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_median_HRF_T2DM_Sub = [nsc_avg_parameters_median_HRF_T2DM_Sub mean(total_median_HRF_parameters_T2DM_Sub(par,11:end),2)];
        nsc_std_parameters_median_HRF_T2DM_Sub = [nsc_std_parameters_median_HRF_T2DM_Sub std(total_median_HRF_parameters_T2DM_Sub(par,11:end),0,2)];
    
    else                                                                                            
        
        % Indexes for the corresponding median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Sub condition in
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_median_HRF_T2DM_Sub = [nsc_median_iqr_indexes_median_HRF_T2DM_Sub par];       

        % Median + interquartile range of the median HRF's parameters with
        % non-normal distribution of T2DM subjects in the Sub condition in 
        % the negative signal change ROIs  
        nsc_median_parameters_median_HRF_T2DM_Sub = [nsc_median_parameters_median_HRF_T2DM_Sub median(total_median_HRF_parameters_T2DM_Sub(par,11:end),2)];
        nsc_iqr_parameters_median_HRF_T2DM_Sub = [nsc_iqr_parameters_median_HRF_T2DM_Sub iqr(total_median_HRF_parameters_T2DM_Sub(par,11:end),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups    
    if sw_nsc_median_HRF_CNT_Sub(1,par)==0 && sw_nsc_median_HRF_T2DM_Sub(1,par)==0                  
        nsc_ttest_indexes_median_HRF_sub = [nsc_ttest_indexes_median_HRF_sub par];                      % indexes for the corresponding HRF parameters from the median HRF with normal distribution in the Sub condition in the negative signal change ROIs
        [h,p,~,stats] = ttest2(total_median_HRF_parameters_CNT_Sub(par,11:end)',total_median_HRF_parameters_T2DM_Sub(par,11:end)');
        nsc_ttest_median_HRF_sub = [nsc_ttest_median_HRF_sub [h;p;stats.tstat]];                        % T-test results of the HRF parameters from the median HRF in the Sub condition in the negative signal change ROIs
    else                                                                                            
        nsc_wilcoxon_indexes_median_HRF_sub = [nsc_wilcoxon_indexes_median_HRF_sub par];                % indexes for the corresponding HRF parameters from the median HRF with non-normal distribution in the Sub condition in the negative signal change ROIs
        [p,h,stats] = ranksum(total_median_HRF_parameters_CNT_Sub(par,11:end)',total_median_HRF_parameters_T2DM_Sub(par,11:end)');
        nsc_wilcoxon_median_HRF_sub = [nsc_wilcoxon_median_HRF_sub [h;p;stats.zval]];                   % Wilcoxon test results of the HRF parameters from the median HRF in the Sub condition in the negative signal change ROIs
    end
   
    

    % ----------------------- AVERAGE CONDITION -------------------------

    
    
    % «««««««««««««««««««««««««« Average HRF »»»»»»»»»»»»»»»»»»»»»»»»»»»»
    
    
    % ******************* Positive Signal Change ROIs *******************

    
    [h,p,stats]=swtest(total_average_HRF_parameters_CNT_Avg(par,1:10),0.05,true);                   % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Avg condition for each ROI of the positive signal change ROIs    
    sw_psc_avg_HRF_CNT_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of CNT subjects in the Avg condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_avg_HRF_CNT_Avg = [psc_avg_std_indexes_avg_HRF_CNT_Avg par];        
        
        % Average + standard-deviation of the average HRF's parameters with 
        % normal distribution of CNT subjects in the Avg condition in the
        % positive signal change ROIs  
        psc_avg_parameters_avg_HRF_CNT_Avg = [psc_avg_parameters_avg_HRF_CNT_Avg mean(total_average_HRF_parameters_CNT_Avg(par,1:10),2)];      
        psc_std_parameters_avg_HRF_CNT_Avg = [psc_std_parameters_avg_HRF_CNT_Avg std(total_average_HRF_parameters_CNT_Avg(par,1:10),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding average HRF's parameters with
        % non-normal distribution of CNT subjects in the Avg condition in
        % the positive signal change ROIs          
        psc_median_iqr_indexes_avg_HRF_CNT_Avg = [psc_median_iqr_indexes_avg_HRF_CNT_Avg par];

        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of CNT subjects in the Avg condition in
        % the positive signal change ROIs 
        psc_median_parameters_avg_HRF_CNT_Avg = [psc_median_parameters_avg_HRF_CNT_Avg median(total_average_HRF_parameters_CNT_Avg(par,1:10),2)];                    
        psc_iqr_parameters_avg_HRF_CNT_Avg = [psc_iqr_parameters_avg_HRF_CNT_Avg iqr(total_average_HRF_parameters_CNT_Avg(par,1:10),2)];

    end
    
    
    [h,p,stats] = swtest(total_average_HRF_parameters_T2DM_Avg(par,1:10), 0.05, true);              % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of T2DM subjects in the Avg condition for each ROI of the positive signal change ROIs                                                                         
    sw_psc_avg_HRF_T2DM_Avg(:,par) = [h;p;stats];    
    
    
    if h==0                                                                                         

        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of T2DM subjects in the Avg condition in the
        % positive signal change ROIs
        psc_avg_std_indexes_avg_HRF_T2DM_Avg = [psc_avg_std_indexes_avg_HRF_T2DM_Avg par];
        
        % Average + standard-deviation of the average HRF's parameters with 
        % normal distribution of T2DM subjects in the Avg condition in the
        % positive signal change ROIs  
        psc_avg_parameters_avg_HRF_T2DM_Avg = [psc_avg_parameters_avg_HRF_T2DM_Avg mean(total_average_HRF_parameters_T2DM_Avg(par,1:10),2)];    
        psc_std_parameters_avg_HRF_T2DM_Avg = [psc_std_parameters_avg_HRF_T2DM_Avg std(total_average_HRF_parameters_T2DM_Avg(par,1:10),0,2)];

    else                                                                                            

        % Indexes for the corresponding average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the positive signal change ROIs  
        psc_median_iqr_indexes_avg_HRF_T2DM_Avg = [psc_median_iqr_indexes_avg_HRF_T2DM_Avg par];

        % Median + interquartile range of the average HRF's parameters with
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the positive signal change ROIs        
        psc_median_parameters_avg_HRF_T2DM_Avg = [psc_median_parameters_avg_HRF_T2DM_Avg median(total_average_HRF_parameters_T2DM_Avg(par,1:10),2)];
        psc_iqr_parameters_avg_HRF_T2DM_Avg = [psc_iqr_parameters_avg_HRF_T2DM_Avg iqr(total_average_HRF_parameters_T2DM_Avg(par,1:10),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_psc_avg_HRF_CNT_Avg(1,par)==0 && sw_psc_avg_HRF_T2DM_Avg(1,par)==0                        
        psc_ttest_indexes_avg_HRF_avg = [psc_ttest_indexes_avg_HRF_avg par];                            % indexes for the corresponding HRF parameters from the average HRF with normal distribution in the Avg condition in the positive signal change ROIs
        [h,p,~,stats] = ttest2(total_average_HRF_parameters_CNT_Avg(par,1:10)',total_average_HRF_parameters_T2DM_Avg(par,1:10)');
        psc_ttest_avg_HRF_avg = [psc_ttest_avg_HRF_avg [h;p;stats.tstat]];                              % T-test results of the HRF parameters from the average HRF in the Avg condition in the positive signal change ROIs
    else                                                                                            
        psc_wilcoxon_indexes_avg_HRF_avg = [psc_wilcoxon_indexes_avg_HRF_avg par];                      % indexes for the corresponding HRF parameters from the average HRF with non-normal distribution in the Avg condition in the positive signal change ROIs
        [p,h,stats] = ranksum(total_average_HRF_parameters_CNT_Avg(par,1:10)',total_average_HRF_parameters_T2DM_Avg(par,1:10)');
        psc_wilcoxon_avg_HRF_avg = [psc_wilcoxon_avg_HRF_avg [h;p;stats.zval]];                         % Wilcoxon test results of the HRF parameters from the average HRF in the Avg condition in the positive signal change ROIs
    end
    
    
    
    % ****************** Negative Signal Change ROIs ********************

   
    [h,p,stats]=swtest(total_average_HRF_parameters_CNT_Avg(par,11:end),0.05,true);                 % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of CNT subjects in the Avg condition for each ROI of the negative signal change ROIs
    sw_nsc_avg_HRF_CNT_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of CNT subjects in the Avg condition in the
        % negative signal change ROIs
        nsc_avg_std_indexes_avg_HRF_CNT_Avg = [nsc_avg_std_indexes_avg_HRF_CNT_Avg par];
        
        % Average + standard-deviation of the average HRF's parameters with  
        % normal distribution of CNT subjects in the Avg condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_avg_HRF_CNT_Avg = [nsc_avg_parameters_avg_HRF_CNT_Avg mean(total_average_HRF_parameters_CNT_Avg(par,11:end),2)];                
        nsc_std_parameters_avg_HRF_CNT_Avg = [nsc_std_parameters_avg_HRF_CNT_Avg std(total_average_HRF_parameters_CNT_Avg(par,11:end),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of CNT subjects in the Avg condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_avg_HRF_CNT_Avg = [nsc_median_iqr_indexes_avg_HRF_CNT_Avg par];
        
        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of CNT subjects in the Avg condition in 
        % the negative signal change ROIs    
        nsc_median_parameters_avg_HRF_CNT_Avg = [nsc_median_parameters_avg_HRF_CNT_Avg median(total_average_HRF_parameters_CNT_Avg(par,11:end),2)];
        nsc_iqr_parameters_avg_HRF_CNT_Avg = [nsc_iqr_parameters_avg_HRF_CNT_Avg iqr(total_average_HRF_parameters_CNT_Avg(par,11:end),2)];
        
    end
    
    
    [h,p,stats] = swtest(total_average_HRF_parameters_T2DM_Avg(par,11:end), 0.05, true);            % estimates h, p and stats of the Shapiro-Wilk test for each average HRF parameter of T2DM subjects in the Avg condition for each ROI of the negative signal change ROIs                              
    sw_nsc_avg_HRF_T2DM_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding average HRF's parameters with  
        % normal distribution of T2DM subjects in the Avg condition in the
        % negative signal change ROIs        
        nsc_avg_std_indexes_avg_HRF_T2DM_Avg = [nsc_avg_std_indexes_avg_HRF_T2DM_Avg par];
        
        % Average + standard-deviation of the average HRF's parameters with 
        % normal distribution of T2DM subjects in the Avg condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_avg_HRF_T2DM_Avg = [nsc_avg_parameters_avg_HRF_T2DM_Avg mean(total_average_HRF_parameters_T2DM_Avg(par,11:end),2)];
        nsc_std_parameters_avg_HRF_T2DM_Avg = [nsc_std_parameters_avg_HRF_T2DM_Avg std(total_average_HRF_parameters_T2DM_Avg(par,11:end),0,2)];
    
    else                                                                                            
        
        % Indexes for the corresponding average HRF's parameters with  
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_avg_HRF_T2DM_Avg = [nsc_median_iqr_indexes_avg_HRF_T2DM_Avg par];       
        
        % Median + interquartile range of the average HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the negative signal change ROIs  
        nsc_median_parameters_avg_HRF_T2DM_Avg = [nsc_median_parameters_avg_HRF_T2DM_Avg median(total_average_HRF_parameters_T2DM_Avg(par,11:end),2)];
        nsc_iqr_parameters_avg_HRF_T2DM_Avg = [nsc_iqr_parameters_avg_HRF_T2DM_Avg iqr(total_average_HRF_parameters_T2DM_Avg(par,11:end),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_nsc_avg_HRF_CNT_Avg(1,par)==0 && sw_nsc_avg_HRF_T2DM_Avg(1,par)==0                        
        nsc_ttest_indexes_avg_HRF_avg = [nsc_ttest_indexes_avg_HRF_avg par];                            % indexes for the corresponding HRF parameters from the average HRF with normal distribution in the Avg condition in the negative signal change ROIs
        [h,p,~,stats] = ttest2(total_average_HRF_parameters_CNT_Avg(par,11:end)',total_average_HRF_parameters_T2DM_Avg(par,11:end)');
        nsc_ttest_avg_HRF_avg = [nsc_ttest_avg_HRF_avg [h;p;stats.tstat]];                              % T-test results of the HRF parameters from the average HRF in the Avg condition in the negative signal change ROIs
    else                                                                                            
        nsc_wilcoxon_indexes_avg_HRF_avg = [nsc_wilcoxon_indexes_avg_HRF_avg par];                      % indexes for the corresponding HRF parameters from the average HRF with non-normal distribution in the Avg condition in the negative signal change ROIs
        [p,h,stats] = ranksum(total_average_HRF_parameters_CNT_Avg(par,11:end)',total_average_HRF_parameters_T2DM_Avg(par,11:end)');
        nsc_wilcoxon_avg_HRF_avg = [nsc_wilcoxon_avg_HRF_avg [h;p;stats.zval]];                         % Wilcoxon test results of the HRF parameters from the average HRF in the Avg condition in the negative signal change ROIs
    end
   
    
        
    % «««««««««««««««««««««««««« Median HRF »»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    
        
    % ******************* Positive Signal Change ROIs *******************

    
    [h,p,stats]=swtest(total_median_HRF_parameters_CNT_Avg(par,1:10),0.05,true);                    % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of CNT subjects in the Avg condition for each ROI of the positive signal change ROIs
    sw_psc_median_HRF_CNT_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of CNT subjects in the Avg condition in the positive
        % signal change ROIs
        psc_avg_std_indexes_median_HRF_CNT_Avg = [psc_avg_std_indexes_median_HRF_CNT_Avg par];
        
        % Average + standard-deviation of the median HRF's parameters with 
        % normal distribution of CNT subjects in the Avg condition in the
        % positive signal change ROIs  
        psc_avg_parameters_median_HRF_CNT_Avg = [psc_avg_parameters_median_HRF_CNT_Avg mean(total_median_HRF_parameters_CNT_Avg(par,1:10),2)];
        psc_std_parameters_median_HRF_CNT_Avg = [psc_std_parameters_median_HRF_CNT_Avg std(total_median_HRF_parameters_CNT_Avg(par,1:10),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Avg condition in 
        % the positive signal change ROIs          
        psc_median_iqr_indexes_median_HRF_CNT_Avg = [psc_median_iqr_indexes_median_HRF_CNT_Avg par];

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Avg condition in 
        % the positive signal change ROIs 
        psc_median_parameters_median_HRF_CNT_Avg = [psc_median_parameters_median_HRF_CNT_Avg median(total_median_HRF_parameters_CNT_Avg(par,1:10),2)];
        psc_iqr_parameters_median_HRF_CNT_Avg = [psc_iqr_parameters_median_HRF_CNT_Avg iqr(total_median_HRF_parameters_CNT_Avg(par,1:10),2)];

    end
    
    
    [h,p,stats] = swtest(total_median_HRF_parameters_T2DM_Avg(par,1:10), 0.05, true);               % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of T2DM subjects in the Avg condition for each ROI of the positive signal change ROIs                                                                         
    sw_psc_median_HRF_T2DM_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         

        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of T2DM subjects in the Avg condition in the 
        % positive signal change ROIs
        psc_avg_std_indexes_median_HRF_T2DM_Avg = [psc_avg_std_indexes_median_HRF_T2DM_Avg par];
        
        % Average + standard-deviation of the median HRF's parameters with 
        % normal distribution of T2DM subjects in the Avg condition in the
        % positive signal change ROIs  
        psc_avg_parameters_median_HRF_T2DM_Avg = [psc_avg_parameters_median_HRF_T2DM_Avg mean(total_median_HRF_parameters_T2DM_Avg(par,1:10),2)];
        psc_std_parameters_median_HRF_T2DM_Avg = [psc_std_parameters_median_HRF_T2DM_Avg std(total_median_HRF_parameters_T2DM_Avg(par,1:10),0,2)];

    else                                                                                            

        % Indexes for the corresponding median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the positive signal change ROIs  
        psc_median_iqr_indexes_median_HRF_T2DM_Avg = [psc_median_iqr_indexes_median_HRF_T2DM_Avg par];

        % Median + interquartile range of the median HRF's parameters with
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the positive signal change ROIs        
        psc_median_parameters_median_HRF_T2DM_Avg = [psc_median_parameters_median_HRF_T2DM_Avg median(total_median_HRF_parameters_T2DM_Avg(par,1:10),2)];
        psc_iqr_parameters_median_HRF_T2DM_Avg = [psc_iqr_parameters_median_HRF_T2DM_Avg iqr(total_median_HRF_parameters_T2DM_Avg(par,1:10),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups 
    if sw_psc_median_HRF_CNT_Avg(1,par)==0 && sw_psc_median_HRF_T2DM_Avg(1,par)==0                  
        psc_ttest_indexes_median_HRF_avg = [psc_ttest_indexes_median_HRF_avg par];                      % indexes for the corresponding HRF parameters from the median HRF with normal distribution in the Avg condition in the positive signal change ROIs
        [h,p,~,stats] = ttest2(total_median_HRF_parameters_CNT_Avg(par,1:10)',total_median_HRF_parameters_T2DM_Avg(par,1:10)');
        psc_ttest_median_HRF_avg = [psc_ttest_median_HRF_avg [h;p;stats.tstat]];                        % T-test results of the HRF parameters from the median HRF in the Avg condition in the positive signal change ROIs
    else                                                                                            
        psc_wilcoxon_indexes_median_HRF_avg = [psc_wilcoxon_indexes_median_HRF_avg par];                % indexes for the corresponding HRF parameters from the median HRF with non-normal distribution in the Avg condition in the positive signal change ROIs
        [p,h,stats] = ranksum(total_median_HRF_parameters_CNT_Avg(par,1:10)',total_median_HRF_parameters_T2DM_Avg(par,1:10)');
        psc_wilcoxon_median_HRF_avg = [psc_wilcoxon_median_HRF_avg [h;p;stats.zval]];                   % Wilcoxon test results of the HRF parameters from the median HRF in the Avg condition in the positive signal change ROIs
    end

    
    
    % ****************** Negative Signal Change ROIs ********************


    [h,p,stats]=swtest(total_median_HRF_parameters_CNT_Avg(par,11:end),0.05,true);                  % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of CNT subjects in the Avg condition for each ROI of the negative signal change ROIs
    sw_nsc_median_HRF_CNT_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of CNT subjects in the Avg condition in the negative
        % signal change ROIs
        nsc_avg_std_indexes_median_HRF_CNT_Avg = [nsc_avg_std_indexes_median_HRF_CNT_Avg par];

        % Average + standard-deviation of the median HRF's parameters with  
        % normal distribution of CNT subjects in the Avg condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_median_HRF_CNT_Avg = [nsc_avg_parameters_median_HRF_CNT_Avg mean(total_median_HRF_parameters_CNT_Avg(par,11:end),2)];
        nsc_std_parameters_median_HRF_CNT_Avg = [nsc_std_parameters_median_HRF_CNT_Avg std(total_median_HRF_parameters_CNT_Avg(par,11:end),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding median HRF's parameters with  
        % non-normal distribution of CNT subjects in the Avg condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_median_HRF_CNT_Avg = [nsc_median_iqr_indexes_median_HRF_CNT_Avg par];
        
        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of CNT subjects in the Avg condition in
        % the negative signal change ROIs    
        nsc_median_parameters_median_HRF_CNT_Avg = [nsc_median_parameters_median_HRF_CNT_Avg median(total_median_HRF_parameters_CNT_Avg(par,11:end),2)];
        nsc_iqr_parameters_median_HRF_CNT_Avg = [nsc_iqr_parameters_median_HRF_CNT_Avg iqr(total_median_HRF_parameters_CNT_Avg(par,11:end),2)];
  
    end
    
    
    [h,p,stats] = swtest(total_median_HRF_parameters_T2DM_Avg(par,11:end), 0.05, true);             % estimates h, p and stats of the Shapiro-Wilk test for each median HRF parameter of T2DM subjects in the Avg condition for each ROI of the negative signal change ROIs                               
    sw_nsc_median_HRF_T2DM_Avg(:,par) = [h;p;stats];
    
    
    if h==0                                                                                         
        
        % Indexes for the corresponding median HRF's parameters with normal 
        % distribution of T2DM subjects in the Avg condition in the 
        % negative signal change ROIs        
        nsc_avg_std_indexes_median_HRF_T2DM_Avg = [nsc_avg_std_indexes_median_HRF_T2DM_Avg par];

        % Average + standard-deviation of the median HRF's parameters with 
        % normal distribution of T2DM subjects in the Avg condition in the
        % negative signal change ROIs          
        nsc_avg_parameters_median_HRF_T2DM_Avg = [nsc_avg_parameters_median_HRF_T2DM_Avg mean(total_median_HRF_parameters_T2DM_Avg(par,11:end),2)];
        nsc_std_parameters_median_HRF_T2DM_Avg = [nsc_std_parameters_median_HRF_T2DM_Avg std(total_median_HRF_parameters_T2DM_Avg(par,11:end),0,2)];

    else                                                                                            
        
        % Indexes for the corresponding median HRF's parameters with  
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the negative signal change ROIs  
        nsc_median_iqr_indexes_median_HRF_T2DM_Avg = [nsc_median_iqr_indexes_median_HRF_T2DM_Avg par];       

        % Median + interquartile range of the median HRF's parameters with 
        % non-normal distribution of T2DM subjects in the Avg condition in 
        % the negative signal change ROIs  
        nsc_median_parameters_median_HRF_T2DM_Avg = [nsc_median_parameters_median_HRF_T2DM_Avg median(total_median_HRF_parameters_T2DM_Avg(par,11:end),2)];
        nsc_iqr_parameters_median_HRF_T2DM_Avg = [nsc_iqr_parameters_median_HRF_T2DM_Avg iqr(total_median_HRF_parameters_T2DM_Avg(par,11:end),2)];

    end
    
    
    % Tests if there are significant differences for the same HRF 
    % parameters between groups
    if sw_nsc_median_HRF_CNT_Avg(1,par)==0 && sw_nsc_median_HRF_T2DM_Avg(1,par)==0                  
        nsc_ttest_indexes_median_HRF_avg = [nsc_ttest_indexes_median_HRF_avg par];                      % indexes for the corresponding HRF parameters from the median HRF with normal distribution in the Avg condition in the negative signal change ROIs
        [h,p,~,stats] = ttest2(total_median_HRF_parameters_CNT_Avg(par,11:end)',total_median_HRF_parameters_T2DM_Avg(par,11:end)');
        nsc_ttest_median_HRF_avg = [nsc_ttest_median_HRF_avg [h;p;stats.tstat]];                        % T-test results of the HRF parameters from the median HRF in the Avg condition in the negative signal change ROIs
    else                                                                                            
        nsc_wilcoxon_indexes_median_HRF_avg = [nsc_wilcoxon_indexes_median_HRF_avg par];                % indexes for the corresponding HRF parameters from the median HRF with non-normal distribution in the Avg condition in the negative signal change ROIs
        [p,h,stats] = ranksum(total_median_HRF_parameters_CNT_Avg(par,11:end)',total_median_HRF_parameters_T2DM_Avg(par,11:end)');
        nsc_wilcoxon_median_HRF_avg = [nsc_wilcoxon_median_HRF_avg [h;p;stats.zval]];                   % Wilcoxon test results of the HRF parameters from the median HRF in the Avg condition in the negative signal change ROIs
    end
    
    
end


%% ADJUSTED P-VALUES

% In this section, we estimate the adjusted p-values according to the 
% Benjamini & Hochberg (1995) approach for the Two-sample T-test and the
% Wilcoxon ranksum test, using a False Discovery Rate (FDR) of 0.10

% The results are stored in the following order:
%                   1. Two-sample T-test Thr
%                   2. Two-sample T-test Sub
%                   3. Two-sample T-test Avg
%                   4. Wilcoxon ranksum test Thr
%                   5. Wilcoxon ranksum test Sub
%                   6. Wilcoxon ranksum test Avg



% ---------------------------- Average HRF -------------------------------


if (isempty(psc_ttest_avg_HRF_thr)==0 || isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 || isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 || isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                        % if one of the statistical tests takes place in all conditions
    if (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                    % if both statistical tests take place in all conditions       
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');         % 'yes' in order to see the report   
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in both Thr and Sub conditions, but only Two-sample T-test takes place in the Avg condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in both Thr and Avg conditions, but only Two-sample T-test takes place in the Sub condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in both Sub and Avg conditions
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Avg condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                             																
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
    elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Sub condition  
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');        																									   
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in both Sub and Avg conditions, but only Two-sample T-test takes place in the Thr condition     
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in both Thr and Avg conditions
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Avg condition    
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in both Thr and Sub conditions     
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in all conditions
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if only Two-sample T-test takes place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Sub condition 
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_avg_HRF_thr)==0 && isempty(psc_wilcoxon_avg_HRF_thr)==1) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if only Two-sample T-test takes place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in both Sub and Avg conditions, but only Wilcoxon ranksum test takes place in the Thr condition     
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Thr condition   
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Sub condition, but only Wilcoxon ranksum test takes place in both Thr and Avg conditions    
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Thr condition     
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in both Sub and Avg conditions and Wilcoxon ranksum test takes place in the Thr condition   
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_sub(2,:) psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==0 && isempty(psc_wilcoxon_avg_HRF_sub)==1) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in both Thr and Avg conditions 
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Wilcoxon ranksum test takes place in both Thr and Sub conditions     
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         	
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==0 && isempty(psc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in both Thr and Sub conditions   
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_ttest_avg_HRF_avg(2,:) psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_avg_HRF_thr)==1 && isempty(psc_wilcoxon_avg_HRF_thr)==0) && (isempty(psc_ttest_avg_HRF_sub)==1 && isempty(psc_wilcoxon_avg_HRF_sub)==0) && (isempty(psc_ttest_avg_HRF_avg)==1 && isempty(psc_wilcoxon_avg_HRF_avg)==0)												% if only Wilcoxon ranksum test takes place in all conditions 
        [psc_avg_HRF_bh_h,psc_avg_HRF_crit_p,psc_avg_HRF_adj_ci,psc_avg_HRF_adj_p]=fdr_bh([psc_wilcoxon_avg_HRF_thr(2,:) psc_wilcoxon_avg_HRF_sub(2,:) psc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
    end
end


if (isempty(nsc_ttest_avg_HRF_thr)==0 || isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 || isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 || isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                        % if one of the statistical tests takes place in all conditions
    if (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                    % if both statistical tests take place in all conditions       
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');         % 'yes' in order to see the report   
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in both Thr and Sub conditions, but only Two-sample T-test takes place in the Avg condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in both Thr and Avg conditions, but only Two-sample T-test takes place in the Sub condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in both Sub and Avg conditions
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Avg condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                             																
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
    elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Sub condition  
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');        																									   
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in both Sub and Avg conditions, but only Two-sample T-test takes place in the Thr condition     
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in both Thr and Avg conditions
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Avg condition    
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in both Thr and Sub conditions     
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in all conditions
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if only Two-sample T-test takes place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Sub condition 
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_avg_HRF_thr)==0 && isempty(nsc_wilcoxon_avg_HRF_thr)==1) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if only Two-sample T-test takes place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in both Sub and Avg conditions, but only Wilcoxon ranksum test takes place in the Thr condition     
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Thr condition   
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if both statistical tests take place in the Sub condition, but only Wilcoxon ranksum test takes place in both Thr and Avg conditions    
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Thr condition     
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in both Sub and Avg conditions and Wilcoxon ranksum test takes place in the Thr condition   
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_sub(2,:) nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==0 && isempty(nsc_wilcoxon_avg_HRF_sub)==1) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in both Thr and Avg conditions 
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)                                                % if both statistical tests take place in the Avg condition, but only Wilcoxon ranksum test takes place in both Thr and Sub conditions     
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');                                         	
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==0 && isempty(nsc_wilcoxon_avg_HRF_avg)==1)												% if only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in both Thr and Sub conditions   
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_ttest_avg_HRF_avg(2,:) nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_avg_HRF_thr)==1 && isempty(nsc_wilcoxon_avg_HRF_thr)==0) && (isempty(nsc_ttest_avg_HRF_sub)==1 && isempty(nsc_wilcoxon_avg_HRF_sub)==0) && (isempty(nsc_ttest_avg_HRF_avg)==1 && isempty(nsc_wilcoxon_avg_HRF_avg)==0)												% if only Wilcoxon ranksum test takes place in all conditions 
        [nsc_avg_HRF_bh_h,nsc_avg_HRF_crit_p,nsc_avg_HRF_adj_ci,nsc_avg_HRF_adj_p]=fdr_bh([nsc_wilcoxon_avg_HRF_thr(2,:) nsc_wilcoxon_avg_HRF_sub(2,:) nsc_wilcoxon_avg_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
    end
end


% ----------------------------- Median HRF ------------------------------
 

if (isempty(psc_ttest_median_HRF_thr)==0 || isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 || isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 || isempty(psc_wilcoxon_median_HRF_avg)==0)                                                      			% if one of the statistical tests takes place in all conditions
    if (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                                				% if both statistical tests take place in all conditions       
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');       % 'yes' in order to see the report   
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in both Thr and Sub conditions, but only Two-sample T-test takes place in the Avg condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in both Thr and Avg conditions, but only Two-sample T-test takes place in the Sub condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in both Sub and Avg conditions
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Avg condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                             																
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
    elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Sub condition  
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');        																									   
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in both Sub and Avg conditions, but only Two-sample T-test takes place in the Thr condition     
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in both Thr and Avg conditions
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Avg condition    
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in both Thr and Sub conditions     
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in all conditions
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if only Two-sample T-test takes place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Sub condition 
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_median_HRF_thr)==0 && isempty(psc_wilcoxon_median_HRF_thr)==1) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if only Two-sample T-test takes place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in both Sub and Avg conditions, but only Wilcoxon ranksum test takes place in the Thr condition     
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Thr condition   
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Sub condition, but only Wilcoxon ranksum test takes place in both Thr and Avg conditions    
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Thr condition     
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in both Sub and Avg conditions and Wilcoxon ranksum test takes place in the Thr condition   
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_sub(2,:) psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==0 && isempty(psc_wilcoxon_median_HRF_sub)==1) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in both Thr and Avg conditions 
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Wilcoxon ranksum test takes place in both Thr and Sub conditions     
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         	
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==0 && isempty(psc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in both Thr and Sub conditions   
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_ttest_median_HRF_avg(2,:) psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(psc_ttest_median_HRF_thr)==1 && isempty(psc_wilcoxon_median_HRF_thr)==0) && (isempty(psc_ttest_median_HRF_sub)==1 && isempty(psc_wilcoxon_median_HRF_sub)==0) && (isempty(psc_ttest_median_HRF_avg)==1 && isempty(psc_wilcoxon_median_HRF_avg)==0)															% if only Wilcoxon ranksum test takes place in all conditions 
        [psc_median_HRF_bh_h,psc_median_HRF_crit_p,psc_median_HRF_adj_ci,psc_median_HRF_adj_p]=fdr_bh([psc_wilcoxon_median_HRF_thr(2,:) psc_wilcoxon_median_HRF_sub(2,:) psc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
    end
end


if (isempty(nsc_ttest_median_HRF_thr)==0 || isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 || isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 || isempty(nsc_wilcoxon_median_HRF_avg)==0)                                                      			% if one of the statistical tests takes place in all conditions
    if (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                                				% if both statistical tests take place in all conditions       
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');       % 'yes' in order to see the report   
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in both Thr and Sub conditions, but only Two-sample T-test takes place in the Avg condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in both Thr and Avg conditions, but only Two-sample T-test takes place in the Sub condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in both Sub and Avg conditions
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Avg condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                             																
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
    elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Thr condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Sub condition  
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');        																									   
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in both Sub and Avg conditions, but only Two-sample T-test takes place in the Thr condition     
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in both Thr and Avg conditions
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Avg condition    
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in both Thr and Sub conditions     
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in all conditions
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if only Two-sample T-test takes place in both Thr and Sub conditions, but only Wilcoxon ranksum test takes place in the Avg condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Thr condition and Wilcoxon ranksum test takes place in the Sub condition 
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in both Thr and Avg conditions, but only Wilcoxon ranksum test takes place in the Sub condition
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_median_HRF_thr)==0 && isempty(nsc_wilcoxon_median_HRF_thr)==1) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if only Two-sample T-test takes place in the Thr condition, but only Wilcoxon ranksum test takes place in both Sub and Avg conditions
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in both Sub and Avg conditions, but only Wilcoxon ranksum test takes place in the Thr condition     
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if both statistical tests take place in the Sub condition, but only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in the Thr condition   
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');		
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if both statistical tests take place in the Sub condition, but only Wilcoxon ranksum test takes place in both Thr and Avg conditions    
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in the Thr condition     
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in both Sub and Avg conditions and Wilcoxon ranksum test takes place in the Thr condition   
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_sub(2,:) nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==0 && isempty(nsc_wilcoxon_median_HRF_sub)==1) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if only Two-sample T-test takes place in the Sub condition and Wilcoxon ranksum test takes place in both Thr and Avg conditions 
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																											
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==0)                                              			% if both statistical tests take place in the Avg condition, but only Wilcoxon ranksum test takes place in both Thr and Sub conditions     
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');                                         	
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==0 && isempty(nsc_wilcoxon_median_HRF_avg)==1)															% if only Two-sample T-test takes place in the Avg condition and Wilcoxon ranksum test takes place in both Thr and Sub conditions   
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_ttest_median_HRF_avg(2,:) nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:)]',0.10,'pdep','yes');	
	elseif (isempty(nsc_ttest_median_HRF_thr)==1 && isempty(nsc_wilcoxon_median_HRF_thr)==0) && (isempty(nsc_ttest_median_HRF_sub)==1 && isempty(nsc_wilcoxon_median_HRF_sub)==0) && (isempty(nsc_ttest_median_HRF_avg)==1 && isempty(nsc_wilcoxon_median_HRF_avg)==0)															% if only Wilcoxon ranksum test takes place in all conditions 
        [nsc_median_HRF_bh_h,nsc_median_HRF_crit_p,nsc_median_HRF_adj_ci,nsc_median_HRF_adj_p]=fdr_bh([nsc_wilcoxon_median_HRF_thr(2,:) nsc_wilcoxon_median_HRF_sub(2,:) nsc_wilcoxon_median_HRF_avg(2,:)]',0.10,'pdep','yes');  																										
    end
end



%% STATISTIC STORAGING:

% In this section, we create tables regarding the HRF parameters of
% interest, as well as the Shapiro-Wilk test, Two-sample T-test and 
% Wilcoxon ranksum test results for them and the FDR correction for the
% last two mentioned tests.


% Column with the statistical test's parameter names to identify its 
% corresponding results    
sw_test_values = array2table({'Hypothesis test result', 'p-value','W-value'}','VariableNames', {'Value'});  
ttest_test_values = array2table({'Hypothesis test result', 'p-value','t-value'}','VariableNames', {'Value'});  
wilcoxon_test_values = array2table({'Hypothesis test result', 'p-value','z-value'}','VariableNames', {'Value'});  


for par = 1:n_parameters
    HRF_parameters{par} = strrep(strip(HRF_parameters{par}),' ','_');              % Replaces spaces by '_'
end


% «««««««««««««««««««««««« Shapiro-Wilk tables »»»»»»»»»»»»»»»»»»»»»»»»»»» 


% Forms tables with the Shapiro-Wilk test results for the average and 
% median HRF's parameters data from each set of ROI in each group and  
% condition as well as its headers 
sw_psc_avg_HRF_CNT_Thr_table = array2table(sw_psc_avg_HRF_CNT_Thr,'VariableNames', HRF_parameters);               
sw_psc_avg_HRF_CNT_Sub_table = array2table(sw_psc_avg_HRF_CNT_Sub,'VariableNames', HRF_parameters);
sw_psc_avg_HRF_CNT_Avg_table = array2table(sw_psc_avg_HRF_CNT_Avg,'VariableNames', HRF_parameters);
sw_psc_avg_HRF_T2DM_Thr_table = array2table(sw_psc_avg_HRF_T2DM_Thr,'VariableNames', HRF_parameters);      
sw_psc_avg_HRF_T2DM_Sub_table = array2table(sw_psc_avg_HRF_T2DM_Sub,'VariableNames', HRF_parameters);  
sw_psc_avg_HRF_T2DM_Avg_table = array2table(sw_psc_avg_HRF_T2DM_Avg,'VariableNames', HRF_parameters); 

sw_psc_median_HRF_CNT_Thr_table = array2table(sw_psc_median_HRF_CNT_Thr,'VariableNames', HRF_parameters);               
sw_psc_median_HRF_CNT_Sub_table = array2table(sw_psc_median_HRF_CNT_Sub,'VariableNames', HRF_parameters); 
sw_psc_median_HRF_CNT_Avg_table = array2table(sw_psc_median_HRF_CNT_Avg,'VariableNames', HRF_parameters);       
sw_psc_median_HRF_T2DM_Thr_table = array2table(sw_psc_median_HRF_T2DM_Thr,'VariableNames', HRF_parameters);      
sw_psc_median_HRF_T2DM_Sub_table = array2table(sw_psc_median_HRF_T2DM_Sub,'VariableNames', HRF_parameters);  
sw_psc_median_HRF_T2DM_Avg_table = array2table(sw_psc_median_HRF_T2DM_Avg,'VariableNames', HRF_parameters);
                
sw_nsc_avg_HRF_CNT_Thr_table = array2table(sw_nsc_avg_HRF_CNT_Thr,'VariableNames', HRF_parameters);               
sw_nsc_avg_HRF_CNT_Sub_table = array2table(sw_nsc_avg_HRF_CNT_Sub,'VariableNames', HRF_parameters);       
sw_nsc_avg_HRF_CNT_Avg_table = array2table(sw_nsc_avg_HRF_CNT_Avg,'VariableNames', HRF_parameters);       
sw_nsc_avg_HRF_T2DM_Thr_table = array2table(sw_nsc_avg_HRF_T2DM_Thr,'VariableNames', HRF_parameters);      
sw_nsc_avg_HRF_T2DM_Sub_table = array2table(sw_nsc_avg_HRF_T2DM_Sub,'VariableNames', HRF_parameters); 
sw_nsc_avg_HRF_T2DM_Avg_table = array2table(sw_nsc_avg_HRF_T2DM_Avg,'VariableNames', HRF_parameters); 

sw_nsc_median_HRF_CNT_Thr_table = array2table(sw_nsc_median_HRF_CNT_Thr,'VariableNames', HRF_parameters);               
sw_nsc_median_HRF_CNT_Sub_table = array2table(sw_nsc_median_HRF_CNT_Sub,'VariableNames', HRF_parameters);       
sw_nsc_median_HRF_CNT_Avg_table = array2table(sw_nsc_median_HRF_CNT_Avg,'VariableNames', HRF_parameters);       
sw_nsc_median_HRF_T2DM_Thr_table = array2table(sw_nsc_median_HRF_T2DM_Thr,'VariableNames', HRF_parameters);      
sw_nsc_median_HRF_T2DM_Sub_table = array2table(sw_nsc_median_HRF_T2DM_Sub,'VariableNames', HRF_parameters); 
sw_nsc_median_HRF_T2DM_Avg_table = array2table(sw_nsc_median_HRF_T2DM_Avg,'VariableNames', HRF_parameters); 


% Merges tables
sw_psc_avg_HRF_CNT_Thr_table = [sw_test_values sw_psc_avg_HRF_CNT_Thr_table];
sw_psc_avg_HRF_CNT_Sub_table = [sw_test_values sw_psc_avg_HRF_CNT_Sub_table];
sw_psc_avg_HRF_CNT_Avg_table = [sw_test_values sw_psc_avg_HRF_CNT_Avg_table];
sw_psc_avg_HRF_T2DM_Thr_table = [sw_test_values sw_psc_avg_HRF_T2DM_Thr_table];
sw_psc_avg_HRF_T2DM_Sub_table = [sw_test_values sw_psc_avg_HRF_T2DM_Sub_table];
sw_psc_avg_HRF_T2DM_Avg_table = [sw_test_values sw_psc_avg_HRF_T2DM_Avg_table];

sw_psc_median_HRF_CNT_Thr_table = [sw_test_values sw_psc_median_HRF_CNT_Thr_table];
sw_psc_median_HRF_CNT_Sub_table = [sw_test_values sw_psc_median_HRF_CNT_Sub_table];
sw_psc_median_HRF_CNT_Avg_table = [sw_test_values sw_psc_median_HRF_CNT_Avg_table];
sw_psc_median_HRF_T2DM_Thr_table = [sw_test_values sw_psc_median_HRF_T2DM_Thr_table];
sw_psc_median_HRF_T2DM_Sub_table = [sw_test_values sw_psc_median_HRF_T2DM_Sub_table];
sw_psc_median_HRF_T2DM_Avg_table = [sw_test_values sw_psc_median_HRF_T2DM_Avg_table];

sw_nsc_avg_HRF_CNT_Thr_table = [sw_test_values sw_nsc_avg_HRF_CNT_Thr_table];
sw_nsc_avg_HRF_CNT_Sub_table = [sw_test_values sw_nsc_avg_HRF_CNT_Sub_table];
sw_nsc_avg_HRF_CNT_Avg_table = [sw_test_values sw_nsc_avg_HRF_CNT_Avg_table];
sw_nsc_avg_HRF_T2DM_Thr_table = [sw_test_values sw_nsc_avg_HRF_T2DM_Thr_table];
sw_nsc_avg_HRF_T2DM_Sub_table = [sw_test_values sw_nsc_avg_HRF_T2DM_Sub_table];
sw_nsc_avg_HRF_T2DM_Avg_table = [sw_test_values sw_nsc_avg_HRF_T2DM_Avg_table];

sw_nsc_median_HRF_CNT_Thr_table = [sw_test_values sw_nsc_median_HRF_CNT_Thr_table];
sw_nsc_median_HRF_CNT_Sub_table = [sw_test_values sw_nsc_median_HRF_CNT_Sub_table];
sw_nsc_median_HRF_CNT_Avg_table = [sw_test_values sw_nsc_median_HRF_CNT_Avg_table];
sw_nsc_median_HRF_T2DM_Thr_table = [sw_test_values sw_nsc_median_HRF_T2DM_Thr_table];
sw_nsc_median_HRF_T2DM_Sub_table = [sw_test_values sw_nsc_median_HRF_T2DM_Sub_table];
sw_nsc_median_HRF_T2DM_Avg_table = [sw_test_values sw_nsc_median_HRF_T2DM_Avg_table];



% ********* Two sample T-test and Wilcoxon ranksum test tables ***********
    

% ««««««««««««««««««««««««« Two sample T-test »»»»»»»»»»»»»»»»»»»»»»»»»»» 


% ............................ Average HRF ..............................


if isempty(psc_ttest_avg_HRF_thr)==0                                                                    % when an average HRF parameter from different groups in the Thr condition has normal distribution  
        
    psc_ttest_avg_HRF_thr_headers = cell(1,size(psc_ttest_avg_HRF_thr,2));                              % headers for the table

    for l = 1:length(psc_ttest_indexes_avg_HRF_thr)
        psc_ttest_avg_HRF_thr_headers(l) = HRF_parameters(psc_ttest_indexes_avg_HRF_thr(l));            % gets the headers of the normal Thr parameters
    end

    % Forms tables with the Two sample T-test results for the average HRF's
    % parameters with normal distribution from different groups in the Thr 
    % condition in the positive signal change ROIs as well as its headers 
    psc_ttest_avg_HRF_thr_table = array2table(psc_ttest_avg_HRF_thr,'VariableNames', psc_ttest_avg_HRF_thr_headers);            

    % Merges the previous table with the statistical test's parameter names 
    psc_ttest_avg_HRF_thr_table = [ttest_test_values psc_ttest_avg_HRF_thr_table];

end


if isempty(nsc_ttest_avg_HRF_thr)==0                                                                    % when an average HRF parameter from different groups in the Thr condition has normal distribution  
        
    nsc_ttest_avg_HRF_thr_headers = cell(1,size(nsc_ttest_avg_HRF_thr,2));                              % headers for the table

    for l = 1:length(nsc_ttest_indexes_avg_HRF_thr)
        nsc_ttest_avg_HRF_thr_headers(l) = HRF_parameters(nsc_ttest_indexes_avg_HRF_thr(l));            % gets the headers of the normal Thr parameters
    end

    % Forms tables with the Two sample T-test results for the average HRF's
    % parameters with normal distribution from different groups in the Thr 
    % condition in the negative signal change ROIs as well as its headers 
    nsc_ttest_avg_HRF_thr_table = array2table(nsc_ttest_avg_HRF_thr,'VariableNames', nsc_ttest_avg_HRF_thr_headers);            

    % Merges the previous table with the statistical test's parameter names 
    nsc_ttest_avg_HRF_thr_table = [ttest_test_values nsc_ttest_avg_HRF_thr_table];

end
 
    
if isempty(psc_ttest_avg_HRF_sub)==0                                                                    % when an average HRF parameter from different groups in the Sub condition has normal distribution       

    psc_ttest_avg_HRF_sub_headers = cell(1,size(psc_ttest_avg_HRF_sub,2));                                          

    for l = 1:length(psc_ttest_indexes_avg_HRF_sub)
        psc_ttest_avg_HRF_sub_headers(l) = HRF_parameters(psc_ttest_indexes_avg_HRF_sub(l));            % gets the headers of the normal Sub parameters
    end
    
    % Forms tables with the Two sample T-test results for the average HRF's
    % parameters with normal distribution from different groups in the Sub 
    % condition in the positive signal change ROIs as well as its headers
    psc_ttest_avg_HRF_sub_table = array2table(psc_ttest_avg_HRF_sub,'VariableNames', psc_ttest_avg_HRF_sub_headers);            

    % Merges the previous table with the statistical test's parameter names 
    psc_ttest_avg_HRF_sub_table = [ttest_test_values psc_ttest_avg_HRF_sub_table];

end


if isempty(nsc_ttest_avg_HRF_sub)==0                                                                    % when an average HRF parameter from different groups in the Sub condition has normal distribution       

    nsc_ttest_avg_HRF_sub_headers = cell(1,size(nsc_ttest_avg_HRF_sub,2));                                          

    for l = 1:length(nsc_ttest_indexes_avg_HRF_sub)
        nsc_ttest_avg_HRF_sub_headers(l) = HRF_parameters(nsc_ttest_indexes_avg_HRF_sub(l));            % gets the headers of the normal Sub parameters
    end
    
    % Forms tables with the Two sample T-test results for the average HRF's
    % parameters with normal distribution from different groups in the Sub 
    % condition in the negative signal change ROIs as well as its headers
    nsc_ttest_avg_HRF_sub_table = array2table(nsc_ttest_avg_HRF_sub,'VariableNames', nsc_ttest_avg_HRF_sub_headers);            

    % Merges the previous table with the statistical test's parameter names 
    nsc_ttest_avg_HRF_sub_table = [ttest_test_values nsc_ttest_avg_HRF_sub_table];

end


if isempty(psc_ttest_avg_HRF_avg)==0                                                                    % when an average HRF parameter from different groups in the Avg condition has normal distribution       

    psc_ttest_avg_HRF_avg_headers = cell(1,size(psc_ttest_avg_HRF_avg,2));                                          

    for l = 1:length(psc_ttest_indexes_avg_HRF_avg)
        psc_ttest_avg_HRF_avg_headers(l) = HRF_parameters(psc_ttest_indexes_avg_HRF_avg(l));            % gets the headers of the normal Avg parameters
    end
    
    % Forms tables with the Two sample T-test results for the average HRF's
    % parameters with normal distribution from different groups in the Avg 
    % condition in the positive signal change ROIs  as well as its headers
    psc_ttest_avg_HRF_avg_table = array2table(psc_ttest_avg_HRF_avg,'VariableNames', psc_ttest_avg_HRF_avg_headers);            

    % Merges the previous table with the statistical test's parameter names 
    psc_ttest_avg_HRF_avg_table = [ttest_test_values psc_ttest_avg_HRF_avg_table];

end


if isempty(nsc_ttest_avg_HRF_avg)==0                                                                    % when an average HRF parameter from different groups in the Avg condition has normal distribution       

    nsc_ttest_avg_HRF_avg_headers = cell(1,size(nsc_ttest_avg_HRF_avg,2));                                          

    for l = 1:length(nsc_ttest_indexes_avg_HRF_avg)
        nsc_ttest_avg_HRF_avg_headers(l) = HRF_parameters(nsc_ttest_indexes_avg_HRF_avg(l));            % gets the headers of the normal Avg parameters
    end
    
    % Forms tables with the Two sample T-test results for the average HRF's
    % parameters with normal distribution from different groups in the Avg 
    % condition in the negative signal change ROIs as well as its headers
    nsc_ttest_avg_HRF_avg_table = array2table(nsc_ttest_avg_HRF_avg,'VariableNames', nsc_ttest_avg_HRF_avg_headers);            

    % Merges the previous table with the statistical test's parameter names 
    nsc_ttest_avg_HRF_avg_table = [ttest_test_values nsc_ttest_avg_HRF_avg_table];

end


% ............................ Median HRF ..............................


if isempty(psc_ttest_median_HRF_thr)==0                                                                 % when a median HRF parameter from different groups in the Thr condition has normal distribution  
        
    psc_ttest_median_HRF_thr_headers = cell(1,size(psc_ttest_median_HRF_thr,2));                        % headers for the table

    for l = 1:length(psc_ttest_indexes_median_HRF_thr)
        psc_ttest_median_HRF_thr_headers(l) = HRF_parameters(psc_ttest_indexes_median_HRF_thr(l));      % gets the headers of the normal Thr parameters
    end

    % Forms tables with the Two sample T-test results for the median HRF's
    % parameters with normal distribution from different groups in the Thr 
    % condition in the positive signal change ROIs as well as its headers 
    psc_ttest_median_HRF_thr_table = array2table(psc_ttest_median_HRF_thr,'VariableNames', psc_ttest_median_HRF_thr_headers);            

    % Merges the previous table with the statistical test's parameter names 
    psc_ttest_median_HRF_thr_table = [ttest_test_values psc_ttest_median_HRF_thr_table];
end


if isempty(nsc_ttest_median_HRF_thr)==0                                                                 % when a median HRF parameter from different groups in the Thr condition has normal distribution  
        
    nsc_ttest_median_HRF_thr_headers = cell(1,size(nsc_ttest_median_HRF_thr,2));                        % headers for the table

    for l = 1:length(nsc_ttest_indexes_median_HRF_thr)
        nsc_ttest_median_HRF_thr_headers(l) = HRF_parameters(nsc_ttest_indexes_median_HRF_thr(l));      % gets the headers of the normal Thr parameters
    end

    % Forms tables with the Two sample T-test results for the median HRF's
    % parameters with normal distribution from different groups in the Thr 
    % condition in the negative signal change ROIs as well as its headers 
    nsc_ttest_median_HRF_thr_table = array2table(nsc_ttest_median_HRF_thr,'VariableNames', nsc_ttest_median_HRF_thr_headers);            

    % Merges the previous table with the statistical test's parameter names 
    nsc_ttest_median_HRF_thr_table = [ttest_test_values nsc_ttest_median_HRF_thr_table];

end
 
    
if isempty(psc_ttest_median_HRF_sub)==0                                                                 % when a median HRF parameter from different groups in the Sub condition has normal distribution       

    psc_ttest_median_HRF_sub_headers = cell(1,size(psc_ttest_median_HRF_sub,2));                                          

    for l = 1:length(psc_ttest_indexes_median_HRF_sub)
        psc_ttest_median_HRF_sub_headers(l) = HRF_parameters(psc_ttest_indexes_median_HRF_sub(l));      % gets the headers of the normal Sub parameters
    end
    
    % Forms tables with the Two sample T-test results for the median HRF's
    % parameters with normal distribution from different groups in the Sub 
    % condition in the positive signal change ROIs as well as its headers
    psc_ttest_median_HRF_sub_table = array2table(psc_ttest_median_HRF_sub,'VariableNames', psc_ttest_median_HRF_sub_headers);            

    % Merges the previous table with the statistical test's parameter names 
    psc_ttest_median_HRF_sub_table = [ttest_test_values psc_ttest_median_HRF_sub_table];

end


if isempty(nsc_ttest_median_HRF_sub)==0                                                                 % when a median HRF parameter from different groups in the Sub condition has normal distribution       

    nsc_ttest_median_HRF_sub_headers = cell(1,size(nsc_ttest_median_HRF_sub,2));                                          

    for l = 1:length(nsc_ttest_indexes_median_HRF_sub)
        nsc_ttest_median_HRF_sub_headers(l) = HRF_parameters(nsc_ttest_indexes_median_HRF_sub(l));      % gets the headers of the normal Sub parameters
    end
    
    % Forms tables with the Two sample T-test results for the median HRF's
    % parameters with normal distribution from different groups in the Sub 
    % condition in the negative signal change ROIs as well as its headers
    nsc_ttest_median_HRF_sub_table = array2table(nsc_ttest_median_HRF_sub,'VariableNames', nsc_ttest_median_HRF_sub_headers);            

    % Merges the previous table with the statistical test's parameter names 
    nsc_ttest_median_HRF_sub_table = [ttest_test_values nsc_ttest_median_HRF_sub_table];

end


if isempty(psc_ttest_median_HRF_avg)==0                                                                    % when a median HRF parameter from different groups in the Avg condition has normal distribution       

    psc_ttest_median_HRF_avg_headers = cell(1,size(psc_ttest_median_HRF_avg,2));                                          

    for l = 1:length(psc_ttest_indexes_median_HRF_avg)
        psc_ttest_median_HRF_avg_headers(l) = HRF_parameters(psc_ttest_indexes_median_HRF_avg(l));         % gets the headers of the normal Avg parameters
    end
    
    % Forms tables with the Two sample T-test results for the median HRF's
    % parameters with normal distribution from different groups in the Avg 
    % condition in the positive signal change ROIs as well as its headers
    psc_ttest_median_HRF_avg_table = array2table(psc_ttest_median_HRF_avg,'VariableNames', psc_ttest_median_HRF_avg_headers);            

    % Merges the previous table with the statistical test's parameter names 
    psc_ttest_median_HRF_avg_table = [ttest_test_values psc_ttest_median_HRF_avg_table];

end


if isempty(nsc_ttest_median_HRF_avg)==0                                                                     % when a median HRF parameter from different groups in the Avg condition has normal distribution       

    nsc_ttest_median_HRF_avg_headers = cell(1,size(nsc_ttest_median_HRF_avg,2));                                          

    for l = 1:length(nsc_ttest_indexes_median_HRF_avg)
        nsc_ttest_median_HRF_avg_headers(l) = HRF_parameters(nsc_ttest_indexes_median_HRF_avg(l));          % gets the headers of the normal Avg parameters
    end
    
    % Forms tables with the Two sample T-test results for the median HRF's
    % parameters with normal distribution from different groups in the Avg 
    % condition in the negative signal change ROIs as well as its headers
    nsc_ttest_median_HRF_avg_table = array2table(nsc_ttest_median_HRF_avg,'VariableNames', nsc_ttest_median_HRF_avg_headers);            

    % Merges the previous table with the statistical test's parameter names 
    nsc_ttest_median_HRF_avg_table = [ttest_test_values nsc_ttest_median_HRF_avg_table];

end


% ««««««««««««««««««««««« Wilcoxon ranksum test »»»»»»»»»»»»»»»»»»»»»»»»»» 


% ............................ Average HRF ..............................


if isempty(psc_wilcoxon_avg_HRF_thr)==0                                                                     % when an average HRF parameter from different groups in the Thr condition has non-normal distribution 

    psc_wilcoxon_avg_HRF_thr_headers = cell(1,size(psc_wilcoxon_avg_HRF_thr,2));                                       

    for l = 1:length(psc_wilcoxon_indexes_avg_HRF_thr)
        psc_wilcoxon_avg_HRF_thr_headers(l) = HRF_parameters(psc_wilcoxon_indexes_avg_HRF_thr(l));          % gets the headers of non-normal Thr HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the average
    % HRF's parameters with non-normal distribution from different groups 
    % in the Thr condition in the positive signal change ROIs as well as 
    % its headers
    psc_wilcoxon_avg_HRF_thr_table = array2table(psc_wilcoxon_avg_HRF_thr,'VariableNames', psc_wilcoxon_avg_HRF_thr_headers);       

    % Merges the previous table with the statistical test's parameter names 
    psc_wilcoxon_avg_HRF_thr_table = [wilcoxon_test_values psc_wilcoxon_avg_HRF_thr_table];

end


if isempty(nsc_wilcoxon_avg_HRF_thr)==0                                                                     % when an average HRF parameter from different groups in the Thr condition has non-normal distribution 

    nsc_wilcoxon_avg_HRF_thr_headers = cell(1,size(nsc_wilcoxon_avg_HRF_thr,2));                                       

    for l = 1:length(nsc_wilcoxon_indexes_avg_HRF_thr)
        nsc_wilcoxon_avg_HRF_thr_headers(l) = HRF_parameters(nsc_wilcoxon_indexes_avg_HRF_thr(l));          % gets the headers of non-normal Thr HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the average
    % HRF's parameters with non-normal distribution from different groups 
    % in the Thr condition in the negative signal change ROIs as well as 
    % its headers
    nsc_wilcoxon_avg_HRF_thr_table = array2table(nsc_wilcoxon_avg_HRF_thr,'VariableNames', nsc_wilcoxon_avg_HRF_thr_headers);       

    % Merges the previous table with the statistical test's parameter names 
    nsc_wilcoxon_avg_HRF_thr_table = [wilcoxon_test_values nsc_wilcoxon_avg_HRF_thr_table];

end


if isempty(psc_wilcoxon_avg_HRF_sub)==0                                                                     % when an average HRF parameter from different groups in the Sub condition has non-normal distribution  

    psc_wilcoxon_avg_HRF_sub_headers = cell(1,size(psc_wilcoxon_avg_HRF_sub,2));                                        

    for l = 1:length(psc_wilcoxon_indexes_avg_HRF_sub)
        psc_wilcoxon_avg_HRF_sub_headers(l) = HRF_parameters(psc_wilcoxon_indexes_avg_HRF_sub(l));          % gets the headers of non-normal Sub HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the average
    % HRF's parameters with non-normal distribution from different groups 
    % in the Sub condition in the positive signal change ROIs as well as 
    % its headers
    psc_wilcoxon_avg_HRF_sub_table = array2table(psc_wilcoxon_avg_HRF_sub,'VariableNames', psc_wilcoxon_avg_HRF_sub_headers);       

    % Merges the previous table with the statistical test's parameter names 
    psc_wilcoxon_avg_HRF_sub_table = [wilcoxon_test_values psc_wilcoxon_avg_HRF_sub_table];

end


if isempty(nsc_wilcoxon_avg_HRF_sub)==0                                                                     % when an average HRF parameter from different groups in the Sub condition has non-normal distribution  

    nsc_wilcoxon_avg_HRF_sub_headers = cell(1,size(nsc_wilcoxon_avg_HRF_sub,2));                                        

    for l = 1:length(nsc_wilcoxon_indexes_avg_HRF_sub)
        nsc_wilcoxon_avg_HRF_sub_headers(l) = HRF_parameters(nsc_wilcoxon_indexes_avg_HRF_sub(l));          % gets the headers of non-normal Sub HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the average
    % HRF's parameters with non-normal distribution from different groups 
    % in the Sub condition in the negative signal change ROIs as well as 
    % its headers
    nsc_wilcoxon_avg_HRF_sub_table = array2table(nsc_wilcoxon_avg_HRF_sub,'VariableNames', nsc_wilcoxon_avg_HRF_sub_headers);       

    % Merges the previous table with the statistical test's parameter names 
    nsc_wilcoxon_avg_HRF_sub_table = [wilcoxon_test_values nsc_wilcoxon_avg_HRF_sub_table];

end


if isempty(psc_wilcoxon_avg_HRF_avg)==0                                                                     % when an average HRF parameter from different groups in the Avg condition has non-normal distribution  

    psc_wilcoxon_avg_HRF_avg_headers = cell(1,size(psc_wilcoxon_avg_HRF_avg,2));                                        

    for l = 1:length(psc_wilcoxon_indexes_avg_HRF_avg)
        psc_wilcoxon_avg_HRF_avg_headers(l) = HRF_parameters(psc_wilcoxon_indexes_avg_HRF_avg(l));          % gets the headers of non-normal Avg HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the average
    % HRF's parameters with non-normal distribution from different groups 
    % in the Avg condition in the positive signal change ROIs as well as 
    % its headers
    psc_wilcoxon_avg_HRF_avg_table = array2table(psc_wilcoxon_avg_HRF_avg,'VariableNames', psc_wilcoxon_avg_HRF_avg_headers);       

    % Merges the previous table with the statistical test's parameter names 
    psc_wilcoxon_avg_HRF_avg_table = [wilcoxon_test_values psc_wilcoxon_avg_HRF_avg_table];

end


if isempty(nsc_wilcoxon_avg_HRF_avg)==0                                                                     % when an average HRF parameter from different groups in the Avg condition has non-normal distribution  

    nsc_wilcoxon_avg_HRF_avg_headers = cell(1,size(nsc_wilcoxon_avg_HRF_avg,2));                                        

    for l = 1:length(nsc_wilcoxon_indexes_avg_HRF_avg)
        nsc_wilcoxon_avg_HRF_avg_headers(l) = HRF_parameters(nsc_wilcoxon_indexes_avg_HRF_avg(l));          % gets the headers of non-normal Avg HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the average 
    % HRF's parameters with non-normal distribution from different groups 
    % in the Avg condition in the negative signal change ROIs as well as 
    % its headers
    nsc_wilcoxon_avg_HRF_avg_table = array2table(nsc_wilcoxon_avg_HRF_avg,'VariableNames', nsc_wilcoxon_avg_HRF_avg_headers);       

    % Merges the previous table with the statistical test's parameter names 
    nsc_wilcoxon_avg_HRF_avg_table = [wilcoxon_test_values nsc_wilcoxon_avg_HRF_avg_table];

end



% ............................ Median HRF ..............................


if isempty(psc_wilcoxon_median_HRF_thr)==0                                                                  % when a median HRF parameter from different groups in the Thr condition has non-normal distribution 

    psc_wilcoxon_median_HRF_thr_headers = cell(1,size(psc_wilcoxon_median_HRF_thr,2));                                       

    for l = 1:length(psc_wilcoxon_indexes_median_HRF_thr)
        psc_wilcoxon_median_HRF_thr_headers(l) = HRF_parameters(psc_wilcoxon_indexes_median_HRF_thr(l));    % gets the headers of non-normal Thr HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the median
    % HRF's parameters with non-normal distribution from different groups 
    % in the Thr condition in the positive signal change ROIs as well as 
    % its headers
    psc_wilcoxon_median_HRF_thr_table = array2table(psc_wilcoxon_median_HRF_thr,'VariableNames', psc_wilcoxon_median_HRF_thr_headers);       

    % Merges the previous table with the statistical test's parameter names 
    psc_wilcoxon_median_HRF_thr_table = [wilcoxon_test_values psc_wilcoxon_median_HRF_thr_table];

end


if isempty(nsc_wilcoxon_median_HRF_thr)==0                                                                  % when a median HRF parameter from different groups in the Thr condition has non-normal distribution 

    nsc_wilcoxon_median_HRF_thr_headers = cell(1,size(nsc_wilcoxon_median_HRF_thr,2));                                       

    for l = 1:length(nsc_wilcoxon_indexes_median_HRF_thr)
        nsc_wilcoxon_median_HRF_thr_headers(l) = HRF_parameters(nsc_wilcoxon_indexes_median_HRF_thr(l));    % gets the headers of non-normal Thr HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the median
    % HRF's parameters with non-normal distribution from different groups 
    % in the Thr condition in the negative signal change ROIs as well as 
    % its headers
    nsc_wilcoxon_median_HRF_thr_table = array2table(nsc_wilcoxon_median_HRF_thr,'VariableNames', nsc_wilcoxon_median_HRF_thr_headers);       

    % Merges the previous table with the statistical test's parameter names 
    nsc_wilcoxon_median_HRF_thr_table = [wilcoxon_test_values nsc_wilcoxon_median_HRF_thr_table];

end


if isempty(psc_wilcoxon_median_HRF_sub)==0                                                                  % when a median HRF parameter from different groups in the Sub condition has non-normal distribution  

    psc_wilcoxon_median_HRF_sub_headers = cell(1,size(psc_wilcoxon_median_HRF_sub,2));                                        

    for l = 1:length(psc_wilcoxon_indexes_median_HRF_sub)
        psc_wilcoxon_median_HRF_sub_headers(l) = HRF_parameters(psc_wilcoxon_indexes_median_HRF_sub(l));    % gets the headers of non-normal Sub HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the median 
    % HRF's parameters with non-normal distribution from different groups 
    % in the Sub condition in the positive signal change ROIs as well as 
    % its headers
    psc_wilcoxon_median_HRF_sub_table = array2table(psc_wilcoxon_median_HRF_sub,'VariableNames', psc_wilcoxon_median_HRF_sub_headers);       

    % Merges the previous table with the statistical test's parameter names 
    psc_wilcoxon_median_HRF_sub_table = [wilcoxon_test_values psc_wilcoxon_median_HRF_sub_table];

end


if isempty(nsc_wilcoxon_median_HRF_sub)==0                                                                  % when a median HRF parameter from different groups in the Sub condition has non-normal distribution  

    nsc_wilcoxon_median_HRF_sub_headers = cell(1,size(nsc_wilcoxon_median_HRF_sub,2));                                        

    for l = 1:length(nsc_wilcoxon_indexes_median_HRF_sub)
        nsc_wilcoxon_median_HRF_sub_headers(l) = HRF_parameters(nsc_wilcoxon_indexes_median_HRF_sub(l));    % gets the headers of non-normal Sub HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the median
    % HRF's parameters with non-normal distribution from different groups 
    % in the Sub condition in the negative signal change ROIs as well as 
    % its headers
    nsc_wilcoxon_median_HRF_sub_table = array2table(nsc_wilcoxon_median_HRF_sub,'VariableNames', nsc_wilcoxon_median_HRF_sub_headers);       

    % Merges the previous table with the statistical test's parameter names 
    nsc_wilcoxon_median_HRF_sub_table = [wilcoxon_test_values nsc_wilcoxon_median_HRF_sub_table];

end


if isempty(psc_wilcoxon_median_HRF_avg)==0                                                                  % when a median HRF parameter from different groups in the Avg condition has non-normal distribution  

    psc_wilcoxon_median_HRF_avg_headers = cell(1,size(psc_wilcoxon_median_HRF_avg,2));                                        

    for l = 1:length(psc_wilcoxon_indexes_avg_HRF_avg)
        psc_wilcoxon_median_HRF_avg_headers(l) = HRF_parameters(psc_wilcoxon_indexes_median_HRF_avg(l));    % gets the headers of non-normal Avg HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the median 
    % HRF's parameters with non-normal distribution from different groups 
    % in the Avg condition in the positive signal change ROIs as well as 
    % its headers
    psc_wilcoxon_median_HRF_avg_table = array2table(psc_wilcoxon_median_HRF_avg,'VariableNames', psc_wilcoxon_median_HRF_avg_headers);       

    % Merges the previous table with the statistical test's parameter names 
    psc_wilcoxon_median_HRF_avg_table = [wilcoxon_test_values psc_wilcoxon_median_HRF_avg_table];

end


if isempty(nsc_wilcoxon_median_HRF_avg)==0                                                                  % when a median HRF parameter from different groups in the Avg condition has non-normal distribution  

    nsc_wilcoxon_median_HRF_avg_headers = cell(1,size(nsc_wilcoxon_median_HRF_avg,2));                                        

    for l = 1:length(nsc_wilcoxon_indexes_median_HRF_avg)
        nsc_wilcoxon_median_HRF_avg_headers(l) = HRF_parameters(nsc_wilcoxon_indexes_median_HRF_avg(l));    % gets the headers of non-normal Avg HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the median
    % HRF's parameters with non-normal distribution from different groups 
    % in the Avg condition in the negative signal change ROIs as well as 
    % its headers
    nsc_wilcoxon_median_HRF_avg_table = array2table(nsc_wilcoxon_median_HRF_avg,'VariableNames', nsc_wilcoxon_median_HRF_avg_headers);       

    % Merges the previous table with the statistical test's parameter names 
    nsc_wilcoxon_median_HRF_avg_table = [wilcoxon_test_values nsc_wilcoxon_median_HRF_avg_table];

end



% ************************** HRF Parameters ******************************
    

% «««««««««««««««««««« Average and standard deviation »»»»»»»»»»»»»»»»»»»»


% ---------------------------- Average HRF ------------------------------


% ................................ CNT ..................................


if isempty(psc_avg_parameters_avg_HRF_CNT_Thr)==0                                                                   % when an average HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each set of ROI

    psc_avg_headers_avg_HRF_CNT_Thr = cell(1,size(psc_avg_parameters_avg_HRF_CNT_Thr,2));                           % headers for the table

    for l = 1:length(psc_avg_std_indexes_avg_HRF_CNT_Thr)
        psc_avg_headers_avg_HRF_CNT_Thr(l) = HRF_parameters(psc_avg_std_indexes_avg_HRF_CNT_Thr(l));                % gets the headers of the CNT Thr average HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Thr
    % average HRF parameters with normal distribution as well as its 
    % headers
    psc_avg_parameters_avg_HRF_CNT_Thr_table = array2table(psc_avg_parameters_avg_HRF_CNT_Thr,'VariableNames', psc_avg_headers_avg_HRF_CNT_Thr);       
    psc_std_parameters_avg_HRF_CNT_Thr_table = array2table(psc_std_parameters_avg_HRF_CNT_Thr,'VariableNames', psc_avg_headers_avg_HRF_CNT_Thr);              

end


if isempty(nsc_avg_parameters_avg_HRF_CNT_Thr)==0                                                                   

    nsc_avg_headers_avg_HRF_CNT_Thr = cell(1,size(nsc_avg_parameters_avg_HRF_CNT_Thr,2));                           

    for l = 1:length(nsc_avg_std_indexes_avg_HRF_CNT_Thr)
        nsc_avg_headers_avg_HRF_CNT_Thr(l) = HRF_parameters(nsc_avg_std_indexes_avg_HRF_CNT_Thr(l));                % gets the headers of the CNT Thr average HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Thr
    % average HRF parameters with normal distribution as well as its 
    % headers
    nsc_avg_parameters_avg_HRF_CNT_Thr_table = array2table(nsc_avg_parameters_avg_HRF_CNT_Thr,'VariableNames', nsc_avg_headers_avg_HRF_CNT_Thr);       
    nsc_std_parameters_avg_HRF_CNT_Thr_table = array2table(nsc_std_parameters_avg_HRF_CNT_Thr,'VariableNames', nsc_avg_headers_avg_HRF_CNT_Thr);              

end


if isempty(psc_avg_parameters_avg_HRF_CNT_Sub)==0                                                                         

    psc_avg_headers_avg_HRF_CNT_Sub = cell(1,size(psc_avg_parameters_avg_HRF_CNT_Sub,2));                                         

    for l = 1:length(psc_avg_std_indexes_avg_HRF_CNT_Sub)
        psc_avg_headers_avg_HRF_CNT_Sub(l) = HRF_parameters(psc_avg_std_indexes_avg_HRF_CNT_Sub(l));                % gets the headers of the CNT Sub average HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Sub
    % average HRF parameters with normal distribution as well as its 
    % headers
    psc_avg_parameters_avg_HRF_CNT_Sub_table = array2table(psc_avg_parameters_avg_HRF_CNT_Sub,'VariableNames', psc_avg_headers_avg_HRF_CNT_Sub);       
    psc_std_parameters_avg_HRF_CNT_Sub_table = array2table(psc_std_parameters_avg_HRF_CNT_Sub,'VariableNames', psc_avg_headers_avg_HRF_CNT_Sub);         

end
 

if isempty(nsc_avg_parameters_avg_HRF_CNT_Sub)==0                                                                         

    nsc_avg_headers_avg_HRF_CNT_Sub = cell(1,size(nsc_avg_parameters_avg_HRF_CNT_Sub,2));                                         

    for l = 1:length(nsc_avg_std_indexes_avg_HRF_CNT_Sub)
        nsc_avg_headers_avg_HRF_CNT_Sub(l) = HRF_parameters(nsc_avg_std_indexes_avg_HRF_CNT_Sub(l));                % gets the headers of the CNT Sub average HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Sub
    % average HRF parameters with normal distribution as well as its 
    % headers
    nsc_avg_parameters_avg_HRF_CNT_Sub_table = array2table(nsc_avg_parameters_avg_HRF_CNT_Sub,'VariableNames', nsc_avg_headers_avg_HRF_CNT_Sub);       
    nsc_std_parameters_avg_HRF_CNT_Sub_table = array2table(nsc_std_parameters_avg_HRF_CNT_Sub,'VariableNames', nsc_avg_headers_avg_HRF_CNT_Sub);         

end


if isempty(psc_avg_parameters_avg_HRF_CNT_Avg)==0                                                                         

    psc_avg_headers_avg_HRF_CNT_Avg = cell(1,size(psc_avg_parameters_avg_HRF_CNT_Avg,2));                                         

    for l = 1:length(psc_avg_std_indexes_avg_HRF_CNT_Avg)
        psc_avg_headers_avg_HRF_CNT_Avg(l) = HRF_parameters(psc_avg_std_indexes_avg_HRF_CNT_Avg(l));                % gets the headers of the CNT Avg average HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Avg
    % average HRF parameters with normal distribution as well as its 
    % headers
    psc_avg_parameters_avg_HRF_CNT_Avg_table = array2table(psc_avg_parameters_avg_HRF_CNT_Avg,'VariableNames', psc_avg_headers_avg_HRF_CNT_Avg);       
    psc_std_parameters_avg_HRF_CNT_Avg_table = array2table(psc_std_parameters_avg_HRF_CNT_Avg,'VariableNames', psc_avg_headers_avg_HRF_CNT_Avg);         

end
 

if isempty(nsc_avg_parameters_avg_HRF_CNT_Avg)==0                                                                         

    nsc_avg_headers_avg_HRF_CNT_Avg = cell(1,size(nsc_avg_parameters_avg_HRF_CNT_Avg,2));                                         

    for l = 1:length(nsc_avg_std_indexes_avg_HRF_CNT_Avg)
        nsc_avg_headers_avg_HRF_CNT_Avg(l) = HRF_parameters(nsc_avg_std_indexes_avg_HRF_CNT_Avg(l));                % gets the headers of the CNT Avg average HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Avg
    % average HRF parameters with normal distribution as well as its 
    % headers
    nsc_avg_parameters_avg_HRF_CNT_Avg_table = array2table(nsc_avg_parameters_avg_HRF_CNT_Avg,'VariableNames', nsc_avg_headers_avg_HRF_CNT_Avg);       
    nsc_std_parameters_avg_HRF_CNT_Avg_table = array2table(nsc_std_parameters_avg_HRF_CNT_Avg,'VariableNames', nsc_avg_headers_avg_HRF_CNT_Avg);         

end



% ............................... T2DM ..................................


if isempty(psc_avg_parameters_avg_HRF_T2DM_Thr)==0                                                                         

    psc_avg_headers_avg_HRF_T2DM_Thr = cell(1,size(psc_avg_parameters_avg_HRF_T2DM_Thr,2));                                       

    for l = 1:length(psc_avg_std_indexes_avg_HRF_T2DM_Thr)
        psc_avg_headers_avg_HRF_T2DM_Thr(l) = HRF_parameters(psc_avg_std_indexes_avg_HRF_T2DM_Thr(l));              % gets the headers of the T2DM Thr average HRF parameters with normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the average and standard deviation for the T2DM Thr
    % average HRF parameters with normal distribution as well as its 
    % headers
    psc_avg_parameters_avg_HRF_T2DM_Thr_table = array2table(psc_avg_parameters_avg_HRF_T2DM_Thr,'VariableNames', psc_avg_headers_avg_HRF_T2DM_Thr);    
    psc_std_parameters_avg_HRF_T2DM_Thr_table = array2table(psc_std_parameters_avg_HRF_T2DM_Thr,'VariableNames', psc_avg_headers_avg_HRF_T2DM_Thr);      

end


if isempty(nsc_avg_parameters_avg_HRF_T2DM_Thr)==0                                                                         

    nsc_avg_headers_avg_HRF_T2DM_Thr = cell(1,size(nsc_avg_parameters_avg_HRF_T2DM_Thr,2));                                       

    for l = 1:length(nsc_avg_std_indexes_avg_HRF_T2DM_Thr)
        nsc_avg_headers_avg_HRF_T2DM_Thr(l) = HRF_parameters(nsc_avg_std_indexes_avg_HRF_T2DM_Thr(l));              % gets the headers of the T2DM Thr average HRF parameters with normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the average and standard deviation for the T2DM Thr
    % average HRF parameters with normal distribution as well as its 
    % headers
    nsc_avg_parameters_avg_HRF_T2DM_Thr_table = array2table(nsc_avg_parameters_avg_HRF_T2DM_Thr,'VariableNames', nsc_avg_headers_avg_HRF_T2DM_Thr);    
    nsc_std_parameters_avg_HRF_T2DM_Thr_table = array2table(nsc_std_parameters_avg_HRF_T2DM_Thr,'VariableNames', nsc_avg_headers_avg_HRF_T2DM_Thr);      

end


if isempty(psc_avg_parameters_avg_HRF_T2DM_Sub)==0                                                                        

    psc_avg_headers_avg_HRF_T2DM_Sub = cell(1,size(psc_avg_parameters_avg_HRF_T2DM_Sub,2));                                       

    for l = 1:length(psc_avg_std_indexes_avg_HRF_T2DM_Sub)
        psc_avg_headers_avg_HRF_T2DM_Sub(l) = HRF_parameters(psc_avg_std_indexes_avg_HRF_T2DM_Sub(l));               % gets the headers of the T2DM Sub average HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Sub
    % average HRF parameters with normal distribution as well as its 
    % headers
    psc_avg_parameters_avg_HRF_T2DM_Sub_table = array2table(psc_avg_parameters_avg_HRF_T2DM_Sub,'VariableNames', psc_avg_headers_avg_HRF_T2DM_Sub);    
    psc_std_parameters_avg_HRF_T2DM_Sub_table = array2table(psc_std_parameters_avg_HRF_T2DM_Sub,'VariableNames', psc_avg_headers_avg_HRF_T2DM_Sub);      

end


if isempty(nsc_avg_parameters_avg_HRF_T2DM_Sub)==0                                                                        

    nsc_avg_headers_avg_HRF_T2DM_Sub = cell(1,size(nsc_avg_parameters_avg_HRF_T2DM_Sub,2));                                       

    for l = 1:length(nsc_avg_std_indexes_avg_HRF_T2DM_Sub)
        nsc_avg_headers_avg_HRF_T2DM_Sub(l) = HRF_parameters(nsc_avg_std_indexes_avg_HRF_T2DM_Sub(l));              % gets the headers of normal T2DM Sub average HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Sub
    % average HRF parameters with normal distribution as well as its 
    % headers
    nsc_avg_parameters_avg_HRF_T2DM_Sub_table = array2table(nsc_avg_parameters_avg_HRF_T2DM_Sub,'VariableNames', nsc_avg_headers_avg_HRF_T2DM_Sub);    
    nsc_std_parameters_avg_HRF_T2DM_Sub_table = array2table(nsc_std_parameters_avg_HRF_T2DM_Sub,'VariableNames', nsc_avg_headers_avg_HRF_T2DM_Sub);      

end


if isempty(psc_avg_parameters_avg_HRF_T2DM_Avg)==0                                                                        

    psc_avg_headers_avg_HRF_T2DM_Avg = cell(1,size(psc_avg_parameters_avg_HRF_T2DM_Avg,2));                                       

    for l = 1:length(psc_avg_std_indexes_avg_HRF_T2DM_Avg)
        psc_avg_headers_avg_HRF_T2DM_Avg(l) = HRF_parameters(psc_avg_std_indexes_avg_HRF_T2DM_Avg(l));               % gets the headers of the T2DM Avg average HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Avg
    % average HRF parameters with normal distribution as well as its 
    % headers
    psc_avg_parameters_avg_HRF_T2DM_Avg_table = array2table(psc_avg_parameters_avg_HRF_T2DM_Avg,'VariableNames', psc_avg_headers_avg_HRF_T2DM_Avg);    
    psc_std_parameters_avg_HRF_T2DM_Avg_table = array2table(psc_std_parameters_avg_HRF_T2DM_Avg,'VariableNames', psc_avg_headers_avg_HRF_T2DM_Avg);      

end


if isempty(nsc_avg_parameters_avg_HRF_T2DM_Avg)==0                                                                        

    nsc_avg_headers_avg_HRF_T2DM_Avg = cell(1,size(nsc_avg_parameters_avg_HRF_T2DM_Avg,2));                                       

    for l = 1:length(nsc_avg_std_indexes_avg_HRF_T2DM_Avg)
        nsc_avg_headers_avg_HRF_T2DM_Avg(l) = HRF_parameters(nsc_avg_std_indexes_avg_HRF_T2DM_Avg(l));              % gets the headers of normal T2DM Avg average HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Avg
    % average HRF parameters with normal distribution as well as its 
    % headers
    nsc_avg_parameters_avg_HRF_T2DM_Avg_table = array2table(nsc_avg_parameters_avg_HRF_T2DM_Avg,'VariableNames', nsc_avg_headers_avg_HRF_T2DM_Avg);    
    nsc_std_parameters_avg_HRF_T2DM_Avg_table = array2table(nsc_std_parameters_avg_HRF_T2DM_Avg,'VariableNames', nsc_avg_headers_avg_HRF_T2DM_Avg);      

end



% ---------------------------- Median HRF -------------------------------


% ................................ CNT ..................................


if isempty(psc_avg_parameters_median_HRF_CNT_Thr)==0                                                                % when a median HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each set of ROI

    psc_avg_headers_median_HRF_CNT_Thr = cell(1,size(psc_avg_parameters_median_HRF_CNT_Thr,2));                     % headers for the table

    for l = 1:length(psc_avg_std_indexes_median_HRF_CNT_Thr)
        psc_avg_headers_median_HRF_CNT_Thr(l) = HRF_parameters(psc_avg_std_indexes_median_HRF_CNT_Thr(l));          % gets the headers of the CNT Thr median HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Thr
    % median HRF parameters with normal distribution as well as its headers
    psc_avg_parameters_median_HRF_CNT_Thr_table = array2table(psc_avg_parameters_median_HRF_CNT_Thr,'VariableNames', psc_avg_headers_median_HRF_CNT_Thr);       
    psc_std_parameters_median_HRF_CNT_Thr_table = array2table(psc_std_parameters_median_HRF_CNT_Thr,'VariableNames', psc_avg_headers_median_HRF_CNT_Thr);              

end


if isempty(nsc_avg_parameters_median_HRF_CNT_Thr)==0                                                                   

    nsc_avg_headers_median_HRF_CNT_Thr = cell(1,size(nsc_avg_parameters_median_HRF_CNT_Thr,2));                           

    for l = 1:length(nsc_avg_std_indexes_median_HRF_CNT_Thr)
        nsc_avg_headers_median_HRF_CNT_Thr(l) = HRF_parameters(nsc_avg_std_indexes_median_HRF_CNT_Thr(l));          % gets the headers of the CNT Thr median HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Thr
    % median HRF parameters with normal distribution as well as its headers
    nsc_avg_parameters_median_HRF_CNT_Thr_table = array2table(nsc_avg_parameters_median_HRF_CNT_Thr,'VariableNames', nsc_avg_headers_median_HRF_CNT_Thr);       
    nsc_std_parameters_median_HRF_CNT_Thr_table = array2table(nsc_std_parameters_median_HRF_CNT_Thr,'VariableNames', nsc_avg_headers_median_HRF_CNT_Thr);              

end


if isempty(psc_avg_parameters_median_HRF_CNT_Sub)==0                                                                         

    psc_avg_headers_median_HRF_CNT_Sub = cell(1,size(psc_avg_parameters_median_HRF_CNT_Sub,2));                                         

    for l = 1:length(psc_avg_std_indexes_median_HRF_CNT_Sub)
        psc_avg_headers_median_HRF_CNT_Sub(l) = HRF_parameters(psc_avg_std_indexes_median_HRF_CNT_Sub(l));          % gets the headers of the CNT Sub median HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Sub
    % median HRF parameters with normal distribution as well as its headers
    psc_avg_parameters_median_HRF_CNT_Sub_table = array2table(psc_avg_parameters_median_HRF_CNT_Sub,'VariableNames', psc_avg_headers_median_HRF_CNT_Sub);       
    psc_std_parameters_median_HRF_CNT_Sub_table = array2table(psc_std_parameters_median_HRF_CNT_Sub,'VariableNames', psc_avg_headers_median_HRF_CNT_Sub);         

end
 

if isempty(nsc_avg_parameters_median_HRF_CNT_Sub)==0                                                                         

    nsc_avg_headers_median_HRF_CNT_Sub = cell(1,size(nsc_avg_parameters_median_HRF_CNT_Sub,2));                                         

    for l = 1:length(nsc_avg_std_indexes_median_HRF_CNT_Sub)
        nsc_avg_headers_median_HRF_CNT_Sub(l) = HRF_parameters(nsc_avg_std_indexes_median_HRF_CNT_Sub(l));          % gets the headers of the CNT Sub median HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Sub
    % median HRF parameters with normal distribution as well as its headers 
    nsc_avg_parameters_median_HRF_CNT_Sub_table = array2table(nsc_avg_parameters_median_HRF_CNT_Sub,'VariableNames', nsc_avg_headers_median_HRF_CNT_Sub);       
    nsc_std_parameters_median_HRF_CNT_Sub_table = array2table(nsc_std_parameters_median_HRF_CNT_Sub,'VariableNames', nsc_avg_headers_median_HRF_CNT_Sub);         

end


if isempty(psc_avg_parameters_median_HRF_CNT_Avg)==0                                                                         

    psc_avg_headers_median_HRF_CNT_Avg = cell(1,size(psc_avg_parameters_median_HRF_CNT_Avg,2));                                         

    for l = 1:length(psc_avg_std_indexes_median_HRF_CNT_Avg)
        psc_avg_headers_median_HRF_CNT_Avg(l) = HRF_parameters(psc_avg_std_indexes_median_HRF_CNT_Avg(l));          % gets the headers of the CNT Avg median HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Avg
    % median HRF parameters with normal distribution as well as its headers 
    psc_avg_parameters_median_HRF_CNT_Avg_table = array2table(psc_avg_parameters_median_HRF_CNT_Avg,'VariableNames', psc_avg_headers_median_HRF_CNT_Avg);       
    psc_std_parameters_median_HRF_CNT_Avg_table = array2table(psc_std_parameters_median_HRF_CNT_Avg,'VariableNames', psc_avg_headers_median_HRF_CNT_Avg);         

end
 

if isempty(nsc_avg_parameters_median_HRF_CNT_Avg)==0                                                                         

    nsc_avg_headers_median_HRF_CNT_Avg = cell(1,size(nsc_avg_parameters_median_HRF_CNT_Avg,2));                                         

    for l = 1:length(nsc_avg_std_indexes_median_HRF_CNT_Avg)
        nsc_avg_headers_median_HRF_CNT_Avg(l) = HRF_parameters(nsc_avg_std_indexes_median_HRF_CNT_Avg(l));          % gets the headers of the CNT Avg median HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the CNT Avg
    % median HRF parameters with normal distribution as well as its headers
    nsc_avg_parameters_median_HRF_CNT_Avg_table = array2table(nsc_avg_parameters_median_HRF_CNT_Avg,'VariableNames', nsc_avg_headers_median_HRF_CNT_Avg);       
    nsc_std_parameters_median_HRF_CNT_Avg_table = array2table(nsc_std_parameters_median_HRF_CNT_Avg,'VariableNames', nsc_avg_headers_median_HRF_CNT_Avg);         

end



% ............................... T2DM ..................................


if isempty(psc_avg_parameters_median_HRF_T2DM_Thr)==0                                                                         

    psc_avg_headers_median_HRF_T2DM_Thr = cell(1,size(psc_avg_parameters_median_HRF_T2DM_Thr,2));                                       

    for l = 1:length(psc_avg_std_indexes_median_HRF_T2DM_Thr)
        psc_avg_headers_median_HRF_T2DM_Thr(l) = HRF_parameters(psc_avg_std_indexes_median_HRF_T2DM_Thr(l));        % gets the headers of the T2DM Thr median HRF parameters with normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the average and standard deviation for the T2DM Thr
    % median HRF parameters with normal distribution as well as its headers
    psc_avg_parameters_median_HRF_T2DM_Thr_table = array2table(psc_avg_parameters_median_HRF_T2DM_Thr,'VariableNames', psc_avg_headers_median_HRF_T2DM_Thr);    
    psc_std_parameters_median_HRF_T2DM_Thr_table = array2table(psc_std_parameters_median_HRF_T2DM_Thr,'VariableNames', psc_avg_headers_median_HRF_T2DM_Thr);      

end


if isempty(nsc_avg_parameters_median_HRF_T2DM_Thr)==0                                                                         

    nsc_avg_headers_median_HRF_T2DM_Thr = cell(1,size(nsc_avg_parameters_median_HRF_T2DM_Thr,2));                                       

    for l = 1:length(nsc_avg_std_indexes_median_HRF_T2DM_Thr)
        nsc_avg_headers_median_HRF_T2DM_Thr(l) = HRF_parameters(nsc_avg_std_indexes_median_HRF_T2DM_Thr(l));        % gets the headers of the T2DM Thr median HRF parameters with normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the average and standard deviation for the T2DM Thr
    % median HRF parameters with normal distribution as well as its headers
    nsc_avg_parameters_median_HRF_T2DM_Thr_table = array2table(nsc_avg_parameters_median_HRF_T2DM_Thr,'VariableNames', nsc_avg_headers_median_HRF_T2DM_Thr);    
    nsc_std_parameters_median_HRF_T2DM_Thr_table = array2table(nsc_std_parameters_median_HRF_T2DM_Thr,'VariableNames', nsc_avg_headers_median_HRF_T2DM_Thr);      

end


if isempty(psc_avg_parameters_median_HRF_T2DM_Sub)==0                                                                        

    psc_avg_headers_median_HRF_T2DM_Sub = cell(1,size(psc_avg_parameters_median_HRF_T2DM_Sub,2));                                       

    for l = 1:length(psc_avg_std_indexes_median_HRF_T2DM_Sub)
        psc_avg_headers_median_HRF_T2DM_Sub(l) = HRF_parameters(psc_avg_std_indexes_median_HRF_T2DM_Sub(l));        % gets the headers of the T2DM Sub median HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Sub
    % median HRF parameters with normal distribution as well as its headers
    psc_avg_parameters_median_HRF_T2DM_Sub_table = array2table(psc_avg_parameters_median_HRF_T2DM_Sub,'VariableNames', psc_avg_headers_median_HRF_T2DM_Sub);    
    psc_std_parameters_median_HRF_T2DM_Sub_table = array2table(psc_std_parameters_median_HRF_T2DM_Sub,'VariableNames', psc_avg_headers_median_HRF_T2DM_Sub);      

end


if isempty(nsc_avg_parameters_median_HRF_T2DM_Sub)==0                                                                        

    nsc_avg_headers_median_HRF_T2DM_Sub = cell(1,size(nsc_avg_parameters_median_HRF_T2DM_Sub,2));                                       

    for l = 1:length(nsc_avg_std_indexes_median_HRF_T2DM_Sub)
        nsc_avg_headers_median_HRF_T2DM_Sub(l) = HRF_parameters(nsc_avg_std_indexes_median_HRF_T2DM_Sub(l));        % gets the headers of normal T2DM Sub median HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Sub
    % median HRF parameters with normal distribution as well as its headers
    nsc_avg_parameters_median_HRF_T2DM_Sub_table = array2table(nsc_avg_parameters_median_HRF_T2DM_Sub,'VariableNames', nsc_avg_headers_median_HRF_T2DM_Sub);    
    nsc_std_parameters_median_HRF_T2DM_Sub_table = array2table(nsc_std_parameters_median_HRF_T2DM_Sub,'VariableNames', nsc_avg_headers_median_HRF_T2DM_Sub);      

end


if isempty(psc_avg_parameters_median_HRF_T2DM_Avg)==0                                                                        

    psc_avg_headers_median_HRF_T2DM_Avg = cell(1,size(psc_avg_parameters_median_HRF_T2DM_Avg,2));                                       

    for l = 1:length(psc_avg_std_indexes_median_HRF_T2DM_Avg)
        psc_avg_headers_median_HRF_T2DM_Avg(l) = HRF_parameters(psc_avg_std_indexes_median_HRF_T2DM_Avg(l));        % gets the headers of the T2DM Avg median HRF parameters with normal distribution in the positive signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Avg
    % median HRF parameters with normal distribution as well as its headers
    psc_avg_parameters_median_HRF_T2DM_Avg_table = array2table(psc_avg_parameters_median_HRF_T2DM_Avg,'VariableNames', psc_avg_headers_median_HRF_T2DM_Avg);    
    psc_std_parameters_median_HRF_T2DM_Avg_table = array2table(psc_std_parameters_median_HRF_T2DM_Avg,'VariableNames', psc_avg_headers_median_HRF_T2DM_Avg);      

end


if isempty(nsc_avg_parameters_median_HRF_T2DM_Avg)==0                                                                        

    nsc_avg_headers_median_HRF_T2DM_Avg = cell(1,size(nsc_avg_parameters_median_HRF_T2DM_Avg,2));                                       

    for l = 1:length(nsc_avg_std_indexes_median_HRF_T2DM_Avg)
        nsc_avg_headers_median_HRF_T2DM_Avg(l) = HRF_parameters(nsc_avg_std_indexes_median_HRF_T2DM_Avg(l));        % gets the headers of normal T2DM Avg median HRF parameters with normal distribution in the negative signal change ROIs
    end

    % Forms tables with the average and standard deviation for the T2DM Avg
    % median HRF parameters with normal distribution as well as its headers
    nsc_avg_parameters_median_HRF_T2DM_Avg_table = array2table(nsc_avg_parameters_median_HRF_T2DM_Avg,'VariableNames', nsc_avg_headers_median_HRF_T2DM_Avg);    
    nsc_std_parameters_median_HRF_T2DM_Avg_table = array2table(nsc_std_parameters_median_HRF_T2DM_Avg,'VariableNames', nsc_avg_headers_median_HRF_T2DM_Avg);      

end



% ««««««««««««««««««« Median and inter-quartile range »»»»»»»»»»»»»»»»»»»»


% ---------------------------- Average HRF ------------------------------


% ................................ CNT ..................................


if isempty(psc_median_parameters_avg_HRF_CNT_Thr)==0                                                                            % when an average HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each set of ROI    

    psc_median_headers_avg_HRF_CNT_Thr = cell(1,size(psc_median_parameters_avg_HRF_CNT_Thr,2));                                 % headers for the table

    for l = 1:length(psc_median_iqr_indexes_avg_HRF_CNT_Thr)
        psc_median_headers_avg_HRF_CNT_Thr(l) = HRF_parameters(psc_median_iqr_indexes_avg_HRF_CNT_Thr(l));                      % gets the headers of the CNT Thr average HRF parameters with non-normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Thr
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_avg_HRF_CNT_Thr_table = array2table(psc_median_parameters_avg_HRF_CNT_Thr,'VariableNames', psc_median_headers_avg_HRF_CNT_Thr);          
    psc_iqr_parameters_avg_HRF_CNT_Thr_table = array2table(psc_iqr_parameters_avg_HRF_CNT_Thr,'VariableNames', psc_median_headers_avg_HRF_CNT_Thr);                            

end


if isempty(nsc_median_parameters_avg_HRF_CNT_Thr)==0                                                                                

    nsc_median_headers_avg_HRF_CNT_Thr = cell(1,size(nsc_median_parameters_avg_HRF_CNT_Thr,2));                                             

    for l = 1:length(nsc_median_iqr_indexes_avg_HRF_CNT_Thr)
        nsc_median_headers_avg_HRF_CNT_Thr(l) = HRF_parameters(nsc_median_iqr_indexes_avg_HRF_CNT_Thr(l));                      % gets the headers of the CNT Thr average HRF parameters with non-normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Thr
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_avg_HRF_CNT_Thr_table = array2table(nsc_median_parameters_avg_HRF_CNT_Thr,'VariableNames', nsc_median_headers_avg_HRF_CNT_Thr);          
    nsc_iqr_parameters_avg_HRF_CNT_Thr_table = array2table(nsc_iqr_parameters_avg_HRF_CNT_Thr,'VariableNames', nsc_median_headers_avg_HRF_CNT_Thr);                            

end


if isempty(psc_median_parameters_avg_HRF_CNT_Sub)==0                                                                                          

    psc_median_headers_avg_HRF_CNT_Sub = cell(1,size(psc_median_parameters_avg_HRF_CNT_Sub,2));                                            

    for l = 1:length(psc_median_iqr_indexes_avg_HRF_CNT_Sub)
        psc_median_headers_avg_HRF_CNT_Sub(l) = HRF_parameters(psc_median_iqr_indexes_avg_HRF_CNT_Sub(l));                      % gets the headers of the CNT Sub average HRF parameters with non-normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Sub
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_avg_HRF_CNT_Sub_table = array2table(psc_median_parameters_avg_HRF_CNT_Sub,'VariableNames', psc_median_headers_avg_HRF_CNT_Sub);          
    psc_iqr_parameters_avg_HRF_CNT_Sub_table = array2table(psc_iqr_parameters_avg_HRF_CNT_Sub,'VariableNames', psc_median_headers_avg_HRF_CNT_Sub);               

end


if isempty(nsc_median_parameters_avg_HRF_CNT_Sub)==0                                                                                          

    nsc_median_headers_avg_HRF_CNT_Sub = cell(1,size(nsc_median_parameters_avg_HRF_CNT_Sub,2));                                            

    for l = 1:length(nsc_median_iqr_indexes_avg_HRF_CNT_Sub)
        nsc_median_headers_avg_HRF_CNT_Sub(l) = HRF_parameters(nsc_median_iqr_indexes_avg_HRF_CNT_Sub(l));                      % gets the headers of the CNT Sub average HRF parameters with non-normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Sub
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_avg_HRF_CNT_Sub_table = array2table(nsc_median_parameters_avg_HRF_CNT_Sub, 'VariableNames', nsc_median_headers_avg_HRF_CNT_Sub);          
    nsc_iqr_parameters_avg_HRF_CNT_Sub_table = array2table(nsc_iqr_parameters_avg_HRF_CNT_Sub, 'VariableNames', nsc_median_headers_avg_HRF_CNT_Sub);               

end


if isempty(psc_median_parameters_avg_HRF_CNT_Avg)==0                                                                                          

    psc_median_headers_avg_HRF_CNT_Avg = cell(1,size(psc_median_parameters_avg_HRF_CNT_Avg,2));                                            

    for l = 1:length(psc_median_iqr_indexes_avg_HRF_CNT_Avg)
        psc_median_headers_avg_HRF_CNT_Avg(l) = HRF_parameters(psc_median_iqr_indexes_avg_HRF_CNT_Avg(l));                      % gets the headers of the CNT Avg average HRF parameters with non-normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Avg
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_avg_HRF_CNT_Avg_table = array2table(psc_median_parameters_avg_HRF_CNT_Avg,'VariableNames', psc_median_headers_avg_HRF_CNT_Avg);          
    psc_iqr_parameters_avg_HRF_CNT_Avg_table = array2table(psc_iqr_parameters_avg_HRF_CNT_Avg,'VariableNames', psc_median_headers_avg_HRF_CNT_Avg);               

end


if isempty(nsc_median_parameters_avg_HRF_CNT_Avg)==0                                                                                          

    nsc_median_headers_avg_HRF_CNT_Avg = cell(1,size(nsc_median_parameters_avg_HRF_CNT_Avg,2));                                            

    for l = 1:length(nsc_median_iqr_indexes_avg_HRF_CNT_Avg)
        nsc_median_headers_avg_HRF_CNT_Avg(l) = HRF_parameters(nsc_median_iqr_indexes_avg_HRF_CNT_Avg(l));                      % gets the headers of the CNT Avg average HRF parameters with non-normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Avg
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_avg_HRF_CNT_Avg_table = array2table(nsc_median_parameters_avg_HRF_CNT_Avg, 'VariableNames', nsc_median_headers_avg_HRF_CNT_Avg);          
    nsc_iqr_parameters_avg_HRF_CNT_Avg_table = array2table(nsc_iqr_parameters_avg_HRF_CNT_Avg, 'VariableNames', nsc_median_headers_avg_HRF_CNT_Avg);               

end


% ............................... T2DM ..................................


if isempty(psc_median_parameters_avg_HRF_T2DM_Thr)==0                                                                                   

    psc_median_headers_avg_HRF_T2DM_Thr = cell(1,size(psc_median_parameters_avg_HRF_T2DM_Thr,2));                                          

    for l = 1:length(psc_median_iqr_indexes_avg_HRF_T2DM_Thr)
        psc_median_headers_avg_HRF_T2DM_Thr(l) = HRF_parameters(psc_median_iqr_indexes_avg_HRF_T2DM_Thr(l));                    % gets the headers of T2DM Thr average HRF parameters with non-normal distribution in the positive signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Thr
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_avg_HRF_T2DM_Thr_table = array2table(psc_median_parameters_avg_HRF_T2DM_Thr, 'VariableNames', psc_median_headers_avg_HRF_T2DM_Thr);       
    psc_iqr_parameters_avg_HRF_T2DM_Thr_table = array2table(psc_iqr_parameters_avg_HRF_T2DM_Thr, 'VariableNames', psc_median_headers_avg_HRF_T2DM_Thr);       

end


if isempty(nsc_median_parameters_avg_HRF_T2DM_Thr)==0                                                                                   

    nsc_median_headers_avg_HRF_T2DM_Thr = cell(1,size(nsc_median_parameters_avg_HRF_T2DM_Thr,2));                                          

    for l = 1:length(nsc_median_iqr_indexes_avg_HRF_T2DM_Thr)
        nsc_median_headers_avg_HRF_T2DM_Thr(l) = HRF_parameters(nsc_median_iqr_indexes_avg_HRF_T2DM_Thr(l));                    % gets the headers of T2DM Thr average HRF parameters with non-normal distribution in the negative signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Thr
    % average HRF parameters with non-normal distribution as well as its
    % headers
    nsc_median_parameters_avg_HRF_T2DM_Thr_table = array2table(nsc_median_parameters_avg_HRF_T2DM_Thr, 'VariableNames', nsc_median_headers_avg_HRF_T2DM_Thr);       
    nsc_iqr_parameters_avg_HRF_T2DM_Thr_table = array2table(nsc_iqr_parameters_avg_HRF_T2DM_Thr, 'VariableNames', nsc_median_headers_avg_HRF_T2DM_Thr);       

end


if isempty(psc_median_parameters_avg_HRF_T2DM_Sub)==0                                                                                   

    psc_median_headers_avg_HRF_T2DM_Sub = cell(1,size(psc_median_parameters_avg_HRF_T2DM_Sub,2));                                          

    for l = 1:length(psc_median_iqr_indexes_avg_HRF_T2DM_Sub)
        psc_median_headers_avg_HRF_T2DM_Sub(l) = HRF_parameters(psc_median_iqr_indexes_avg_HRF_T2DM_Sub(l));                    % gets the headers of the T2DM Sub average HRF parameters with non-normal distribution in the positive signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Sub
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_avg_HRF_T2DM_Sub_table = array2table(psc_median_parameters_avg_HRF_T2DM_Sub, 'VariableNames', psc_median_headers_avg_HRF_T2DM_Sub);       
    psc_iqr_parameters_avg_HRF_T2DM_Sub_table = array2table(psc_iqr_parameters_avg_HRF_T2DM_Sub, 'VariableNames', psc_median_headers_avg_HRF_T2DM_Sub);       

end


if isempty(nsc_median_parameters_avg_HRF_T2DM_Sub)==0                                                                                   

    nsc_median_headers_avg_HRF_T2DM_Sub = cell(1,size(nsc_median_parameters_avg_HRF_T2DM_Sub,2));                                          

    for l = 1:length(nsc_median_iqr_indexes_avg_HRF_T2DM_Sub)
        nsc_median_headers_avg_HRF_T2DM_Sub(l) = HRF_parameters(nsc_median_iqr_indexes_avg_HRF_T2DM_Sub(l));                    % gets the headers of the T2DM Sub average HRF parameters with non-normal distribution in the negative signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Sub
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_avg_HRF_T2DM_Sub_table = array2table(nsc_median_parameters_avg_HRF_T2DM_Sub, 'VariableNames', nsc_median_headers_avg_HRF_T2DM_Sub);       
    nsc_iqr_parameters_avg_HRF_T2DM_Sub_table = array2table(nsc_iqr_parameters_avg_HRF_T2DM_Sub, 'VariableNames', nsc_median_headers_avg_HRF_T2DM_Sub);       

end


if isempty(psc_median_parameters_avg_HRF_T2DM_Avg)==0                                                                                   

    psc_median_headers_avg_HRF_T2DM_Avg = cell(1,size(psc_median_parameters_avg_HRF_T2DM_Avg,2));                                          

    for l = 1:length(psc_median_iqr_indexes_avg_HRF_T2DM_Avg)
        psc_median_headers_avg_HRF_T2DM_Avg(l) = HRF_parameters(psc_median_iqr_indexes_avg_HRF_T2DM_Avg(l));                    % gets the headers of the T2DM Avg average HRF parameters with non-normal distribution in the positive signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Avg
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_avg_HRF_T2DM_Avg_table = array2table(psc_median_parameters_avg_HRF_T2DM_Avg, 'VariableNames', psc_median_headers_avg_HRF_T2DM_Avg);       
    psc_iqr_parameters_avg_HRF_T2DM_Avg_table = array2table(psc_iqr_parameters_avg_HRF_T2DM_Avg, 'VariableNames', psc_median_headers_avg_HRF_T2DM_Avg);       

end


if isempty(nsc_median_parameters_avg_HRF_T2DM_Avg)==0                                                                                   

    nsc_median_headers_avg_HRF_T2DM_Avg = cell(1,size(nsc_median_parameters_avg_HRF_T2DM_Avg,2));                                          

    for l = 1:length(nsc_median_iqr_indexes_avg_HRF_T2DM_Avg)
        nsc_median_headers_avg_HRF_T2DM_Avg(l) = HRF_parameters(nsc_median_iqr_indexes_avg_HRF_T2DM_Avg(l));                    % gets the headers of the T2DM Avg average HRF parameters with non-normal distribution in the negative signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Avg
    % average HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_avg_HRF_T2DM_Avg_table = array2table(nsc_median_parameters_avg_HRF_T2DM_Avg, 'VariableNames', nsc_median_headers_avg_HRF_T2DM_Avg);       
    nsc_iqr_parameters_avg_HRF_T2DM_Avg_table = array2table(nsc_iqr_parameters_avg_HRF_T2DM_Avg, 'VariableNames', nsc_median_headers_avg_HRF_T2DM_Avg);       

end



% ---------------------------- Median HRF -------------------------------


% ................................ CNT ..................................


if isempty(psc_median_parameters_median_HRF_CNT_Thr)==0                                                                         % when a median HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each set of ROI    

    psc_median_headers_median_HRF_CNT_Thr = cell(1,size(psc_median_parameters_median_HRF_CNT_Thr,2));                           % headers for the table

    for l = 1:length(psc_median_iqr_indexes_median_HRF_CNT_Thr)
        psc_median_headers_median_HRF_CNT_Thr(l) = HRF_parameters(psc_median_iqr_indexes_median_HRF_CNT_Thr(l));                % gets the headers of the CNT Thr median HRF parameters with non-normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Thr
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_median_HRF_CNT_Thr_table = array2table(psc_median_parameters_median_HRF_CNT_Thr,'VariableNames', psc_median_headers_median_HRF_CNT_Thr);          
    psc_iqr_parameters_median_HRF_CNT_Thr_table = array2table(psc_iqr_parameters_median_HRF_CNT_Thr,'VariableNames', psc_median_headers_median_HRF_CNT_Thr);                            

end


if isempty(nsc_median_parameters_median_HRF_CNT_Thr)==0                                                                                

    nsc_median_headers_median_HRF_CNT_Thr = cell(1,size(nsc_median_parameters_median_HRF_CNT_Thr,2));                                             

    for l = 1:length(nsc_median_iqr_indexes_median_HRF_CNT_Thr)
        nsc_median_headers_median_HRF_CNT_Thr(l) = HRF_parameters(nsc_median_iqr_indexes_median_HRF_CNT_Thr(l));                % gets the headers of the CNT Thr median HRF parameters with non-normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Thr
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_median_HRF_CNT_Thr_table = array2table(nsc_median_parameters_median_HRF_CNT_Thr,'VariableNames', nsc_median_headers_median_HRF_CNT_Thr);          
    nsc_iqr_parameters_median_HRF_CNT_Thr_table = array2table(nsc_iqr_parameters_median_HRF_CNT_Thr,'VariableNames', nsc_median_headers_median_HRF_CNT_Thr);                            

end


if isempty(psc_median_parameters_median_HRF_CNT_Sub)==0                                                                                          

    psc_median_headers_median_HRF_CNT_Sub = cell(1,size(psc_median_parameters_median_HRF_CNT_Sub,2));                                            

    for l = 1:length(psc_median_iqr_indexes_median_HRF_CNT_Sub)
        psc_median_headers_median_HRF_CNT_Sub(l) = HRF_parameters(psc_median_iqr_indexes_median_HRF_CNT_Sub(l));                % gets the headers of the CNT Sub median HRF parameters with non-normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Sub
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_median_HRF_CNT_Sub_table = array2table(psc_median_parameters_median_HRF_CNT_Sub,'VariableNames', psc_median_headers_median_HRF_CNT_Sub);          
    psc_iqr_parameters_median_HRF_CNT_Sub_table = array2table(psc_iqr_parameters_median_HRF_CNT_Sub,'VariableNames', psc_median_headers_median_HRF_CNT_Sub);               

end


if isempty(nsc_median_parameters_median_HRF_CNT_Sub)==0                                                                                          

    nsc_median_headers_median_HRF_CNT_Sub = cell(1,size(nsc_median_parameters_median_HRF_CNT_Sub,2));                                            

    for l = 1:length(nsc_median_iqr_indexes_median_HRF_CNT_Sub)
        nsc_median_headers_median_HRF_CNT_Sub(l) = HRF_parameters(nsc_median_iqr_indexes_median_HRF_CNT_Sub(l));                % gets the headers of the CNT Sub median HRF parameters with non-normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Sub
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_median_HRF_CNT_Sub_table = array2table(nsc_median_parameters_median_HRF_CNT_Sub, 'VariableNames', nsc_median_headers_median_HRF_CNT_Sub);          
    nsc_iqr_parameters_median_HRF_CNT_Sub_table = array2table(nsc_iqr_parameters_median_HRF_CNT_Sub, 'VariableNames', nsc_median_headers_median_HRF_CNT_Sub);               

end


if isempty(psc_median_parameters_median_HRF_CNT_Avg)==0                                                                                          

    psc_median_headers_median_HRF_CNT_Avg = cell(1,size(psc_median_parameters_median_HRF_CNT_Avg,2));                                            

    for l = 1:length(psc_median_iqr_indexes_median_HRF_CNT_Avg)
        psc_median_headers_median_HRF_CNT_Avg(l) = HRF_parameters(psc_median_iqr_indexes_median_HRF_CNT_Avg(l));                % gets the headers of the CNT Avg median HRF parameters with non-normal distribution in the positive signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Avg
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_median_HRF_CNT_Avg_table = array2table(psc_median_parameters_median_HRF_CNT_Avg,'VariableNames', psc_median_headers_median_HRF_CNT_Avg);          
    psc_iqr_parameters_median_HRF_CNT_Avg_table = array2table(psc_iqr_parameters_median_HRF_CNT_Avg,'VariableNames', psc_median_headers_median_HRF_CNT_Avg);               

end


if isempty(nsc_median_parameters_median_HRF_CNT_Avg)==0                                                                                          

    nsc_median_headers_median_HRF_CNT_Avg = cell(1,size(nsc_median_parameters_median_HRF_CNT_Avg,2));                                            

    for l = 1:length(nsc_median_iqr_indexes_median_HRF_CNT_Avg)
        nsc_median_headers_median_HRF_CNT_Avg(l) = HRF_parameters(nsc_median_iqr_indexes_median_HRF_CNT_Avg(l));                % gets the headers of the CNT Avg median HRF parameters with non-normal distribution in the negative signal change ROIs
    end
    
    % Forms tables with the median and interquartile range for the CNT Avg
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_median_HRF_CNT_Avg_table = array2table(nsc_median_parameters_median_HRF_CNT_Avg, 'VariableNames', nsc_median_headers_median_HRF_CNT_Avg);          
    nsc_iqr_parameters_median_HRF_CNT_Avg_table = array2table(nsc_iqr_parameters_median_HRF_CNT_Avg, 'VariableNames', nsc_median_headers_median_HRF_CNT_Avg);               

end


% ............................... T2DM ..................................


if isempty(psc_median_parameters_median_HRF_T2DM_Thr)==0                                                                                   

    psc_median_headers_median_HRF_T2DM_Thr = cell(1,size(psc_median_parameters_median_HRF_T2DM_Thr,2));                                          

    for l = 1:length(psc_median_iqr_indexes_median_HRF_T2DM_Thr)
        psc_median_headers_median_HRF_T2DM_Thr(l) = HRF_parameters(psc_median_iqr_indexes_median_HRF_T2DM_Thr(l));              % gets the headers of T2DM Thr median HRF parameters with non-normal distribution in the positive signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Thr
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_median_HRF_T2DM_Thr_table = array2table(psc_median_parameters_median_HRF_T2DM_Thr, 'VariableNames', psc_median_headers_median_HRF_T2DM_Thr);       
    psc_iqr_parameters_median_HRF_T2DM_Thr_table = array2table(psc_iqr_parameters_median_HRF_T2DM_Thr, 'VariableNames', psc_median_headers_median_HRF_T2DM_Thr);       

end


if isempty(nsc_median_parameters_median_HRF_T2DM_Thr)==0                                                                                   

    nsc_median_headers_median_HRF_T2DM_Thr = cell(1,size(nsc_median_parameters_median_HRF_T2DM_Thr,2));                                          

    for l = 1:length(nsc_median_iqr_indexes_median_HRF_T2DM_Thr)
        nsc_median_headers_median_HRF_T2DM_Thr(l) = HRF_parameters(nsc_median_iqr_indexes_median_HRF_T2DM_Thr(l));              % gets the headers of T2DM Thr median HRF parameters with non-normal distribution in the negative signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Thr
    % median HRF parameters with non-normal distribution as well as its
    % headers
    nsc_median_parameters_median_HRF_T2DM_Thr_table = array2table(nsc_median_parameters_median_HRF_T2DM_Thr, 'VariableNames', nsc_median_headers_median_HRF_T2DM_Thr);       
    nsc_iqr_parameters_median_HRF_T2DM_Thr_table = array2table(nsc_iqr_parameters_median_HRF_T2DM_Thr, 'VariableNames', nsc_median_headers_median_HRF_T2DM_Thr);       

end


if isempty(psc_median_parameters_median_HRF_T2DM_Sub)==0                                                                                   

    psc_median_headers_median_HRF_T2DM_Sub = cell(1,size(psc_median_parameters_median_HRF_T2DM_Sub,2));                                          

    for l = 1:length(psc_median_iqr_indexes_median_HRF_T2DM_Sub)
        psc_median_headers_median_HRF_T2DM_Sub(l) = HRF_parameters(psc_median_iqr_indexes_median_HRF_T2DM_Sub(l));              % gets the headers of the T2DM Sub median HRF parameters with non-normal distribution in the positive signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Sub
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_median_HRF_T2DM_Sub_table = array2table(psc_median_parameters_median_HRF_T2DM_Sub, 'VariableNames', psc_median_headers_median_HRF_T2DM_Sub);       
    psc_iqr_parameters_median_HRF_T2DM_Sub_table = array2table(psc_iqr_parameters_median_HRF_T2DM_Sub, 'VariableNames', psc_median_headers_median_HRF_T2DM_Sub);       

end


if isempty(nsc_median_parameters_median_HRF_T2DM_Sub)==0                                                                                   

    nsc_median_headers_median_HRF_T2DM_Sub = cell(1,size(nsc_median_parameters_median_HRF_T2DM_Sub,2));                                          

    for l = 1:length(nsc_median_iqr_indexes_median_HRF_T2DM_Sub)
        nsc_median_headers_median_HRF_T2DM_Sub(l) = HRF_parameters(nsc_median_iqr_indexes_median_HRF_T2DM_Sub(l));              % gets the headers of the T2DM Sub median HRF parameters with non-normal distribution in the negative signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Sub
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_median_HRF_T2DM_Sub_table = array2table(nsc_median_parameters_median_HRF_T2DM_Sub, 'VariableNames', nsc_median_headers_median_HRF_T2DM_Sub);       
    nsc_iqr_parameters_median_HRF_T2DM_Sub_table = array2table(nsc_iqr_parameters_median_HRF_T2DM_Sub, 'VariableNames', nsc_median_headers_median_HRF_T2DM_Sub);       

end


if isempty(psc_median_parameters_median_HRF_T2DM_Avg)==0                                                                                   

    psc_median_headers_median_HRF_T2DM_Avg = cell(1,size(psc_median_parameters_median_HRF_T2DM_Avg,2));                                          

    for l = 1:length(psc_median_iqr_indexes_median_HRF_T2DM_Avg)
        psc_median_headers_median_HRF_T2DM_Avg(l) = HRF_parameters(psc_median_iqr_indexes_median_HRF_T2DM_Avg(l));              % gets the headers of the T2DM Avg median HRF parameters with non-normal distribution in the positive signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Avg
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    psc_median_parameters_median_HRF_T2DM_Avg_table = array2table(psc_median_parameters_median_HRF_T2DM_Avg, 'VariableNames', psc_median_headers_median_HRF_T2DM_Avg);       
    psc_iqr_parameters_median_HRF_T2DM_Avg_table = array2table(psc_iqr_parameters_median_HRF_T2DM_Avg, 'VariableNames', psc_median_headers_median_HRF_T2DM_Avg);       

end


if isempty(nsc_median_parameters_median_HRF_T2DM_Avg)==0                                                                                   

    nsc_median_headers_median_HRF_T2DM_Avg = cell(1,size(nsc_median_parameters_median_HRF_T2DM_Avg,2));                                          

    for l = 1:length(nsc_median_iqr_indexes_median_HRF_T2DM_Avg)
        nsc_median_headers_median_HRF_T2DM_Avg(l) = HRF_parameters(nsc_median_iqr_indexes_median_HRF_T2DM_Avg(l));              % gets the headers of the T2DM Avg median HRF parameters with non-normal distribution in the negative signal change ROIs
    end

    % Forms tables with the median and interquartile range for the T2DM Avg
    % median HRF parameters with non-normal distribution as well as its 
    % headers
    nsc_median_parameters_median_HRF_T2DM_Avg_table = array2table(nsc_median_parameters_median_HRF_T2DM_Avg, 'VariableNames', nsc_median_headers_median_HRF_T2DM_Avg);       
    nsc_iqr_parameters_median_HRF_T2DM_Avg_table = array2table(nsc_iqr_parameters_median_HRF_T2DM_Avg, 'VariableNames', nsc_median_headers_median_HRF_T2DM_Avg);       

end


% *********************** FDR Correction tables *************************
   

% Column with the statistical test's parameter names to identify its 
% corresponding results    
bh_test_values = array2table({'Hypothesis test result', 'Adjusted p-value'}','VariableNames', {'Value'});  


% ---------------------------- Average HRF ------------------------------


% ««««««««««««««««««««««««« Two sample T-test »»»»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(psc_ttest_avg_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the average HRF parameters in the  
    % Thr condition as well as its headers
    psc_bh_ttest_avg_HRF_thr_table = array2table([psc_avg_HRF_bh_h(1:size(psc_ttest_avg_HRF_thr,2),1) psc_avg_HRF_adj_p(1:size(psc_ttest_avg_HRF_thr,2),1)]','VariableNames', psc_ttest_avg_HRF_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_ttest_avg_HRF_thr_table = [bh_test_values psc_bh_ttest_avg_HRF_thr_table];
    
end


if isempty(nsc_ttest_avg_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the average HRF parameters in the 
    % Thr condition as well as its headers
    nsc_bh_ttest_avg_HRF_thr_table = array2table([nsc_avg_HRF_bh_h(1:size(nsc_ttest_avg_HRF_thr,2),1) nsc_avg_HRF_adj_p(1:size(nsc_ttest_avg_HRF_thr,2),1)]','VariableNames', nsc_ttest_avg_HRF_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_ttest_avg_HRF_thr_table = [bh_test_values nsc_bh_ttest_avg_HRF_thr_table];
    
end


if isempty(psc_ttest_avg_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the average HRF parameters in the 
    % Sub condition as well as its headers
    psc_bh_ttest_avg_HRF_sub_table = array2table([psc_avg_HRF_bh_h(1+size(psc_ttest_avg_HRF_thr,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2),1) psc_avg_HRF_adj_p(1+size(psc_ttest_avg_HRF_thr,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2),1)]','VariableNames', psc_ttest_avg_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_ttest_avg_HRF_sub_table = [bh_test_values psc_bh_ttest_avg_HRF_sub_table];

end


if isempty(nsc_ttest_avg_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the average HRF parameters in the 
    % Sub condition as well as its headers
    nsc_bh_ttest_avg_HRF_sub_table = array2table([nsc_avg_HRF_bh_h(1+size(nsc_ttest_avg_HRF_thr,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2),1) nsc_avg_HRF_adj_p(1+size(nsc_ttest_avg_HRF_thr,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2),1)]','VariableNames', nsc_ttest_avg_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_ttest_avg_HRF_sub_table = [bh_test_values nsc_bh_ttest_avg_HRF_sub_table];

end


if isempty(psc_ttest_avg_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the average HRF parameters in the 
    % Avg condition as well as its headers
    psc_bh_ttest_avg_HRF_avg_table = array2table([psc_avg_HRF_bh_h(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2),1) psc_avg_HRF_adj_p(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2),1)]','VariableNames', psc_ttest_avg_HRF_avg_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_ttest_avg_HRF_avg_table = [bh_test_values psc_bh_ttest_avg_HRF_avg_table];

end


if isempty(nsc_ttest_avg_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the average HRF parameters in the 
    % Avg condition as well as its headers    
    nsc_bh_ttest_avg_HRF_avg_table = array2table([nsc_avg_HRF_bh_h(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2),1) nsc_avg_HRF_adj_p(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2),1)]','VariableNames', nsc_ttest_avg_HRF_avg_headers);

    % Merges the previous table with the statistical test's parameter names
    nsc_bh_ttest_avg_HRF_avg_table = [bh_test_values nsc_bh_ttest_avg_HRF_avg_table];

end



% «««««««««««««««««««««««« Wilcoxon ranksum test »»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(psc_wilcoxon_avg_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the average HRF 
    % parameters in the Thr condition as well as its headers
    psc_bh_wilcoxon_avg_HRF_thr_table = array2table([psc_avg_HRF_bh_h(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2),1) psc_avg_HRF_adj_p(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2),1)]','VariableNames', psc_wilcoxon_avg_HRF_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_wilcoxon_avg_HRF_thr_table = [bh_test_values psc_bh_wilcoxon_avg_HRF_thr_table];

end


if isempty(nsc_wilcoxon_avg_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the average HRF 
    % parameters in the Thr condition as well as its headers
    nsc_bh_wilcoxon_avg_HRF_thr_table = array2table([nsc_avg_HRF_bh_h(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2),1) nsc_avg_HRF_adj_p(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2),1)]','VariableNames', nsc_wilcoxon_avg_HRF_thr_headers);
  
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_wilcoxon_avg_HRF_thr_table = [bh_test_values nsc_bh_wilcoxon_avg_HRF_thr_table];

end


if isempty(psc_wilcoxon_avg_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the average HRF 
    % parameters in the Sub condition as well as its headers
    psc_bh_wilcoxon_avg_HRF_sub_table = array2table([psc_avg_HRF_bh_h(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2)+size(psc_wilcoxon_avg_HRF_sub,2),1) psc_avg_HRF_adj_p(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2):size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2)+size(psc_wilcoxon_avg_HRF_sub,2),1)]','VariableNames', psc_wilcoxon_avg_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_wilcoxon_avg_HRF_sub_table = [bh_test_values psc_bh_wilcoxon_avg_HRF_sub_table];

end


if isempty(nsc_wilcoxon_avg_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the average HRF 
    % parameters in the Sub condition as well as its headers
    nsc_bh_wilcoxon_avg_HRF_sub_table = array2table([nsc_avg_HRF_bh_h(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2)+size(nsc_wilcoxon_avg_HRF_sub,2),1) nsc_avg_HRF_adj_p(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2):size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2)+size(nsc_wilcoxon_avg_HRF_sub,2),1)]','VariableNames', nsc_wilcoxon_avg_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_wilcoxon_avg_HRF_sub_table = [bh_test_values nsc_bh_wilcoxon_avg_HRF_sub_table];
    
end


if isempty(psc_wilcoxon_avg_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the average HRF 
    % parameters in the Avg condition as well as its headers
    psc_bh_wilcoxon_avg_HRF_avg_table = array2table([psc_avg_HRF_bh_h(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2)+size(psc_wilcoxon_avg_HRF_sub,2):end,1) psc_avg_HRF_adj_p(1+size(psc_ttest_avg_HRF_thr,2)+size(psc_ttest_avg_HRF_sub,2)+size(psc_ttest_avg_HRF_avg,2)+size(psc_wilcoxon_avg_HRF_thr,2)+size(psc_wilcoxon_avg_HRF_sub,2):end,1)]','VariableNames', psc_wilcoxon_avg_HRF_avg_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_wilcoxon_avg_HRF_avg_table = [bh_test_values psc_bh_wilcoxon_avg_HRF_avg_table];

end


if isempty(nsc_wilcoxon_avg_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the average HRF 
    % parameters in the Avg condition as well as its headers
    nsc_bh_wilcoxon_avg_HRF_avg_table = array2table([nsc_avg_HRF_bh_h(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2)+size(nsc_wilcoxon_avg_HRF_sub,2):end,1) nsc_avg_HRF_adj_p(1+size(nsc_ttest_avg_HRF_thr,2)+size(nsc_ttest_avg_HRF_sub,2)+size(nsc_ttest_avg_HRF_avg,2)+size(nsc_wilcoxon_avg_HRF_thr,2)+size(nsc_wilcoxon_avg_HRF_sub,2):end,1)]','VariableNames', nsc_wilcoxon_avg_HRF_avg_headers);   
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_wilcoxon_avg_HRF_avg_table = [bh_test_values nsc_bh_wilcoxon_avg_HRF_avg_table];
    
end



% ---------------------------- Median HRF -------------------------------


% ««««««««««««««««««««««««« Two sample T-test »»»»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(psc_ttest_median_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the median HRF parameters in the  
    % Thr condition as well as its headers
    psc_bh_ttest_median_HRF_thr_table = array2table([psc_median_HRF_bh_h(1:size(psc_ttest_median_HRF_thr,2),1) psc_median_HRF_adj_p(1:size(psc_ttest_median_HRF_thr,2),1)]','VariableNames', psc_ttest_median_HRF_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_ttest_median_HRF_thr_table = [bh_test_values psc_bh_ttest_median_HRF_thr_table];
    
end


if isempty(nsc_ttest_median_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the median HRF parameters in the 
    % Thr condition as well as its headers
    nsc_bh_ttest_median_HRF_thr_table = array2table([nsc_median_HRF_bh_h(1:size(nsc_ttest_median_HRF_thr,2),1) nsc_median_HRF_adj_p(1:size(nsc_ttest_median_HRF_thr,2),1)]','VariableNames', nsc_ttest_median_HRF_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_ttest_median_HRF_thr_table = [bh_test_values nsc_bh_ttest_median_HRF_thr_table];
    
end


if isempty(psc_ttest_median_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the median HRF parameters in the 
    % Sub condition as well as its headers
    psc_bh_ttest_median_HRF_sub_table = array2table([psc_median_HRF_bh_h(1+size(psc_ttest_median_HRF_thr,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2),1) psc_median_HRF_adj_p(1+size(psc_ttest_median_HRF_thr,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2),1)]','VariableNames', psc_ttest_median_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_ttest_median_HRF_sub_table = [bh_test_values psc_bh_ttest_median_HRF_sub_table];

end


if isempty(nsc_ttest_median_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the median HRF parameters in the 
    % Sub condition as well as its headers
    nsc_bh_ttest_median_HRF_sub_table = array2table([nsc_median_HRF_bh_h(1+size(nsc_ttest_median_HRF_thr,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2),1) nsc_median_HRF_adj_p(1+size(nsc_ttest_median_HRF_thr,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2),1)]','VariableNames', nsc_ttest_median_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_ttest_median_HRF_sub_table = [bh_test_values nsc_bh_ttest_median_HRF_sub_table];

end


if isempty(psc_ttest_median_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the median HRF parameters in the 
    % Avg condition as well as its headers
    psc_bh_ttest_median_HRF_avg_table = array2table([psc_median_HRF_bh_h(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2),1) psc_median_HRF_adj_p(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2),1)]','VariableNames', psc_ttest_median_HRF_avg_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_ttest_median_HRF_avg_table = [bh_test_values psc_bh_ttest_median_HRF_avg_table];

end


if isempty(nsc_ttest_median_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Two-sample t-test for the median HRF parameters in the 
    % Avg condition as well as its headers    
    nsc_bh_ttest_median_HRF_avg_table = array2table([nsc_median_HRF_bh_h(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2),1) nsc_median_HRF_adj_p(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2),1)]','VariableNames', nsc_ttest_median_HRF_avg_headers);

    % Merges the previous table with the statistical test's parameter names
    nsc_bh_ttest_median_HRF_avg_table = [bh_test_values nsc_bh_ttest_median_HRF_avg_table];

end



% «««««««««««««««««««««««« Wilcoxon ranksum test »»»»»»»»»»»»»»»»»»»»»»»»» 


if isempty(psc_wilcoxon_median_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the median HRF 
    % parameters in the Thr condition as well as its headers
    psc_bh_wilcoxon_median_HRF_thr_table = array2table([psc_median_HRF_bh_h(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2),1) psc_median_HRF_adj_p(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2),1)]','VariableNames', psc_wilcoxon_median_HRF_thr_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_wilcoxon_median_HRF_thr_table = [bh_test_values psc_bh_wilcoxon_median_HRF_thr_table];

end


if isempty(nsc_wilcoxon_median_HRF_thr)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the median HRF 
    % parameters in the Thr condition as well as its headers
    nsc_bh_wilcoxon_median_HRF_thr_table = array2table([nsc_median_HRF_bh_h(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2),1) nsc_median_HRF_adj_p(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2),1)]','VariableNames', nsc_wilcoxon_median_HRF_thr_headers);
  
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_wilcoxon_median_HRF_thr_table = [bh_test_values nsc_bh_wilcoxon_median_HRF_thr_table];

end


if isempty(psc_wilcoxon_median_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the median HRF 
    % parameters in the Sub condition as well as its headers
    psc_bh_wilcoxon_median_HRF_sub_table = array2table([psc_median_HRF_bh_h(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2)+size(psc_wilcoxon_median_HRF_sub,2),1) psc_median_HRF_adj_p(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2):size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2)+size(psc_wilcoxon_median_HRF_sub,2),1)]','VariableNames', psc_wilcoxon_median_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_wilcoxon_median_HRF_sub_table = [bh_test_values psc_bh_wilcoxon_median_HRF_sub_table];

end


if isempty(nsc_wilcoxon_median_HRF_sub)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the median HRF  
    % parameters in the Sub condition as well as its headers
    nsc_bh_wilcoxon_median_HRF_sub_table = array2table([nsc_median_HRF_bh_h(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2)+size(nsc_wilcoxon_median_HRF_sub,2),1) nsc_median_HRF_adj_p(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2):size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2)+size(nsc_wilcoxon_median_HRF_sub,2),1)]','VariableNames', nsc_wilcoxon_median_HRF_sub_headers);
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_wilcoxon_median_HRF_sub_table = [bh_test_values nsc_bh_wilcoxon_median_HRF_sub_table];

end


if isempty(psc_wilcoxon_median_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the median HRF 
    % parameters in the Avg condition as well as its headers
    psc_bh_wilcoxon_median_HRF_avg_table = array2table([psc_median_HRF_bh_h(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2)+size(psc_wilcoxon_median_HRF_sub,2):end,1) psc_median_HRF_adj_p(1+size(psc_ttest_median_HRF_thr,2)+size(psc_ttest_median_HRF_sub,2)+size(psc_ttest_median_HRF_avg,2)+size(psc_wilcoxon_median_HRF_thr,2)+size(psc_wilcoxon_median_HRF_sub,2):end,1)]','VariableNames', psc_wilcoxon_median_HRF_avg_headers);
    
    % Merges the previous table with the statistical test's parameter names
    psc_bh_wilcoxon_median_HRF_avg_table = [bh_test_values psc_bh_wilcoxon_median_HRF_avg_table];

end


if isempty(nsc_wilcoxon_median_HRF_avg)==0
    
    % Forms tables with the Benjamini & Hochberg approach to control the
    % FDR of the Wilcoxon ranksum test results for the median HRF  
    % parameters in the Avg condition as well as its headers
    nsc_bh_wilcoxon_median_HRF_avg_table = array2table([nsc_median_HRF_bh_h(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2)+size(nsc_wilcoxon_median_HRF_sub,2):end,1) nsc_median_HRF_adj_p(1+size(nsc_ttest_median_HRF_thr,2)+size(nsc_ttest_median_HRF_sub,2)+size(nsc_ttest_median_HRF_avg,2)+size(nsc_wilcoxon_median_HRF_thr,2)+size(nsc_wilcoxon_median_HRF_sub,2):end,1)]','VariableNames', nsc_wilcoxon_median_HRF_avg_headers);   
    
    % Merges the previous table with the statistical test's parameter names
    nsc_bh_wilcoxon_median_HRF_avg_table = [bh_test_values nsc_bh_wilcoxon_median_HRF_avg_table];  

end
