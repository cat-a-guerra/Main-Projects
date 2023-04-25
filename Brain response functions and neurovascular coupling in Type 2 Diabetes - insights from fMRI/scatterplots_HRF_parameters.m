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
load('headers.mat');
load('covariates_HRFdata.mat');
load('avg_covariates_HRFdata.mat');


%% CONSTANTS:

n_parameters = 11;
n_covariates = (2:width(covariates_HRFdata));               % covariate columns used to the plots

labels = {'CNT','T2DM'};   
colors = [0 0.7 0.93; 0.96 0.30 0.33];                      % Thr (1st, blue) and Sub (2nd, red) condition colors

titles = cell(1,n_parameters*size(ROI_regions,2));

% Plot titles: ROI name + HRF parameter of interest
for h = 1:length(headers)
    titles(1,h) = regexprep(headers(h), '_', ' ');
end


%% PLOTTING:

% In this section, we produce descriptive univariate scatter plots of the 
% HRF parameters of interest. These are organized per row and group 
% conditions per line as it follows: 
% 1: Thr CNT | 2: Sub CNT | 3: Thr T2DM | 4: Sub T2DM

  
% *********************** Individual data plots ***************************

r = 1;

for c = 1:length(n_covariates)

    covariate_data = covariates_HRFdata(:,[1,n_covariates(c)]);             % we count the first column since it distinguishes data on the scatterplot

    if ismember(c,1:n_parameters:width(covariates_HRFdata))==1              % peak amplitude (1st parameter out of 11)
        figure('Name', sprintf('%s',ROI_regions{r}))                        % creates a new figure for each ROI
        subplot(3,4,1)
    end
    if ismember(c,2:n_parameters:width(covariates_HRFdata))==1              % peak latency (2nd parameter out of 11) 
        subplot(3,4,2)
    end
    if ismember(c,3:n_parameters:width(covariates_HRFdata))==1              % relative slope to peak (3rd parameter out of 11)
        subplot(3,4,3)
    end
    if ismember(c,4:n_parameters:width(covariates_HRFdata))==1              % area under the curve (4th parameter out of 11)
        subplot(3,4,4)
    end
    if ismember(c,5:n_parameters:width(covariates_HRFdata))==1              % area of the positive curve sections (5th parameter out of 11)
        subplot(3,4,5)
    end
    if ismember(c,6:n_parameters:width(covariates_HRFdata))==1              % first positive curve section area (6th parameter out of 11)
        subplot(3,4,6)
    end
    if ismember(c,7:n_parameters:width(covariates_HRFdata))==1              % second positive curve section area (7th parameter out of 11)
        subplot(3,4,7)
    end
    if ismember(c,8:n_parameters:width(covariates_HRFdata))==1              % third positive curve section area (8th parameter out of 11)
        subplot(3,4,8)
    end
    if ismember(c,9:n_parameters:width(covariates_HRFdata))==1              % area of the negative curve sections (9th parameter out of 11)
        subplot(3,4,9)
    end
    if ismember(c,10:n_parameters:width(covariates_HRFdata))==1             % initial dip area (10th parameter out of 11)
        subplot(3,4,10)
    end
    if ismember(c,11:n_parameters:width(covariates_HRFdata))==1             % undershoot dip area (11th parameter out of 11)
        subplot(3,4,11)
    end

    
    % Scatterplot which compares HRF parameters data from CNT and T2DM 
    % subjects
    [xPoints,yPoints,Label,RangeCut] = UnivarScatter([covariate_data(1:2:height(covariates_HRFdata),:);covariate_data(2:2:height(covariates_HRFdata),:)],'Label',{labels,labels});
    hold on
   
    % Selecting data regarding each HRF parameter in CNT subjects
    xPoints_CNT = xPoints(:,1);
    xPoints_CNT = xPoints_CNT(~isnan(xPoints_CNT));
    yPoints_CNT = yPoints(:,1);
    yPoints_CNT = yPoints_CNT(~isnan(yPoints_CNT));

    % Coloring different condition data regarding each HRF parameter in CNT 
    % subjects
    n_conditions = 0;
    for condition_start = 1:length(xPoints_CNT)/2:length(xPoints_CNT)                   % for each condition (step = length/2 --> number of CNT subjects)
        n_conditions = n_conditions+1;
        condition_subjects = condition_start+(0:(length(xPoints_CNT)/2)-1);             % CNT subjects of each condition
        scatter(xPoints_CNT(condition_subjects),yPoints_CNT(condition_subjects),20,'MarkerEdgeColor','k','MarkerFaceColor',colors(n_conditions,:),'LineWidth',0.25)
    end
   
    
    % Selecting data regarding each HRF parameter in T2DM subjects
    xPoints_T2DM = xPoints(:,2);
    xPoints_T2DM = xPoints_T2DM(~isnan(xPoints_T2DM));
    yPoints_T2DM = yPoints(:,2);
    yPoints_T2DM = yPoints_T2DM(~isnan(yPoints_T2DM));

    % Coloring different condition data regarding each HRF parameter in 
    % T2DM subjects
    n_conditions = 0;
    for condition_start = 1:length(xPoints_T2DM)/2:length(xPoints_T2DM)                 % for each condition (step = length/2 --> number of T2DM subjects)
        n_conditions = n_conditions+1;
        condition_subjects = condition_start+(0:(length(xPoints_T2DM)/2)-1);            % T2DM subjects of each condition
        scatter(xPoints_T2DM(condition_subjects),yPoints_T2DM(condition_subjects),20,'MarkerEdgeColor','k','MarkerFaceColor',colors(n_conditions,:),'LineWidth',0.25)
    end

    
    % Further plot arrangements
    title(titles{(c)})
    
    if ismember(c,1:n_parameters:width(covariates_HRFdata))==1              % peak amplitude
        ylabel('Beta values')
    end
    if ismember(c,2:n_parameters:width(covariates_HRFdata))==1              % peak latency
        ylabel('Time (s)')
    end
    if ismember(c,3:n_parameters:width(covariates_HRFdata))==1              % relative slope to peak
        ylabel('Beta values/Time (s^{-1})')
    end
    if ismember(c,4:n_parameters:width(covariates_HRFdata))==1              % area under the curve
        ylabel('Area units (a.u.)')
    end
    if ismember(c,5:n_parameters:width(covariates_HRFdata))==1              % area of the positive curve sections
        ylabel('Area units (a.u.)')
    end
    if ismember(c,6:n_parameters:width(covariates_HRFdata))==1              % first positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,7:n_parameters:width(covariates_HRFdata))==1              % second positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,8:n_parameters:width(covariates_HRFdata))==1              % third positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,9:n_parameters:width(covariates_HRFdata))==1              % area of the negative curve sections
        ylabel('Area units (a.u.)')
    end
    if ismember(c,10:n_parameters:width(covariates_HRFdata))==1             % initial dip area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,11:n_parameters:width(covariates_HRFdata))==1             % undershoot dip area
        ylabel('Area units (a.u.)')
        r = r+1;                                                            % changes ROI after the plotting each ROI's last parameter
    end
    pbaspect([6,4,1])
end



% ************************ Average data plots ***************************

r = 1;

for c = 1:length(n_covariates)

    avg_covariate_data = avg_covariates_HRFdata(:,[1,n_covariates(c)]);     % we count the first column since it distinguishes data on the scatterplot

    if ismember(c,1:n_parameters:width(avg_covariates_HRFdata))==1          % peak amplitude (1st parameter out of 11)
        figure('Name', sprintf('%s',ROI_regions{r}))                        % creates a new figure for each ROI
        subplot(3,4,1)
    end
    if ismember(c,2:n_parameters:width(avg_covariates_HRFdata))==1          % peak latency (2nd parameter out of 11) 
        subplot(3,4,2)
    end
    if ismember(c,3:n_parameters:width(avg_covariates_HRFdata))==1          % relative slope to peak (3rd parameter out of 11)
        subplot(3,4,3)
    end
    if ismember(c,4:n_parameters:width(avg_covariates_HRFdata))==1          % area under the curve (4th parameter out of 11)
        subplot(3,4,4)
    end
    if ismember(c,5:n_parameters:width(avg_covariates_HRFdata))==1          % area of the positive curve sections (5th parameter out of 11)
        subplot(3,4,5)
    end
    if ismember(c,6:n_parameters:width(avg_covariates_HRFdata))==1          % first positive curve section area (6th parameter out of 11)
        subplot(3,4,6)
    end
    if ismember(c,7:n_parameters:width(avg_covariates_HRFdata))==1          % second positive curve section area (7th parameter out of 11)
        subplot(3,4,7)
    end
    if ismember(c,8:n_parameters:width(avg_covariates_HRFdata))==1          % third positive curve section area (8th parameter out of 11)
        subplot(3,4,8)
    end
    if ismember(c,9:n_parameters:width(avg_covariates_HRFdata))==1          % area of the negative curve sections (9th parameter out of 11)
        subplot(3,4,9)
    end
    if ismember(c,10:n_parameters:width(avg_covariates_HRFdata))==1         % initial dip area (10th parameter out of 11)
        subplot(3,4,10)
    end
    if ismember(c,11:n_parameters:width(avg_covariates_HRFdata))==1         % undershoot dip area (11th parameter out of 11)
        subplot(3,4,11)
    end

    % Scatterplot which compares HRF parameters data from CNT and T2DM 
    % subjects
    [xPoints,yPoints,Label,RangeCut] = UnivarScatter([avg_covariate_data(1:2:height(avg_covariates_HRFdata),:);avg_covariate_data(2:2:height(avg_covariates_HRFdata),:)],'Label',{labels,labels});
    hold on

    % Selecting data regarding each HRF parameter in CNT subjects
    xPoints_CNT = xPoints(:,1);
    xPoints_CNT = xPoints_CNT(~isnan(xPoints_CNT));
    yPoints_CNT = yPoints(:,1);
    yPoints_CNT = yPoints_CNT(~isnan(yPoints_CNT));

    % Coloring different condition data regarding each HRF parameter in CNT 
    % subjects
    n_conditions = 0;
    for condition_start = 1:length(xPoints_CNT)/2:length(xPoints_CNT)                   % for each condition (step = length/2 --> number of CNT subjects)
        n_conditions = n_conditions+1;
        condition_subjects = condition_start+(0:(length(xPoints_CNT)/2)-1);             % CNT subjects of each condition
        scatter(xPoints_CNT(condition_subjects),yPoints_CNT(condition_subjects),20,'MarkerEdgeColor','k','MarkerFaceColor',colors(n_conditions,:),'LineWidth',0.25)
    end
    
    
    % Selecting data regarding each HRF parameter in T2DM subjects
    xPoints_T2DM = xPoints(:,2);
    xPoints_T2DM = xPoints_T2DM(~isnan(xPoints_T2DM));
    yPoints_T2DM = yPoints(:,2);
    yPoints_T2DM = yPoints_T2DM(~isnan(yPoints_T2DM));

    % Coloring different condition data regarding each HRF parameter in 
    % T2DM subjects
    n_conditions = 0;
    for condition_start = 1:length(xPoints_T2DM)/2:length(xPoints_T2DM)                 % for each condition (step = length/2 --> number of T2DM subjects)
        n_conditions = n_conditions+1;
        condition_subjects = condition_start+(0:(length(xPoints_T2DM)/2)-1);            % T2DM subjects of each condition
        scatter(xPoints_T2DM(condition_subjects),yPoints_T2DM(condition_subjects),20,'MarkerEdgeColor','k','MarkerFaceColor',colors(n_conditions,:),'LineWidth',0.25)
    end
    
    
    % Further plot arrangements
    title(titles{(c)})

    if ismember(c,1:n_parameters:width(avg_covariates_HRFdata))==1          % peak amplitude
        ylabel('Beta values')
    end
    if ismember(c,2:n_parameters:width(avg_covariates_HRFdata))==1          % peak latency
        ylabel('Time (s)')
    end
    if ismember(c,3:n_parameters:width(avg_covariates_HRFdata))==1          % relative slope to peak
        ylabel('Beta values/Time (s^{-1})')
    end
    if ismember(c,4:n_parameters:width(avg_covariates_HRFdata))==1          % area under the curve
        ylabel('Area units (a.u.)')
    end
    if ismember(c,5:n_parameters:width(avg_covariates_HRFdata))==1          % area of the positive curve sections
        ylabel('Area units (a.u.)')
    end
    if ismember(c,6:n_parameters:width(avg_covariates_HRFdata))==1          % first positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,7:n_parameters:width(avg_covariates_HRFdata))==1          % second positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,8:n_parameters:width(avg_covariates_HRFdata))==1          % third positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,9:n_parameters:width(avg_covariates_HRFdata))==1          % area of the negative curve sections
        ylabel('Area units (a.u.)')
    end
    if ismember(c,10:n_parameters:width(avg_covariates_HRFdata))==1         % initial dip area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,11:n_parameters:width(avg_covariates_HRFdata))==1         % undershoot dip area
        ylabel('Area units (a.u.)')
        r = r+1;                                                            % changes ROI after the plotting each ROI's last parameter
    end
    pbaspect([6,4,1])
end
   

%% BOXPLOTS 

% In this section, we draw boxplots that show how the values of the HRF 
% parameters of interest were distributed per group and condition in each 
% ROI

% Setting the boxplot's properties
boxplot_positions = [1.35 1.50 1.70 1.85];                                  % position of each group + condition boxplot
boxplot_colors = [[1 0.37 0.37],[0.86 0 0],[0.47 0.86 1],[0 0.54 0.74]];    % colors of each group + condition boxplot - ordered back to front
boxplot_labels = cell(282,1);                                               % data labels of each subject's group + condition

% Identifies each subject's group and condition
for h=1:height(covariates_HRFdata)
    if h<n_CNT*2+1 && rem(h,2)~=0                                           % odd rows up to the 77th subject --> CNT Thr
        boxplot_labels{h} = 'CNT Thr';
    elseif h<n_CNT*2+1 && rem(h,2)==0                                       % even rows up to the 77th subject --> CNT Sub
        boxplot_labels{h} = 'CNT Sub';
    elseif h>=n_CNT*2+1 && h<n_CNT*2+1+n_T2DM && rem(h,2)~=0                % odd rows from the 77th subject forward --> T2DM Thr
        boxplot_labels{h} = 'T2DM Thr';
    else                                                                    % even rows up to the 77th subject forward --> T2DM Sub
        boxplot_labels{h} ='T2DM Sub';
    end
end


r=1;

for c = 1:length(n_covariates)
    
    % Selecting the data from the HRF parameters of interest according to 
    % each group and condition
    CNT_Thr_data = table2array(covariates_HRFdata(1:2:154,c+1));                    
    CNT_Sub_data = table2array(covariates_HRFdata(2:2:154,c+1));                    
    T2DM_Thr_data = table2array(covariates_HRFdata(155:2:end,c+1));                 
    T2DM_Sub_data = table2array(covariates_HRFdata(156:2:end,c+1));                 

    % Rearranging the order of each set of parameters from each group and 
    % condition
    parameter_data = [CNT_Thr_data; CNT_Sub_data; T2DM_Thr_data; T2DM_Sub_data];    

  
    if ismember(c,1:n_parameters:width(covariates_HRFdata))==1              % peak amplitude (1st parameter out of 11)
        figure('Name', sprintf('%s',ROI_regions{r}))                        % creates a new figure for each ROI
        subplot(3,4,1)
    end
    if ismember(c,2:n_parameters:width(covariates_HRFdata))==1              % peak latency (2nd parameter out of 11)  
        subplot(3,4,2)
    end
    if ismember(c,3:n_parameters:width(covariates_HRFdata))==1              % relative slope to peak (3rd parameter out of 11) 
        subplot(3,4,3)
    end
    if ismember(c,4:n_parameters:width(covariates_HRFdata))==1              % area under the curve (4th parameter out of 11)
        subplot(3,4,4)
    end
    if ismember(c,5:n_parameters:width(covariates_HRFdata))==1              % area of the positive curve sections (5th parameter out of 11)
        subplot(3,4,5)
    end
    if ismember(c,6:n_parameters:width(covariates_HRFdata))==1              % first positive curve section area (6th parameter out of 11)
        subplot(3,4,6)
    end
    if ismember(c,7:n_parameters:width(covariates_HRFdata))==1              % second positive curve section area (7th parameter out of 11)
        subplot(3,4,7)
    end
    if ismember(c,8:n_parameters:width(covariates_HRFdata))==1              % third positive curve section area (8th parameter out of 11)
        subplot(3,4,8)
    end
    if ismember(c,9:n_parameters:width(covariates_HRFdata))==1              % area of the negative curve sections (9th parameter out of 11)
        subplot(3,4,9)
    end
    if ismember(c,10:n_parameters:width(covariates_HRFdata))==1             % initial dip area (10th parameter out of 11)
        subplot(3,4,10)
    end
    if ismember(c,11:n_parameters:width(covariates_HRFdata))==1             % undershoot dip area (11th parameter out of 11)
        subplot(3,4,11)
    end
    
    
    % Boxplot which compares HRF parameters data from CNT and T2DM 
    % subjects per condition in each ROI
    boxplot(parameter_data,boxplot_labels,'positions',boxplot_positions,'Colors','k','Symbol',' ')
    xlabel('Group Condition')
    title(titles{(c)})

    % Patching each group + condition boxplot with its corresponding colors 
    box = findobj(gca,'Tag','Box');                                                                         % finds the boxes for each group and condition
    for b=1:length(box)
        patch(get(box(b),'XData'),get(box(b),'YData'),boxplot_colors(1,(b-1)*3+1:b*3),'FaceAlpha',.5);
    end
    
    
    % Further plot arrangements
    if ismember(c,1:n_parameters:width(covariates_HRFdata))==1              % peak amplitude
        ylabel('Beta weights')
    end
    if ismember(c,2:n_parameters:width(covariates_HRFdata))==1              % peak latency 
        ylabel('Time (s)')
    end
    if ismember(c,3:n_parameters:width(covariates_HRFdata))==1              % relative slope to peak 
        ylabel('Beta weights / time (s^{-1})')
    end
    if ismember(c,4:n_parameters:width(covariates_HRFdata))==1              % area under the curve
        ylabel('Area units (a.u.)')
    end
    if ismember(c,5:n_parameters:width(covariates_HRFdata))==1              % area of the positive curve sections
        ylabel('Area units (a.u.)')
    end
    if ismember(c,6:n_parameters:width(covariates_HRFdata))==1              % first positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,7:n_parameters:width(covariates_HRFdata))==1              % second positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,8:n_parameters:width(covariates_HRFdata))==1              % third positive curve section area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,9:n_parameters:width(covariates_HRFdata))==1              % area of the negative curve sections 
        ylabel('Area units (a.u.)')
    end
    if ismember(c,10:n_parameters:width(covariates_HRFdata))==1             % initial dip area
        ylabel('Area units (a.u.)')
    end
    if ismember(c,11:n_parameters:width(covariates_HRFdata))==1             % undershoot dip area 
        ylabel('Area units (a.u.)')
        r=r+1;                                                              % changes ROI after the plotting each ROI's last parameter                                              
    end

    % Renewing the set of parameters from each group and condition for the
    % following ROI
    parameter_data=[];
end
