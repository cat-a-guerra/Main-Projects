
% SCHIZOPHRENIA MORPHOMETRY
% João Duarte | January 2019

% Descriptive univariate scatter plots

% Data is organized in 'covariates' with variables per row and particiapnts
% per line (1-20: BPD | 21-40: CNT | 41-60: SCZ)

clear
close all
clc

% Select covariates to plot
covariateNum = [2,3,7,1,11,17];
% 1 - Age (years)
% 2 - TIV (ml)
% 3 - GM (ml)
% 4 - WM (ml)
% 5 - CSF (ml)
% 6 - WMH (%)
% 7 - GM globus pallidus (ml)
% 8 - Gender
% 9 - Schizobipolar scale
% 10 - Disease onset (years)
% 11 - Disease duration (years)
% 12 - Number of hospitalisations
% 13 - Early stage disease
% 14 - Psychosis history
% 15 - Suicidary hystory
% 16 - Substance abuse history
% 17 - CPZE - clorpromazine
% 18 - DDD - define daily dose
% 19 - BPRS - general psychopathology
% 20 - PANSS - SCZ symptoms}
% 21 - PANSS-P - positive symptoms
% 22 - PANSS-N - negative symptoms
% 23 - PANSS-G - general symptoms
% 24 - PSP - functioning scale
% 25 - ITAQ - insight scale

% Load 'covariates' with data
load('covariates.mat')
labels = {'BPD','CNT','SCZ'};
titles = {'Age (years)', 'TIV (ml)', 'Total GM volume (ml)', 'WM (ml)', 'CSF (ml)', 'WMH (%)', ...
    'Volume GM globus pallidus (ml)', 'Gender',...
    'Schizobipolar scale', 'Disease onset (years)', 'Disease duration (years)',...
    'Number of hospitalisations', 'Early stage disease', 'Psychosis history',...
    'Suicidary hystory', 'Substance abuse history', 'CPZE - clorpromazine',...
    'DDD - define daily dose', 'BPRS - general psychopathology', 'PANSS - SCZ symptoms',...
    'PANSS-P - positive symptoms', 'PANSS-N - negative symptoms', 'PANSS-G - general symptoms',...
    'PSP - functioning scale', 'ITAQ - insight scale'
    };

% % Choose colors for groups or keep default colors
% answer = questdlg('Do you want to choose the colors from the palette?','','Yes','Use Default','Use Default');
% if strcmp('Yes',answer)
%     Colors = ColorCoder(3);
% else
Colors = [0.93 0.69 0.13; 0.47 0.67 0.19; 0.85 0.33 0.10];
% end

% Generate plots for each selected covariate
for c = 1:length(covariateNum)
    
    covariateData = covariates(:,[1,covariateNum(c)+1]);
    Colors = [0.93 0.69 0.13; 0.47 0.67 0.19; 0.85 0.33 0.10];
    
    if ismember(covariateNum(c),[9:25])
        labels = {'BPD','SCZ'};
        covariateData = covariates([1:20,41:60],[1,covariateNum(c)+1]);
        Colors = [0.93 0.69 0.13; 0.85 0.33 0.10];
    end
    
    if ismember(covariateNum(c),[20,21,22,23])
        labels = {'SCZ'};
        covariateData = covariates([41:60],[1,covariateNum(c)+1]);
        Colors = [0.85 0.33 0.10];
    end
       
    
    %     subplot(2,3,c)
    figure
    [~,~,~,RangeCut] = UnivarScatter(covariateData,'Label',labels,'MarkerFaceColor',Colors,'Width',0.4,'Compression',5);
    %         UnivarScatter(covariateData,'Label',labels,'MarkerFaceColor',Colors,'RangeCut',RangeCut+3);
    ylabel(titles{covariateNum(c)})
    %     title(titles{covariateNum(c)})
    pbaspect([6,4,1])
    %     hold on
    
    %     % With light/dark grey boxes for SEM/std
    %     figure
    %     UnivarScatter(covariateData,'Label',labels,'MarkerFaceColor',Colors,'Width',0.4,'Compression',5,'RangeCut',25);
    %     title(titles{covariateNum(c)})
    %     pbaspect([6,4,1])
    
    %     % % With same-group-color boxes for SEM/std
    %     figure
    %     [~,~,~,RangeCut] = UnivarScatter(covariateData,'Label',labels,'MarkerFaceColor',Colors,'SEMColor',Colors/1.5,'StdColor',Colors/2);
    %     title(titles{covariateNum(c)})
    %     % pbaspect([6,4,1])
    
end






