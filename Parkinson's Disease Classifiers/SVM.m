%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                       BIOSTATISTICS FINAL PROJECT 
% 
%                       Catarina Guerra | 2015240209
%                       University of Coimbra, 2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all
clc 
 
 
%% EXTRACT INFORMATION

% Load csv data
database = readtable('parkinsons_database.csv');

% Extract feature names 
feature_names = database.Properties.VariableNames;
feature_names = feature_names([2:17 19:end]);           % excludes subject (1) and classification (18) columns

% Extract class column
class = table2array(database(:,18));                    % 0 = healthy control; 1 = Parkinson's Disease

% Rearrange dataset by also excluding subject and classification columns
data = table2array(database(:,[2:17 19:end])); 

% Normalize data
norm_data = zscore(data);


%% INITIALIZATION

% Define number of folds for the K-fold cross-validation used on SVM
k_fold = 5; 

% Perform SVM and get its performance metrics
[n_correct, n_incorrect, accuracy, recall, precision, specificity, f1_score] = performance(norm_data, class, k_fold);

% Get current number of features
n_features = size(norm_data,2);

% To save the removed feature's F1-Scores and names, as well as its
% remaining feature data in each loop iteration. Each index will correspond 
% to the (number of features+1) removed from the dataset
removed_f1_scores = zeros(1, n_features+1);
removed_feature_names = cell(1, n_features+1);
norm_data_cell_list = cell(1,n_features+1);

% Initialize the first instance of the previous arrays with the original 
% dataset's data, F1-Score, and removed feature name 
removed_f1_scores(1,1) = f1_score;                      
removed_feature_names(1,1) = {'None'};              % the original dataset won't have any removed variable
norm_data_cell_list(1,1) = {norm_data(:,:)};

% To increment throughout the while loop for indexing purposes
c = 1;


%% EXECUTION

while n_features > 0
    
    % Initialize in each iteration since its size will change
    f1_score_without_feature_i = zeros(1,n_features);
    
    % Iteratively remove one feature from the dataset and estimate the 
    % F1-Score of the remaining dataset
    for i = 1:n_features  
      
        if i == 1
            norm_data_without_feature_i = norm_data(:, i+1:end);
        else
            norm_data_without_feature_i = norm_data(:, [1:i-1, i+1:end]);
        end
        
        [~,~,~,~,~,~,new_f1_score] = performance(norm_data_without_feature_i, class, k_fold); 
        f1_score_without_feature_i(1,i) = new_f1_score;
  
    end
    
    % Determine the feature index that should be removed in order to have
    % the highest F1-Score
    removed_feature_index = find(f1_score_without_feature_i == max(f1_score_without_feature_i));

    % Add the correspondent F1-Score to an array in which its indices will 
    % represent the (number of features - 1) removed from the dataset 
    removed_f1_scores(1,c+1) = f1_score_without_feature_i(removed_feature_index(1,1));
     
    % Add the correspondent feature name to a cell in which its indices 
    % will represent the features removed from the dataset     
    removed_feature_names(1,c+1) = feature_names(removed_feature_index(1,1));
    
    % Renew the feature names' cell in each loop iteration
    if removed_feature_index(1,1) == 1
        feature_names = feature_names(:, removed_feature_index(1,1)+1:end);
    else
        feature_names = feature_names(:,[1:removed_feature_index(1,1)-1, removed_feature_index(1,1)+1:end]);
    end
    
    % Remove the feature data that allows us to get a higher F1-Score
    if c < (size(removed_feature_names,2)-2)            % in this case 21 so that removing columns within the for loop can work properly
       if removed_feature_index(1,1) == 1
            norm_data = norm_data(:, removed_feature_index(1,1)+1:end);
        else
            norm_data = norm_data(:, [1:removed_feature_index(1,1)-1, removed_feature_index(1,1)+1:end]);
       end
    end
   
    % Save the dataset with the kept features
    norm_data_cell_list(1,c+1) = {norm_data(:,:)};
  
    % Count the number of iterations from 1 onwards
    c = c+1;

    % Decrement the number of features per iteration
    n_features = n_features - 1;
  
end

%% FINAL CALCULATIONS

% Retrieve the maximum F1-Score and the index of its occurrence
max_f1 = max(removed_f1_scores);
index_max_f1 = find(removed_f1_scores == max(removed_f1_scores));
index_max_f1 = index_max_f1(1);

% Get the features that cause the highest F1-Score
best_features = removed_feature_names(index_max_f1+1:end);

% Printing statements
fprintf('Max F1-Score: %.2f \n%', max_f1)
fprintf('Number of variables used: %d \n', size(removed_feature_names(index_max_f1+1:end),2))
fprintf('Selected variables: \n')
fprintf(1, '%s\n', best_features{:})