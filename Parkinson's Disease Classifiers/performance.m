function [n_correct, n_incorrect, accuracy, recall, precision, specificity, f1_score] = performance(x, y, k)

	% Randomly partioning data into k training and test sets
	partioned_data = cvpartition(y, 'Kfold', k);     

	% Variable assignment
	n_cases = size(x,1);
	TP = 0;
	TN = 0;
	FP = 0;
	FN = 0;

    for i = 1:k
		
		% Find the indices for training and test data
        training_idx = find(partioned_data.training(i) == 1);    
        testing_idx = find(partioned_data.test(i) == 1);       
        
		% Get the training data and its original classification
        training_data = x(training_idx,:);          
        training_classification = y(training_idx,:);       
        
		% Get test data and its original classification
		test_data = x(testing_idx,:);           
		test_classification = y(testing_idx,:);  
		
	    % Fit SVM hyperplane into the training data
		SVM_model = fitcsvm(training_data, training_classification);   
		
		% Classify test data according to SVM model
		SVM_classification = predict(SVM_model, test_data);   

        % Get occurrences of all positive and negative cases in the test set
		positive_test_classification = (test_classification==1);
		negative_test_classification = (test_classification==0);
    
		% Get occurrences of all positive and negative cases output by the SVM model
		positive_SVM_classification = (SVM_classification==1);
		negative_SVM_classification = (SVM_classification==0);
		
		% Get TP, TN, FP, FN in each loop
		tp = sum(positive_test_classification .* positive_SVM_classification);
		fn = sum(positive_test_classification) - tp;
		tn = sum(negative_test_classification .* negative_SVM_classification);
		fp = sum(negative_test_classification) - tn;
		
		% Cumulative sum of TP, TN, FP, FN 
		TP = TP + tp;  
		FN = FN + fn;
		TN = TN + tn;              
		FP = FP + fp;
		
    end
       
	
	% Get number of correct and incorrect classifications
	n_correct = TP+TN;
	n_incorrect = FP+FN;
	
	% Get classification report metrics: Accuracy, Precision, Recall, Specificity, F1-Score
	accuracy = 100*(TP+TN)/n_cases;
	recall = 100*TP/(TP+FN);
	precision = 100*TP/(TP+FP);
	specificity = 100*TN/(TN+FP);
	f1_score = 100*(2*TP)/(2*TP+FN+FP);
  
end