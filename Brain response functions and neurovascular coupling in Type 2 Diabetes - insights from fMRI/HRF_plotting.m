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

% Side note: part of this code was retrieved from http://www.math.mcgill.ca/keith/BICstat/ 

%% CONSTANTS:

%hrf_parameters = [5.4 5.2 10.8 7.35 0.35 0];
hrf_parameters = [5.4 5.2 10.5 7.35 0.35 0];
time = (0:175)/10;                                                                  % time of interest

hrf = fmridesign(time,0,[1 0],[],hrf_parameters);                                   % designs the HRF according to a set of parameters

%% PLOTS:

% In this section, we plot the canonical and the initial dip HRF models


% ------------------------ Canonical HRF model ---------------------------

% Canonical HRF plot
figure
plot(time,hrf,'LineWidth',2)
xlabel('Time (seconds)')
ylabel('Intensity')
title("Glover's Canonical HRF model")


% ----------------------- Initial dip HRF model --------------------------

% Here, we manually create an initial dip according to the canonical HRF we
% previously designed and then we insert it the original HRF. I know it is 
% not the most elegant solution, but at the time it was the only one I 
% could see working. 

% The initial dip resembles to a positive quadratic function shape. Thus,
% to manually create an initial dip, we get the coefficients of a quadratic 
% function which could fit well into our original HRF. 
% According to the shape of our original HRF, we would like to plot a
% positive quadratic function between 0 and 2.1 s, with a vortex at 0.9 s.
% Hence, by using the quadratic function equation (a*x^2 + b*x + c = 0) 
% and the parabola vortex formula (h = -b/(2*a)), and considering the 
% intended requirements, we implemented the following linear equation
% system to solve a, b and c (given we know the x's):

    % equation 1 ==> a*(0^2) + b*0 + c = f(0)
    % equation 2 ==> a*(2.1^2) + b*2.1 + c = f(2.1)
    % equation 3 ==> a*(2*0.9) + b + 0*c = 0 --> from the parabola vortex
% formula: h=0.9, and therefore, 0.9 = -b/(2*a) <==> (2*0.9)*a + b = 0

% Solving the linear equation system
syms a b c
eq1 = a*(time(1)^2) + b*time(1) + c == hrf(1);
eq2 = a*(time(22)^2) + b*time(22) + c == hrf(22);
eq3 = a*(2*0.9) + 1*b + 0 == 0;
[A,B] = equationsToMatrix([eq1, eq2, eq3], [a, b, c]);                              % A - coeficients of a,b,c; B - 2nd term of the equation
result = double(linsolve(A,B));                                                     % solution (a,b,c)

x = (0:21)/10;                                                                      % 21 because we want an initial dip from 0 to 2.1 s (21/10=2.1)
initial_dip_function = (((result(1)).*x.^2+(result(2)).*x+result(3)))';             % quadratic function for the initial dip
id_hrf = [initial_dip_function; hrf(23:end)];                                       % the initial dip HRF will be composed by the initial dip section (0 to 2.1 s) and the original HRF from 2.2 s onwards


% Adjusting the quadratic function and the original HRF slopes so that the 
% resulting HRF appears more homogeneous by using a linear function. 
coeffs = polyfit(time(26:40)', hrf(26:40), 1);                                      % fitting the HRF values from 2.5 to 3.9 s into a linear function
id_hrf(18:25)=(time(18:25)*coeffs(1) + coeffs(2))';                                 % using the linear function to replace the initial dip HRF values from 1.7 to 2.5 s according to the linear function


% Initial dip HRF plot
figure
plot(time,id_hrf,'LineWidth',2)
xlabel('Time (seconds)')
ylabel('Intensity')
title("Initial Dip HRF model")


%% SEGMENTS OF INTEREST - INITIAL DIP HRF:

% In this section, we plot the initial dip HRF and patch each of its 
% segments of interest. To do so, we need to define a parameter for each 
% segment.

figure
total_area = area(time, id_hrf, min(id_hrf), 'FaceColor', [0.4 0.78 0.94]);         % patching the area under the initial dip HRF until its minimum

roots = find_intersection(time,id_hrf);                                             % finds the intersections between the initial dip HRF and y=0


% ------------------------- Negative Sections ----------------------------

neg = id_hrf<0;                                                                     % checks for initial dip HRF negative values (1 - negative; 0 - positive)        
negative_time = time'.*neg;                                                         % returns the time elements in which the initial dip HRF is lower than 0, the others are null
negative_hrf = id_hrf.*neg;                                                         % returns the HRF elements in which the initial dip HRF is lower than 0, the others are null

negative_time = negative_time(negative_time>0);                                     % gets every time point corresponding to negative values of the initial dip HRF
negative_hrf = negative_hrf(negative_hrf<0);                                        % gets every negative value of the initial dip HRF

[peak,volume] = max(id_hrf);                                                        % gets the peak amplitude and and its time volume
peak_latency = time(volume);                                                        % converts time volume into latency (in seconds)
        

% «««««««««««««««««««««««««««« Initial dip »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

% In order to patch this HRF segment, we defined a perimeter comprised of a 
% closed loop which starts and ends at the time point where the initial dip 
% begins and contains the values of HRF and y=0 at the time window 
% corresponding to the initial dip.

% Perimeter part which contains the HRF values and time values from 
% ]beginning of initial dip; end of initial dip] s
if isequal(time(negative_time<=peak_latency),negative_time(negative_time<=peak_latency))==1     % if the HRF is always negative before the initial dip
    initial_dip_time = [negative_time(negative_time<=peak_latency); roots(1)];                  % gets the time points of the negative HRF values prior to the peak latency and of the initial dip -> rise peak transition zero - it is the 1st intersection of the HRF with y=0 if the HRF is always negative before the end of the initial dip
else                                                                                            
    initial_dip_time = [negative_time(negative_time<=peak_latency); roots(2)];                  % gets the time points of the negative HRF values prior to the peak latency and of the initial dip -> rise peak transition zero - it is the 2nd intersection of the HRF with y=0 if the HRF has a zero before the end of the initial dip
end
initial_dip_hrf = [negative_hrf(negative_time<=peak_latency); 0];                               % gets the negative HRF values that were prior to the peak latency and the zero from the initial dip -> rise peak transition

% Perimeter part which contains the y=0 values and time values from 
% [beginning of initial dip; end of initial dip[ s
if isequal(time(negative_time<=peak_latency),negative_time(negative_time<=peak_latency))==1     
    initial_dip_time_y0 = flipud(negative_time(initial_dip_time<roots(1)));                     % retrieves again the time points of the negative HRF values that were prior to the peak latency and flips them up to down
else
    initial_dip_time_y0 = flipud([roots(1); negative_time(initial_dip_time<roots(2))]);         % retrieves again the time points of the negative HRF values that were prior to the peak latency and the time point of the first intersection between the HRF and y=0 and flips them up to down                                
end
initial_dip_hrf_y0 = zeros(length(initial_dip_time_y0),1);                                      % values of y=0 for each of these time points


% We close the (beginning of initial dip -> end of initial dip -> beginning 
% of initial dip) loop in time and values and patch it
initial_dip = patch([initial_dip_time; initial_dip_time_y0],[initial_dip_hrf; initial_dip_hrf_y0],'white');
initial_dip.EdgeColor = 'none';
hatch(initial_dip,20,[1 0.08 0.08],'--',3,1) 


% «««««««««««««««««««««««««««« Undershoot »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

% In order to patch this HRF segment, we defined a perimeter comprised of a 
% closed loop which starts and ends at the peak fall -> undershoot
% transition zero instant and contains the values of HRF and y=0 at the 
% time window corresponding to the undershoot.

% Perimeter part which contains the HRF values and time values from 
% [beginning of undershoot; end of undershoot[ s
if isequal(time(negative_time<=peak_latency),negative_time(negative_time<=peak_latency))==1     
    undershoot_time = [roots(2); negative_time(negative_time>roots(2))];                        % gets the peak fall -> undershoot transition zero - it is the 2nd intersection of the HRF with y=0 if the HRF is always negative before the end of the initial dip - and the time points of the negative HRF values posterior to that instant
    undershoot_hrf = [0; negative_hrf(negative_time>roots(2))];                                 % gets the zero from the peak fall -> undershoot transition and the negative HRF values that were posterior to that instant
else
    undershoot_time = [roots(3); negative_time(negative_time>roots(3))];                        % gets the peak fall -> undershoot transition zero - it is the 3rd intersection of the HRF with y=0 if the HRF has a zero before the end of the initial dip - and the time points of the negative HRF values posterior to that instant
    undershoot_hrf = [0; negative_hrf(negative_time>roots(3))];
end

% Perimeter part which contains the y=0 values and time values from 
% ]beginning of undershoot; end of undershoot] s
if isequal(time(negative_time<=peak_latency),negative_time(negative_time<=peak_latency))==1     
    undershoot_time_y0 = flipud(negative_time(negative_time>roots(2)));                         % retrieves again the time points of the negative HRF values that were posterior to the rise fall -> undershoot transition zero and flips them up to down
else
    undershoot_time_y0 = flipud(negative_time(negative_time>roots(3)));                         
end
undershoot_hrf_y0 = zeros(length(undershoot_time_y0),1);                                        % values of y=0 for each of these time points

% We close the (beginning of undershoot -> end of undershoot -> beginning 
% of undershoot) loop in time and values and patch it
undershoot = patch([undershoot_time; undershoot_time_y0],[undershoot_hrf; undershoot_hrf_y0],'white');
undershoot.EdgeColor = 'none';
hatch(undershoot,27,[0.1 0.23 0.47],'--',3,1)



% ------------------------- Positive Sections ----------------------------
 
pos = id_hrf>0;                                                                     % checks for initial dip HRF positive values (1 - positive; 0 - negative)        
positive_time = time'.*pos;                                                         % returns the time elements in which the initial dip HRF is higher than 0, the others are null
positive_hrf = id_hrf.*pos;                                                         % returns the HRF elements in which the initial dip HRF is higher than 0, the others are null

positive_time = positive_time(positive_time>0);                                     % gets every time point corresponding to positive values of the initial dip HRF
positive_hrf = positive_hrf(positive_hrf>0);                                        % gets every positive value of the initial dip HRF
        

% «««««««««««««««««« First positive HRF section (0-5s) »»»»»»»»»»»»»»»»»»»

% In order to patch this HRF segment, we first have to check if the HRF
% has positive values before the initial dip. If so, this segment is then
% composed of two parts - the positive beginning of the HRF and the rise
% peak - which will be patched separately.

first_positive_hrf = positive_hrf(positive_time<=5);                                % gets every time value corresponding to the values of the first positive HRF section
first_positive_time = positive_time(positive_time<=5);                              % gets every value of the first positive HRF section


% *************** Positive HRF beginning (hypothetical) ******************

% In case there is a positive HRF part before the initial dip, in order to
% patch this HRF segment, we have to define a parameter comprised of a 
% closed loop which starts and ends at t=0s and contains the values of HRF 
% and y=0 at the time window corresponding to that positive HRF beginning.

if isempty(time(positive_time<roots(1)))==0                                                                 % if there are positive HRF values before the initial dip
    
    % Perimeter part which contains the HRF values and time values from 
    % ]positive HRF beginning; first zero-value of HRF] s
    positive_beginning_time = [first_positive_time(first_positive_time<roots(1)); roots(1)];                % gets the time points of the positive HRF values prior to the initial dip and of the positive beginning -> initial dip transition zero - it is the 1st intersection of the HRF with y=0 if the HRF is positive before the initial dip
    positive_beginning_hrf = [first_positive_hrf(first_positive_time<roots(1)); 0];                         % gets the positive HRF values that were prior to the initial dip and the positive HRF beginning -> initial dip transition zero
    
    % Perimeter part which contains the y=0 values and time values from 
    % [positive HRF beginning; first zero-value of HRF[ s
    positive_beginning_time_y0 = flipud([roots(1); first_positive_time(first_positive_time<roots(1))]);     % retrieves again the time points of the positive HRF values that were prior to the initial dip but including the instant of the first intersection between the HRF and y=0 (the positive HRF beginning -> initial dip transition) and flips them up to down       
    positive_beginning_hrf_y0 = zeros(length(positive_beginning_time_y0),1);                                % values of y=0 for each of these time points
    
    % We close the (positive HRF beginning -> first zero-value of HRF -> 
    % positive HRF beginning) loop in time and values and patch it
    positive_beginning = patch([positive_beginning_time; positive_beginning_time_y0],[positive_beginning_hrf; positive_beginning_hrf_y0],[0.4 0.78 0.94]);
    positive_beginning.EdgeColor = 'none';
    hatch(positive_beginning,0,'w',':',6,2)
end


% ******** End of initial dip -> rise to peak transition - 5 s **********

% In order to patch this HRF segment, we defined a perimeter comprised of a 
% closed loop which starts and ends at the end of initial dip -> rise to
% peak transition zero instant and contains the values of HRF and y=0 from
% that instant up to 5 s. 

% Perimeter part which contains the HRF values and time values from 
% [end of initial dip -> rise to peak transition; 5] s
if isempty(time(positive_time<roots(1)))==1 || isempty(time(negative_time<roots(1)))==0         % if there is not a positive HRF beginning or if the HRF is not always negative before the initial dip   
    if roots(2)<5                                                                               % if the initial dip ends before the time point of 5 s
        up_to_five_time = [roots(2); first_positive_time(first_positive_time>roots(2))];        % gets the end of initial dip -> rise to peak transition zero - it is the 2nd intersection of the HRF with y=0 in this case - and the time points of the positive HRF values from that instant to 5 s
        up_to_five_hrf = [0; first_positive_hrf(first_positive_time>roots(2))];                 % gets the zero from the end of initial dip -> rise to peak transition and the positive HRF values from that instant to 5 s
    end
else
    if roots(1)<5
        up_to_five_time = [roots(1); first_positive_time(first_positive_time>roots(1))];        % gets the end of initial dip -> rise to peak transition zero - it is the 1st intersection of the HRF with y=0 if there is a positive HRF beginning or an always negative HRF before the initial dip - and the time points of the positive HRF values from that instant to 5 s
        up_to_five_hrf = [0; first_positive_hrf(first_positive_time>roots(1))];
    end
end

% Perimeter part which contains the y=0 values and time values from 
% ]end of initial dip -> rise to peak transition; 5] s 
if isempty(time(positive_time<roots(1)))==1 || isempty(time(negative_time<roots(1)))==0         
    if roots(2)<5                                                                               
        up_to_five_time_y0 = flipud(first_positive_time(first_positive_time>roots(2)));         % retrieves again the time points of the positive HRF values that were after the initial dip and before 5 s and flips them up to down       
    end
else
    if roots(1)<5
        up_to_five_time_y0 = flipud(first_positive_time(first_positive_time>roots(1)));       
    end
end
up_to_five_hrf_y0 = zeros(length(up_to_five_time_y0),1);                                        % values of y=0 for each of these time points

% We close the (end of initial dip // rise to peak transition -> 5 s -> 
% end of initial dip // rise to peak transition) loop in time and values 
% and patch it
up_to_five = patch([up_to_five_time; up_to_five_time_y0],[up_to_five_hrf; up_to_five_hrf_y0],[0.4 0.78 0.94]);
up_to_five.EdgeColor = 'none';
hatch(up_to_five,0,'w',':',6,2)
      

% ««««««««««««««« Second positive HRF section (5-10 s) »»»»»»»»»»»»»»»»»»

% In order to patch this HRF segment, we defined a perimeter comprised of a 
% closed loop which starts and ends at 5 s, and contains the values of HRF 
% and y=0 at the time window from 5 to 10 s.

second_positive_hrf = positive_hrf(positive_time>=5 & positive_time<10);            % gets every time value corresponding to the values of the second positive HRF section
second_positive_time = positive_time(positive_time>=5 & positive_time<10);          % gets every value of the second positive HRF section


% Perimeter part which contains the HRF values and time values from 
% [5; 10] s
if isempty(time(positive_time<roots(1)))==1 || isempty(time(negative_time<roots(1)))==0         
    if roots(3)<10                                                                              % if the peak fall ends before the time point of 10 s
        five_to_ten_time = [second_positive_time(second_positive_time<roots(3)); roots(3)];     % gets the time points of the positive HRF values from 5 s to the end of peak fall -> undershoot transition zero - it is the 3rd intersection of the HRF with y=0 in this case 
        five_to_ten_hrf = [second_positive_hrf(second_positive_time<roots(3)); 0];              % gets the positive HRF values from 5 s to the end of peak fall -> undershoot transition zero              
    else                                                                                        % if the peak fall ends after the time point of 10 s
        five_to_ten_time = [second_positive_time(second_positive_time<10); 10];                 % gets the time points of the positive HRF values from 5 s to 10 s 
        five_to_ten_hrf = [second_positive_hrf(second_positive_time<10); 0];                    % gets the positive HRF values from 5 s to 10 s
    end
else
    if roots(2)<10
        five_to_ten_time = [second_positive_time(second_positive_time<roots(2)); roots(2)];
        five_to_ten_hrf = [second_positive_hrf(second_positive_time<roots(2)); 0];
     else
        five_to_ten_time = [second_positive_time(second_positive_time<10); 10];
        five_to_ten_hrf = [second_positive_hrf(second_positive_time<10); 0];
    end
end

% Perimeter part which contains the y=0 values and time values from 
% [5; 10[ s 
if isempty(time(positive_time<roots(1)))==1 || isempty(time(negative_time<roots(1)))==0         
    if roots(3)<10                                                                              
        five_to_ten_time_y0 = flipud(second_positive_time(second_positive_time<roots(3)));      % retrieves again the time points of the positive HRF values from 5 s to the end of peak fall -> undershoot transition zero - it is the 3rd intersection of the HRF with y=0 in this case and flips them up to down  
    else                                                                                        
        five_to_ten_time_y0 = flipud(second_positive_time(second_positive_time<10));            % retrieves again the time points of the positive HRF values from 5 s to 10 s and flips them up to down  
    end
else
    if roots(2)<10
        five_to_ten_time_y0 = flipud(second_positive_time(second_positive_time<roots(2)));
     else
        five_to_ten_time_y0 = flipud(second_positive_time(second_positive_time<10));
    end
end
five_to_ten_hrf_y0 = zeros(length(five_to_ten_time_y0),1);                                      % values of y=0 for each of these time points

% We close the (5 s -> 10 s -> 5 s) loop in time and values and patch it
five_to_ten = patch([five_to_ten_time; five_to_ten_time_y0],[five_to_ten_hrf; five_to_ten_hrf_y0],[0.4 0.78 0.94]);
five_to_ten.EdgeColor = 'none';
hatch(five_to_ten,0,[1 0.93 0],':',6,2)
hold on


% ««««««««««««««« Third positive HRF section (10-17.5 s) »»»»»»»»»»»»»»»»

% If there are positive values after the time point of 10 s, we need to
% patch this extra HRF segment. To do so, we defined a perimeter comprised 
% of a closed loop which starts and ends at 10 s, and contains the values 
% of HRF and y=0 at the time window from 10 to 17.5 s (the end of our HRF).

third_positive_hrf = positive_hrf(positive_time>=10 & positive_time<17.5);          % gets every time value corresponding to the values of the third positive HRF section     
third_positive_time = positive_time(positive_time>=10 & positive_time<17.5);        % gets every value of the third positive HRF section


if isempty(third_positive_hrf)==0 && isempty(third_positive_time)==0                % if from 10 to 17.5 s there are positive HRF values (the peak fall -> undershoot transition zero takes place in this time window)
    
    % Perimeter part which contains the HRF values and time values from 
    % [10; 17.5[ s
    if isempty(time(positive_time<roots(1)))==1 || isempty(time(negative_time<roots(1)))==0         
        if roots(3)>10 
            ten_to_end_time = [third_positive_time(third_positive_time<roots(3)); roots(3)];        % gets the time points of the positive HRF values from 10 s to the end of peak fall -> undershoot transition zero - it is the 3rd intersection of the HRF with y=0 in this case 
            ten_to_end_hrf = [third_positive_hrf(third_positive_time<roots(3)); 0];                 % gets the positive HRF values from 10 s to the end of peak fall -> undershoot transition zero              
        end
    else
        if roots(2)>10
             ten_to_end_time = [third_positive_time(third_positive_time<roots(2)); roots(2)];       % gets the time points of the positive HRF values from 10 s to the end of peak fall -> undershoot transition zero - it is the 2rd intersection of the HRF with y=0 in this case 
             ten_to_end_hrf = [third_positive_hrf(third_positive_time<roots(2)); 0];                % gets the positive HRF values from 10 s to the end of peak fall -> undershoot transition zero              
        end
    end

    % Perimeter part which contains the y=0 values and time values from 
    % [10; 17.5[ s 
    if isempty(time(positive_time<roots(1)))==1 || isempty(time(negative_time<roots(1)))==0         
        if roots(3)>10                                                                              
            ten_to_end_time_y0 = flipud(third_positive_time(third_positive_time<roots(3)));         % retrieves again the time points of the positive HRF values from 10 s to the end of peak fall -> undershoot transition zero - it is the 3rd intersection of the HRF with y=0 in this case and flips them up to down  
        end
    else
        if roots(2)>10
            ten_to_end_time_y0 = flipud(third_positive_time(third_positive_time<roots(2)));
        end
    end
    ten_to_end_hrf_y0 = zeros(length(ten_to_end_time_y0),1);                                        % values of y=0 for each of these time points


    % We close the (10 s -> 17.5 s -> 10 s) loop in time and values and 
    % patch it
    ten_to_end = patch([ten_to_end_time; ten_to_end_time_y0],[ten_to_end_hrf; ten_to_end_hrf_y0],[0.4 0.78 0.94]);
    ten_to_end.EdgeColor = 'none';
    hatch(ten_to_end,0,[1 0.48 0.13],':',6,2)
    hold on
end

% Final plot arrangements
plot(time, zeros(1,length(time)), '--k')
hold on 
plot(time,id_hrf,'LineWidth',2, 'MarkerFaceColor', [0 0.4470 0.7410])                   % replotting the initial dip HRF curve
hold on 
set(gca,'XTick',0:2.5:17.5);                                                            % time range: [0; 17.5] s with a 2.5 s step
xlabel('Time (seconds)')
ylabel('Intensity')
title("Initial Dip HRF model - areas of interest")

%% BLOCK STIMULUS CONVOLUTION

% In this section we plot a block stimulus to then be convolved with the
% HRF

% Building the block stimulus
block = zeros(length(time),1);
block(time>=4 & time<8) = 0.2;          
block(time>=12 & time<16) = 0.2;
% We just want to have non-zero values in the time ranges [4-8]s and 
% [12-16]s, the others maintain its original value of zero. 

% Plotting the block estimulus
figure
plot(time,block,'r','LineWidth',2)
xlabel('Time (seconds)')
ylabel('Intensity')
title('Block stimulus')


% Convolving the block stimulus with the HRF
convolution = conv(hrf,block);
convolution_time =(0:(length(convolution)-1))/10;

% Plotting the convolved function
figure
plot(convolution_time,convolution,'g','LineWidth',2);
xlabel('Time (seconds)')
ylabel('Intensity')
title('Block stimulus * HRF Convolution')