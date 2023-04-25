function area_neg = negative_area(x,y)
 
        yline = 0;
        zero_crossing = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                        % Returns Approximate Zero-Crossing Indices Of Argument Vector(y);
        idx = zero_crossing(y - yline);                                                       % Indices Near ‘y=0’, 
        
        if isempty(idx)                                                                       % No intersection with baseline
           x_intersect=[];
        else
            if idx(end)==length(y)
                for j = 1:numel(idx)-1
                    x_intersect(j) = interp1(y(idx(j):idx(j)+1),x(idx(j):idx(j)+1),yline,'linear','extrap');  % Find ‘x’ At ‘y=0’
                end                                                                               
            else
                for j = 1:numel(idx)
                    x_intersect(j) = interp1(y(idx(j):idx(j)+1),x(idx(j):idx(j)+1),yline,'linear','extrap');  % Find ‘x’ At ‘y=0’
                end
            end
        end
        
        neg=y<0;                                    % Checks y negative values (1 - negative; 0 - positive)        
        x_neg=x'.*neg;                              % Returns the x elements in which y<0, the others are null
        y_neg=y.*neg;                               % Returns the y elements in which y<0, the others are null
        
        x_neg=x_neg(x_neg>0);                       % Vector with every x which corresponds to negative y values
        y_neg=y_neg(y_neg<0);                       % Vector with every negative y value

        if x(1)==0 && neg(1)==1                     % When x=0 is negative
            x_neg=[0; x_neg];
        end
        
        x_total=sort([x_neg',x_intersect]);         % Sorts the set of points: negative and intersection points
        [~,ia,~] = intersect(x_total,x_intersect);  % Gets the intersection indexes between all x points and intersection points relative to all x_points (x_total)
        y_total=ones(1,length(x_total(1,:)));       % Similar to a logical vector
        y_total(ia)=0;                              % Intersection points turn to 0

        if x(1)==0 && neg(1)==0                     % When x=0 is positive
            x(x==0)=[];
            [~,ia,ib] = intersect(x_total,x);       % Gets the intersection indexes between all x points and x values (ia => x_total, ib => x)
            if isempty(ib)==1 && isempty(ia)==1     % No intersection indexes
                y_total(ia)=abs(y(ib));
            else
                y_total(ia)=abs(y(ib+1));           % Corresponding vectores by indexes
            end
        else
            [~,ia,ib] = intersect(x_total,x);       % Gets the intersection indexes between all x points and x values (ia => x_total, ib => x)
            y_total(ia)=abs(y(ib));                 % Corresponding vectores by indexes
        end
        
        area_neg = trapz(x_total,y_total);
        
end