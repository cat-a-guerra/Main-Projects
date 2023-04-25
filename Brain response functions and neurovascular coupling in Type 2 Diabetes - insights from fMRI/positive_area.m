function area_pos = positive_area(x,y)
 
        yline = 0;
        zero_crossing = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                        % Returns Approximate Zero-Crossing Indices Of Argument Vector(y);
        idx = zero_crossing(y - yline);                                                       % Indices Near ‘y=1’, % 0 is the function
                
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
        
        pos=y>0;                        % Checks y positive values (1 - positive; 0 - negative)
        x_pos=x'.*pos;                  % Returns the x elements in which y>0, the others are null
        y_pos=y.*pos;                   % Returns the y elements in which y>0, the others are null

        x_pos=x_pos(x_pos>0);           % Vector with every x which corresponds to positive y values
        y_pos=y_pos(y_pos>0);           % Vector with every positive y value
        
        if x(1)==0 && pos(1)==1         % When x=0 is positive 
            x_pos=[0; x_pos];
        end

        x_total=sort([x_pos',x_intersect]);         % Sorts the set of points: positive and intersection points
        [~,ia,~] = intersect(x_total,x_intersect);  % Gets the intersection indexes between all x points and intersection points relative to all x_points (x_total)
        y_total=ones(1,length(x_total(1,:)));       % Similar to a logical vector
        y_total(ia)=0;                              % Intersection points turn to 0
        
        if x(1)==0 && pos(1)==0                     % When x=0 is negative
            x(x==0)=[];
            [~,ia,ib] = intersect(x_total,x);       % Gets the intersection indexes between all x points and x values (ia => x_total, ib => x)
            if isempty(ib)==1 && isempty(ia)==1     % No intersection indexes
                y_total(ia)=y(ib);
            else
                y_total(ia)=y(ib+1);                % Corresponding vectores by indexes
            end
        else
            [~,ia,ib] = intersect(x_total,x);       % Gets the intersection indexes between all x points and x values (ia => x_total, ib => x)
            y_total(ia)=y(ib);                      % Corresponding vectores by indexes
        end
        
        area_pos = trapz(x_total,y_total);

end