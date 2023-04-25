function roots = find_intersection(x,y)

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

    roots = x_intersect;

end