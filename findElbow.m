function R = findElbow(singular_values)
    % Converts the singular values into a curve
    numPoints = length(singular_values);
    curve = [1:numPoints; singular_values']';
    
    % Line vector from first to last point
    lineVec = curve(end,:) - curve(1,:);
    lineVecN = lineVec / norm(lineVec);
    
    % Vector between all points and the first point
    vecFromFirst = bsxfun(@minus, curve, curve(1,:));
    
    % Project vecFromFirst onto lineVecN
    scalarProduct = dot(vecFromFirst, repmat(lineVecN,numPoints,1), 2);
    vecFromFirstParallel = scalarProduct * lineVecN;
    
    % Compute the distance from the line to get the perpendicular distance
    vecToLine = vecFromFirst - vecFromFirstParallel;
    distToLine = sqrt(sum(vecToLine.^2, 2));
    
    % Find the index of the point with maximum distance to the line
    [~, R] = max(distToLine);
end
