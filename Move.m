function [newRow newCol newIdx] = Move(inDs,numOfCreatures,inBoardId,inBSize) 
%Move creatures randomly one step. If two creatires will move to the
%same cell, both will stay put.
    [l u r d] = GetSurrounding(inDs,inBoardId,inBSize,Defs.HAS_DIRECTION);
    possibleDirections = ones(5,numOfCreatures);
    possibleDirections(1,:) = l;
    possibleDirections(2,:) = u * 2;
    possibleDirections(3,:) = r * 3;
    possibleDirections(4,:) = d * 4;
    possibleDirections(5,:) = 1 * 5;

    %Generate random moves (according to legal moves):  
    [rowCount,colCount] = size(possibleDirections);
    nonZeroCount = sum(possibleDirections ~= 0);
    index = round(rand(1,colCount) .* nonZeroCount +0.5);
    [nonZeroIndices,~] = find(possibleDirections);
    index(2:end) = index(2:end) + cumsum(nonZeroCount(1:end-1));
    randomMoves = possibleDirections(nonZeroIndices(index)+(0:rowCount:(rowCount*colCount-1))');        

    %Move to new position (values in new vector):
    rowChange = [0 -1 0 1 0];
    colChange = [-1 0 1 0 0];
    newRow = inDs.Row + rowChange(randomMoves);
    newCol = inDs.Col + colChange(randomMoves);
    newIdx = sub2ind([inBSize inBSize],newRow,newCol);

    %Check for intersections - move only the non intersecting creatures:
    sv = sort(newIdx);
    idx = sv(2:end) == sv(1:end-1);
    uniqueVector = ~ismember(newIdx,sv(idx));  
    
    newRow = inDs.Row + rowChange(randomMoves) .* uniqueVector;
    newCol = inDs.Col + colChange(randomMoves) .* uniqueVector;
    newIdx = sub2ind([inBSize inBSize],newRow,newCol);        
end