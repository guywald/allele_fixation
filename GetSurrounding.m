function [left up right down] = GetSurrounding(inDs, inBoardId, bSize, dataType)
%GETSURROUNDING Obtain information about surrounding neighboring cells.
    % dataType: 1 for free positions
    %           2 for neighboring creatures by indexes
    %           3 for neighboring by ids       
    %   1) Get raw indexes of neighbors from a direction
    %   2) Calculate legal neighbors indx values according to border
    %   3) Merge legal binary vector (by multiplying) with neighboring
    %      indxs
    %   4) Update the legal move with more elimination from neighbors:
    %       Find Neighbors by getting from legal indexes on the board places
    %       where a neighbor doesn't exist legal<Direction>
    %   5) Merge with legal border moves
        

        currentLocationIndex = sub2ind([bSize bSize],inDs.Row,inDs.Col);
        
        %Mathematical value of neighboring cells:
        leftIndex = currentLocationIndex - bSize;
        upIndex = currentLocationIndex - 1; 
        rightIndex = currentLocationIndex + bSize;
        downIndex = currentLocationIndex + 1;
        
        %Check if neighboring cells are within the board:
        hasLeft = currentLocationIndex > bSize;
        hasUp = mod(currentLocationIndex,bSize) ~= 1;
        hasRight =  rightIndex <= bSize*bSize;
        hasDown = mod(currentLocationIndex,bSize) ~= 0;
        
        %Obtain data of content of neighboring cells:
        left = GetDataOfDirection(leftIndex, hasLeft, inDs, inBoardId, dataType);
        up = GetDataOfDirection(upIndex, hasUp, inDs, inBoardId, dataType);
        right = GetDataOfDirection(rightIndex, hasRight, inDs, inBoardId, dataType);
        down = GetDataOfDirection(downIndex, hasDown, inDs, inBoardId, dataType);
end

function [directionData] = GetDataOfDirection(inRawIndex, inBoarderHasDirection, inDs, inBoardId, dataType)
%GetDataOfDirection Obtain the data of neighboring cells of a certain
%direction.
    
    %We change the illegal cells which have 0 to 1 so when searching in the
    %creatures data structure the index won't be out of bound. This has no
    %effect since we know with "hasDirection" if that direction is even
    %legal.
    legalDirection = inBoarderHasDirection .* inRawIndex;
    legalDirection(legalDirection==0) = 1;

    %IDs vector:
    directionIds = inBoarderHasDirection .* inBoardId(legalDirection); 
    if (dataType == 1) %# ID
        directionData = directionIds;
        return;
    end

    %SEX vector:
    legalDirectionIds = directionIds;
    legalDirectionIds(legalDirectionIds==0) = 1;
    if (dataType == 2) %# SEX
        directionSex = inDs.Sex(legalDirectionIds) .* (directionIds ~= 0);
        directionData = directionSex;
        return;
    end

    %AGE vector:
    if (dataType == 3) %# AGE
        directionAge = inDs.Age(legalDirectionIds) .* (directionIds ~= 0);
        directionData = directionAge;
        return;
    end

    %OCCUPANCY & TURNS movement ability vector:
    %If moving creature has no more turns than all neighboring cells are
    %not available for him to move.
    hasDirection = inBoarderHasDirection(:)' .* (inBoardId(legalDirection)==0) .* (inDs.Turns > 0); 
    if (dataType == 4) %# HAS_DIRECTION
        directionData = hasDirection;
        return;
    end

    %DIRECTION INDICES (Exclude non relevants with hasDirection):
    directionIndex = hasDirection .* inRawIndex;                 
    directionData = directionIndex;
end