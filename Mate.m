

function [offsprings numOfNewChildren] = Mate(inDs, inBoardId, bSize, numOfCreatures, inReproductionAge)      

        offsprings = 0; %Initiaz value

        [l u r d] = GetSurrounding(inDs, inBoardId, bSize,Defs.ID); %ID
        [ls us rs ds] = GetSurrounding(inDs,inBoardId,bSize,Defs.SEX); %SEX
        [la ua ra da] = GetSurrounding(inDs,inBoardId,bSize,Defs.AGE); %AGE

        %Extract only the neighbors of the females:
        possibleMaleNeighbors = ones(4,numOfCreatures);
        possibleMaleNeighbors(1,:) = CalcPossible(ls,la,l,inDs,inReproductionAge);
        possibleMaleNeighbors(2,:) = CalcPossible(us,ua,u,inDs,inReproductionAge);
        possibleMaleNeighbors(3,:) = CalcPossible(rs,ra,r,inDs,inReproductionAge);
        possibleMaleNeighbors(4,:) = CalcPossible(ds,da,d,inDs,inReproductionAge);

        %For each female - pick a random male neighbor
        [r,c] = size(possibleMaleNeighbors);
        [~, idx] = max(rand(r, c).*(possibleMaleNeighbors>0), [], 1);
        randomMaleCreatureForFemales = possibleMaleNeighbors(idx+(0:c-1)*r);           

        %Check for non unique choice:
        sv = sort(randomMaleCreatureForFemales);
        idx = sv(2:end) == sv(1:end-1);

        %The value is the male, the column is the famle:
        malesChosenToMateWith = ~ismember(randomMaleCreatureForFemales,sv(idx)) .* randomMaleCreatureForFemales;               

        matingFemales = inDs.Id .* (malesChosenToMateWith ~= 0); 
        minimizedMatingFemalesList = find(matingFemales);
        minimizedMatingMalesList = malesChosenToMateWith(minimizedMatingFemalesList);
        [~, couplesCount] = size(minimizedMatingFemalesList);
        numOfNewChildren = couplesCount;
        %MATE:
        if (couplesCount > 0)
            fatherAlleles = GetAlleles(inDs,minimizedMatingMalesList);
            motherAlleles = GetAlleles(inDs,minimizedMatingFemalesList);
            
            [allele1Child1 allele2Child1] = MendelianInheritance(fatherAlleles,motherAlleles);
            [allele1Child2 allele2Child2] = MendelianInheritance(fatherAlleles,motherAlleles);
            tmpOffSpringsSex = randi(2,couplesCount,1);
            offsprings = struct('Id',[minimizedMatingFemalesList minimizedMatingMalesList],...
                              'Allele1',[allele1Child1 allele1Child2],...
                              'Allele2',[allele2Child1 allele2Child2],...
                              'Sex',[tmpOffSpringsSex (3-tmpOffSpringsSex)]);
        end        

end
    
function [fatherAllele motherAllele] = MendelianInheritance(inFatherAlleles, inMotherAlleles)
    fatherAllele = ChoseAlleles(inFatherAlleles);
    motherAllele = ChoseAlleles(inMotherAlleles);
end
    
function chosenAllele = ChoseAlleles(inAlleles) 
    %Function chooses from a set of allele pairs of a parent one of
    %it's allels randomly.
    rndAllele = randi(2,size(inAlleles.Allele1,2),1); %Choose allele from a random column (either all1 or all2).
    allelesIdx = sub2ind([size(inAlleles.Allele1,2) 2], 1:size(inAlleles.Allele1,2),rndAllele');         
    
    vec = [inAlleles.Allele1; inAlleles.Allele2];
    chosenAllele = vec(allelesIdx)';
end

function allelesStruct = GetAlleles(inDs, listOfIds)
    allelesStruct = struct('Allele1',inDs.Allele1(listOfIds),'Allele2',inDs.Allele2(listOfIds));
end

function res = CalcPossible(inSex, inAge, inID, inDs, inReproductionAge)
    res = (inSex == Defs.MALE) .* (inAge == inReproductionAge) .* (inDs.Sex==Defs.FEMALE).* (inDs.Age == inReproductionAge) .* inID;
end
