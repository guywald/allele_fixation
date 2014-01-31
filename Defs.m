classdef (Sealed) Defs
    %DEFINITIONS Global and constant definitions of the game
    
    properties(Constant)
        ID = 1;
        SEX = 2;
        AGE = 3;
        HAS_DIRECTION = 4;
        INDEX = 5;
        
        MALE = 1;
        FEMALE = 2;
    end
    
    methods (Access = private)
        function out = Defs
        end
    end
end

