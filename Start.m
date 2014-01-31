function Start()
%% --- General Definitions -------
    clear;
    pause on;
    
    global ds; %# Data structure which holds the creatures
    global board; %# Board of the sim game
    global run; %# Game advance mode (continues or steps)
    global iteration; %# Holds the current sim iteration
    global turnsCalc; %# Holds the Alleles to Turn converter.
    global prevFreq; %# Holder for previous turns allele frequency
    
    %# Game close function
    global closeAll; 
    global isInLoop;

    %# Init of game run mode:
    run=0; 
    isInLoop = 0;
    closeAll = 0;
    iteration = 0;    
     
    %# Initial sim conditions:
    config = struct('BoardDim',25,...
                    'NumOfCreatures',100,...
                    'Initial_AA_Phenotype',25,...
                    'Initial_Aa_Phenotype',50,...
                    'Initial_aa_Phenotype',25,...
                    'Iterations',500,...
                    'AA_Phenotype',30,...
                    'Aa_Phenotype',20,...
                    'aa_Phenotype',10,...
                    'MinBoardSize',3,...
                    'MaxBoardSize',100,...
                    'MinCreatures',2,...
                    'MaxCreatures',9999,...
                    'MaxIterations',1000000,...
                    'MinSteps',0,...
                    'MaxSteps',100,...
                    'HasChanged',0,...
                    'TurnsToDraw',10,...
                    'TurnsToDrawMax',1000,...
                    'ReproductionAge',30,...
                    'ReproductionAgeMax',1000,...
                    'DrawBoard',1);
    
	%# Initial analytical model conditions:
    analytical = struct('A_Frequency',0.5,...
                        'Generations',500,...
                        'MinGenerations',2,...
                        'MaxGenerations',1000000,...
                        'AA_Fitness',1.01,...
                        'Aa_Fitness',1,...
                        'aa_Fitness',0.99);
                    
    BGColor = [0.3922 0.4745 0.6353];                
    PColor = [0.4000 0.6000 1.0000];
                
%% --- General Fucntion -----------------
%# Input legality
function res = NearestNumInRange(inNum, inMinVal, inMaxVal)
%NearestNumInRange returns the nearest number in the rage. If the number is
%in range the same num will be returned. If Out of range, nearest border.
    if (inNum < inMinVal)
        res = inMinVal;
        return;
    end
    
    if (inNum > inMaxVal)
        res = inMaxVal;
        return;
    end
    
    res = inNum;
end
function res = NearestPossibleNumberOfCreatures(inNumberOfCreatures, inBoardDimensions)
    res = min(inBoardDimensions*inBoardDimensions,inNumberOfCreatures) - ((inBoardDimensions*inBoardDimensions) <= inNumberOfCreatures);
end
function res = NearestPossibleBoardDimensions(inNumberOfCreatures, inBoardDimensions)
    res = max(ceil(sqrt(inNumberOfCreatures)),inBoardDimensions) + (inBoardDimensions <= ceil(sqrt(inNumberOfCreatures))) * (ceil(sqrt(inNumberOfCreatures)) ~= config.MaxBoardSize);
end
function res = GetLegalNumeric(inInput, inAlternative)
%GetLegalNumeric calculates if input is number. If not, returns inAlternative.    
    if (isnan(inInput) == 0)
        res = inInput;
    else
       res = inAlternative; 
    end
end
function ChangeFlag(inPrevVal, inNewVal)
%ChangeFlag returns 1 if changes were made, 0 otherwise.    
    if (inPrevVal ~= inNewVal)
        config.HasChanged = 1;
        Reset();
    end
end

function Reset()
%Reset resets many values which are needed for a new run of the simulation.

        run = 0;
        iteration  = 0;
        set(iterNumLabel,'string',num2str(iteration));
        RuntimeEnabled('on');
        
        if (config.Initial_AA_Phenotype + config.Initial_Aa_Phenotype + config.Initial_aa_Phenotype ~= config.NumOfCreatures)
            config.Initial_AA_Phenotype = 0;
            config.Initial_Aa_Phenotype = config.NumOfCreatures;
            config.Initial_aa_Phenotype = 0;
            
             set(initial_AAInput.TextBoxHandler,'String',num2str(config.Initial_AA_Phenotype));    
             set(initial_AaInput.TextBoxHandler,'String',num2str(config.Initial_Aa_Phenotype));    
             set(initial_aaInput.TextBoxHandler,'String',num2str(config.Initial_aa_Phenotype)); 
             InitialAllelesTextBoxCallback();
        end        
        
        [ds board] = Init(config);

        DrawBoard(ds,boardAxes); 
        cla(simFrequencyAxes);
        set(simFrequencyAxes,'YLim',[0 1]);
        prevFreq = FrequencyCount();
        
end

function figure_CloseRequestFcn(hObject, eventdata, handles)
%figure_CloseRequestFcn overrides the close function of the 'X' button. 
    if isequal (get( f, 'waitstatus' ), 'inactive')
        close force;
        closeAll = 1;
        return;
    end

    closeAll = 1;
    uiresume(f);
    try
        if (isInLoop == 0)
            close all;
            return;
        elseif (run == 0)
            uiresume(f);
            return;
        end
    catch
        disp('Some error occured. Brute force stop initiated');
        close all;
        delete(f);
        return;
    end
    clf;
    delete(f);
end

function frequencyStruct = FrequencyCount()
    frequencyStruct = struct('A',(sum(ds.Allele1 == 1) + sum(ds.Allele2 == 1))/(2 * config.NumOfCreatures),...
                             'AA',sum(ds.Allele1 == 1 & ds.Allele2 == 1)/config.NumOfCreatures,...
                             'Aa',sum((ds.Allele1 == 1 & ds.Allele2 == 2) | (ds.Allele1 == 2 & ds.Allele2 == 1))/config.NumOfCreatures,...
                             'aa',sum(ds.Allele1 == 2 & ds.Allele2 == 2)/config.NumOfCreatures);
end

function res = AnalyticalAlleleFrequency(t, x)
    fA = 1-x;
    wAA = analytical.AA_Fitness;
    wAa = analytical.Aa_Fitness;
    waa = analytical.aa_Fitness;
    res = (fA*(1-fA)*(fA*(wAa - wAA) + (1-fA)*(waa - wAa)))/(fA*fA*wAA + 2*fA*(1-fA)*wAa + (1-fA)*(1-fA)*waa);
end

function UpdateChildren(offsprings)
%UpdateChildren - Updates new children/offsprings after mating step.
    ds.Allele1(offsprings.Id) = offsprings.Allele1;
    ds.Allele2(offsprings.Id) = offsprings.Allele2;
    ds.Sex(offsprings.Id) = offsprings.Sex;
    ds.Turns(offsprings.Id) = 0;
    ds.Age(offsprings.Id) = 1;
end

function previousFrequency = UpdatePlot(inAxes, inFrequencyMemory, prevLastFreq, inCurrentIteration)
    currIterMod = mod(inCurrentIteration,config.TurnsToDraw);
    if (currIterMod ~= 0) %If the plot is updated before reaching the interval (when allele fixated before interval).
        span = (inCurrentIteration-currIterMod+1):inCurrentIteration;
        try
        inFrequencyMemory.A = inFrequencyMemory.A(1:length(span));
        inFrequencyMemory.AA = inFrequencyMemory.AA(1:length(span));
        inFrequencyMemory.Aa = inFrequencyMemory.Aa(1:length(span));
        inFrequencyMemory.aa = inFrequencyMemory.aa(1:length(span));
        catch
            disp('Hi');
        end
    else
        span = (inCurrentIteration-config.TurnsToDraw+1):inCurrentIteration;
    end    
        freqSpan = 1:length(span);
        plot(inAxes,[span(1)-1 span],[prevLastFreq.A inFrequencyMemory.A(freqSpan)],'r-',...
                    [span(1)-1 span],1-[prevLastFreq.A inFrequencyMemory.A(freqSpan)],'b-',...
                    [span(1)-1 span],[prevLastFreq.AA inFrequencyMemory.AA(freqSpan)],'g:',...
                    [span(1)-1 span],[prevLastFreq.Aa inFrequencyMemory.Aa(freqSpan)],'k:',...
                    [span(1)-1 span],[prevLastFreq.aa inFrequencyMemory.aa(freqSpan)],'m:'); %Data
        SetSimPlotAxis(inAxes); 
        previousFrequency = struct('A',inFrequencyMemory.A(length(span)),...
                                   'AA',inFrequencyMemory.AA(length(span)),...
                                   'Aa',inFrequencyMemory.Aa(length(span)),...
                                   'aa',inFrequencyMemory.aa(length(span)));
end

function DoPause()
            currIter = iteration;
            DrawBoard(ds,boardAxes);
            UpdatePlot(simFrequencyAxes,frequencyMem,prevFreq,iteration);
            set(iterNumLabel,'string',num2str(iteration));
            uiwait(f);
            %If a reset was done while in step mode, puse again (to allow
            %graphics refreash.
            if (iteration == 0 & currIter ~= 0)
                DoPause();
            end
end

function res = AdvanceIteration()
%AdvanceIteration - If all creatures moved all their turns, advance the
%iteration. If fixation was achieved, stop. If end of iterations reached,
%stop.
    if (sum(sum(ds.Turns)) == 0 | (iteration == 0))
        %Handle Iteration counting:

        set(iterNumLabel,'string',num2str(iteration));
        if (closeAll == 1)
            res = 1;
            return;
        end
        if (iteration == config.Iterations) %Reached end of iterations
            UpdatePlot(simFrequencyAxes,frequencyMem,prevFreq,iteration);
            uiwait(f);
            if (closeAll ~= 1)
                cla(simFrequencyAxes);
                Reset();
            end
            res = 1;
            return;
        end

        iterMod = mod(iteration,config.TurnsToDraw);
        if (iterMod==0 & iteration~=0)
            prevFreq = UpdatePlot(simFrequencyAxes,frequencyMem,prevFreq,iteration);
            frequencyMem = struct('A',zeros(1,config.TurnsToDraw),...
                                  'AA',zeros(1,config.TurnsToDraw),...
                                  'Aa',zeros(1,config.TurnsToDraw),...
                                  'aa',zeros(1,config.TurnsToDraw));
        end
        
        iteration = iteration + 1; 
        currentFreq = FrequencyCount();
        frequencyMem.A(iterMod+1) = currentFreq.A;
        frequencyMem.AA(iterMod+1) = currentFreq.AA;        
        frequencyMem.Aa(iterMod+1) = currentFreq.Aa;
        frequencyMem.aa(iterMod+1) = currentFreq.aa;
        
        set(iterNumLabel,'string',num2str(iteration));
        %#If we've reached fixation:
        if (currentFreq.A == 1 | 1-currentFreq.A == 1)
            run = 0;
            UpdatePlot(simFrequencyAxes,frequencyMem,prevFreq,iteration);
            uiwait(f);
            cla(simFrequencyAxes);
            Reset();
            res = 1;
            return;
        end
        ds.Turns = turnsCalc((ds.Allele1+ds.Allele2)-1); %Update turns according to genetics
    end
    ds.Age = ds.Age + ((ds.Age < config.ReproductionAge) .* (ds.Turns > 0)); %Increment age
    res = 0;
end

function MainLoop()
    
    while(closeAll == 0)
            if (run==0) %Pause if needed:
                DoPause();
            end
            if (closeAll == 1) break; end;
            %When no creature has more iterations then a new turn begins:
            shouldSkipCurrentIteration = AdvanceIteration();
            if (shouldSkipCurrentIteration == 1)
                continue;
            end
            ds.Turns = ds.Turns - (ds.Turns > 0); %Decrese turns by 1 to positive turns
            [ds.Row ds.Col newIdx] = Move(ds,config.NumOfCreatures,board,config.BoardDim);
            board = zeros(config.BoardDim,config.BoardDim);
            board(newIdx) = ds.Id;
            [offsprings numOfNewChildren] = Mate(ds,board,config.BoardDim, config.NumOfCreatures, config.ReproductionAge);
            if (numOfNewChildren ~= 0)
                UpdateChildren(offsprings);
            end
      if (mod(iteration,config.TurnsToDraw)==0)   
        if (closeAll == 1) break; end;
        DrawBoard(ds,boardAxes);
        pause(0.00001);
      end
    end
end




%% --- Graphical Functions --------------
function [bh th] = LegendAndText(inHandle, inColor, inText, inLeftPos, inBottomPos)
    bh = uipanel(inHandle,'Units','characters','BorderType','line','BorderWidth',2,...
        'BackgroundColor',inColor,'Position',[inLeftPos inBottomPos+0.15 4 1.5]);
    th = uicontrol(inHandle,'Style','Text','string',inText,'Units','characters',...
        'Position',[inLeftPos+4 inBottomPos+0.5 8 1]);
end
    
function SetSimPlotAxis(inAxes)
    set(inAxes,'XLim',[0 min(floor(iteration/100)*100+100,config.Iterations)]); 
end

function hs = InputLabelTextSlider(inParent ,inText, inInitialVal, inMinVal, inMaxVal, inTextBoxCallback, inSliderCallback, inLeftPos, inBottomPos)
    hs = struct('GroupHandler','','LabelHandler',0,'TextBoxHandler',0,'SliderHandler',0);

    textLength = length(inText)+2;
    textBoxLength = ceil(log10(inMaxVal))+5;
    sliderLength = 15;
    extraLength = 2;
    totalLength = textLength + textBoxLength + sliderLength + extraLength;

    hs.GroupHandler = uibuttongroup(inParent,...
        'Units','characters',...
        'Bordertype','none',...
        'Position',[inLeftPos inBottomPos totalLength 2.1]);

    hs.LabelHandler = uicontrol(hs.GroupHandler,...
                'Style','Text',...
                'String',inText,...
                'Units','characters',...
                'Position',[0 -0.4 textLength 2],...
                'HorizontalAlignment','left');

    tmp = get(hs.LabelHandler,'Position');         
    hs.TextBoxHandler = uicontrol(hs.GroupHandler,...
                'Style','Edit',...
                'String',inInitialVal,...
                'Units','characters',...
                'Position',[(tmp(1)+tmp(3)+1) 0 textBoxLength 2],...
                'CallBack', inTextBoxCallback);

    tmp = get(hs.TextBoxHandler,'Position');   
    sliderStep = 1 / inMaxVal;
    hs.SliderHandler = uicontrol(hs.GroupHandler,...
                'Style','slider',...
                'Value',inInitialVal,...
                'Units','characters',...
                'Position', [(tmp(1)+tmp(3)+1) 0 sliderLength 2],...
                'Min',inMinVal,...
                'Max',inMaxVal,...
                'SliderStep', [sliderStep sliderStep],...
                'CallBack', inSliderCallback);        
end

function hs = InputLabelText(inParent ,inText, inInitialVal, inMaxVal, inTextBoxCallback, inLeftPos, inBottomPos)
    hs = struct('GroupHandler','','LabelHandler',0,'TextBoxHandler',0);
    
    labelLength = length(inText)+2;
    textBoxLength = ceil(log10(inMaxVal))+7;
    extraLength = 8;
    totalLength = labelLength + textBoxLength + extraLength;    
    
    hs.GroupHandler = uibuttongroup(inParent,...
        'Units','characters',...
        'Bordertype','none',...
        'Position',[inLeftPos inBottomPos totalLength 2.1]);
    
    hs.LabelHandler = uicontrol(hs.GroupHandler,...
                'Style','Text',...
                'Units','characters',...
                'String',inText,...
                'Units','characters',...
                'Position',[0 -0.05 labelLength 2],...
                'HorizontalAlignment','left');

    tmp = get(hs.LabelHandler,'Position');         
    hs.TextBoxHandler = uicontrol(hs.GroupHandler,...
                'Style','Edit',...
                'String',inInitialVal,...
                'Units','characters',...
                'Position',[(tmp(1)+tmp(3)+1) 0 textBoxLength 2],...
                'CallBack', inTextBoxCallback);    
end

function btn = PushButton(inParent, inText, inCallback, inLeftPos, inBottomPos)
    btn = uicontrol(inParent,'Style','PushButton',...
        'String',inText,...
        'Units','characters',...
        'Position',[inLeftPos,inBottomPos,8,2],...
        'CallBack', inCallback);
end

function btn = PushButtonSmall(inParent, inText, inCallback, inLeftPos, inBottomPos)
    btn = PushButton(inParent, inText, inCallback, inLeftPos, inBottomPos);
    set(btn,'Position',[inLeftPos,inBottomPos,8,3.5]);
end

function DrawBoard(inDs,inAxes)
    if (config.DrawBoard == 0)
        return;
    end
    s = config.BoardDim;
    tmpBoard = zeros(s, s);
    tempIdx = sub2ind([s s],inDs.Row,inDs.Col);
    tmpBoard(tempIdx) = ColorConverter(inDs.Age,inDs.Sex);                    

    imagesc(tmpBoard,'Parent',inAxes);
    caxis(boardAxes,[1 25]);

    turnsCalc = [config.AA_Phenotype config.Aa_Phenotype config.aa_Phenotype];
    function [colorValue] = ColorConverter(inAge, inSex)
        colorValue = (inSex * 10) + (inAge==config.ReproductionAge)*5;
        %10 or 20 for youngs
        %11 or 21 for adult
    end
end

%% --- Callback Functions --------------- 
function NumOfCreaturesTextBoxCallback(varargin)
    num = str2double(get(numOfCreaturesBtns.TextBoxHandler,'String'));
    oldVal = config.NumOfCreatures;
    num = GetLegalNumeric(num,oldVal);
    config.NumOfCreatures = NearestNumInRange(num,config.MinCreatures,config.MaxCreatures);
    config.NumOfCreatures = NearestPossibleNumberOfCreatures(config.NumOfCreatures,config.BoardDim);
    set(numOfCreaturesBtns.TextBoxHandler,'String',num2str(config.NumOfCreatures));
    set(numOfCreaturesBtns.SliderHandler,'Value',config.NumOfCreatures);  
    ChangeFlag(config.NumOfCreatures,oldVal);
end

function NumOfCreaturesSliderCallback(varargin)
    num = get(numOfCreaturesBtns.SliderHandler,'Value');
    num = round(num);
    oldVal = config.NumOfCreatures;
    config.NumOfCreatures = num;
    config.NumOfCreatures = NearestPossibleNumberOfCreatures(config.NumOfCreatures,config.BoardDim);
    set(numOfCreaturesBtns.TextBoxHandler,'String',num2str(config.NumOfCreatures));
    set(numOfCreaturesBtns.SliderHandler,'Value',config.NumOfCreatures);  
    ChangeFlag(config.NumOfCreatures,oldVal);
end

function BoardSizeTextBoxCallback(varargin)
    num = str2double(get(boardSizeBtns.TextBoxHandler,'String'));
    oldVal = config.BoardDim;
    num = GetLegalNumeric(num,oldVal);
    config.BoardDim = NearestNumInRange(num,config.MinBoardSize,config.MaxBoardSize);    
    config.BoardDim = NearestPossibleBoardDimensions(config.NumOfCreatures,config.BoardDim);
    set(boardSizeBtns.TextBoxHandler,'String',num2str(config.BoardDim));
    set(boardSizeBtns.SliderHandler,'Value',config.BoardDim);
    ChangeFlag(config.BoardDim,oldVal);
end

function BoardSizeSliderCallback(varargin)    
    num = get(boardSizeBtns.SliderHandler,'Value');
    oldVal = config.BoardDim;    
    num = round(num);
    config.BoardDim = num;
    config.BoardDim = NearestPossibleBoardDimensions(config.NumOfCreatures,config.BoardDim);
    set(boardSizeBtns.TextBoxHandler,'String',num2str(config.BoardDim));
    set(boardSizeBtns.SliderHandler,'Value',config.BoardDim);
    ChangeFlag(config.BoardDim,oldVal);
end

function IterationsTextBoxCallback(varargin)
    num = str2double(get(iterationsBtns.TextBoxHandler,'String'));
    oldVal = config.Iterations;
    num = GetLegalNumeric(num,oldVal);
    config.Iterations = NearestNumInRange(num,1,config.MaxIterations);
    set(iterationsBtns.TextBoxHandler,'String',num2str(config.Iterations));
    set(iterationsBtns.SliderHandler,'Value',config.Iterations);
    ChangeFlag(config.Iterations,oldVal);   
end

function IterationsSliderCallback(varargin)    
    num = get(iterationsBtns.SliderHandler,'Value');
    oldVal = config.Iterations;
    num = GetLegalNumeric(num,oldVal);
    num = round(num);
    config.Iterations = NearestNumInRange(num,1,config.MaxIterations);
    %config.Iterations = num;
    set(iterationsBtns.TextBoxHandler,'String',num2str(config.Iterations));
    ChangeFlag(config.Iterations,oldVal);
    
end

function res = CalculateGoodChange(inOriginalValue, inTextBoxHandler, inMinValue, inMaxValue)
%CalculateGoodChange - Checks the legality of the input, sets a calculated
%legal value to the textbox and updates the slider. The function returns
%the new value so the current calue can be updated. There is min/max value
%as input so the validity of the input is checked.
    num = str2double(get(inTextBoxHandler,'String'));
    oldVal = inOriginalValue;
    num = GetLegalNumeric(num,oldVal);
    num = NearestNumInRange(num,inMinValue,inMaxValue);
    set(inTextBoxHandler,'String',num2str(num)); 
    res = num;
end

function StepsTextBoxCallback(varargin)
    oldVal_AA = config.AA_Phenotype;
    oldVal_Aa = config.Aa_Phenotype;
    oldVal_aa = config.aa_Phenotype;
    
    config.AA_Phenotype = CalculateGoodChange(config.AA_Phenotype,AAInput.TextBoxHandler,config.MinSteps,config.MaxSteps);
    config.Aa_Phenotype = CalculateGoodChange(config.Aa_Phenotype,AaInput.TextBoxHandler,config.MinSteps,config.MaxSteps);
    config.aa_Phenotype = CalculateGoodChange(config.aa_Phenotype,aaInput.TextBoxHandler,config.MinSteps,config.MaxSteps);
    
    ChangeFlag(oldVal_AA,config.AA_Phenotype);
    ChangeFlag(oldVal_Aa,config.Aa_Phenotype);
    ChangeFlag(oldVal_aa,config.aa_Phenotype);      
end   

function InitialAllelesTextBoxCallback(varargin)   
    config.Initial_AA_Phenotype = CalculateGoodChange(config.Initial_AA_Phenotype,initial_AAInput.TextBoxHandler,0,config.NumOfCreatures);
    config.Initial_Aa_Phenotype = CalculateGoodChange(config.Initial_Aa_Phenotype,initial_AaInput.TextBoxHandler,0,config.NumOfCreatures);
    config.Initial_aa_Phenotype = CalculateGoodChange(config.Initial_aa_Phenotype,initial_aaInput.TextBoxHandler,0,config.NumOfCreatures);
    
    if (config.Initial_AA_Phenotype + config.Initial_Aa_Phenotype + config.Initial_aa_Phenotype ~= config.NumOfCreatures)
            set(runBtn,'Enable','off');
            set(stepBtn,'Enable','off');
    else
            set(runBtn,'Enable','on');
            set(stepBtn,'Enable','on');        
    end
    
end  

function AFreqTextBoxCallback(varargin)
    num = str2double(get(freqTextBox,'String'));
    oldVal = analytical.A_Frequency;
    num = GetLegalNumeric(num,oldVal);
    analytical.A_Frequency = NearestNumInRange(num,0,1);
    set(freqTextBox,'String',num2str(analytical.A_Frequency));
    set(freqSlider,'Value',analytical.A_Frequency);  
end

function AFreqSliderCallback(varargin)
    num = get(freqSlider,'Value');
    analytical.A_Frequency = num;
    set(freqTextBox,'String',num2str(analytical.A_Frequency));
end

function FitnessTextBoxCallback(varargin)   
    analytical.AA_Fitness = CalculateGoodChange(analytical.AA_Fitness,fitnessAAInput.TextBoxHandler,0,100);
    analytical.Aa_Fitness = CalculateGoodChange(analytical.Aa_Fitness,fitnessAaInput.TextBoxHandler,0,100);
    analytical.aa_Fitness = CalculateGoodChange(analytical.aa_Fitness,fitnessaaInput.TextBoxHandler,0,100);    
end   

function GenerationsTextBoxCallback(varargin)
    num = str2double(get(generationsBtns.TextBoxHandler,'String'));
    oldVal = analytical.Generations;
    num = GetLegalNumeric(num,oldVal);   
    analytical.Generations = num;
    analytical.Generations = NearestNumInRange(analytical.Generations,analytical.MinGenerations,analytical.MaxGenerations)
    set(generationsBtns.TextBoxHandler,'String',num2str(analytical.Generations));   
    set(generationsBtns.SliderHandler,'Value',analytical.Generations);
    set(analyticalFrequencyAxes,'XLim',[0 analytical.Generations]);
end   

function GenerationsSliderCallback(varargin)    
    num = get(generationsBtns.SliderHandler,'Value');
    num = round(num);
    analytical.Generations = NearestNumInRange(num,1,analytical.MaxGenerations);
    set(generationsBtns.TextBoxHandler,'String',num2str(analytical.Generations));
end

function ResetBtnCallback(h, eventdata)
    Reset();
    uiresume(f);
end

function RunBtnCallback(h, eventdata)
    if (iteration == 0)
        Reset();
    end
    run = 1;
    RuntimeEnabled('off');
    uiresume(f);
end

function StepBtnCallback(h, eventdata)
    run = 0;
    RuntimeEnabled('off');
    uiresume(f);
end

function RuntimeEnabled(inState)
    set(turnsToDrawTextBox,'Enable',inState);
    
    set(reproductionAgeTextBox,'Enable',inState);   
    
    set(initial_AAInput.TextBoxHandler,'Enable',inState);
    set(initial_AaInput.TextBoxHandler,'Enable',inState);
    set(initial_aaInput.TextBoxHandler,'Enable',inState);
    
    set(AAInput.TextBoxHandler,'Enable',inState);
    set(AaInput.TextBoxHandler,'Enable',inState);
    set(aaInput.TextBoxHandler,'Enable',inState);
    
    set(iterationsBtns.SliderHandler,'Enable',inState);
    set(boardSizeBtns.SliderHandler,'Enable',inState);
    set(numOfCreaturesBtns.SliderHandler,'Enable',inState);
    
    set(iterationsBtns.TextBoxHandler,'Enable',inState);
    set(boardSizeBtns.TextBoxHandler,'Enable',inState);
    set(numOfCreaturesBtns.TextBoxHandler,'Enable',inState);    
end

function AnalyzeBtnCallback(h, eventdata)
    [tv f1] = ode23(@AnalyticalAlleleFrequency,[1 analytical.Generations],[analytical.A_Frequency]);
    analyticalPlot = plot(analyticalFrequencyAxes,tv,f1,'b-');      
end

function HWBtnCallback(h, eventdata)
    [tv f1] = ode23(@AnalyticalAlleleFrequency,[1 analytical.Generations],[analytical.A_Frequency]);
    analyticalPlot = plot(analyticalFrequencyAxes,tv,((1-f1).*(1-f1)),'g-',tv,(2*(f1.*(1-f1))),'k-',tv,(f1).*(f1),'m-');      
end

function AnalyzeResetBtnCallback(h, eventdata)
    cla(analyticalFrequencyAxes);
end

function DrawBoardCallback(h, eventdata)
%DrawBoardCallback sets if to draw the board or not.
    config.DrawBoard = get(drawBoardCheckbox,'Value');
    if (config.DrawBoard == 1)
        DrawBoard(ds,boardAxes);
    end
end
 

function TurnsToDrawTextBoxCallback(varargin)
    num = str2double(get(turnsToDrawTextBox,'String'));
    oldVal = config.TurnsToDraw;
    num = GetLegalNumeric(num,oldVal);
    config.TurnsToDraw = NearestNumInRange(num,1,config.TurnsToDrawMax);
    set(turnsToDrawTextBox,'String',num2str(config.TurnsToDraw));
end

function ReproductionAgeTextBoxCallback(varargin)
    num = str2double(get(reproductionAgeTextBox,'String'));
    oldVal = config.ReproductionAge;
    num = GetLegalNumeric(num,oldVal);
    config.ReproductionAge = NearestNumInRange(num,3,config.ReproductionAgeMax);
    set(reproductionAgeTextBox,'String',num2str(config.ReproductionAge));
end

%% --- GUI Layout --
%% Figure
f = figure('Units','pixels',...
        'Position',[100 100 925 850],...
        'HandleVisibility','callback',...
        'CloseRequestFcn',@figure_CloseRequestFcn,...
        'IntegerHandle','off',...
        'Renderer','painters',...
        'Toolbar','none',...
        'MenuBar','none',...
        'NumberTitle','off',...
        'Resize','off',...
        'Name','Allele Fixation - Mini Project');
         

topPanel = uipanel('Units','pixels',...
    'Position', [5 452.5 917 395],...
    'bordertype','none',...
    'BackgroundColor',BGColor,...
    'Parent',f);
    
bottomPanel = uipanel('Units','pixels',...
    'Position', [5 5 917 445],...
    'BackgroundColor',BGColor,...
    'bordertype','none',...
    'Parent',f);    

%% --- Top panel ----------------------
topPanelDim = get(topPanel,'Position');

%% Control Panel
simControlPanel = uipanel('Units','pixels',...
    'bordertype','none',...
    'Position', [5 5 450 topPanelDim(4)-10],...
    'BackgroundColor',PColor,... 
    'Parent',topPanel);    
    
uicontrol(simControlPanel,'Style','Text',...
    'String','Allele Fixation - Guy Wald',...
    'FontSize',18,...
    'Units','characters',...
    'HorizontalAlignment','Center',...
    'Position',[1 26 88 3]);
            
simCPGroup = uipanel(simControlPanel,'Units','characters',...
    'BorderType','line',...
    'BorderWidth',2,...
    'Position',[1 10 88 14],...
    'Title','Simulation control');

simLegend = uipanel(simControlPanel,'Units','characters',...
    'bordertype','line',...
    'BorderWidth',2,...
    'Position',[23 23.6 66 2.2]);
uicontrol(simLegend,'units','characters','string','Legend:',...
    'style','text','Position',[0 0.5 9 1]);
[bh1 th1] = LegendAndText(simLegend, [0.52 0 0],'F Adult',10,0);
[bh2 th2] = LegendAndText(simLegend, [0.875 1 0.1254],'M Adult',24,0);
[bh3 th3] = LegendAndText(simLegend, [1.0 0.314 0],'F Child',38,0);
[bh4 th4] = LegendAndText(simLegend, [0 1 1],'M Child',52,0);


%# Draw every number of turns
uicontrol(simCPGroup,'units','characters','string','Draw every:',...
'style','text','Position',[1 0.5 10 2.3],'HorizontalAlignment','left');
turnsToDrawTextBox = uicontrol(simCPGroup,'String',num2str(config.TurnsToDraw),'units','characters','Position',[9 0.7 10 1.5],...
    'style','edit','callback',@TurnsToDrawTextBoxCallback,'Enable','on');

%# Draw board?
uicontrol(simCPGroup,'units','characters','string','Draw Board:',...
'style','text','Position',[20 1 15 1]);
drawBoardCheckbox = uicontrol(simCPGroup,'Value',1,'units','characters','Position',[34.2 1 3 1],...
    'style','checkbox','callback',@DrawBoardCallback);

%# Reproduction age:
uicontrol(simCPGroup,'units','characters','string','Rep. Age:',...
'style','text','Position',[37 0 15 2]);
reproductionAgeTextBox = uicontrol(simCPGroup,'String',num2str(config.ReproductionAge),'units','characters','Position',[50 0.7 10 1.5],...
    'style','edit','callback',@ReproductionAgeTextBoxCallback,'Enable','on');

resetBtn = PushButton(simCPGroup,'Reset',@ResetBtnCallback,61.8,0.2);
runBtn = PushButton(simCPGroup,'Run',@RunBtnCallback,70,0.2);
stepBtn = PushButton(simCPGroup,'Step',@StepBtnCallback,78.2,0.2);  

numOfCreaturesBtns = InputLabelTextSlider(simCPGroup,...
    'Number of creatures: ',...
    config.NumOfCreatures,config.MinCreatures,config.MaxCreatures,...
    @NumOfCreaturesTextBoxCallback,...
    @NumOfCreaturesSliderCallback,1,8);

boardSizeBtns = InputLabelTextSlider(simCPGroup,...
    'Dimension of board:  ',...
    config.BoardDim,config.MinBoardSize,config.MaxBoardSize,...
    @BoardSizeTextBoxCallback,...
    @BoardSizeSliderCallback,1,5.5);

iterationsBtns = InputLabelTextSlider(simCPGroup,...
    'Number of iterations:',...
    config.Iterations,1,config.MaxIterations,...
    @IterationsTextBoxCallback,...
    @IterationsSliderCallback,1,3);                         

stepsLabel = uicontrol(simCPGroup,...
  'Style','Text',...
  'String','Number of steps (Phenotype):',...
  'Units','characters',...
  'Position',[1 10 25 2.5],...
  'HorizontalAlignment','left');

AAInput = InputLabelText(simCPGroup,...
    'AA:',config.AA_Phenotype,config.MaxSteps,...
    @StepsTextBoxCallback,19,10.5);                      
AaInput = InputLabelText(simCPGroup,...
    'Aa:',config.Aa_Phenotype,config.MaxSteps,...
    @StepsTextBoxCallback,35,10.5);                          
aaInput = InputLabelText(simCPGroup,...
    'aa:',config.aa_Phenotype,config.MaxSteps,...
    @StepsTextBoxCallback,52,10.5);      

initialFreqLabel = uicontrol(simCPGroup,...
  'Style','Text',...
  'String','Initial  Genotypes:',...
  'Units','characters',...
  'Position',[55 6 15 2.5],...
  'HorizontalAlignment','left');

initial_AAInput = InputLabelText(simCPGroup,...
    'AA x',config.Initial_AA_Phenotype,config.MaxCreatures,...
    @InitialAllelesTextBoxCallback,68.5,8);                      
initial_AaInput = InputLabelText(simCPGroup,...
    'Aa x',config.Initial_Aa_Phenotype,config.MaxCreatures,...
    @InitialAllelesTextBoxCallback,68.5,5.5);                          
initial_aaInput = InputLabelText(simCPGroup,...
    'aa x',config.Initial_aa_Phenotype,config.MaxCreatures,...
    @InitialAllelesTextBoxCallback,68.5,3);      

analyticalCPGroup = uipanel(simControlPanel,'Units','characters',...
    'BorderType','line',...
    'BorderWidth',1,...
    'Position',[1 0.5 88 9],...
    'Title','Analytical Model control');
    
    analyzeResetBtn = PushButtonSmall(analyticalCPGroup,'Reset',@AnalyzeResetBtnCallback,61.8,0.2);
    analyzeBtn = PushButtonSmall(analyticalCPGroup,'Allele',@AnalyzeBtnCallback,70,0.2);
    HWBtn = PushButtonSmall(analyticalCPGroup,'H-W',@HWBtnCallback,78.2,0.2);

    stepsLabel = uicontrol(analyticalCPGroup,...
                          'Style','Text',...
                          'String','Fitness:',...
                          'Units','characters',...
                          'Position',[1 5.1 10 2.5],...
                          'HorizontalAlignment','left');
                      
    fitnessAAInput = InputLabelText(analyticalCPGroup,...
        'AA:',analytical.AA_Fitness,config.MaxSteps,...
        @FitnessTextBoxCallback,1,4);                      
    fitnessAaInput = InputLabelText(analyticalCPGroup,...
        'Aa:',analytical.Aa_Fitness,config.MaxSteps,...
        @FitnessTextBoxCallback,18,4);                          
    fitnessaaInput = InputLabelText(analyticalCPGroup,...
        'aa:',analytical.aa_Fitness,config.MaxSteps,...
        @FitnessTextBoxCallback,35,4);                          

    
   freqLabel = uicontrol(analyticalCPGroup,'Style','Text',...
        'String','''a'' frequency:',...
        'Units','characters','Position',[53 6.7 15 1],...
        'HorizontalAlignment','left');
    
   freqTextBox = uicontrol(analyticalCPGroup,'Style','Edit',...
       'Units','characters','Position',[55 4 15 2],...
       'String',analytical.A_Frequency,...
       'Callback',@AFreqTextBoxCallback); 
   
   freqSlider = uicontrol(analyticalCPGroup,...
                'Style','slider',...
                'Value',analytical.A_Frequency,...
                'Units','characters',...
                'Position', [71 4 15 2],...
                'Min',0,...
                'Max',1,...
                'SliderStep', [0.01 0.01],...
                'CallBack', @AFreqSliderCallback);   
            
    generationsBtns = InputLabelTextSlider(analyticalCPGroup,...
        'Number of generations: ',...
        analytical.Generations,analytical.MinGenerations,analytical.MaxGenerations,...
        @GenerationsTextBoxCallback,...
        @GenerationsSliderCallback,1,0.2);
            
            
    
    
%% Simulation board
simBoardPanel = uipanel('Units','pixels',...
    'bordertype','none',...
    'Position', [460 5 450 topPanelDim(4)-10],...
    'BackgroundColor',PColor,... 
    'Parent',topPanel);   

iterNumLabel = uicontrol(simBoardPanel,'style','text', ...
       'string','1', ...
       'fontsize',12, ...
       'position',[0,0,100,20]);

%% Bottom panel
bottomPanelDim = get(bottomPanel,'Position');

analyticalModelPanel = uipanel('Units','pixels',...
    'bordertype','none',...
    'Position', [5 5 450 bottomPanelDim(4)-10],...
    'BackgroundColor',PColor,... 
    'Parent',bottomPanel); 
%% Simulation Frequency Plot:
simFrequencyPanel = uipanel('Units','pixels',...
    'bordertype','none',...
    'Position', [460 5 450 bottomPanelDim(4)-10],...
    'BackgroundColor',PColor,... 
    'Parent',bottomPanel); 
simFrequencyAxes = axes('parent',simFrequencyPanel,...
    'Visible','on', 'Units','pixels','Position',[55 75 375 340]);


SetSimPlotAxis(simFrequencyAxes);
simPlotHandle = plot(simFrequencyAxes,0,0,'r-',0,0,'b-',0,0,'g:',0,0,'k:',0,0,'m:');
simPlotLegend = legend(simPlotHandle,'A','a','AA','Aa','aa');
set(simPlotLegend,'Units','pixels','Position',[55 12 385 20]);

set(simPlotLegend,'Orientation','horizontal');
set(simFrequencyAxes,'NextPlot','add');
xlabel(simFrequencyAxes,'Generation');
ylabel(simFrequencyAxes,'Frequency'); 

%% Analytical Frequency plot:
analyticalFrequencyAxes = axes('parent',analyticalModelPanel,...
    'Visible','on','NextPlot','add','Units','pixels','Position',[55 75 375 340]);
InitAnalyticalFrequencyAxes();

function InitAnalyticalFrequencyAxes()
    caxis(analyticalFrequencyAxes,[1 25]);
    set(analyticalFrequencyAxes,'YLim',[0 1]);
    set(analyticalFrequencyAxes,'XLim',[0 analytical.Generations]);
    xlabel(analyticalFrequencyAxes,'Generation');
    ylabel(analyticalFrequencyAxes,'Frequency');  
    analyticalPlot = plot(analyticalFrequencyAxes,0,0,'b-',0,0,'g-',0,0,'k-',0,0,'m-');
    analyticalPlotLegend = legend(analyticalPlot,'a','AA','Aa','aa');
    %set(analyticalPlotLegend,'Location','SouthOutside');
    set(analyticalPlotLegend,'Orientation','horizontal');  
    set(analyticalPlotLegend,'Units','pixels','Position',[50 12 385 20]);
end



%% Main Run:

[ds board] = Init(config);
boardAxes = axes('parent',simBoardPanel,...
'Visible','off');

DrawBoard(ds,boardAxes)
turnsCalc = [config.AA_Phenotype config.Aa_Phenotype config.aa_Phenotype];
iteration = 0;
isInLoop = 1;
frequencyMem = struct('A',zeros(1,config.TurnsToDraw),...
                      'AA',zeros(1,config.TurnsToDraw),...
                      'Aa',zeros(1,config.TurnsToDraw),...
                      'aa',zeros(1,config.TurnsToDraw));
prevFreq = FrequencyCount();


MainLoop();
close force;

    
end

