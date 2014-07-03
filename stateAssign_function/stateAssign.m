function [newData eventData eventSequences] = stateAssign(data, settings)
% Allows user to make state assignments for single molecule fluorescent
% observations over multiple color channels and multiple assignment
% classes.
%
% Tested in Matlab R2012a
%
% INPUTS:
%
% -data: Array. where each row is an observation frame and each column is fluorescence intensity from each fluorophore in that frame.
% First column is time of each frame (assumed ot be in seconds with each frame representing the same ammount of time).
% Following columns are groups of fluorophores grouped by columns.
% For example:
%
% data =
% 1 0 0 0 0
% 2 1 0 0 0
% 3 0 1 0 0
% 4 0 0 1 0
%
% is a single trace, taken over four seconds at a 1 second framerate, with
% four colors. the first, second, and third colors fluoresced at 2, 3, and
% 4 second respectively. The fourth color never fluoresced.
% 
% -settings: Structure. Summarizes all of the settings for the
% function. The fields are:
%     -traceColors: RGB matrix. Each row has RGB specifiers for how each color should be plotted
%     -resumeFlag: 0 or 1. If 0 function starts without loading any previous state assignmets. If 1 function will load previous assignments stored in resumeFile
%     -numClasses: Integer. Specifies number of different assignment classes
%     -numColors: Integer. Specifies number of different colors per trace
%     -classColors: RGB matrix. Each row has RGB specifiers for how each class assignment should be plotted
%     -resumeFile: String. Name of file where temporary assignments are made when function is stopped before all traces are analyzed. 
%
% OUTPUTS:
%
% -newData: Array. Version of data input which reflects any changes to background levels made by user. 
% 
% -eventData: Cell. For every assginment class, sumarizes every assignment made. Includes, among other things, 
% what  events preceed and follow it, event start and end times, as well as fluorescence values of all channels during the event.  
%
% -eventSequences: Cell. For every assignment class summarizes the event sequence of each trace.             
%
% INSTRUCTIONS:
%
% 1) press 'q' to set COLOR used for assignments
% 2) press 'w' to set assignemnt CLASS
% 3) press 'e' to set STATE being assigned
% ASSIGNMENT: click mouse in three places: right of the peak, left of the peak, and at the signal cutoff (so that enything above is considered signal)
% BACKGROUND SETTING:
%       All traces: press 'b', and select area to become background  all
%       traces will shift so that region is centered on 0 intensity.
%       trace-by-trace: press 'b', and select area to become background of
%       whatever color is selected. That color's trace will shift.
% MOVING: Move from one trace to the next by pressing ',' or '.'
% ZOOMING: Press 'z' and define the zoom boundaries with mouse, zoom-out by pressing 'u'
% STOPPING TO RESUME LATER: Press 's'. Program will stop and save a
% temporary file that can be used to restart assignments. To restart set
% resumeFlag value in settings to 1, and make sure the resumeFile is in the
% working directory

%% MAIN FUNCTION

%%% INITIALIZATION BLOCK

% skip initialization if restarting from previous session
if settings.resumeFlag == 1
    
    load(settings.resumeFile);

else
    
    % make sure settins file makes sense
    settingsChar(settings)
    
    % characterize data file
    [t ttotal framerate nFrames traces] = dataChar(data, settings);
    
    % generate temporary assignment structures:
    [HMMseqP STATES sMean stateNum] = stateMaker(settings.numClasses, nFrames, traces);
    
    % split data into easily plotted rank 3 array, separated by color
    colorTraces = colorTracesMaker(ttotal, nFrames, traces, settings.numColors);
    
    % to plot state assignments generate a rank 3 array, separated by
    % classes
    HMM = zeros(nFrames, traces, settings.numClasses);
    
    % variables initialized for plotting
    zooming = 0;
    zoomRegion = [min(t) max(t)];
    n=1;
    colorPick = 1;
    newState = 1;
    classPick = 1;
    
    % These variables generate the color strings to help identify each
    % color and class colors on the plot
    rainbowTitle1 = rainbowTitleMaker('Color ', settings.numColors, settings.traceColors);
    rainbowTitle2 = rainbowTitleMaker('Class ', settings.numClasses, settings.classColors);
    
end

%%% FIGURE PLOTTING LOOP

while n <= traces
    %%% PLOT BLOCK
    
    % generate subaxes and plotted traces
    ax(1) = subplot(2,1,1);
    ax(2) = subplot(2,1,2);
    currentTraces =  reshape(colorTraces(:,n,:),[],settings.numColors);
    currentHMM = reshape(HMM(:,n,:),[],settings.numClasses);
    
    % plot 2 windows and set their handles
    ph1 = plot(ax(1), t, currentTraces);
    ph2 = plot(ax(2), t, currentHMM);
    
    % changes x-axis values based on whether the zoom command is active
    if zooming
        if zoomRegion(1) > zoomRegion(2)
            xlim([zoomRegion(2) zoomRegion(1)])
        else
            xlim([zoomRegion(1) zoomRegion(2)])
        end
    else
        set(ax(:), 'xlim', [min(t) max(t)])
    end
    
    % set the title, axes, and trace colors for the top window
    set(ax(1),'Xcolor',[1 1 1],'YColor',[1 1 1],'Color',[0 0 0]);
    title(ax(1),['Trace ' num2str(n), ' ---> ' num2str(traces - n) ' remaining;     Color used to assign (q): ' num2str(colorPick) ...
        10 rainbowTitle1], ...
        'FontSize', 24)
    ylabel(ax(1),'Fluorescence intensity (AU)', 'FontSize', 20)
    set(ax(1),'FontSize', 16)
    
    for i=1:settings.numColors
        set(ph1(i),'Color',settings.traceColors(i,:));
    end
    grid(ax(1),'on');
    
    % set the title, axes, and trace colors for the bottom window
    set(ax(2),'Xcolor',[1 1 1],'YColor',[1 1 1],'Color',[0 0 0]);
    title(ax(2),['Class of state being assigned (w): ' num2str(classPick) ';     State being assigned (e): ' num2str(newState) ...
        10 rainbowTitle2], ...
        'FontSize', 24)
    xlabel(ax(2),'time (s)', 'FontSize', 20)
    ylabel(ax(2),'State assignment', 'FontSize', 20)
    set(ax(2),'FontSize', 16)
    ylim(ax(2),[-0.05 1.1])
    
    grid(ax(2),'on');
    for i=1:settings.numClasses
        set(ph2(i),'Color',settings.classColors(i,:));
    end
    
    grid on
    linkaxes(ax,'x')
    
    %%% FIGURE MANIPULATION BLOCK.
    check1 = waitforbuttonpress;
    
    % assign states when mouse is clicked
    if check1 == 0
        
        [setLevelX,setLevelY] = ginput(3);
        assignmentTrace = currentTraces(:,colorPick);
        changedIndices = assignmentTrace >= setLevelY(3);
        
        % determines assignment start and end points. Corrects them if the
        % outside of the plot is clicked
        start = round((setLevelX(1))/framerate);
        if start < 1
            start = 1;
        end
        
        finish = round((setLevelX(2))/framerate);
        if finish > length(t)
            finish = length(t);
        end
        
        onMask = zeros(length(t),1);
        onMask(start:finish) = 1;
        changedIndices = changedIndices & onMask == 1;
        
        % put the state assignment in the correct class array:
        HMMseqP{classPick}(changedIndices,n) = sMean{classPick}(newState);
        currentHMM(changedIndices,classPick) = sMean{classPick}(newState);
        HMM(:,n,:) = currentHMM;
        
        refresh(gcf);
        
    % set the COLOR to be used for assignments when 'q' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'q') == 1
        
        check2 = waitforbuttonpress;
        colorPick = str2double(get(gcf,'CurrentCharacter'));
        
        while isempty(colorPick) || isnan(colorPick) ||colorPick == 0 || colorPick > settings.numColors;
            beep
            check2 = waitforbuttonpress;
            colorPick = str2double(get(gcf,'CurrentCharacter'));
        end
        
    % set the ASSIGNMENT CLASS when 'w' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'w') == 1
        
        check3 = waitforbuttonpress;
        classPick = str2double(get(gcf,'CurrentCharacter'));
        
        while isempty(classPick) || isnan(classPick) || classPick == 0 || classPick > settings.numClasses;
            beep
            check3 = waitforbuttonpress;
            classPick = str2double(get(gcf,'CurrentCharacter'));
        end
        
    % set the STATE to be assigned when 'e' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'e') == 1
        
        check4 = waitforbuttonpress;
        newState = str2double(get(gcf,'CurrentCharacter'));
        while isempty(newState) || isnan(newState) ||newState == 0 || newState > length(sMean{classPick});
            beep
            check4 = waitforbuttonpress;
            newState = str2double(get(gcf,'CurrentCharacter'));
        end
        
    % advance to next molecule if '.' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'.') == 1
        
        n = n+1;
        zooming = 0;
        
   % return to previous molecule if ',' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),',') == 1 && n>1
        
        n = n-1;
        zooming = 0;
        
    % set up zoom axes and commands plot block above to
    % create zoom version when 'z' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'z')
        
        zoomRegion = ginput(2);
        zooming = 1;
        
    % return figure to unzoomed axes when 'u' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'u')
        
        zooming = 0;
        
    % set all assignments in current assignment class to 1 when 'c' is typed
    elseif strcmp(get(gcf,'CurrentCharacter'),'c')
        
        HMMseqP{classPick}(:,n) = sMean{classPick}(1);
        currentHMM(:,classPick) = sMean{classPick}(1);
        HMM(:,n,:) = currentHMM;
        
    % sets the area selected color as background (intensity 0) for all traces when 'b' is typed and area selected
    elseif strcmp(get(gcf,'CurrentCharacter'),'b')
        backRegion = ginput(2);
        backRegion(1) = round(backRegion(1)/framerate);
        backRegion(2) = round(backRegion(2)/framerate);

        colorTraces = totalBacker(backRegion, currentTraces, colorTraces, NaN, n, t);
        
            % sets the area selected color as background (intensity 0) for all traces when 'b' is typed and area selected
    elseif strcmp(get(gcf,'CurrentCharacter'),'v')
        backRegion = ginput(2);
        backRegion(1) = round(backRegion(1)/framerate);
        backRegion(2) = round(backRegion(2)/framerate);
         
        colorTraces = totalBacker(backRegion, currentTraces, colorTraces, colorPick, n, t);

        
    % create temporary file so that program can stop and resume later.
    elseif strcmp(get(gcf,'CurrentCharacter'),'s')
        save(settings.resumeFile);
        break
        
    else
        beep
    end
end

%%% SUMMARY GENERATION BLOCK

% transfer assignments to STATE variable for downstream analysis
for class = 1:settings.numClasses
    for state = 1:stateNum
        STATES{class}{state} = HMMseqP{class}==sMean{class}(state);
    end
end

% generate new_data variable that takes into account any changes
% that might have been made (i.e. background levels)
newData = [];
for  i= 1:traces
    newData = [newData reshape(colorTraces(:,i,:),[],settings.numColors)];
end
newData = [t newData];

% generate the event summary structures
[eventData eventSequences] = singleMolSummarizer(data, settings, STATES);

% save a version of all available variables in case troubleshooting is
% necessary after assignments are made
save stateAssignComplete

close all;

end


%% ACCESSORY FUNCTIONS

function settingsChar(settings)
% Makes sure fields in settings variable make sense

assert(settings.numColors <= size(settings.traceColors,1), 'traceColor array is too small')
assert(settings.numClasses <= size(settings.classColors,1), 'classColors array is too small')
if settings.resumeFlag == 1
    assert(~ismepty(settings.resumeFile), 'resumeFile cannot be empty is resumeFlag is set to 1')
end
end

function [t ttotal framerate nFrames traces] = dataChar(data, settings)
% Basic characterization of data variable

ttotal = data;
t = ttotal(:,1);
ttotal(:,1) = [];
framerate = t(1);
[nFrames, nCol] = size(ttotal);
traces = nCol/settings.numColors;

assert(round(std(diff(t))) <=1e-3, 'Time column in data is not equally spaced')
assert(mod(nCol,settings.numColors) == 0, 'numColor setting is wrong')
end

function [HMMseqP STATES sMean stateNum] = stateMaker(numClasses, nFrames, traces)
% Generates temporary state assignment structures
%
% INPUTS:
%
% -numClasses: Integer. Specifies number of different assignment classes
% -nFrames: Integer. Number of frames in data for which the structures are
% being made
% -traces: Number of traces that will be assigned
%
% OUTPUTS:
%
% -HMMseqP: Cell. Keeps track of the state assingment value for each trace
% in each class.
% -STATES: Cell. Keeps track of the state assignment index for each trace
% in each class.
% -sMean: Cell. Determines values graphed for each state within assignement
% class.
% -stateNum: Integer. Number of states to be assigned per class (9 is max
% and default, one state for each digit with 1 as NO STATE).

stateNum = 9;
STATES = cell(1,numClasses);
HMMseqP = cell(1, numClasses);
sMean = cell(1, numClasses);

% A basic value array for state assignments, assumed to be the same for all assingment classes
basicStateMean = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

% populate variables from above
for class = 1:numClasses
    
    HMMseqP{class} = zeros (nFrames, traces);
    sMean{class} = basicStateMean;
    
    for state = 1:stateNum 
        STATES{class}{state}= zeros(nFrames, traces);
    end
end
end

function colorTraces = colorTracesMaker(ttotal, nFrames, traces, numColors)
% Generate an easily manipulated rank 3 array, separated by color.
%
% INPUTS:
%
% -ttotal: Array. Derived from dataChar function. Standard data format
% except with t vector prefix removed.
% -nFrames: Integer. Number of frames per frame.
% -traces: Integer. Number of traces in data.
% -numColors: Integer. Number of colors per trace. 
%
% OUTPUT:
%
% -colorTraces: Array. Rank 3 array, separated by color of all traces.

% initalization variables
colorTraces = zeros(nFrames, traces, numColors);
nTraces = 1:traces;

% loop through ttotal to make the colorTraces array
for i=1:numColors
    colorTraces(:,:,i) = ttotal(:,nTraces * numColors - numColors + i);
end
end

function rainbowTitle = rainbowTitleMaker(string, numColors, RGBcolor)
% Generates color-coordinated title strings, based on RGB color arrays
% provided. (This is necessary because matlab bug doesn't allow for legends
% in ginput function is used.)
%
% INPUTS:
%
% -string: String. Specifies the string to be repeated and colored.
% -numColors: Integer. Number of times the string should be repeated in
% different colors.
% -RGBcolor: RGB array. Each row is an RGB specifier for what color each
% repeated string should be. 
%
% OUTPUT:
%
% -rainbowTitle: String. Formated so that when given to a matlab plot text
% command it will be interpreted as a multi-colored string

rainbowTitle = '';
% Add color suffixes to string
namesColor = nameMaker(string, numColors);

% Generate string
for n = 1:length(namesColor)
    rainbowTitle = [rainbowTitle, '\color[rgb]{', ...
        num2str(RGBcolor(n,1)), ' ', num2str(RGBcolor(n,2)), ' ', num2str(RGBcolor(n,3)), ...
        '}', namesColor{n}, '      '];
end
end

function nameCell = nameMaker(string, num)
% Generates a name cell wiht string prefixing a number iterator. Useful for rainbowTitleMaker
% function
%
% INPUTS:
%
% -string: String. Determines the string that precedes each number
% -num: number of iterations
%
% OUTPUT:
% 
% -nameCell: Cell. Each location holds unique iterated string.

nameCell = cell(1, num);

for n = 1:num
    nameCell{n} = [string num2str(n)];
end
end

function colorTraces = totalBacker(backRegion, currentTraces, colorTraces, colorPick, n, t)
% Makes selected area the background level of all traces. 
% 
% INPUTS: 
%
% -backRegion: Array. Sets limits for area to be used as background. 
% -currentTraces: Array. Traces to be modified.
% -colorTraces: Array: Source of traces to be modified, and ultimately
% where the changes need to be stored. 
% -colorPick: Integer or NaN. If set to NaN, the function will apply changes to all traces in currentTraces, 
% otherwise changes will only be applied to the current colorPick trace.
% -n: Reflects trace being analyzed in plotting loop 
% -t: time vector of data. 
%
% OUTPUT:
%
% -colorTraces: Array. Original modified to reflect changes made in
% background.

        % make sure backregion limits make sense
        if backRegion(1) > backRegion(2)
            
            if backRegion(1) > length(t)
                backRegion(1) = length(t);
            elseif backRegion(2) < 1
                backRegion(2) = 1;
            end
            
            if isnan(colorPick)
                back = currentTraces(backRegion(2):backRegion(1), :);
            else
                back = currentTraces(backRegion(2):backRegion(1), colorPick);
            end
            
        else
            
            if backRegion(2) > length(t)
                backRegion(2) = length(t);
            elseif backRegion(1) < 1
                backRegion(1) = 1;
            end
            
            if isnan(colorPick)
                back = currentTraces(backRegion(1):backRegion(2), :);
            else
                back = currentTraces(backRegion(1):backRegion(2), colorPick);
            end
        end
        
        % substracts background from traces
        backMean = mean(back);
        
        if isnan(colorPick)
            for color = 1: length(backMean)
                currentTraces(:,color) = currentTraces(:,color) - backMean (color);
            end
        else
           currentTraces(:,colorPick) = currentTraces(:,colorPick) - backMean; 
        end
        
        colorTraces(:,n,:) = currentTraces;
        
end

function [eventData eventSequences] = singleMolSummarizer(data, settings, STATES)
% Generates a summary of the state assigments created by stateAssign
% function.
%
% INPUTS:
% -data: Array. Same as used for stateAssign
% -settings: Structure. Same as used for state Assign.
% - STATES: Cell. Keeps track of the state assignment index for each trace
% in each class.
%
% OUTPUTS: 
%
%-eventData: Cell. For every assginment class, sumarizes every assignment made. Includes, among other things, 
% what  events preceed and follow it, event start and end times, as well as fluorescence values of all channels during the event.  
%-eventSequences: Cell. For every assignment class summarizes the event sequence of each trace.   

% characterize data
[t ttotal framerate nFrames traces] = dataChar(data, settings);

% generate rank-3 array to extract the trace information in the event
% summaries below
colorTraces = colorTracesMaker(ttotal, nFrames, traces, settings.numColors);

% generate names for color fields in eventData structure
namesColor = nameMaker('color_', settings.numColors);

%%% loop to generate eventData structure and eventSequences cell

for class = 1:settings.numClasses
    
    for trace=1:traces
        
        % construct a vector containing STATES for the trace
        traceStates = zeros(nFrames,1);
        
        for state = 1:size(STATES{class},2)
            
            traceStates = traceStates + state*STATES{class}{state}(:,trace);
        
        end
        
        currentTraces =  reshape(colorTraces(:,trace,:),[],settings.numColors);
        
        %characterize states in trace 
        transitions = diff(traceStates);
        transitions = find(transitions ~= 0);
        eventsInMol = [traceStates(1) traceStates(transitions+1)'];
        
        for event = 1:length(eventsInMol)
            
            eventData{class}{event,trace}= recordEvent(event, transitions, eventsInMol, currentTraces ,t, namesColor);
        
        end
        
        eventSequences{class}{trace} = eventsInMol;
        
    end
end
end

function record = recordEvent(event, transitions, eventsInMol, currentTraces ,t, namesColor)
% Records event analysis variables into record structure.
%
% INPUTS: All inputs derived from singleMolSummarizer function, except for:
%
% -t: Vector: times at which each frame happened.
% - namesColor: 
%
% OUTPUT:
%
% -record. Structure. Records event characteristics.

record.state = eventsInMol(event);

% determine event start time
if (event == 1)
    record.stateBefore = 100;
    start = 1;
    record.start = t(start);
else
    record.stateBefore = eventsInMol(event-1);
    start = transitions(event-1)+1;
    record.start = t(start);
end

% determine event finish time
if (event == length(eventsInMol))
    record.stateAfter = 100;
    finish = length(currentTraces);
    record.finish = t(finish);
else
    record.stateAfter = eventsInMol(event+1);
    finish = transitions(event);
    record.finish = t(finish);
end

% determine additional parameters based on event start and finish times
record.lifetime = t(1)*(finish - start + 1);
record.numEvents = length(find(eventsInMol==record.state));
record.position = length(find(eventsInMol(1:event)==record.state));
record.relativePosition = event;
record.beforeSequence = eventsInMol(1:event);
record.afterSequence = eventsInMol(event:length(eventsInMol));

% saves traces during the event
for color = 1:length(namesColor)
    record = setfield(record, namesColor{color}, currentTraces(start:finish,color));
end
end
