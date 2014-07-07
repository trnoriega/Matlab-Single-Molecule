function newData = fretSelector(data, settings)
% Selects of traces from a data file that are likely to show
% anti-correlated fluorescence resonace energy transfer (FRET) behavior.
%
% Tested in Matlab R2012a
%
% INPUTS:
%
% -data: Array. where each row is an observation frame and each column is fluorescence intensity from each fluorophore in that frame.
% First column is time of each frame (assumed ot be in seconds with each frame representing the same ammount of time).
% Following columns are groups of fluorophores grouped by columns. One
% column should represent a FRET DONOR signal, another column should
% specify FRET ACCEPTOR SIGNAL. Which is which is specified in settings.
%
% For example:
%
% data =
% 1 1 0 0 0
% 2 1 0 0 0
% 3 0 1 0 0
% 4 1 0 1 0
%
% is a single trace, taken over four seconds at a 1 second framerate, with
% four colors. the first is the FRET DONOR, the second the FRET ACCEPTOR, the third  and fourth are additional channels
% The FRET event happened at second 3. The third color fluoresced at second 4. The fourth color never fluoresced.
%
% -settings: Structure. Summarizes all of the settings for the
% function. The fields are:
%     -numColors: Integer. Specifies number of different colors per trace
%     -donorSD, acceptorSD, fretSD: Float. Spefies number of standard deviations above
%     background that should be considered as a donor, acceptor, or FRET event. Determined
%     empirically based on experimental conditions.
%     -fretPair: 1X2 vector. number positions of the donor and acceptor
%     signals in the data file. In the example above the fretPair vector
%     would be: [1 2]
%     -smoothFlag: 1 or 0. Indicates whether user wants data to be
%     smoothed using a discrete transfer function. 1 indicates yes, 0
%     indicates no. 
%     -smoothFrames: Integer from 1-5. If smoothFlag = 1, number of frames overwhich the 
%     smoothing should be done.
%
% OUTPUT:
%
% -newData: Array. Version of data with only the traces that only includes the traces that might have anti-correlated FRET events.

% characterize data file
ttotal = data;
t = ttotal(:,1);
ttotal(:,1) = [];
nCol = size(ttotal,2);
traces = nCol/settings.numColors;
assert(round(std(diff(t))) <=1e-3, 'Time column in data is not equally spaced')
assert(mod(nCol,settings.numColors) == 0, 'numColor setting is wrong')

% Smooths the data if smoothFlag is set to 1
if settings.smoothFlag == 1
    
    ttotal = double(ttotal);
    
    for n = 1:size(ttotal,2)
        ttotal(:,n) = filtfilt(ones(1,settings.smoothFrames),settings.smoothFrames,ttotal(:,n));
    end
    
end


% extract FRET, donor, and acceptor traces
nTraces = 1:traces;
donor = ttotal(:,(nTraces * settings.numColors - (settings.numColors-settings.fretPair(1))));
acceptor = ttotal(:,(nTraces * settings.numColors -(settings.numColors-settings.fretPair(2))));
fret = acceptor./(donor+acceptor);

% make donor, acceptor and FRET trace-by-trace frame-by-frame variation
% matrices
differenceDonor = diff(donor);
differenceAcceptor = diff(acceptor);
differenceFret = diff(fret);

% determine summary statistics
meanDonor = mean(differenceDonor);
stdDonor = std(differenceDonor);
meanAcceptor = mean(differenceAcceptor);
stdAcceptor = std(differenceAcceptor);
meanFret = mean(differenceFret);
stdFret = std(differenceFret);

% determine cutoffs based on settings file parameters
cutDonorUp = meanDonor + settings.donorSD * stdDonor;
cutDonorDown = meanDonor - settings.donorSD * stdDonor;
cutAcceptorUp = meanAcceptor + settings.acceptorSD * stdAcceptor;
cutAcceptorDown = meanAcceptor - settings.acceptorSD * stdAcceptor;
cutFretUp = meanFret + settings.fretSD * stdFret;
cutFretDown = meanFret - settings.fretSD * stdFret;

% find all points in traces where all three cutoffs (FRET, donor, and
% acceptor) are met in rising (up) FRET events. Point-by-point comparisons
% is important because it assures temporal anti-correlation.

% initialize summary with tally of rising (up) FRET events
testUp = zeros(1,traces);

for i = nTraces
    
    indexFretUp = find(differenceFret(:,i) > cutFretUp(i));
    indexDonorDown = find(differenceDonor(:,i) < cutDonorDown(i));
    indexAcceptorUp = find(differenceAcceptor(:,i) > cutAcceptorUp(i));
    finalOneUp = indexFretUp(ismember(indexFretUp,indexDonorDown));
    finalTwoUp = finalOneUp(ismember(finalOneUp,indexAcceptorUp));
    if ~isempty(finalTwoUp)
        testUp(i) = 1;
    end
    
end

% find all points in traces where all three cutoffs (FRET, donor, and
% acceptor) are met in falling (down) FRET events. Point-by-point comparisons
% is important because it assures temporal anti-correlation.

% initialize summary with tally of falling (down) FRET events
testDown = zeros(1,traces);

for i = nTraces
    
    indexFretDown = find(differenceFret(:,i) < cutFretDown(i));
    indexDonorUp = find(differenceDonor(:,i) > cutDonorUp(i));
    indexAcceptorDown = find(differenceAcceptor(:,i) < cutAcceptorDown(i));
    finalOneDown = indexFretDown(ismember(indexFretDown,indexDonorUp));
    finalTwoDown = finalOneDown(ismember(finalOneDown,indexAcceptorDown));
    if ~isempty(finalTwoDown)
        testDown(i) = 1;
    end
    
end

% make a new data array with only traces that have anti-correlated donor
% and acceptor signal (either up or down)
tracesFinal = testUp + testDown;
tracesFinal = tracesFinal > 0;
indexFinal = zeros(1, traces*settings.numColors);

for i = 0:settings.numColors-1
   indexFinal(:, nTraces*settings.numColors - i) = tracesFinal;
end

newData = ttotal(:, logical(indexFinal));
newData = [t newData];

end