function events = eventDataExtractor(eventData, targetClass, targetState, parameter)
% Extracts particular fields from the eventData summary structure created
% by the assignState function.
%
% INPUTS:
%
% -eventData: Cell of structures. Summary structure created by stateAssign
% function. Contains data to be extracted.
% -targetClass: Integer. Class of assignment to be extracted. 
% -targetState: Integer. State of assignment to be extracted.
% -parameter: String. Fieldname specifying which property of the assigned
% event to be extracter. For example 'start', 'finish', 'color_1' for start
% times, end times or trace values of color 1 for the event assingments
% specified
%
% OUTPUT:
%
% events: Array or Cell. Array or cell organized by event number in rows
% and trace number in columns. Contains extracted values. Will be a cell if
% values are not single numbers. Will be a, more easily manipulated, array otherwise. 

% characterize eventData and initialize output
assert(isfield(eventData{targetClass}{1}, parameter), 'parameter string does not specify a fieldname in eventData')
sizeEventData = size(eventData{targetClass});
events = cell(sizeEventData);
sizes = zeros(sizeEventData);

% Loop fills in events cell
for trace=1:sizeEventData(2)
    traceEventIndex = 1;
    for event = 1:sizeEventData(1)
        if ~isempty(eventData{targetClass}{event, trace}) &&  eventData{targetClass}{event, trace}.state == targetState
            events{traceEventIndex, trace} = getfield(eventData{targetClass}{event, trace}, parameter);
            sizes(traceEventIndex, trace) = length(events{traceEventIndex, trace});
            traceEventIndex = traceEventIndex + 1;
        end
    end
end

% Determines if events cell contains single numbers
sizeTest = sizes > 1;
sizeTest = sum(sum(sizeTest));

% if events has only single numbers converts is to a more easily
% manipulated array. 
if sizeTest == 0
    eventsTempo = zeros(size(events));
    for trace = 1:size(events,2)
        for event = 1:size(events,1)
            if isempty(events{event, trace})
                eventsTempo(event, trace) = 0;
            else
                eventsTempo(event, trace) = events{event, trace};
            end
        end
    end
    events = eventsTempo;
end

end