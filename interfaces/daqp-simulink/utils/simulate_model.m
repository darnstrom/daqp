function data = simulate_model(modelName,params)

assert(isstruct(params),'params needs to be a struct')

numExperiments = length(params);
paramNames = fieldnames(params);
lenParams = length(paramNames);

input(numExperiments) = Simulink.SimulationInput(modelName);
if numExperiments > 1
    input(1:numExperiments-1) = input(numExperiments);
end
for idxExperiment = 1:numExperiments
    for idxParam = 1:lenParams
        paramName = paramNames{idxParam};
        paramValue = params(idxExperiment).(paramName);
        input(idxExperiment) = input(idxExperiment).setVariable(paramName,paramValue);
    end
    input(idxExperiment) = input(idxExperiment).setModelParameter('SaveFormat', 'Structure');
end


%% Change the save format
% open_system(modelName);
% mode = get_param(modelName,'SaveFormat');
% set_param(modelName, 'SaveFormat', 'Structure');


%% Run the simulation
data = struct();
% try
    % Attempt to run the simulation for each experiment.
    % This block is used to catch any errors that may occur during the simulation process,
    % such as issues with model parameters or simulation settings.
    % tic
    for idxExperiment = 1:numExperiments
        simout = sim(input(idxExperiment));
        
        if idxExperiment == 1
            data = simout2data(simout);
        else
            data(idxExperiment) = simout2data(simout);
        end
    end
% catch x
%     warning(x.message)
%     if ~exist('data','var')
%         data = struct();
%     end
% end


%% Reset the save format
% set_param(modelName, 'SaveFormat', mode);
% close_system(modelName);

end

function data = simout2data(simout)
    youtLength = length(simout.yout.signals);
    youtSignals = cell(1,youtLength);
    timeLen = length(simout.tout);
    
    % Get the signal names from the yout signals
    for i = 1:youtLength
        if ~isempty(simout.yout.signals(i).label)
            youtSignals{i} = simout.yout.signals(i).label;
        else
            signalParts = strsplit(simout.yout.signals(i).blockName,'/');
            youtSignals{i} = signalParts{end};
        end
    end
    
    % Get the signal names from the logsout signals
    if isempty(simout.logsout)
        logsOutLength = 0;
        logsOutSignals = cell(1,0);
    else
        logsOutLength = simout.logsout.numElements;
        logsOutSignals = cell(1,logsOutLength);
        for i = 1:logsOutLength
            logsOutSignals{i} = simout.logsout{i}.Name;
        end
    end
    
    % Create the data structure
    % The first row contains the signal names, the second row contains the signal values
    signalLength = 1 + youtLength + logsOutLength;
    dataCell = cell(2,signalLength);
    dataCell{1} = simout.tout;
    for idxYOut = 1:youtLength
        dataCell{1,idxYOut} = youtSignals{idxYOut};
        dataCell{2,idxYOut} = permuteSignalMatrix(simout.yout.signals(idxYOut).values,timeLen);
    end

    % Add the signals from the logsout
    idxSignal = idxYOut + 1;
    for idxLogsOut = 1:logsOutLength
        % Check if the signal is already in the data structure
        if ismember(logsOutSignals{idxLogsOut},youtSignals)
            warning('Duplicate signal "%s" in the simulation output',logsOutSignals{idxLogsOut})
            continue
        end
        % Signal name
        dataCell{1,idxSignal} = logsOutSignals{idxLogsOut};
        
        % Signal values
        dataCell{2,idxSignal} = permuteSignalMatrix(simout.logsout{idxLogsOut}.Values.Data,timeLen);
        idxSignal = idxSignal + 1;
    end
    
    % Convert the cell array to a struct
    data = struct(dataCell{:,1:idxSignal-1});
end
    
function mat = permuteSignalMatrix(dataMat, timeLen)
    % Permute the signal matrix so that the time dimension is the last dimension
    % 
    % 2D-input:
    % - the time dimension will always be the last dimension
    %
    %  3D-input: 
    % - if one of the data dimensions is 1, then it will result in a 2D matrix 
    %       with the time dimension as the last dimension
    % - if none of the data dimensions is 1, then the time dimension will be 
    %       the last dimension

    dims = size(dataMat);
    idxTime = find(dims == timeLen,1);
    notTime = find(dims ~= timeLen);

    if isempty(idxTime)
        error('The time dimension is not found in the signal matrix')
    end

    if all(dims == 1)
        mat = dataMat;
    elseif length(notTime) == 1
        mat = permute(dataMat,[notTime,idxTime]);
    elseif all(dims(notTime) == 1)
        mat = permute(dataMat,[idxTime,notTime])';
    elseif any(dims(notTime) == 1)
        notTimeAndNotOne = find(dims ~= timeLen & dims ~= 1);
        notTimeAndOne = find(dims ~= timeLen & dims == 1);
        mat = permute(dataMat,[notTimeAndNotOne,idxTime,notTimeAndOne]);
    else
        mat = permute(dataMat,[notTime,idxTime]);
    end
end