folderName = 'SimDataCSVs';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

for i = 1:length(results)
    columnNames = cell(1, length(results(i).DataInfo));
    for j = 1:length(results(i).DataInfo)
        if isfield(results(i).DataInfo{j}, 'Compartment')
            compartment = results(i).DataInfo{j}.Compartment;
        else
            compartment = 'UnknownCompartment';
        end
        type = results(i).DataInfo{j}.Type;
        name = results(i).DataInfo{j}.Name;
        columnNames{j} = matlab.lang.makeUniqueStrings(strjoin({compartment, type, name}, '_'));
    end
    data = results(i).Data;
    % Extract the time vector for the current SimData object
    timeVector = results(i).Time; 
    % Check if time is already in minutes; otherwise convert it
    
    % Ensure dataWithTime has the same number of rows as data
    if length(timeVector) == size(data, 1)
        dataWithTime = [timeVector, data];
        columnNames = ['Time', columnNames];
        T = array2table(dataWithTime, 'VariableNames', columnNames);
        filename = fullfile(folderName, sprintf('SimDataResults_%d.csv', i));
        writetable(T, filename);
    else
        error('Time vector and data have different number of rows for result index %d.', i);
    end
end
