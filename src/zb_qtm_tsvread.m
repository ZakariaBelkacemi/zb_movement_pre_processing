
function[tsvdata, tsvheader] = zb_qtm_tsvread(tsvfilename,varargin)

header_length = 12;
ignore_data = 0;

for i = 1:2:numel(varargin)
        switch(varargin{i})
            case 'header_length'
                header_length = varargin{i+1};
            case 'ignore_data'
                ignore_data = varargin{i+1};
            otherwise
        end
end

% Open the TSV file for reading
% tsvfilename = 'C:\Users\ARNGDC\Desktop\tsv.tsv';
fileID = fopen(tsvfilename, 'r');

% Initialize the variable to store the data
tsvheader = {};

% Ici nous considérons que le header correspond aux 12 premières lignes
i = header_length;
while ~feof(fileID) && i > 0
    line = fgetl(fileID); % Read one line at a time
    %         disp(line)
    data = strsplit(line, '\t'); % Split the line by tabs
    % Process the data if needed (here, we store the data as a cell array)
    data_str_array = [];
    for ii = 1:numel(data)
        data_str_array = [data_str_array, convertCharsToStrings(data{ii})];
    end
    tsvheader{end+1} = {data_str_array}; % Append the data to the variable
    
    i = i - 1;
end

tsvheader = tsvheader';

tsvdata = [];

if ignore_data == 0

    while ~feof(fileID)
        line = fgetl(fileID); % Read one line at a time
        data = strsplit((line), '\t'); % Split the line by tabs
        % Process the data if needed (here, we store the data as a cell array)
        data_array = [];
        for ii = 1:numel(data)
            data_array(end+1) = str2num(data{ii});
        end
        tsvdata(end+1,:) = data_array; % Append the data to the variable
    end

    % Close the file
    fclose(fileID);
end

end