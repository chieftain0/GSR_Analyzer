GSR_fear_file = 'GSR_FEAR.csv';
GSR_baseline_file = 'GSR_Baseline.csv';
fear_matrix = read_gsr_data(GSR_fear_file);
baseline_matrix = read_gsr_data(GSR_baseline_file);
disp('GSR data read successfully');

%Determine Sampling Frequency
dt = seconds(mean(diff(fear_matrix.TimeStamp)));
Fs = round(1/dt);

subplot(6, 2, 1);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Original GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 2);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Original GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

fear_matrix.GSRValue = medfilt1(fear_matrix.GSRValue);
baseline_matrix.GSRValue = medfilt1(baseline_matrix.GSRValue);
disp('Median filter applied successfully');

subplot(6, 2, 3);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Median filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 4);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Median filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;


fear_matrix.GSRValue = lowpass(fear_matrix.GSRValue, 20, Fs);
baseline_matrix.GSRValue = lowpass(baseline_matrix.GSRValue, 20, Fs);
disp('Low pass filter applied successfully');

subplot(6, 2, 5);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Low Pass filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 6);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Low Pass filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

fear_matrix.GSRValue = ten_point_average_filter(fear_matrix.GSRValue);
baseline_matrix.GSRValue = ten_point_average_filter(baseline_matrix.GSRValue);
disp('10 point average filter applied successfully');

subplot(6, 2, 7);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('10 point average filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 8);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('10 point average filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

%Due to moving average filter, we need to cut out a few seconds of data
fear_matrix = cut(fear_matrix, 2);
baseline_matrix = cut(baseline_matrix, 2);
disp('Cuts applied successfully');


[fear_maxima, fear_minima] = find_local_extrema(fear_matrix);
[baseline_maxima, baseline_minima] = find_local_extrema(baseline_matrix);

% Plot Fear data with local maxima and minima
subplot(6, 2, 9);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
hold on;
scatter(fear_matrix.TimeStamp(fear_maxima), fear_matrix.GSRValue(fear_maxima), 'r', 'filled');
scatter(fear_matrix.TimeStamp(fear_minima), fear_matrix.GSRValue(fear_minima), 'g', 'filled');
hold off;
title('Processed GSR Values with Local Maxima and Minima (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

% Plot Baseline data with local maxima and minima
subplot(6, 2, 10);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
hold on;
scatter(baseline_matrix.TimeStamp(baseline_maxima), baseline_matrix.GSRValue(baseline_maxima), 'r', 'filled');
scatter(baseline_matrix.TimeStamp(baseline_minima), baseline_matrix.GSRValue(baseline_minima), 'g', 'filled');
hold off;
title('Processed GSR Values with Local Maxima and Minima (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

















function data_table = read_gsr_data(file_name)
% Reads GSR data from a CSV file and returns a table
% Input:
%   - file_name: The name of the CSV file
%   - time_format: The format of the time stamp in the file

% Read data from the file
timeFormat = 'mm:ss:SSS';
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["TimeStamp", "GSRValue"];
opts.VariableTypes = ["string", "double"];

% Combine file path and name
file_path = fullfile(pwd, file_name);

% Check if the file exists
if exist(file_path, 'file') ~= 2
    error('File not found: %s', file_path);
end

% Read the data
data_table = readtable(file_path, opts);

% Convert the TimeStamp column to datetime
data_table.TimeStamp = datetime(data_table.TimeStamp, 'Format', timeFormat);
end

function filtered_array = ten_point_average_filter(input_array)
b = ones(1, 10) * 1/10;
a = 1;
filtered_array = filter(b, a, input_array);
end


function [fear_maxima, fear_minima] = find_local_extrema(fear_matrix)
minProminence = 1.5;
fear_maxima = islocalmax(fear_matrix.GSRValue, 'MinProminence', minProminence);
fear_minima = islocalmin(fear_matrix.GSRValue, 'MinProminence', minProminence);

end

function [cut_fear_matrix] = cut(fear_matrix, time)
% Determine the duration of the first two seconds
duration_to_cut = seconds(time);

% Find the indices of the data points to keep (after the first two seconds)
indices_to_keep = fear_matrix.TimeStamp >= fear_matrix.TimeStamp(1) + duration_to_cut;

% Apply the indices to both fear_matrix and baseline_matrix
cut_fear_matrix = fear_matrix(indices_to_keep, :);
end

