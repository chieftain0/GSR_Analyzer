GSR_fear_file = 'GSR_FEAR.csv';
GSR_baseline_file = 'GSR_Baseline.csv';
fear_matrix = read_gsr_data(GSR_fear_file);
baseline_matrix = read_gsr_data(GSR_baseline_file);
disp('GSR data read successfully');

%Determine Sampling Frequency
%I know it was given in the assignment definition
%Should've read it thoroughly before going into matlab
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

fear_matrix.GSRValue = third_order_median_filter(fear_matrix.GSRValue);
baseline_matrix.GSRValue = third_order_median_filter(baseline_matrix.GSRValue);
disp('Median filter applied successfully');
%or use these functions like a sane person
%fear_matrix.GSRValue = medfilt1(fear_matrix.GSRValue);
%baseline_matrix.GSRValue = medfilt1(baseline_matrix.GSRValue);

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
disp('Point Average filter applied successfully');

subplot(6, 2, 7);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Point Average filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 8);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Point Average filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

%Due to filters, we need to cut out a few seconds of data
fear_matrix = cut(fear_matrix, 2);
baseline_matrix = cut(baseline_matrix, 2);
disp('Cuts applied successfully');

fear_matrix.GSRValue = map_to_0_100(fear_matrix.GSRValue);
baseline_matrix.GSRValue = map_to_0_100(baseline_matrix.GSRValue);
disp('Data normalized successfully');

subplot(6, 2, 9);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Normalized GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 10);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Normalized GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;



[fear_maxima, fear_minima] = find_local_extrema(fear_matrix);
[baseline_maxima, baseline_minima] = find_local_extrema(baseline_matrix);

% Plot Fear data with local maxima and minima
subplot(6, 2, 11);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
hold on;
scatter(fear_matrix.TimeStamp(fear_maxima), fear_matrix.GSRValue(fear_maxima), 'g', 'filled');
scatter(fear_matrix.TimeStamp(fear_minima), fear_matrix.GSRValue(fear_minima), 'r', 'filled');
hold off;
title('Processed GSR Values with Local Maxima and Minima (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

% Plot Baseline data with local maxima and minima
subplot(6, 2, 12);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
hold on;
scatter(baseline_matrix.TimeStamp(baseline_maxima), baseline_matrix.GSRValue(baseline_maxima), 'g', 'filled');
scatter(baseline_matrix.TimeStamp(baseline_minima), baseline_matrix.GSRValue(baseline_minima), 'r', 'filled');
hold off;
title('Processed GSR Values with Local Maxima and Minima (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;


[F1_Fear, F2_Fear, F3_Fear, F4_Fear, F5_Fear, F6_Fear, F7_Fear, F8_Fear, F9_Fear, F10_Fear] = analyze(fear_matrix, fear_maxima, fear_minima);
[F1_Baseline, F2_Baseline, F3_Baseline, F4_Baseline, F5_Baseline, F6_Baseline, F7_Baseline, F8_Baseline, F9_Baseline, F10_Baseline] = analyze(baseline_matrix, baseline_maxima, baseline_minima);

Fear_Index_Fear = calculate_fear_index(F1_Fear, F2_Fear, F3_Fear, F4_Fear, F5_Fear, F6_Fear, F7_Fear, F8_Fear, F9_Fear, F10_Fear);
Fear_Index_Baseline = calculate_fear_index(F1_Baseline, F2_Baseline, F3_Baseline, F4_Baseline, F5_Baseline, F6_Baseline, F7_Baseline, F8_Baseline, F9_Baseline, F10_Baseline);
format longG;
disp('Fear Index (Fearfull): ');
disp(Fear_Index_Fear);
disp('Fear Index (Baseline): ');
disp(Fear_Index_Baseline);




function data_table = read_gsr_data(file_name)
% Define the time format for 'TimeStamp' column
timeFormat = 'mm:ss:SSS';

% Specify import options for reading the CSV file
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["TimeStamp", "GSRValue"];
opts.VariableTypes = ["string", "double"];

% Construct the full file path
file_path = fullfile(pwd, file_name);

% Check if the file exists
if exist(file_path, 'file') ~= 2
    error('File not found: %s', file_path);
end

% Read the CSV file into a table
data_table = readtable(file_path, opts);

% Convert 'TimeStamp' column to datetime format
data_table.TimeStamp = datetime(data_table.TimeStamp, 'Format', timeFormat);
end


function filtered_array = third_order_median_filter(input_array)
% Get the length of the input array
N = length(input_array);

% Initialize the filtered array
filtered_array = zeros(size(input_array));

% Apply the third-order median filter
for i = 2:N-1
    % Extract the neighborhood around the current element
    neighborhood = input_array(i-1:i+1);
    
    % Calculate the median of the neighborhood and assign it to the filtered array
    filtered_array(i) = median(neighborhood);
end
end

function filtered_array = ten_point_average_filter(input_array)
% Define the coefficients for a 10-point moving average filter
b = ones(1, 10) * 1/10;
a = 1;

% Apply the filter to the input array
filtered_array = filter(b, a, input_array);
end

function normalized_data = map_to_0_100(data)
% Find the minimum and maximum values in the input data
min_val = min(data);
max_val = max(data);

% Normalize the data to the range [0, 100]
normalized_data = 100 * (data - min_val) / (max_val - min_val);
end

function [fear_maxima, fear_minima] = find_local_extrema(fear_matrix)
% Set the minimum prominence for identifying extrema
minProminence = 3;

% Identify local maxima based on GSR values
fear_maxima = islocalmax(fear_matrix.GSRValue, 'MinProminence', minProminence);

% Identify local minima based on GSR values
fear_minima = islocalmin(fear_matrix.GSRValue, 'MinProminence', minProminence);
end

function [cut_fear_matrix] = cut(fear_matrix, time)
% Convert time to duration for comparison
duration_to_cut = seconds(time);

% Determine indices to keep based on the specified time
indices_to_keep = fear_matrix.TimeStamp >= fear_matrix.TimeStamp(1) + duration_to_cut;

% Extract the relevant portion of the fear_matrix
cut_fear_matrix = fear_matrix(indices_to_keep, :);
end

function [mean_value, variance_value, total_amplitude, total_rise_time, total_energy, max_min_difference, time_difference_max_min, num_maxima, power, gsr_bandwidth] = analyze(data, maxima, minima)
gsr_values = data.GSRValue;
% Calculate mean and variance
sum_values = 0;
for i = 1:length(gsr_values)
    sum_values = sum_values + gsr_values(i);
end
mean_value = sum_values / length(gsr_values);

% Initialize variance calculation
sum_squared_diff = 0;

% Continue calculating variance
for i = 1:length(gsr_values)
    squared_diff = (gsr_values(i) - mean_value)^2;
    sum_squared_diff = sum_squared_diff + squared_diff;
end
variance_value = sum_squared_diff / length(gsr_values);


% Peak analysis
max_indices = find(maxima);
min_indices = find(minima);
total_rise_time = 0;
total_amplitude = 0;
total_energy = 0;
max_min_difference = 0;
time_difference_max_min = 0;
num_maxima = length(max_indices);
gsr_bandwidth_numerator = 0;
gsr_bandwidth_denominator = 0;

for i = 1:length(max_indices)
    current_max_index = max_indices(i);
    
    % Find the nearest preceding minimum index
    preceding_min_index = max(min_indices(min_indices < current_max_index));
    
    % Calculate amplitude and add to total
    total_amplitude = total_amplitude + data.GSRValue(current_max_index);
    
    % Calculate rise time and add to total
    rise_time = seconds(data.TimeStamp(current_max_index) - data.TimeStamp(preceding_min_index));
    total_rise_time = total_rise_time + rise_time;
    
    % Calculate energy and add to total
    total_energy = total_energy + 0.5 * data.GSRValue(current_max_index) * rise_time;
    
    % Calculate max-min difference and update if it's the highest so far
    current_difference = data.GSRValue(current_max_index) - data.GSRValue(preceding_min_index);
    if current_difference > max_min_difference
        max_min_difference = current_difference;
        
        % Calculate time difference between highest minimum and maximum with highest amplitude difference
        time_difference_max_min = seconds(data.TimeStamp(current_max_index) - data.TimeStamp(preceding_min_index));
    end
    
    % Calculate numerator and denominator for GSR bandwidth
    if i > 1
        gsr_bandwidth_numerator = gsr_bandwidth_numerator + (data.GSRValue(current_max_index) - data.GSRValue(current_max_index - 1))^2;
    end
    gsr_bandwidth_denominator = gsr_bandwidth_denominator + data.GSRValue(current_max_index)^2;
end

% Calculate power: total energy / time difference between first and last data points
power = total_energy / seconds(data.TimeStamp(end) - data.TimeStamp(1));

% Calculate GSR bandwidth
gsr_bandwidth = (1 / (2 * pi)) * sqrt(gsr_bandwidth_numerator / gsr_bandwidth_denominator);
end

function Fear_Index = calculate_fear_index(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10)
Fear_Index = F1 + 2 * F2 + F3 + 0.5 * F4 + F5 + 2 * F6 + F7 + 5 * F8 + 0.001 * F9 + 0.5 * F10;
end