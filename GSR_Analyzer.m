%Our MATLAB licence does not allow to install Signal Processing Toolbox
%so I had to write my own functions


GSR_fear_file = 'GSR_FEAR.csv';
GSR_baseline_file = 'GSR_Baseline.csv';
fear_matrix = read_gsr_data(GSR_fear_file);
baseline_matrix = read_gsr_data(GSR_baseline_file);
disp('GSR data read successfully');


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

%Apply median filter
fear_matrix.GSRValue = third_order_median_filter(fear_matrix.GSRValue);
baseline_matrix.GSRValue = third_order_median_filter(baseline_matrix.GSRValue);
disp('Median filter applied successfully');

subplot(6, 2, 3);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Median Filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 4);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Median Filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

%Apply low pass filter
fear_matrix.GSRValue = low_pass_filter(fear_matrix.GSRValue, 20/7.5, 15, 3);
baseline_matrix.GSRValue = low_pass_filter(baseline_matrix.GSRValue, 20/7.5, 15, 3);
disp('Low pass filter applied successfully');

subplot(6, 2, 5);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Low Pass Filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 6);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Low Pass Filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

%Apply 10-point moving average filter
fear_matrix.GSRValue = moving_average_filter(fear_matrix.GSRValue, 11);
baseline_matrix.GSRValue = moving_average_filter(baseline_matrix.GSRValue, 11);
disp('Moving average filter applied successfully');

subplot(6, 2, 7);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Moving Average Filtered GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 8);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Moving Average Filtered GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

%Normalize GSR values
fear_matrix.GSRValue = normalize_gsr(fear_matrix.GSRValue);
baseline_matrix.GSRValue = normalize_gsr(baseline_matrix.GSRValue);
disp('GSR values normalized successfully');

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

% Cut out the first two seconds of data
start_time = fear_matrix.TimeStamp(1) + seconds(2);
fear_matrix = fear_matrix(fear_matrix.TimeStamp >= start_time, :);

start_time = baseline_matrix.TimeStamp(1) + seconds(2);
baseline_matrix = baseline_matrix(baseline_matrix.TimeStamp >= start_time, :);

subplot(6, 2, 11);
plot(fear_matrix.TimeStamp, fear_matrix.GSRValue);
title('Processed GSR Values (Fear)');
xlabel('Time');
ylabel('GSR Value');
grid on;

subplot(6, 2, 12);
plot(baseline_matrix.TimeStamp, baseline_matrix.GSRValue);
title('Processed GSR Values (Baseline)');
xlabel('Time');
ylabel('GSR Value');
grid on;

%Calculate mean GSR
mean_gsr_fear = mean(fear_matrix.GSRValue);
mean_gsr_baseline = mean(baseline_matrix.GSRValue);

%Calculate variance GSR
variance_gsr_fear = var(fear_matrix.GSRValue);
variance_gsr_baseline = var(baseline_matrix.GSRValue);


%Calculate sum of peak rise times
[pear_rise_time_sum_fear, peak_count_fear] = calculate_sum_peak_rise_time(fear_matrix.GSRValue, fear_matrix.TimeStamp);
[pear_rise_time_sum_baseline, peak_count_baseline] = calculate_sum_peak_rise_time(baseline_matrix.GSRValue, baseline_matrix.TimeStamp);

% Calculate sum of peak amplitudes
peak_amplitude_sum_fear = calculate_sum_peak_amplitude(fear_matrix.GSRValue);
peak_amplitude_sum_baseline = calculate_sum_peak_amplitude(baseline_matrix.GSRValue);

% Calculate sum peak energy
energy_sum_fear = calculate_peak_energy_sum(fear_matrix.GSRValue, fear_matrix.TimeStamp);
energy_sum_baseline = calculate_peak_energy_sum(baseline_matrix.GSRValue, baseline_matrix.TimeStamp);

[highest_SCR_fear, highest_SCR_rise_time_fear] = find_max_amplitude_and_rise_time(fear_matrix.GSRValue, fear_matrix.TimeStamp);
[highest_SCR_baseline, highest_SCR_rise_time_baseline] = find_max_amplitude_and_rise_time(baseline_matrix.GSRValue, baseline_matrix.TimeStamp);





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


function filtered_data = third_order_median_filter(data)
% Pad the data at the beginning and end with zeros
padded_data = [zeros(2, 1); data; zeros(2, 1)];

% Initialize the filtered data vector
filtered_data = zeros(size(data));

% Apply the 3rd order median filter
for i = 1:length(data)
    window = padded_data(i:i+4);
    filtered_data(i) = median(window);
end
end

function filtered_signal = low_pass_filter(input_signal, passband_freq, sampling_rate, order)
% Normalize the cutoff frequency
normalized_cutoff = passband_freq / (sampling_rate/2);

% Design the filter using a Butterworth filter
[b, a] = butter(order, normalized_cutoff, 'low');

% Apply the filter to the input signal
filtered_signal = filter(b, a, input_signal);
end


function smoothed_signal = moving_average_filter(input_signal, N)
% Apply the N-point moving average filter
b = ones(1, N) / N;
smoothed_signal = conv(input_signal, b, 'same');
end


function normalized_data = normalize_gsr(data)
% Calculate the minimum and maximum values of the data
min_value = min(data);
max_value = max(data);

% Normalize the data to the range [0, 100]
normalized_data = 100 * (data - min_value) / (max_value - min_value);
end


function [sum_rise_time, peak_count] = calculate_sum_peak_rise_time(gsr_values, timestamps)
% Identify peaks
[~, locs] = findpeaks(gsr_values, 'MinPeakProminence', 2);

% Calculate rise times from each peak to the base (3 seconds earlier)
rise_times = zeros(size(locs));

for i = 1:length(locs)
    peak_location = locs(i);
    
    % Find the timestamp 3 seconds earlier
    base_timestamp = max(timestamps(peak_location) - seconds(3), timestamps(1));
    
    % Find the index of the timestamp closest to the base_timestamp
    [~, base_index] = min(abs(timestamps - base_timestamp));
    
    % Calculate the rise time in seconds
    rise_times(i) = seconds(timestamps(peak_location) - timestamps(base_index));
end

% Sum of peak rise times
sum_rise_time = sum(rise_times);

% Count of peaks
peak_count = length(locs);
end

function sum_amplitude = calculate_sum_peak_amplitude(gsr_values)
% Identify peaks
[peaks, ~] = findpeaks(gsr_values,'MinPeakProminence', 2);

% Sum of peak amplitudes
sum_amplitude = sum(peaks);
end

function sum_energy = calculate_peak_energy_sum(gsr_values, timestamps)
% Identify peaks
[peaks, locs] = findpeaks(gsr_values,'MinPeakProminence', 2);

% Calculate rise times
rise_times = diff(timestamps(locs));

% Calculate peak energy sum
sum_energy = sum(0.5 * peaks(1:end-1) .* seconds(rise_times));
end

function [max_amplitude, rise_time] = find_max_amplitude_and_rise_time(gsr_values, timestamps)
% Identify peaks
[peaks, locs] = findpeaks(gsr_values, 'MinPeakProminence', 2);

% Initialize variables to store maximum amplitude and rise time
max_amplitude = 0;
rise_time = 0;

% Loop through all peaks to find the highest amplitude
for i = 1:length(peaks)
    % Find the index of the start and end of the peak
    peak_start_index = locs(i);
    if i < length(locs)
        peak_end_index = locs(i + 1);
    else
        peak_end_index = length(gsr_values);
    end
    
    % Extract the portion of the signal between the peak start and end
    peak_signal = gsr_values(peak_start_index:peak_end_index);
    timestamp_peak = timestamps(peak_start_index:peak_end_index);
    
    % Find the timestamp 3 seconds earlier
    base_timestamp = max(timestamp_peak(1) - seconds(3), timestamps(1));
    
    % Find the index of the timestamp closest to the base_timestamp
    [~, base_index] = min(abs(timestamp_peak - base_timestamp));
    
    % Calculate the amplitude between peak and base
    current_amplitude = peaks(i) - min(peak_signal);
    
    % Calculate the rise time correctly and convert to seconds
    [~, max_index] = max(peak_signal);
    current_rise_time = seconds(abs(timestamp_peak(max_index) - timestamps(base_index)));
    
    % Update max amplitude and rise time if needed
    if current_amplitude > max_amplitude
        max_amplitude = current_amplitude;
        rise_time = current_rise_time;
    end
end
end

