function eeg_data_timewindow = time_window_shifting(eeg_data, window_length, overlap)
% INPUTS:
%          - eeg_data: samples x channel
%          - window_length: length of each time window
%          - overlap: the overlap between windows
% OUTPUTS:
%          - eeg_data_timewindow:  channles x samples x windows 

[sample_num, Ch_num] = size(eeg_data);
num_windows = floor((sample_num  - window_length) / (window_length - overlap))+1; 
eeg_data_timewindow = zeros(Ch_num, window_length, num_windows);

for i = 1:Ch_num
    for j = 1:num_windows
        start_sample = (j-1)*(window_length - overlap)+1;
        end_sample = start_sample+window_length-1;
   
        if end_sample <= sample_num
            window_sample = start_sample:end_sample;
            eeg_data_timewindow(i,:,j) = eeg_data(window_sample,i);
        else
            break;
        end
    end
end
end

