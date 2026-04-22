clear;
clc;
dataPath='D:\Matlab_project\dataset\Epilepsy_dataset\seizure8';

M = 7;
Epilepsy_data(M,1) = struct('id',[],'signal',[]);
Preictal_data(M,1) = struct('id',[],'signal',[]);
Ictal_data(M,1) = struct('id',[],'signal',[]);


for i = 1:M
    file_index = 1;
    filename_pre = sprintf('sz%d_pre_clean.mat', file_index );
    filename_ict = sprintf('sz%d_ict_clean.mat', file_index );
    Data_pre = load(fullfile(dataPath, filename_pre));
    Data_ict = load(fullfile(dataPath, filename_ict));
    eeg_data_pre  = Data_pre.pre_eeg;
    eeg_data_ict  = Data_ict.ict_eeg;
    
    % IEEG data all period
    eeg_data = [eeg_data_pre;eeg_data_ict];
    Epilepsy_data(i).id = i;
    Epilepsy_data(i).signal =  eeg_data;
    
    %IEEG data preictal period
    Preictal_data(i).id = i;
    Preictal_data(i).signal =  eeg_data_pre;
    %IEEG data ictal period
    Ictal_data(i).id = i;
    Ictal_data(i).signal =  eeg_data_ict;
    
end

