%% Loading data and performing HRV analysis

if exist("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Rest_output", 'dir')
    rmdir("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Rest_output", "s");
end
if exist("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Only_video_output", 'dir')
    rmdir("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Only_video_output", "s");
end
if exist("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Video_output", 'dir')
    rmdir("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Video_output", "s");
end
if exist("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\HRV_Output", 'dir')
    rmdir("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\HRV_Output", "s");
end



files = dir("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\spider_database\spider_ecg\**\BitalinoECG.txt");
triggers_files = dir("C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\spider_database\spider_ecg\**\Triggers.txt");
% Settings
fs = 100;
plot_ecg = 0;
run_not_mean = 1;
run_mean = 0;
enable_normalize = 0;


for i = 1:length(files)
    temp = readtable(strcat(files(i).folder, "\", files(i).name));
    ch1 = table2array(temp(:, "Var1")); % ECG signal
    time_stamps = table2array(temp(:, "Var2")); % Channel 2
    N = length(ch1);
    empty_t = [];
    
    % HRV Initialize settings
    sub_temp = strsplit(files(i).folder, "\");
    subid = string(sub_temp(end));
    low_quality_threshold = 0.6;
    peak_threshold = 0.55;
    win_len = 180;
    increment = 30;
    video = ch1(N-1201*fs:N-301*fs);

    
    % Splits data according to timestamps
    triggers = readtable(strcat(triggers_files(i).folder, "\", triggers_files(i).name));

    % Finds which rows contain timestamps
    row_rest = find(strcmp('BIOFEEDBACK-REST',table2array(triggers(:, 1))));
    % Skips if the signal doesn't contain marked rest
    if isempty(row_rest)
        msgbox(strcat(subid, " skipped. Rest period not marked."), "WARNING", "warn")
        continue;
    end
    IndexClip = strfind(table2array(triggers(:, 1)),'CLIP');
    IndexClip = find(not(cellfun('isempty',IndexClip)));
    row_video_start = min(IndexClip) + 1;  % Skips demo clip
    row_video_end = max(IndexClip);

    % Finds timestamps in file
    video_start = table2array(triggers(row_video_start, 2));
    video_end = table2array(triggers(row_video_end, 3));
    rest_start = table2array(triggers(row_rest, 2));
    rest_end = table2array(triggers(row_rest, 3));

    if rest_end > max(time_stamps)
        rest_end = max(time_stamps);
    end
    
    % Converts timestamps into indices
    i_video_start = find((video_start - time_stamps) < 0.1, 1);
    i_video_end = find((video_end - time_stamps) < 0.1, 1);
    i_rest_start = find((rest_start - time_stamps) < 0.1, 1);
    i_rest_end = find((rest_end - time_stamps) < 0.1, 1);

    % Splits raw_ecg into stressed and rest ecg
    video_ecg = ch1(i_video_start:i_video_end);
    rest_ecg = ch1(i_rest_start:i_rest_end);
    

    % Only videos
    clear HRVparams
    HRVparams = InitializeHRVparams('');
    HRVparams.Fs = fs;
    HRVparams.windowlength = win_len;
    HRVparams.increment = increment;
    HRVparams.sqi.LowQualityThreshold = low_quality_threshold;
    HRVparams.PeakDetect.THRES = peak_threshold;
    HRVparams.HRT.on = 0;
    HRVparams.DFA.on = 0;
    HRVparams.MSE.on = 0;
    HRVparams.af.on = 0;
    HRVparams.output.separate = 0;
    HRVparams.writedata = 'Only_video_output';

    [temp_result_only_video, file_temp_only_video] = Main_HRV_Analysis(video_ecg, empty_t, 'ECGWaveform', HRVparams, subid);
    if (length(temp_result_only_video) == 1)
        msgbox(strcat(subid, " skipped, HRV analysis failed. Window size is to small for signal (small signal length)."), "ERROR", "error")
        continue
    end
    if (i == 1)
        only_video_results_mean = mean(temp_result_only_video, 1);
    else
        only_video_results_mean = [only_video_results_mean; mean(temp_result_only_video, 1)];
    end
    
    %{
    % Video
    clear HRVparams
    HRVparams = InitializeHRVparams('');
    HRVparams.Fs = fs;
    HRVparams.windowlength = win_len;
    HRVparams.increment = increment;
    HRVparams.sqi.LowQualityThreshold = low_quality_threshold;
    HRVparams.PeakDetect.THRES = peak_threshold;
    HRVparams.HRT.on = 0;
    HRVparams.DFA.on = 0;
    HRVparams.MSE.on = 0;
    HRVparams.af.on = 0;
    HRVparams.output.separate = 0;
    HRVparams.writedata = 'Video_output';
    
    [temp_result_video, file_temp_video] = Main_HRV_Analysis(video, empty_t, 'ECGWaveform', HRVparams, subid);
    if (i == 1)
        video_results_mean = mean(temp_result_video);
    else
        video_results_mean = [video_results_mean; mean(temp_result_video)];
    end
    %}
    
    if(plot_ecg == 1)
        fh = figure();
        fh.WindowState = 'maximized';
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(video_ecg,HRVparams,subid);
        plot(video_ecg)
        hold on
        plot(jqrs_ann,video_ecg(jqrs_ann),'o') 
        title(strcat("Video ECG ", subid))
        hold off
    end

    
    % Rest
    clear HRVparams
    HRVparams = InitializeHRVparams('');
    HRVparams.Fs = fs;
    HRVparams.windowlength = win_len;
    HRVparams.increment = increment;
    HRVparams.sqi.LowQualityThreshold = low_quality_threshold;
    HRVparams.PeakDetect.THRES = peak_threshold;
    HRVparams.HRT.on = 0;
    HRVparams.DFA.on = 0;
    HRVparams.MSE.on = 0;
    HRVparams.af.on = 0;
    HRVparams.output.separate = 0;
    HRVparams.writedata = 'Rest_output';

    [temp_result_rest, file_temp_rest] = Main_HRV_Analysis(rest_ecg, empty_t, 'ECGWaveform', HRVparams, subid);
    if (length(temp_result_rest) == 1)
        only_video_results_mean(end, :) = []; % Need to remove the row from video for this subject that was added prior
        msgbox(strcat(subid, " skipped, HRV analysis failed. Window size is to small for signal (small signal length)."), "ERROR", "error")
        continue
    end
    if (i == 1)
        rest_results_mean = mean(temp_result_rest, 1);
    else
        rest_results_mean = [rest_results_mean; mean(temp_result_rest, 1)];
    end
         
    if(plot_ecg == 1)
        fh = figure();
        fh.WindowState = 'maximized';
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(rest_ecg,HRVparams,subid);
        plot(rest_ecg)
        hold on
        plot(jqrs_ann,rest_ecg(jqrs_ann),'o') 
        title(strcat("Rest ECG ", subid))
        hold off
    end
end


%% !! Not mean !! 
if (run_not_mean == 1)
%% Plot (not mean)    

    d            = dir('C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Only_video_output\AllPatients_HRV_results_allwindows_*.csv');
    [~, index]   = max([d.datenum]);
    youngestFile_video = fullfile(d(index).folder, d(index).name);  
    only_video_results = readtable(youngestFile_video);

    %{
    d            = dir('C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Video_output\AllPatients_HRV_results_allwindows_*.csv');
    [~, index]   = max([d.datenum]);
    youngestFile_video = fullfile(d(index).folder, d(index).name);  
    video_results = readtable(youngestFile_video);
    %}
    
    d            = dir('C:\Users\Mikael\OneDrive\Desktop\Elektronik\ELA412 Advancerad signalbehandling\Project\Rest_output\AllPatients_HRV_results_allwindows_*.csv');
    [~, index]   = max([d.datenum]);
    youngestFile_rest = fullfile(d(index).folder, d(index).name);
    rest_results = readtable(youngestFile_rest);

    
    if length(table2array(only_video_results(:, 1))) < length(table2array(rest_results(:, 1)))
        arr_len = length(table2array(only_video_results(:, 1)));
    else
        arr_len = length(table2array(rest_results(:, 1)));
    end
    
    
    if (enable_normalize == 1)
        only_video_results(:, 2:end) = normalize(only_video_results(:, 2:end));
        rest_results(:, 2:end) = normalize(rest_results(:, 2:end));
    end
    
    
    % Time
    %{
    nnmean_video = table2array(only_video_results(:, 4));
    nnmean_rest = table2array(rest_results(:, 4));
    
    nnmedian_video = table2array(only_video_results(:, 5));
    nnmedian_rest = table2array(rest_results(:, 5));
    %}
    
    sdnn_video = table2array(only_video_results(:, 11));
    sdnn_rest = table2array(rest_results(:, 11));
    
    rmssd_video = table2array(only_video_results(:, 12));
    rmssd_rest = table2array(rest_results(:, 12));
    
    pnn50_video = table2array(only_video_results(:, 13));
    pnn50_rest = table2array(rest_results(:, 13));
    
    
    %{
    figure
    boxplot([nnmean_rest(1:arr_len), nnmean_video(1:arr_len), nnmedian_rest(1:arr_len), nnmedian_video(1:arr_len)],'Notch','on','Labels', ...
        {'Rest NNmean', 'Video NNmean', 'Rest NNmedian', 'Video NNmedian'})
    title('Box plot of time domain results from all patients, video vs rest (not mean)')
    hold off
    
    figure
    violinplot([nnmean_rest(1:arr_len), nnmean_video(1:arr_len), nnmedian_rest(1:arr_len), nnmedian_video(1:arr_len)], ...
        {'Rest NNmean', 'Video NNmean', 'Rest NNmedian', 'Video NNmedian'});
    title('Violin plot of time domain results from all patients, video vs rest (not mean)')
    hold off
    %}
    
    
    figure
    boxplot([sdnn_rest(1:arr_len), sdnn_video(1:arr_len), rmssd_rest(1:arr_len), rmssd_video(1:arr_len)],'Notch','on','Labels', ...
        {'Rest SDNN', 'Video SDNN', 'Rest RMSSD', 'Video RMSSD'})
    title('Box plot of time domain results from all patients, video vs rest (not mean)')
    hold off
    
    figure
    violinplot([sdnn_rest(1:arr_len), sdnn_video(1:arr_len), rmssd_rest(1:arr_len), rmssd_video(1:arr_len)], ...
        {'Rest SDNN', 'Video SDNN', 'Rest RMSSD', 'Video RMSSD'});
    title('Violin plot of time domain results from all patients, video vs rest (not mean)')
    hold off
    
    
    figure
    boxplot([pnn50_rest(1:arr_len), pnn50_video(1:arr_len)],'Notch','on','Labels', ...
        {'Rest pnn50', 'Video pnn50'})
    title('Box plot of time domain results from all patients, video vs rest (not mean)')
    hold off
    
    figure
    violinplot([pnn50_rest(1:arr_len), pnn50_video(1:arr_len)], ...
        {'Rest pnn50', 'Video pnn50'});
    title('Violin plot of time domain results from all patients, video vs rest (not mean)')
    hold off

    
    
    
    % Frequency
    video_results_lf = table2array(only_video_results(:, 19));
    rest_results_lf = table2array(rest_results(:, 19));

    video_results_hf = table2array(only_video_results(:, 20));
    rest_results_hf = table2array(rest_results(:, 20));

    video_results_lfhf = table2array(only_video_results(:, 21));
    rest_results_lfhf = table2array(rest_results(:, 21));

    
    figure
    boxplot([rest_results_lf(1:arr_len), video_results_lf(1:arr_len), rest_results_hf(1:arr_len), video_results_hf(1:arr_len)],'Notch','on','Labels',{'Rest LF', 'Video LF', 'Rest HF', 'Video HF'})
    title('Box plot of frequency results (LF and HF) from all patients, video vs rest (not mean)')
    hold off

    figure
    violinplot([rest_results_lf(1:arr_len), video_results_lf(1:arr_len), rest_results_hf(1:arr_len), video_results_hf(1:arr_len)], {'Rest LF', 'Video LF', 'Rest HF', 'Video HF'});
    title('Violin plot of frequency results (LF and HF) from all patients, video vs rest (not mean)')
    hold off


    figure
    boxplot([rest_results_lfhf(1:arr_len), video_results_lfhf(1:arr_len)],'Notch','on','Labels',{'Rest LF/HF', 'Video LF/HF'})
    title('Box plot of LF/HF from all patients, video vs rest (not mean)')
    hold off

    figure
    violinplot([rest_results_lfhf(1:arr_len), video_results_lfhf(1:arr_len)], {'Rest LF/HF', 'Video LF/HF'});
    title('Violin plot of LF/HF from all patients, video vs rest (not mean)')
    hold off

%% paired t-test (not mean)

    % Time
    %{
    [test_decision_null_hypothesis_nnmean, p_nnmean, confidence_interval_nnmean, stats_nnmean] = ttest(nnmean_rest(1:arr_len), nnmean_video(1:arr_len));
    [test_decision_null_hypothesis_nnmedian, p_nnmedian, confidence_interval_nnmedian, stats_nnmedian] = ttest(nnmedian_rest(1:arr_len), nnmedian_video(1:arr_len));
    %}
    [test_decision_null_hypothesis_sdnn, p_sdnn, confidence_interval_sdnn, stats_sdnn] = ttest(sdnn_rest(1:arr_len), sdnn_video(1:arr_len));
    [test_decision_null_hypothesis_rmssd, p_rmssd, confidence_interval_rmssd, stats_rmssd] = ttest(rmssd_rest(1:arr_len), rmssd_video(1:arr_len));
    [test_decision_null_hypothesis_pnn50, p_pnn50, confidence_interval_pnn50, stats_pnn50] = ttest(pnn50_rest(1:arr_len), pnn50_video(1:arr_len));

    paired_t_test_time_table = table([test_decision_null_hypothesis_sdnn, p_sdnn], ...
        [test_decision_null_hypothesis_rmssd, p_rmssd], [test_decision_null_hypothesis_pnn50, p_pnn50]);
    paired_t_test_time_table.Properties.VariableNames = ["SDNN", "RMSSD", "pnn50"]

    
    % Frequency
    [test_decision_null_hypothesis_lf, p_lf, confidence_interval_lf, stats_lf] = ttest(rest_results_lf(1:arr_len), video_results_lf(1:arr_len));
    [test_decision_null_hypothesis_hf, p_hf, confidence_interval_hf, stats_hf] = ttest(rest_results_hf(1:arr_len), video_results_hf(1:arr_len));
    [test_decision_null_hypothesis_lfhf, p_lfhf, confidence_interval_lfhf, stats_lfhf] = ttest(rest_results_lfhf(1:arr_len), video_results_lfhf(1:arr_len));

    paired_t_test_freq_table = table([test_decision_null_hypothesis_lf, p_lf], [test_decision_null_hypothesis_hf, p_hf], [test_decision_null_hypothesis_lfhf, p_lfhf]);
    paired_t_test_freq_table.Properties.VariableNames = ["LF", "HF", "LF/HF"]
    
end


%% !! Mean !!
if (run_mean == 1)
%% Plot (mean)

    if (enable_normalize == 1)
        only_video_results_mean = normalize(only_video_results_mean);
        rest_results_mean = normalize(rest_results_mean);
    end
    
    % Time
    %{
    nnmean_video_mean = only_video_results_mean(:, 3);
    nnmean_rest_mean = rest_results_mean(:, 3);
    
    nnmedian_video_mean = only_video_results_mean(:, 4);
    nnmedian_rest_mean = rest_results_mean(:, 4);
    %}
    
    sdnn_video_mean = only_video_results_mean(:, 10);
    sdnn_rest_mean = rest_results_mean(:, 10);
    
    rmssd_video_mean = only_video_results_mean(:, 11);
    rmssd_rest_mean = rest_results_mean(:, 11);
    
    pnn50_video_mean = only_video_results_mean(:, 12);
    pnn50_rest_mean = rest_results_mean(:, 12);
    
    
    %{
    figure
    boxplot([nnmean_rest_mean, nnmean_video_mean, nnmedian_rest_mean, nnmedian_video_mean],'Notch','on','Labels', ...
        {'Rest NNmean', 'Video NNmean', 'Rest NNmedian', 'Video NNmedian'})
    title('Box plot of time domain results from all patients, video vs rest (mean)')
    hold off
    
    figure
    violinplot([nnmean_rest_mean, nnmean_video_mean, nnmedian_rest_mean, nnmedian_video_mean], ...
        {'Rest NNmean', 'Video NNmean', 'Rest NNmedian', 'Video NNmedian'});
    title('Violin plot of time domain results from all patients, video vs rest (mean)')
    hold off
    %}
    
    
    figure
    boxplot([sdnn_rest_mean, sdnn_video_mean, rmssd_rest_mean, rmssd_video_mean],'Notch','on','Labels', ...
        {'Rest SDNN', 'Video SDNN', 'Rest RMSSD', 'Video RMSSD'})
    title('Box plot of time domain results from all patients, video vs rest (mean)')
    hold off
    
    figure
    violinplot([sdnn_rest_mean, sdnn_video_mean, rmssd_rest_mean, rmssd_video_mean], ...
        {'Rest SDNN', 'Video SDNN', 'Rest RMSSD', 'Video RMSSD'});
    title('Violin plot of time domain results from all patients, video vs rest (mean)')
    hold off
    
    
    figure
    boxplot([pnn50_rest_mean, pnn50_video_mean],'Notch','on','Labels', ...
        {'Rest pnn50', 'Video pnn50'})
    title('Box plot of time domain results from all patients, video vs rest (mean)')
    hold off
    
    figure
    violinplot([pnn50_rest_mean, pnn50_video_mean], ...
        {'Rest pnn50', 'Video pnn50'});
    title('Violin plot of time domain results from all patients, video vs rest (mean)')
    hold off
    
    
    
    
    % Frequency
    video_results_mean_lf = only_video_results_mean(:, 18);
    rest_results_mean_lf = rest_results_mean(:, 18);

    video_results_mean_hf = only_video_results_mean(:, 19);
    rest_results_mean_hf = rest_results_mean(:, 19);

    video_results_mean_lfhf = only_video_results_mean(:, 20);
    rest_results_mean_lfhf = rest_results_mean(:, 20);
    
    
    figure
    boxplot([rest_results_mean_lf, video_results_mean_lf, rest_results_mean_hf, video_results_mean_hf],'Notch','on','Labels',{'Rest LF', 'Video LF', 'Rest HF', 'Video HF'})
    title('Box plot of frequency results (LF and HF) from all patients, video vs rest (mean)')
    hold off

    figure
    violinplot([rest_results_mean_lf, video_results_mean_lf, rest_results_mean_hf, video_results_mean_hf], {'Rest LF', 'Video LF', 'Rest HF', 'Video HF'});
    title('Violin plot of frequency results (LF and HF) from all patients, video vs rest (mean)')
    hold off


    figure
    boxplot([rest_results_mean_lfhf, video_results_mean_lfhf],'Notch','on','Labels',{'Rest LF/HF', 'Video LF/HF'})
    title('Box plot of LF/HF from all patients, video vs rest (mean)')
    hold off

    figure
    violinplot([rest_results_mean_lfhf, video_results_mean_lfhf], {'Rest LF/HF', 'Video LF/HF'});
    title('Violin plot of LF/HF from all patients, video vs rest (mean)')
    hold off
    
%% paired t-test (mean)

    % Time
    %{
    [test_decision_null_hypothesis_mean_nnmean, p_mean_nnmean, confidence_interval_mean_nnmean, stats_mean_nnmean] = ttest(nnmean_rest_mean, nnmean_video_mean);
    [test_decision_null_hypothesis_mean_nnmedian, p_mean_nnmedian, confidence_interval_mean_nnmedian, stats_mean_nnmedian] = ttest(nnmedian_rest_mean, nnmedian_video_mean);
    %}
    [test_decision_null_hypothesis_mean_sdnn, p_mean_sdnn, confidence_interval_mean_sdnn, stats_mean_sdnn] = ttest(sdnn_rest_mean, sdnn_video_mean);
    [test_decision_null_hypothesis_mean_rmssd, p_mean_rmssd, confidence_interval_mean_rmssd, stats_mean_rmssd] = ttest(rmssd_rest_mean, rmssd_video_mean);
    [test_decision_null_hypothesis_mean_pnn50, p_mean_pnn50, confidence_interval_mean_pnn50, stats_mean_pnn50] = ttest(pnn50_rest_mean, pnn50_video_mean);

    paired_t_test_time_mean_table = table([test_decision_null_hypothesis_mean_sdnn, p_mean_sdnn], ...
        [test_decision_null_hypothesis_mean_rmssd, p_mean_rmssd], [test_decision_null_hypothesis_mean_pnn50, p_mean_pnn50]);
    paired_t_test_time_mean_table.Properties.VariableNames = ["SDNN", "RMSSD", "pnn50"]

    
    
    % Frequency
    [test_decision_null_hypothesis_mean_lf, p_mean_lf, confidence_interval_mean_lf, stats_mean_lf] = ttest(rest_results_mean_lf, video_results_mean_lf);
    [test_decision_null_hypothesis_mean_hf, p_mean_hf, confidence_interval_mean_hf, stats_mean_hf] = ttest(rest_results_mean_hf, video_results_mean_hf);
    [test_decision_null_hypothesis_mean_lfhf, p_mean_lfhf, confidence_interval_mean_lfhf, stats_mean_lfhf] = ttest(rest_results_mean_lfhf, video_results_mean_lfhf);

    paired_t_test_freq_mean_table = table([test_decision_null_hypothesis_mean_lf, p_mean_lf], [test_decision_null_hypothesis_mean_hf, p_mean_hf], [test_decision_null_hypothesis_mean_lfhf, p_mean_lfhf]);
    paired_t_test_freq_mean_table.Properties.VariableNames = ["LF", "HF", "LF/HF"]

end