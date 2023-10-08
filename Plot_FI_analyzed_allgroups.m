% plot aspects of FI curves from .mat files from Wei's code

%To plot an specific condition uncomment the two groups you will like to
%regroup. 

%WT Animals (All conditions plotted as one)
% FI_DRP = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/DRPthrough062222to063022.mat');
% FI_EVP= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/EVPthrough062222to063022.mat');
% FI_DRP_WT= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/DRNthrough062222to063022.mat');
% FI_EVP_WT= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/EVNthrough062222to063022.mat');


%SCN2A Animals (All conditions plotted as one)
%FI_DRP = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/DRPthrough081022to100522.mat'); 
%FI_EVP = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/EVPthrough081022to100522.mat');
%FI_DRP_WT = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/DRNthrough081022to100522.mat');
%FI_EVP_WT = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/EVNthrough081022to100522.mat');

%SCN2A vs WT
FI_DRP = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/nuevoexp2023/SCN2A_DRP_Aug10_22toMay22_23');
FI_EVP = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/SCN2AmergeEVN_EVP_June22_22toMay22_23.mat');

%oad('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023_h5_files/mergeEVN_EVP_062222to021223.mat');
%FI_DRP_WT =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/olddata_DRP_WT_062222to063022.mat');
FI_DRP_WT =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/WT_DRP_June22_22toMay4_23.mat');
FI_EVP_WT =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/WT_mergeEVN_EVP_Jun22_22toMay4_23.mat');


%FI_DRP_WT = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/olddata_DRP_WT_062222to063022.mat');
%FI_EVP_WT=load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023_h5_files/mergeEVN_EVP_062222to021223.mat');

FI_DRP = FI_DRP.FI_DATA; %DR+
FI_DRP_WT = FI_DRP_WT.FI_DATA; %DR-
FI_EVP = FI_EVP.FI_DATA; %EV+
FI_EVP_WT = FI_EVP_WT.FI_DATA; %EV-



DRP_IFR_ave = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_IFR_ave = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_IFR_ave_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_IFR_ave_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));

DRP_MFR = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_MFR = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_MFR_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_MFR_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));


DRP_IFR_AUC = nan(length(FI_DRP.filename), 1);
EVP_IFR_AUC = nan(length(FI_EVP.filename), 1);
DRP_IFR_AUC_WT = nan(length(FI_DRP_WT.filename), 1);
EVP_IFR_AUC_WT = nan(length(FI_EVP_WT.filename), 1);


DRP_rheobase = nan(length(FI_DRP.filename), 1);
EVP_rheobase = nan(length(FI_EVP.filename), 1);
DRP_rheobase_WT = nan(length(FI_DRP_WT.filename), 1);
EVP_rheobase_WT = nan(length(FI_EVP_WT.filename), 1);

% 
DRP_latency = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_latency = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_latency_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_latency_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));

DRP_width = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_width = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_width_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_width_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));

DRP_Vth = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_Vth = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_Vth_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_Vth_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));


DRP_ADP = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_ADP = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_ADP_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_ADP_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));


DRP_AP_peak_ave = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_AP_peak_ave = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_AP_peak_ave_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_AP_peak_ave_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));

DRP_spike_count = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_spike_count = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_spike_count_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_spike_count_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));

DRP_AP_peak_1st = nan(length(FI_DRP.filename), length(FI_DRP.curr_inj{1,1}));
EVP_AP_peak_1st = nan(length(FI_EVP.filename), length(FI_EVP.curr_inj{1,1}));
DRP_AP_peak_1st_WT = nan(length(FI_DRP_WT.filename), length(FI_DRP_WT.curr_inj{1,1}));
EVP_AP_peak_1st_WT = nan(length(FI_EVP_WT.filename), length(FI_EVP_WT.curr_inj{1,1}));

grayColor = [.5 .5 .5];
burgande= [.70 .130 .80];
blue= [.30 .74 .73];
green=[0.4660 0.6740 0.1880];


grayColor1 = [.80 .80 .80];
burgande1= [.70 .130 .80];
blue1= [.45 .89 .90];
green1=[0.59 0.80 0.31];
black= [0.2 0.2 0.2];
Colors1 = [black; green1; grayColor1; blue1];
Colors = [0 0 0; green; grayColor; blue];

count_DRP = 0;
for h = 1:length(FI_DRP.filename)
    
    %Exclusion criteria for FI
    if mean(FI_DRP.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
        %exclude cells with an input resistance lower 90Mohms and higher
        %than 350Mohms
    elseif mean(FI_DRP.Rin{1,h}) < 90000000 || mean(FI_DRP.Rin{1,h}) > 350000000 
%         continue
%    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
    %younger cells 
    %potentially memebrane voltage
        continue
    end
    
    count_DRP = count_DRP + 1;
    
    DRP_IFR_ave(h, :) = FI_DRP.IFR_ave{1,h}';
    DRP_MFR(h, :) = FI_DRP.MFR{1,h}';
    DRP_IFR_AUC(h) = trapz(FI_DRP.IFR_ave{1,h}');
    DRP_rheobase(h) = FI_DRP.rheobase(h);
    DRP_latency(h, :) = FI_DRP.lat{1,h}';
    DRP_width(h, :) = FI_DRP.width_median{1,h}';
    DRP_Vth(h, :) = FI_DRP.V_th_ave{1,h}';
    DRP_ADP(h, :) = FI_DRP.adp_index{1,h}';
    DRP_AP_peak_ave(h,:)= FI_DRP.AP_peak_ave{1,h}';
    DRP_spike_count(h,:)= FI_DRP.spike_count{1,h}';
    DRP_AP_peak_1st(h,:)= FI_DRP.AP_peak_1st{1,h}';
end

DRP_all_cell_dv_sec= [];

for h = 1:size(FI_DRP.dV_sec,2)
    
    if mean(FI_DRP.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP.Rin{1,h}) < 90000000 || mean(FI_DRP.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
    end
    
    DRP_dv_sec_cell = FI_DRP.dV_sec{1,h};
    
    DRP_all_trace_dv_sec = [];
    for trace_i = 1:length(DRP_dv_sec_cell)
        DRP_all_trace_dv_sec = [DRP_all_trace_dv_sec, DRP_dv_sec_cell{trace_i}'];
    end
    DRP_all_cell_dv_sec(end+1) = mean(DRP_all_trace_dv_sec, 'omitnan');
end

count_EVP = 0;
for c = 1:length(FI_EVP.filename)
    
    if mean(FI_EVP.Ra{1,c}) > 20000000
        continue
    elseif mean(FI_EVP.Rin{1,c}) < 90000000 || mean(FI_EVP.Rin{1,c}) > 350000000
%         continue
%    elseif mean(FI_EVP.Cp{1,c}) < 0.0000000001
         continue
    end
    
    count_EVP = count_EVP + 1;
    
    EVP_IFR_ave(c, :) = FI_EVP.IFR_ave{1,c}';
    EVP_MFR(c, :) = FI_EVP.MFR{1,c}';
    EVP_IFR_AUC(c) = trapz(FI_EVP.IFR_ave{1,c}');
    EVP_rheobase(c) = FI_EVP.rheobase(c);
    EVP_latency(c, :) = FI_EVP.lat{1,c}';
    EVP_width(c, :) = FI_EVP.width_median{1,c}';
    EVP_Vth(c, :) = FI_EVP.V_th_ave{1,c}';
    EVP_ADP(c, :) = FI_EVP.adp_index{1,c}';
    EVP_AP_peak_ave(c, :) = FI_EVP.AP_peak_ave{1,c}';
    EVP_spike_count(c, :) = FI_EVP.spike_count{1,c}';
    EVP_AP_peak_1st(c, :)= FI_EVP.AP_peak_1st{1,c}';
end

EVP_all_cell_dv_sec= [];

for c = 1:size(FI_EVP.dV_sec,2)
    
%     if mean(FI_EVP.Ra{1,c}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
%         continue
% %         %exclude cells with an input resistance lower 90Mohms and higher
% %         %than 350Mohms
%     elseif mean(FI_EVP.Rin{1,c}) < 90000000 || mean(FI_EVP.Rin{1,c}) > 350000000 
%          continue
% % %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
% %     %younger cells 
% %     %potentially memebrane voltage
% %         continue
%     end
    
    EVP_dv_sec_cell = FI_EVP.dV_sec{1,c};
    
    EVP_all_trace_dv_sec = [];
    for trace_k = 1:length(DRP_dv_sec_cell)
        EVP_all_trace_dv_sec = [EVP_all_trace_dv_sec, EVP_dv_sec_cell{trace_k}'];
    end
    EVP_all_cell_dv_sec(end+1) = mean(EVP_all_trace_dv_sec, 'omitnan');
end

count_DRP_WT = 0;
for h = 1:length(FI_DRP_WT.filename)
    
    %Exclusion criteria for FI
    if mean(FI_DRP_WT.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
        %exclude cells with an input resistance lower 90Mohms and higher
        %than 35Mohms
    elseif mean(FI_DRP_WT.Rin{1,h}) < 90000000 || mean(FI_DRP_WT.Rin{1,h}) > 350000000 
%         continue
%    elseif mean(FI_DRP_WT.Cp{1,h}) < 0.0000000001 %capacitance is lower for
    %younger cells 
        continue
    end
    
    count_DRP_WT = count_DRP_WT + 1;
    
    DRP_IFR_ave_WT(h, :) = FI_DRP_WT.IFR_ave{1,h}';
    DRP_MFR_WT(h, :) =  FI_DRP_WT.MFR{1,h}';
    DRP_IFR_AUC_WT(h) = trapz(FI_DRP_WT.IFR_ave{1,h}');
    DRP_rheobase_WT(h) = FI_DRP_WT.rheobase(h);
    DRP_latency_WT(h, :) = FI_DRP_WT.lat{1,h}';
    DRP_width_WT(h, :) = FI_DRP_WT.width_median{1,h}';
    DRP_Vth_WT(h, :) = FI_DRP_WT.V_th_ave{1,h}';
    DRP_ADP_WT(h, :) = FI_DRP_WT.adp_index{1,h}';
    DRP_AP_peak_ave_WT(h, :) = FI_DRP_WT.AP_peak_ave{1,h}';
    DRP_spike_count_WT(h, :) = FI_DRP_WT.spike_count{1,h}';
    DRP_AP_peak_1st_WT(h, :) = FI_DRP_WT.AP_peak_1st{1,h}'; 
end


DRP_all_cell_dv_sec_WT= [];

for h = 1:size(FI_DRP_WT.dV_sec,2)
    
%     if mean(FI_DRP_WT.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
%         continue
% %         %exclude cells with an input resistance lower 90Mohms and higher
% %         %than 350Mohms
%     elseif mean(FI_DRP_WT.Rin{1,h}) < 90000000 || mean(FI_DRP_WT.Rin{1,h}) > 350000000 
%          continue
% % %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
% %     %younger cells 
% %     %potentially memebrane voltage
% %         continue
%     end
    
    DRP_dv_sec_cell_WT = FI_DRP_WT.dV_sec{1,h};
    
    DRP_all_trace_dv_sec_WT = [];
    for trace_i = 1:length(DRP_dv_sec_cell)
        DRP_all_trace_dv_sec_WT = [DRP_all_trace_dv_sec_WT, DRP_dv_sec_cell_WT{trace_i}'];
    end
    DRP_all_cell_dv_sec_WT(end+1) = mean(DRP_all_trace_dv_sec_WT, 'omitnan');
end

count_EVP_WT = 0;
for c = 1:length(FI_EVP_WT.filename)
    
    if mean(FI_EVP_WT.Ra{1,c}) > 20000000
        continue
    elseif mean(FI_EVP_WT.Rin{1,c}) < 90000000 || mean(FI_EVP_WT.Rin{1,c}) > 350000000
        continue
%     elseif mean(FI_EVP_WT.Cp{1,c}) < 0.0000000001
%          continue
    end
    
    count_EVP_WT = count_EVP_WT + 1;
    
    EVP_IFR_ave_WT(c, :) = FI_EVP_WT.IFR_ave{1,c}';
    EVP_MFR_WT(c, :) = FI_EVP_WT.MFR{1,c}';
    EVP_IFR_AUC_WT(c) = trapz(FI_EVP_WT.IFR_ave{1,c}');
    EVP_rheobase_WT(c) = FI_EVP_WT.rheobase(c);
    EVP_latency_WT(c, :) = FI_EVP_WT.lat{1,c}';
    EVP_width_WT(c, :) = FI_EVP_WT.width_median{1,c}';
    EVP_Vth_WT(c, :) = FI_EVP_WT.V_th_ave{1,c}';
    EVP_ADP_WT(c, :) = FI_EVP_WT.adp_index{1,c}';
    EVP_AP_peak_ave_WT(c, :) = FI_EVP_WT.AP_peak_ave{1,c}';
    EVP_spike_count_WT(c, :) = FI_EVP_WT.spike_count{1,c}';
    EVP_AP_peak_1st_WT(c, :) = FI_EVP_WT.AP_peak_1st{1, c}'; 
end

EVP_all_cell_dv_sec_WT= [];

for c = 1:size(FI_EVP_WT.dV_sec,2)
    
%     if mean(FI_EVP_WT.Ra{1,c}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
%         continue
% %         %exclude cells with an input resistance lower 90Mohms and higher
% %         %than 350Mohms
%     elseif mean(FI_EVP_WT.Rin{1,c}) < 90000000 || mean(FI_EVP_WT.Rin{1,c}) > 350000000 
%          continue
% % %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
% %     %younger cells 
% %     %potentially memebrane voltage
% %         continue
%     end
    
    EVP_dv_sec_cell_WT = FI_EVP_WT.dV_sec{1,c};
    
    EVP_all_trace_dv_sec_WT = [];
    for trace_k = 1:length(DRP_dv_sec_cell_WT)
        EVP_all_trace_dv_sec_WT = [EVP_all_trace_dv_sec_WT, EVP_dv_sec_cell_WT{trace_k}'];
    end
    EVP_all_cell_dv_sec_WT(end+1) = mean(EVP_all_trace_dv_sec_WT, 'omitnan');
end


i_inj = FI_DRP_WT.curr_inj{1,1};





%Calculates error bar for IFR
mean_DRP_IFR = nanmean(DRP_IFR_ave, 1);
sem_DRP_IFR = nanstd(DRP_IFR_ave)./sqrt(count_DRP-1);

mean_EVP_IFR = nanmean(EVP_IFR_ave, 1);
sem_EVP_IFR = nanstd(EVP_IFR_ave)./sqrt(count_EVP-1);

mean_DRP_IFR_WT = nanmean(DRP_IFR_ave_WT, 1);
sem_DRP_IFR_WT = nanstd(DRP_IFR_ave_WT)./sqrt(count_DRP_WT-1);

mean_EVP_IFR_WT = nanmean(EVP_IFR_ave_WT, 1);
sem_EVP_IFR_WT = nanstd(EVP_IFR_ave_WT)./sqrt(count_EVP_WT-1);



%Calculate error bars for latency 
mean_DRP_latency = nanmean(DRP_latency, 1);
sem_DRP_latency = nanstd(DRP_latency)./sqrt(count_DRP-1);
mean_EVP_latency = nanmean(EVP_latency, 1);
sem_EVP_latency = nanstd(EVP_latency)./sqrt(count_EVP-1);

mean_DRP_latency_WT = nanmean(DRP_latency_WT, 1);
sem_DRP_latency_WT = nanstd(DRP_latency_WT)./sqrt(count_DRP_WT-1);
mean_EVP_latency_WT = nanmean(EVP_latency_WT, 1);
sem_EVP_latency_WT = nanstd(EVP_latency_WT)./sqrt(count_EVP_WT-1);

%Calculate error bars for mean firing rate
mean_DRP_MFR = nanmean(DRP_MFR, 1);
sem_DRP_MFR = nanstd(DRP_MFR)./sqrt(count_DRP-1);
mean_EVP_MFR = nanmean(EVP_MFR, 1);
sem_EVP_MFR = nanstd(EVP_MFR)./sqrt(count_EVP-1);

mean_DRP_MFR_WT = nanmean(DRP_MFR_WT, 1);
sem_DRP_MFR_WT = nanstd(DRP_MFR_WT)./sqrt(count_DRP_WT-1); %sqrt(count_DRN-1);
mean_EVP_MFR_WT = nanmean(EVP_MFR_WT, 1);
sem_EVP_MFR_WT = nanstd(EVP_MFR_WT)./sqrt(count_EVP_WT-1); %sqrt(count_EVN-1);


%Error bars for the width amplitudes
mean_DRP_width = nanmean(DRP_width, 1);
sem_DRP_width = nanstd(DRP_width)./sqrt(count_DRP-1);
mean_EVP_width = nanmean(EVP_width, 1);
sem_EVP_width = nanstd(EVP_width)./sqrt(count_EVP-1);

mean_DRP_width_WT = nanmean(DRP_width_WT, 1);
sem_DRP_width_WT = nanstd(DRP_width_WT)./sqrt(count_DRP_WT-1);
mean_EVP_width_WT = nanmean(EVP_width_WT, 1);
sem_EVP_width_WT = nanstd(EVP_width_WT)./sqrt(count_EVP_WT-1);

%V threshold 
mean_DRP_Vth = nanmean(DRP_Vth, 1);
sem_DRP_Vth = nanstd(DRP_Vth)./sqrt(count_DRP-1);
mean_EVP_Vth = nanmean(EVP_Vth, 1);
sem_EVP_Vth = nanstd(EVP_Vth)./sqrt(count_EVP-1);

mean_DRP_Vth_WT= nanmean(DRP_Vth_WT, 1);
sem_DRP_Vth_WT = nanstd(DRP_Vth_WT)./sqrt(count_DRP_WT-1);
mean_EVP_Vth_WT = nanmean(EVP_Vth_WT, 1);
sem_EVP_Vth_WT = nanstd(EVP_Vth_WT)./sqrt(count_EVP_WT-1);

% Area  of something lol
mean_DRP_ADP = nanmean(DRP_ADP, 1);
sem_DRP_ADP = nanstd(DRP_ADP)./sqrt(count_DRP-1);
mean_EVP_ADP = nanmean(EVP_ADP, 1);
sem_EVP_ADP = nanstd(EVP_ADP)./sqrt(count_EVP_WT-1);
mean_DRP_ADP_WT = nanmean(DRP_ADP_WT, 1);
sem_DRP_ADP_WT = nanstd(DRP_ADP_WT)./sqrt(count_DRP_WT-1);
mean_EVP_ADP_WT = nanmean(EVP_ADP_WT, 1);
sem_EVP_ADP_WT = nanstd(EVP_ADP_WT)./sqrt(count_EVP_WT-1);

if length(DRP_IFR_AUC) > length(EVP_IFR_AUC)
    add = length(DRP_IFR_AUC) - length(EVP_IFR_AUC);
    EVP_IFR_AUC(end+1:end+add) = NaN;
elseif length(DRP_IFR_AUC) < length(EVP_IFR_AUC)
    add = length(EVP_IFR_AUC) - length(DRP_IFR_AUC);
    DRP_IFR_AUC(end+1:end+add) = NaN;
end

if length(DRP_rheobase) > length(EVP_rheobase)
    add = length(DRP_rheobase) - length(EVP_rheobase);
    EVP_rheobase(end+1:end+add) = NaN;
elseif length(DRP_rheobase) < length(EVP_rheobase)
    add = length(EVP_rheobase) - length(DRP_rheobase);
    DRP_rheobase(end+1:end+add) = NaN;
end

if length(DRP_IFR_AUC_WT) > length(EVP_IFR_AUC_WT)
    add = length(DRP_IFR_AUC_WT) - length(EVP_IFR_AUC_WT);
    EVP_IFR_AUC_WT(end+1:end+add) = NaN;
elseif length(DRP_IFR_AUC_WT) < length(EVP_IFR_AUC_WT)
    add = length(EVP_IFR_AUC_WT) - length(DRP_IFR_AUC_WT);
    DRP_IFR_AUC_WT(end+1:end+add) = NaN;
end

if length(DRP_rheobase_WT) > length(EVP_rheobase_WT)
    add = length(DRP_rheobase_WT) - length(EVP_rheobase_WT);
    EVP_rheobase_WT(end+1:end+add) = NaN;
elseif length(DRP_rheobase_WT) < length(EVP_rheobase_WT)
    add = length(EVP_rheobase_WT) - length(DRP_rheobase_WT);
    DRP_rheobase_WT(end+1:end+add) = NaN;
end

%%
mean_DRP_AP_peak_ave = nanmean(DRP_AP_peak_ave, 1);
sem_DRP_AP_peak_ave = nanstd(DRP_AP_peak_ave)./sqrt(count_DRP-1);

mean_EVP_AP_peak_ave = nanmean(EVP_AP_peak_ave, 1);
sem_EVP_AP_peak_ave = nanstd(EVP_AP_peak_ave)./sqrt(count_EVP-1);

mean_DRP_AP_peak_ave_WT = nanmean(DRP_AP_peak_ave_WT, 1);
sem_DRP_AP_peak_ave_WT = nanstd(DRP_AP_peak_ave_WT)./sqrt(count_DRP_WT-1);

mean_EVP_AP_peak_ave_WT = nanmean(EVP_AP_peak_ave_WT, 1);
sem_EVP_AP_peak_ave_WT = nanstd(EVP_AP_peak_ave_WT)./sqrt(count_EVP_WT-1);

data_toplot_2= [mean_EVP_AP_peak_ave_WT', mean_DRP_AP_peak_ave_WT', mean_EVP_AP_peak_ave',  mean_DRP_AP_peak_ave'];

UnivarScatter(data_toplot_2, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','box','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('AP peak average (mv)','FontSize', 20, 'LineWidth', 15)
ylim([25 45])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);

%%
%dv/dt
%Plots mean of dv/dt 

pair_WT= padmat(EVP_all_cell_dv_sec_WT', DRP_all_cell_dv_sec_WT', 2);
data_pair_het= padmat( pair_WT,EVP_all_cell_dv_sec', 2);
data_toplot_4= padmat(data_pair_het, DRP_all_cell_dv_sec',2);


UnivarScatter(data_toplot_4, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','box','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('dV/dt','FontSize', 20, 'LineWidth', 15)
%ylim([-12 8])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);


%%
%Plots the spike counts
mean_DRP_spike_count = nanmean(DRP_spike_count, 1);
sem_DRP_spike_count = nanstd(DRP_spike_count)./sqrt(count_DRP-1);

mean_EVP_spike_count= nanmean(EVP_spike_count, 1);
sem_EVP_spike_count = nanstd(EVP_spike_count)./sqrt(count_EVP-1);

mean_DRP_spike_count_WT = nanmean(DRP_spike_count_WT, 1);
sem_DRP_spike_count_WT = nanstd(DRP_spike_count_WT)./sqrt(count_DRP_WT-1);

mean_EVP_spike_count_WT = nanmean(EVP_spike_count_WT, 1);
sem_EVP_spike_count_WT = nanstd(EVP_spike_count_WT)./sqrt(count_EVP_WT-1);

data_toplot_5= [mean_EVP_spike_count_WT', mean_DRP_spike_count_WT', mean_EVP_spike_count', mean_DRP_spike_count']; 

UnivarScatter_ATP(data_toplot_5, 'BoxType','SEM',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','lines','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('Average Spike Frequency (count)','FontSize', 20, 'LineWidth', 15)
ylim([0 30])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);



%%

mean_DRP_spike_count = nanmean(DRP_spike_count, 1);
sem_DRP_spike_count = nanstd(DRP_spike_count)./sqrt(count_DRP-1);

mean_EVP_spike_count= nanmean(EVP_spike_count, 1);
sem_EVP_spike_count = nanstd(EVP_spike_count)./sqrt(count_EVP-1);

mean_DRP_spike_count_WT = nanmean(DRP_spike_count_WT, 1);
sem_DRP_spike_count_WT = nanstd(DRP_spike_count_WT)./sqrt(count_DRP_WT-1);

mean_EVP_spike_count_WT = nanmean(EVP_spike_count_WT, 1);
sem_EVP_spike_count_WT = nanstd(EVP_spike_count_WT)./sqrt(count_EVP_WT-1);


figure_coords = [200,800, 250, 400];

%Location in the graph 

deviation = 0.35;
WT_control_spike_count_jitter = 1 + rand(1,length(EVP_spike_count_WT)).*deviation - deviation/2;
WT_DR_spike_count_jitter = 2 + rand(1,length(DRP_spike_count_WT)).*deviation - deviation/2;

SCN2A_control_spike_count_jitter = 3 + rand(1,length(EVP_spike_count)).*deviation - deviation/2;
SCN2A_DR_spike_count_jitter = 4 + rand(1,length(DRP_spike_count)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(WT_control_spike_count_jitter, EVP_spike_count_WT, 'o', 'Color', black);
hold on
plot(WT_DR_spike_count_jitter, DRP_spike_count_WT, 'o', 'Color', green);
hold on 
plot(SCN2A_control_spike_count_jitter, EVP_spike_count, 'o', 'Color', grayColor);
hold on
plot(SCN2A_DR_spike_count_jitter, DRP_spike_count, 'o', 'Color', blue);
hold on
bar(1,mean_EVP_spike_count_WT ,'EdgeColor',black,'LineWidth',4,'FaceAlpha',0.0);
hold on 
bar(2,mean_DRP_spike_count_WT,'EdgeColor',green,'LineWidth',4,'FaceAlpha',0.0);
hold on 
bar(3,mean_EVP_spike_count,'EdgeColor',grayColor,'LineWidth',4,'FaceAlpha',0.0);
hold on
bar(4,mean_DRP_spike_count,'EdgeColor',blue,'LineWidth',4,'FaceAlpha',0.0);
hold on 
errorbar(1, nanmean(WT_control_AP_rise), WT_control_AP_rise_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
hold on 
errorbar(2, nanmean(WT_DR_AP_rise), WT_DR_AP_rise_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
hold on 
errorbar(3, nanmean(SCN2A_control_AP_rise), SCN2A_control_AP_rise_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
hold on 
errorbar(4, nanmean(SCN2A_DR_AP_rise), SCN2A_DR_AP_rise_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
ylabel('1st Action potential peak (amp) ')
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'}, 'FontSize', 20)





%%

grayColor = [.5 .5 .5];
burgande= [.70 .130 .80];
blue= [.30 .74 .73];
green=[0.4660 0.6740 0.1880];

figure()
hold on
errorbar(i_inj, mean_EVP_IFR, sem_EVP_IFR, 'Color', grayColor, 'LineWidth', 4); %sem_EVP_IFR,
errorbar(i_inj, mean_DRP_IFR, sem_DRP_IFR, 'Color', blue, 'LineWidth', 4) %sem_DRP_IFR, 
errorbar(i_inj, mean_EVP_IFR_WT, sem_EVP_IFR_WT, 'Color', 'k', 'LineWidth', 4); %sem_EVP_IFR,
errorbar(i_inj, mean_DRP_IFR_WT, sem_DRP_IFR_WT, 'Color', green, 'LineWidth', 4) %sem_DRP_IFR,
set(gca, 'box', 'off', 'Fontsize', 20)
%legend('WT(WT)', 'WT(DREADDs)', 'SCN2A (WT)', 'SCN2A (DREADDs)')
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)') %'SCN2A (DR+) es el segundo'
xlabel('Current Injection (pA)')
ylabel('IFR')

%% Mean Firing rate 

figure()
hold on
errorbar(i_inj, mean_DRP_MFR, sem_DRP_MFR, 'Color', blue, 'LineWidth', 4)
errorbar(i_inj, mean_EVP_MFR, sem_EVP_MFR, 'Color', grayColor, 'LineWidth', 4); %sem_EVP_IFR,
errorbar(i_inj, mean_DRP_MFR_WT, sem_DRP_MFR_WT, 'Color', green, 'LineWidth', 4)
errorbar(i_inj, mean_EVP_MFR_WT, sem_EVP_MFR_WT, 'Color', 'k', 'LineWidth', 4); %sem_EVP_IFR,
 
set(gca, 'box', 'off', 'Fontsize', 20)
ylim([0 40])
%legend('WT', 'DREADDs', 'Location', 'NorthWest')
%legend('SCN2A +/- (DREADDs)', 'SCN2A +/- (Control)','WT(DREADDs)', 'WT(Control)', 'Location', 'northeastoutside') %'SCN2A (DR+) es el segundo'
xlabel('Current Injection (pA)')
ylabel('Mean Firing Rate (Hz)')

%Separated plots for SCN2A animals 
figure()
hold on
errorbar(i_inj, mean_DRP_MFR, sem_DRP_MFR, 'Color', blue, 'LineWidth', 4)
errorbar(i_inj, mean_EVP_MFR, sem_EVP_MFR, 'Color', grayColor, 'LineWidth', 4); %sem_EVP_IFR,
set(gca, 'box', 'off', 'Fontsize', 20)
ylim([0 40])
xlabel('Current Injection (pA)')
ylabel('Mean Firing Rate (Hz)')
legend ('DREADDs', 'Control', 'Location', 'northwest')
title('SCN2A +/-')


figure()
hold on
errorbar(i_inj, mean_DRP_MFR_WT, sem_DRP_MFR_WT, 'Color', green, 'LineWidth', 4)
errorbar(i_inj, mean_EVP_MFR_WT, sem_EVP_MFR_WT, 'Color', 'k', 'LineWidth', 4);
set(gca, 'box', 'off', 'Fontsize', 20)
ylim([0 40])
xlabel('Current Injection (pA)')
ylabel('Mean Firing Rate (Hz)')
legend ('DREADDs', 'Control', 'Location', 'northwest')
title('WildType')


figure()
hold on
errorbar(i_inj, mean_EVP_MFR, sem_EVP_MFR, 'Color', grayColor, 'LineWidth', 4);
errorbar(i_inj, mean_EVP_MFR_WT, sem_EVP_MFR_WT, 'Color', 'k', 'LineWidth', 4);
set(gca, 'box', 'off', 'Fontsize', 20)
ylim([0 40])
xlabel('Current Injection (pA)')
ylabel('Mean Firing Rate (Hz)')
legend ('SCN2A+/-', 'WildType', 'Location', 'northwest')
title('Baseline Excitability')



%%
%if I want it to be sparse I need to plot the jitter 


% data_toplot = padmat(mean_Vth_WT, mean_Vth_SCN2A, 1);
%data_toplot = padmat(data_toplot,cond3_amp_avgs', 2);

data_toplot= [mean_EVP_Vth_WT', mean_DRP_Vth_WT', mean_EVP_Vth',  mean_DRP_Vth'];

grayColor = [.5 .5 .5];
burgande= [.70 .130 .80];
blue= [.30 .74 .73];
green=[0.4660 0.6740 0.1880];


grayColor1 = [.80 .80 .80];
burgande1= [.70 .130 .80];
blue1= [.45 .89 .90];
green1=[0.59 0.80 0.31];
black= [0.2 0.2 0.2];
Colors1 = [black; green1; grayColor1; blue1];
Colors = [0 0 0; green; grayColor; blue];

UnivarScatter_ATP(data_toplot, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','box','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('Threshold (mV)','FontSize', 20, 'LineWidth', 15)
ylim([-45 -25])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);


%%
%Action potential peak 1st 

mean_DRP_AP_peak_1st = nanmean(DRP_AP_peak_1st, 1);
sem_DRP_AP_peak_1st = nanstd(DRP_AP_peak_1st)./sqrt(count_DRP-1);

mean_EVP_AP_peak_1st = nanmean(EVP_AP_peak_1st, 1);
sem_EVP_AP_peak_1st = nanstd(EVP_AP_peak_1st)./sqrt(count_EVP-1);

mean_DRP_AP_peak_1st_WT = nanmean(DRP_AP_peak_1st_WT, 1);
sem_DRP_AP_peak_1st_WT = nanstd(DRP_AP_peak_1st_WT)./sqrt(count_DRP_WT-1);

mean_EVP_AP_peak_1st_WT = nanmean(EVP_AP_peak_1st_WT, 1);
sem_EVP_AP_peak_1st_WT = nanstd(EVP_AP_peak_1st_WT)./sqrt(count_EVP_WT-1);



data_toplot_3= [mean_EVP_AP_peak_1st_WT', mean_DRP_AP_peak_1st_WT', mean_EVP_AP_peak_1st',  mean_DRP_AP_peak_1st'];

UnivarScatter_ATP(data_toplot_3, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','box','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('1st Action potential (mv)','FontSize', 20, 'LineWidth', 15)
ylim([25 45])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);



%%

figure()
boxplot([mean_EVP_Vth_WT',  mean_DRP_Vth_WT', mean_EVP_Vth',  mean_DRP_Vth'],'BoxStyle', 'outline', 'Colors',... 
[0.5, 0.5, 0.5]);
hold on
% boxplot(mean_DRP_Vth_WT, 'BoxStyle', 'filled', 'Colors', green);
% hold on

for i = 1:length(mean_DRP_Vth_WT)
    plot(1, mean_EVP_Vth_WT(i), 'ko', 'MarkerSize', 10)
    plot(2, mean_DRP_Vth_WT(i), 'go', 'MarkerSize', 10)
    plot(3, mean_EVP_Vth(i), 'ko', 'MarkerSize', 10)
    plot(4, mean_DRP_Vth(i), 'bo', 'MarkerSize', 10)
end
xticks([1 2 3 4])
set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'Control',...
    'DREADDs', 'Control','DREADDs'})
ylim([-42 -26])
ylabel('Average Threeshold (mV)')



figure()
boxplot([mean_EVP_Vth_WT',  mean_DRP_Vth_WT'],'BoxStyle', 'outline', 'Colors',... 
[0.5, 0.5, 0.5]);
hold on
% boxplot(mean_DRP_Vth_WT, 'BoxStyle', 'filled', 'Colors', green);
% hold on

for i = 1:length(mean_DRP_Vth_WT)
    plot(1, mean_EVP_Vth_WT(i), 'ko', 'MarkerSize', 10)
    plot(2, mean_DRP_Vth_WT(i), 'go', 'MarkerSize', 10)
end
xticks([1 2])
set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'Control',...
    'DREADDs'})
ylim([-42 -26])
ylabel('Average Threeshold (mV)')


figure()
boxplot([mean_EVP_Vth',  mean_DRP_Vth'],'BoxStyle', 'outline', 'Colors',... 
[0.5, 0.5, 0.5]);
hold on
% boxplot(mean_DRP_Vth_WT, 'BoxStyle', 'filled', 'Colors', green);
% hold on

for i = 1:length(mean_DRP_Vth)
    plot(1, mean_EVP_Vth(i), 'ko', 'MarkerSize', 10)
    plot(2, mean_DRP_Vth(i), 'bo', 'MarkerSize', 10)
end
xticks([1 2])
set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'Control',...
    'DREADDs'})
ylim([-42 -26])
ylabel('Average Threeshold (mV)')



%%
figure()
boxplot([mean_EVP_AP_peak_ave_WT',  mean_DRP_AP_peak_ave_WT', mean_EVP_AP_peak_ave',  mean_DRP_AP_peak_ave'],'BoxStyle', 'outline', 'Colors',... 
[0.5, 0.5, 0.5]);
hold on
% boxplot(mean_DRP_Vth_WT, 'BoxStyle', 'filled', 'Colors', green);
% hold on

for i = 1:length(mean_DRP_AP_peak_ave_WT)
    plot(1, mean_EVP_AP_peak_ave_WT(i), 'ko', 'MarkerSize', 10)
    plot(2, mean_DRP_AP_peak_ave_WT(i), 'go', 'MarkerSize', 10)
    plot(3, mean_EVP_AP_peak_ave(i), 'ko', 'MarkerSize', 10)
    plot(4, mean_DRP_AP_peak_ave(i), 'bo', 'MarkerSize', 10)
end
xticks([1 2 3 4])
set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'Control',...
    'DREADDs', 'Control','DREADDs'})
ylim([30 45])
ylabel('Action potential peak ave. (amp)')


%%







%%
%plot spike count 


mean_DRP_spike_count = nanmean(DRP_spike_count, 1);
sem_DRP_spike_count = nanstd(DRP_spike_count)./sqrt(count_DRP-1);

mean_EVP_spike_count= nanmean(EVP_spike_count, 1);
sem_EVP_spike_count = nanstd(EVP_spike_count)./sqrt(count_EVP-1);

mean_DRP_spike_count_WT = nanmean(DRP_spike_count_WT, 1);
sem_DRP_spike_count_WT = nanstd(DRP_spike_count_WT)./sqrt(count_DRP_WT-1);

mean_EVP_spike_count_WT = nanmean(EVP_spike_count_WT, 1);
sem_EVP_spike_count_WT = nanstd(EVP_spike_count_WT)./sqrt(count_EVP_WT-1);

% deviation = 0.35;
% EVP_spike_count_jitter = 2 + rand(1,length(EVP_spike_count)).*deviation - deviation/2);
% DRP_spike_count_jitter = 1 + rand(1,length(DRP_spike_count)).*deviation - deviation/2);
% DRP_spike_count_WT_jitter = 4 + rand(1,length(DRP_spike_count_WT)).*deviation - deviation/2);
% EVP_spike_count_WT_jitter = 3 + rand(1,length(EVP_spike_count_WT).*deviation - deviation/2);

grayColor = [.5 .5 .5];
burgande= [.70 .130 .80];
blue= [.30 .74 .73];
green=[0.4660 0.6740 0.1880];

% plot(EVP_spike_count_jitter, EVP_spike_count, 'o', 'Color', Colors(2,:));
% hold on;
% plot(DRP_spike_count_jitter, DRP_spike_count, 'o', 'Color', Colors(1,:));
% plot(DRP_spike_count_WT_jitter , DRP_spike_count_WTa, 'o', 'Color', Colors(3,:));
% plot(EVP_spike_count_WT_jitter, EVP_spike_count_WT, 'o', 'Color', Colors(4,:));
%%
% figure();
% hold on
% bar(nanmean(mean_DRP_spike_count'))
% hold on 
% bar(nanmean(mean_EVP_spike_count'))
% hold on 
% xticks([1 2 3 4])

% 
% bar(nanmean(mean_spike_count_WT'))
% bar(nanmean(mean_spike_count_WT'))
% hold on 
%%
% errorbar(1, nanmean(DRP_spike_count), sem_DRP_spike_count, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
% errorbar(2, nanmean(EVP_spike_count), sem_EVP_spike_count, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
% errorbar(3, nanmean(EVP_spike_count_WT), sem_EVP_spike_count_WT, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
% errorbar(4, nanmean(DRP_spike_count_WT), sem_DRP_spike_count_WT, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
% box off
% ylabel('')
% xticks([1 2 3 4])
% set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'Control',...
%     'DREADDs', 'Control','DREADDs'})


%%



% figure()
% boxplot([mean_EVP_spike_count_WT',  mean_DRP_spike_count_WT', mean_EVP_spike_count',  mean_DRP_spike_count'],'BoxStyle', 'outline', 'Colors',... 
% [0.5, 0.5, 0.5]);
% hold on
% % boxplot(mean_DRP_Vth_WT, 'BoxStyle', 'filled', 'Colors', green);
% % hold on
% 
% for i = 1:length(mean_DRP_spike_count_WT)
%     plot(1, mean_EVP_spike_count_WT(i), 'ko', 'MarkerSize', 10)
%     plot(2, mean_DRP_spike_count_WT(i), 'go', 'MarkerSize', 10)
%     plot(3, mean_EVP_spike_count(i), 'ko', 'MarkerSize', 10)
%     plot(4, mean_DRP_spike_count(i), 'bo', 'MarkerSize', 10)
% end
% xticks([1 2 3 4])
% set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'Control',...
%     'DREADDs', 'Control','DREADDs'})
% ylim([])
% ylabel('Action potential peak ave. (amp)')




%%
% plot latency 

figure()
hold on
errorbar(i_inj, mean_EVP_latency, sem_EVP_latency, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_latency, sem_DRP_latency, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_EVP_latency_WT, sem_EVP_latency_WT, 'Color', 'm', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_latency_WT, sem_DRP_latency_WT, 'Color', 'g', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)') 
xlabel('Current Injection (pA)')
ylabel('Latency')

% plot AP width  (median)

figure()
hold on
errorbar(i_inj, mean_EVP_width, sem_EVP_width, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_width, sem_DRP_width, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_EVP_width_WT, sem_EVP_width_WT, 'Color', 'g', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_width_WT, sem_DRP_width_WT, 'Color', 'm', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)') 
xlabel('Current Injection (pA)')
ylabel('AP Width (median)')

% plot V threshold (average) 

figure()
hold on
errorbar(i_inj, mean_EVP_Vth, sem_EVP_Vth, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_Vth, sem_DRP_Vth, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_EVP_Vth_WT, sem_EVP_Vth_WT, 'Color', 'g', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_Vth_WT, sem_DRP_Vth_WT, 'Color', 'm', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)') 
xlabel('Current Injection (pA)')
ylabel('V threshold (average)')

% plot adaptation index 

figure()
hold on
errorbar(i_inj, mean_EVP_ADP, sem_EVP_ADP, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_ADP, sem_DRP_ADP, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_EVP_ADP_WT, sem_EVP_ADP_WT, 'Color', 'g', 'LineWidth', 3)
errorbar(i_inj, mean_DRP_ADP_WT, sem_DRP_ADP_WT, 'Color', 'm', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)', 'WT (EV+)', 'WT(DR+)')
xlabel('Current Injection (pA)')
ylabel('Adaptation Index')





% plot AUC
% add jitter!
% figure(6)
% boxplot([EVP_IFR_AUC, DRP_IFR_AUC, EVP_IFR_AUC_WT, DRP_IFR_AUC_WT],...
%     'BoxStyle', 'filled', 'Colors', [0.5, 0.5, 0.5]);
% hold on
% 
% 
% for i = 1:length(EVP_IFR_AUC)
%     plot(1, EVP_IFR_AUC(i), 'ko', 'MarkerSize', 10)
%     plot(2, DRP_IFR_AUC(i), 'bo', 'MarkerSize', 10)
%     plot(3, EVP_IFR_AUC_WT(i), 'go', 'MarkerSize', 10)
%     plot(4, DRP_IFR_AUC_WT(i), 'mo', 'MarkerSize', 10)
% end
% set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'SCN2A (EV+)',...
%     'SCN2A (DR+)', 'WT (EV+)', 'WT(DR+)'})

% ylabel('IFR AUC')

% plot rheobase
% add jitter!
% figure()
% boxplot([EVP_rheobase, DRP_rheobase, EVP_rheobase_WT, DRP_rheobase_WT],...
%     'BoxStyle', 'filled', 'Colors', [0.5, 0.5, 0.5])
% hold on
% 
% for i = 1:length(EVP_rheobase)
%     plot(1, EVP_rheobase(i), 'ko', 'MarkerSize', 10)
%     plot(2, DRP_rheobase(i), 'bo', 'MarkerSize', 10)
%     plot(3, EVP_rheobase_WT(i), 'go', 'MarkerSize', 10)
%     plot(4, DRP_rheobase_WT(i), 'mo', 'MarkerSize', 10)
% end
% 
% set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'SCN2A (EV+)',...
%     'SCN2A (DR+)', 'WT (EV+)', 'WT(DR+)'})
% ylabel('Rheobase')

