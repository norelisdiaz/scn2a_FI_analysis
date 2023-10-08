% plot aspects of FI curves from .mat files from Wei's code

%To plot an specific condition uncomment the two groups you will like to
%regroup. 

%WT Animals (All conditions plotted as one)
% FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/DRPthrough062222to063022.mat');
% FI_CNOonly= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/EVPthrough062222to063022.mat');
% FI_Hm4di_Control= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/DRNthrough062222to063022.mat');
% FI_CNOonly_Control= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/EVNthrough062222to063022.mat');


%SCN2A Animals (All conditions plotted as one)
%FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/DRPthrough081022to100522.mat'); %FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/EVPthrough081022to100522.mat');
FI_Hm4di_Control = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/DRNthrough081022to100522.mat');
FI_CNOonly_Control = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/EVNthrough081022to100522.mat');

%Sorted by Age for WT animals
%FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023/sortedbyage/DRP_p27_101722.mat'); 
%FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023/sortedbyage/EVP_p27_101722.mat');

%FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023/sortedbyage/DRP_p27_101722.mat'); 
%FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023/sortedbyage/EVP_p27_101722.mat');

%SCN2A vs WT
%FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/nuevoexp2023/Only_Analyzed/DRPthrou081022to030923.mat'); 
%FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/DRNthrough081022to100522.mat');
%FI_Hm4di_Control= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/DRPthrough062222to063022.mat');
%FI_CNOonly_Control= load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Only_Analyzed/EVPthrough062222to063022.mat');
%FI_Hm4di_Control = load('/Users/norelis/Documents/MATLAB/FI_DATA/DRP_22junthrough30jun.mat');
%FI_CNOonly_Control= load('/Users/norelis/Documents/MATLAB/FI_DATA/DRN_22junthrough30jun.mat');

FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Mat_data_update/EVP_WT_062222to022123_alldata.mat'); 
FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/olddata_DRP_WT_062222to063022.mat');


FI_Hm4di = FI_Hm4di.FI_DATA; %DR+
FI_Hm4di_Control = FI_Hm4di_Control.FI_DATA; %DR-
FI_CNOonly = FI_CNOonly.FI_DATA; %EV+
FI_CNOonly_Control = FI_CNOonly_Control.FI_DATA; %EV-


Hm4di_MFR = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_MFR = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
Hm4di_MFR_Control = nan(length(FI_Hm4di_Control.filename), length(FI_Hm4di_Control.curr_inj{1,1}));
CNOonly_MFR_Control = nan(length(FI_CNOonly_Control.filename), length(FI_CNOonly_Control.curr_inj{1,1}));

Hm4di_IFR_ave = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_IFR_ave = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
Hm4di_IFR_ave_Control = nan(length(FI_Hm4di_Control.filename), length(FI_Hm4di_Control.curr_inj{1,1}));
CNOonly_IFR_ave_Control = nan(length(FI_CNOonly_Control.filename), length(FI_CNOonly_Control.curr_inj{1,1}));

Hm4di_IFR_AUC = nan(length(FI_Hm4di.filename), 1);
CNOonly_IFR_AUC = nan(length(FI_CNOonly.filename), 1);
Hm4di_IFR_AUC_Control = nan(length(FI_Hm4di_Control.filename), 1);
CNOonly_IFR_AUC_Control = nan(length(FI_CNOonly_Control.filename), 1);


Hm4di_rheobase = nan(length(FI_Hm4di.filename), 1);
CNOonly_rheobase = nan(length(FI_CNOonly.filename), 1);
Hm4di_rheobase_Control = nan(length(FI_Hm4di_Control.filename), 1);
CNOonly_rheobase_Control = nan(length(FI_CNOonly_Control.filename), 1);


Hm4di_latency = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_latency = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
Hm4di_latency_Control = nan(length(FI_Hm4di_Control.filename), length(FI_Hm4di_Control.curr_inj{1,1}));
CNOonly_latency_Control = nan(length(FI_CNOonly_Control.filename), length(FI_CNOonly_Control.curr_inj{1,1}));

Hm4di_width = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_width = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
Hm4di_width_Control = nan(length(FI_Hm4di_Control.filename), length(FI_Hm4di_Control.curr_inj{1,1}));
CNOonly_width_Control = nan(length(FI_CNOonly_Control.filename), length(FI_CNOonly_Control.curr_inj{1,1}));

Hm4di_Vth = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_Vth = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
Hm4di_Vth_Control = nan(length(FI_Hm4di_Control.filename), length(FI_Hm4di_Control.curr_inj{1,1}));
CNOonly_Vth_Control = nan(length(FI_CNOonly_Control.filename), length(FI_CNOonly_Control.curr_inj{1,1}));


Hm4di_ADP = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_ADP = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
Hm4di_ADP_Control = nan(length(FI_Hm4di_Control.filename), length(FI_Hm4di_Control.curr_inj{1,1}));
CNOonly_ADP_Control = nan(length(FI_CNOonly_Control.filename), length(FI_CNOonly_Control.curr_inj{1,1}));

count_Hm4di = 0;
for h = 1:length(FI_Hm4di.filename)
    
    %Exclusion criteria for FI
    if mean(FI_Hm4di.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
        %exclude cells with an input resistance lower 90Mohms and higher
        %than 35Mohms
    elseif mean(FI_Hm4di.Rin{1,h}) < 90000000 || mean(FI_Hm4di.Rin{1,h}) > 350000000 
        continue
 %  elseif mean(FI_Hm4di.Cp{1,h}) < 0.0000000001 %capacitance is lower for
    %younger cells 
  %      continue
    end
    
    count_Hm4di = count_Hm4di + 1;
    
    Hm4di_IFR_ave(h, :) = FI_Hm4di.IFR_ave{1,h}';
    Hm4di_MFR(h, :) = FI_Hm4di.MFR{1,h}';
    Hm4di_IFR_AUC(h) = trapz(FI_Hm4di.IFR_ave{1,h}');
    Hm4di_rheobase(h) = FI_Hm4di.rheobase(h);
    Hm4di_latency(h, :) = FI_Hm4di.lat{1,h}';
    Hm4di_width(h, :) = FI_Hm4di.width_median{1,h}';
    Hm4di_Vth(h, :) = FI_Hm4di.V_th_ave{1,h}';
    Hm4di_ADP(h, :) = FI_Hm4di.adp_index{1,h}';
    
end

count_CNOonly = 0;
for c = 1:length(FI_CNOonly.filename)
    
    if mean(FI_CNOonly.Ra{1,c}) > 20000000
        continue
    elseif mean(FI_CNOonly.Rin{1,c}) < 90000000 || mean(FI_CNOonly.Rin{1,c}) > 350000000
        continue
  % elseif mean(FI_CNOonly.Cp{1,c}) < 0.0000000001
   %      continue
    end
    
    count_CNOonly = count_CNOonly + 1;
    
    CNOonly_IFR_ave(c, :) = FI_CNOonly.IFR_ave{1,c}';
    CNOonly_MFR(c, :) = FI_CNOonly.MFR{1,c}';
    CNOonly_IFR_AUC(c) = trapz(FI_CNOonly.IFR_ave{1,c}');
    CNOonly_rheobase(c) = FI_CNOonly.rheobase(c);
    CNOonly_latency(c, :) = FI_CNOonly.lat{1,c}';
    CNOonly_width(c, :) = FI_CNOonly.width_median{1,c}';
    CNOonly_Vth(c, :) = FI_CNOonly.V_th_ave{1,c}';
    CNOonly_ADP(c, :) = FI_CNOonly.adp_index{1,c}';
    
end

count_Hm4di_Control = 0;
for h = 1:length(FI_Hm4di_Control.filename)
    
    %Exclusion criteria for FI
    if mean(FI_Hm4di_Control.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
        %exclude cells with an input resistance lower 90Mohms and higher
        %than 35Mohms
    elseif mean(FI_Hm4di_Control.Rin{1,h}) < 90000000 || mean(FI_Hm4di_Control.Rin{1,h}) > 350000000 
        continue
  % elseif mean(FI_Hm4di_Control.Cp{1,h}) < 0.0000000001 %capacitance is lower for
    %younger cells 
   %     continue
    end
    
    count_Hm4di_Control = count_Hm4di_Control + 1;
    
    Hm4di_IFR_ave_Control(h, :) = FI_Hm4di_Control.IFR_ave{1,h}';
    Hm4di_MFR_Control(h, :) = FI_Hm4di_Control.MFR{1,h}';
    Hm4di_IFR_AUC_Control(h) = trapz(FI_Hm4di_Control.IFR_ave{1,h}');
    Hm4di_rheobase_Control(h) = FI_Hm4di_Control.rheobase(h);
    Hm4di_latency_Control(h, :) = FI_Hm4di_Control.lat{1,h}';
    Hm4di_width_Control(h, :) = FI_Hm4di_Control.width_median{1,h}';
    Hm4di_Vth_Control(h, :) = FI_Hm4di_Control.V_th_ave{1,h}';
    Hm4di_ADP_Control(h, :) = FI_Hm4di_Control.adp_index{1,h}';
    
end 

count_CNOonly_Control = 0;
for c = 1:length(FI_CNOonly_Control.filename)
    
    if mean(FI_CNOonly_Control.Ra{1,c}) > 20000000
        continue
    elseif mean(FI_CNOonly_Control.Rin{1,c}) < 90000000 || mean(FI_CNOonly_Control.Rin{1,c}) > 350000000
        continue
    %elseif mean(FI_CNOonly_Control.Cp{1,c}) < 0.0000000001
      %   continue
    end
    
    count_CNOonly_Control = count_CNOonly_Control + 1;
    
    CNOonly_IFR_ave_Control(c, :) = FI_CNOonly_Control.IFR_ave{1,c}';
    CNOonly_MFR_Control(c, :) = FI_CNOonly_Control.MFR{1,c}';
    CNOonly_IFR_AUC_Control(c) = trapz(FI_CNOonly_Control.IFR_ave{1,c}');
    CNOonly_rheobase_Control(c) = FI_CNOonly_Control.rheobase(c);
    CNOonly_latency_Control(c, :) = FI_CNOonly_Control.lat{1,c}';
    CNOonly_width_Control(c, :) = FI_CNOonly_Control.width_median{1,c}';
    CNOonly_Vth_Control(c, :) = FI_CNOonly_Control.V_th_ave{1,c}';
    CNOonly_ADP_Control(c, :) = FI_CNOonly_Control.adp_index{1,c}';
    
end



i_inj = FI_Hm4di.curr_inj{1,1};

%Calculate error bars for IFR ave 
mean_Hm4di_IFR = nanmean(Hm4di_IFR_ave, 1);
sem_Hm4di_IFR = nanstd(Hm4di_IFR_ave)./sqrt(count_Hm4di-1);
mean_CNOonly_IFR = nanmean(CNOonly_IFR_ave, 1);
sem_CNOonly_IFR = nanstd(CNOonly_IFR_ave)./sqrt(count_CNOonly-1);

mean_Hm4di_IFR_Control = nanmean(Hm4di_IFR_ave_Control, 1);
sem_Hm4di_IFR_Control = nanstd(Hm4di_IFR_ave_Control)./sqrt(count_Hm4di_Control-1);
mean_CNOonly_IFR_Control = nanmean(CNOonly_IFR_ave_Control, 1);
sem_CNOonly_IFR_Control = nanstd(CNOonly_IFR_ave_Control)./sqrt(count_CNOonly_Control-1);

%Calculate error bars for mean FR 
mean_Hm4di_MFR = nanmean(Hm4di_MFR, 1);
sem_Hm4di_MFR = nanstd(Hm4di_MFR)./sqrt(count_Hm4di-1);
mean_CNOonly_MFR = nanmean(CNOonly_MFR, 1);
sem_CNOonly_MFR = nanstd(CNOonly_MFR)./sqrt(count_CNOonly-1);

mean_Hm4di_MFR_Control = nanmean(Hm4di_MFR_Control, 1);
sem_Hm4di_MFR_Control = nanstd(Hm4di_MFR_Control)./sqrt(count_Hm4di_Control-1);
mean_CNOonly_MFR_Control = nanmean(CNOonly_MFR_Control, 1);
sem_CNOonly_MFR_Control = nanstd(CNOonly_MFR_Control)./sqrt(count_CNOonly_Control-1);

%Calculate error bars for latency 
mean_Hm4di_latency = nanmean(Hm4di_latency, 1);
sem_Hm4di_latency = nanstd(Hm4di_latency)./sqrt(count_Hm4di-1);
mean_CNOonly_latency = nanmean(CNOonly_latency, 1);
sem_CNOonly_latency = nanstd(CNOonly_latency)./sqrt(count_CNOonly-1);

mean_Hm4di_latency_Control = nanmean(Hm4di_latency_Control, 1);
sem_Hm4di_latency_Control = nanstd(Hm4di_latency_Control)./sqrt(count_Hm4di_Control-1);
mean_CNOonly_latency_Control = nanmean(CNOonly_latency_Control, 1);
sem_CNOonly_latency_Control = nanstd(CNOonly_latency_Control)./sqrt(count_CNOonly_Control-1);

%Error bars for the width amplitudes
mean_Hm4di_width = nanmean(Hm4di_width, 1);
sem_Hm4di_width = nanstd(Hm4di_width)./sqrt(count_Hm4di-1);
mean_CNOonly_width = nanmean(CNOonly_width, 1);
sem_CNOonly_width = nanstd(CNOonly_width)./sqrt(count_CNOonly-1);


mean_Hm4di_width_Control = nanmean(Hm4di_width, 1);
sem_Hm4di_width_Control = nanstd(Hm4di_width)./sqrt(count_Hm4di-1);
mean_CNOonly_width_Control = nanmean(CNOonly_width, 1);
sem_CNOonly_width_Control = nanstd(CNOonly_width)./sqrt(count_CNOonly-1);

%V threshold 
mean_Hm4di_Vth = nanmean(Hm4di_Vth, 1);
sem_Hm4di_Vth = nanstd(Hm4di_Vth)./sqrt(count_Hm4di-1);
mean_CNOonly_Vth = nanmean(CNOonly_Vth, 1);
sem_CNOonly_Vth = nanstd(CNOonly_Vth)./sqrt(count_CNOonly-1);

mean_Hm4di_Vth_Control= nanmean(Hm4di_Vth_Control, 1);
sem_Hm4di_Vth_Control = nanstd(Hm4di_Vth_Control)./sqrt(count_Hm4di_Control-1);
mean_CNOonly_Vth_Control = nanmean(CNOonly_Vth_Control, 1);
sem_CNOonly_Vth_Control = nanstd(CNOonly_Vth_Control)./sqrt(count_CNOonly_Control-1);

% Area  of something lol
mean_Hm4di_ADP = nanmean(Hm4di_ADP, 1);
sem_Hm4di_ADP = nanstd(Hm4di_ADP)./sqrt(count_Hm4di-1);
mean_CNOonly_ADP = nanmean(CNOonly_ADP, 1);
sem_CNOonly_ADP = nanstd(CNOonly_ADP)./sqrt(count_CNOonly_Control-1);
mean_Hm4di_ADP_Control = nanmean(Hm4di_ADP_Control, 1);
sem_Hm4di_ADP_Control = nanstd(Hm4di_ADP_Control)./sqrt(count_Hm4di_Control-1);
mean_CNOonly_ADP_Control = nanmean(CNOonly_ADP_Control, 1);
sem_CNOonly_ADP_Control = nanstd(CNOonly_ADP_Control)./sqrt(count_CNOonly_Control-1);

if length(Hm4di_IFR_AUC) > length(CNOonly_IFR_AUC)
    add = length(Hm4di_IFR_AUC) - length(CNOonly_IFR_AUC);
    CNOonly_IFR_AUC(end+1:end+add) = NaN;
elseif length(Hm4di_IFR_AUC) < length(CNOonly_IFR_AUC)
    add = length(CNOonly_IFR_AUC) - length(Hm4di_IFR_AUC);
    Hm4di_IFR_AUC(end+1:end+add) = NaN;
end

if length(Hm4di_rheobase) > length(CNOonly_rheobase)
    add = length(Hm4di_rheobase) - length(CNOonly_rheobase);
    CNOonly_rheobase(end+1:end+add) = NaN;
elseif length(Hm4di_rheobase) < length(CNOonly_rheobase)
    add = length(CNOonly_rheobase) - length(Hm4di_rheobase);
    Hm4di_rheobase(end+1:end+add) = NaN;
end

% plot it all!
% ADD Y LIMS!!


% plot FI curve
%grayColor = [.7 .7 .7];


%%
figure(1)
hold on
errorbar(i_inj, mean_CNOonly_IFR, sem_CNOonly_IFR, 'Color', 'k', 'LineWidth', 3); %sem_CNOonly_IFR,
errorbar(i_inj, mean_Hm4di_IFR, sem_Hm4di_IFR, 'Color', 'b', 'LineWidth', 3) %sem_Hm4di_IFR, 
errorbar(i_inj, mean_CNOonly_IFR_Control, sem_CNOonly_IFR_Control, 'Color', 'g', 'LineWidth', 3); %sem_CNOonly_IFR,
errorbar(i_inj, mean_Hm4di_IFR_Control, sem_Hm4di_IFR_Control, 'Color', 'm', 'LineWidth', 3) %sem_Hm4di_IFR,
set(gca, 'box', 'off', 'Fontsize', 20)
legend('newdata', 'olddata','WT(EV+)', 'WT(DR+)') %'SCN2A (DR+)
%legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)') %'SCN2A (DR+) es el segundo'
xlabel('Current Injection (pA)')
ylabel('IFR (Hz)')

% plot mean firing rate 
figure(2)
hold on
errorbar(i_inj, mean_CNOonly_MFR, sem_CNOonly_MFR, 'Color', 'k', 'LineWidth', 3); %sem_CNOonly_IFR,
errorbar(i_inj, mean_Hm4di_MFR, sem_Hm4di_MFR, 'Color', 'b', 'LineWidth', 3) %sem_Hm4di_IFR, 
errorbar(i_inj, mean_CNOonly_MFR_Control, sem_CNOonly_MFR_Control, 'Color', 'g', 'LineWidth', 3); %sem_CNOonly_IFR,
errorbar(i_inj, mean_Hm4di_MFR_Control, sem_Hm4di_MFR_Control, 'Color', 'm', 'LineWidth', 3) %sem_Hm4di_IFR,
set(gca, 'box', 'off', 'Fontsize', 20)
legend('DREADDS', 'Empty Vector','WT(EV+)', 'WT(DR+)')
%legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)') %'SCN2A (DR+) es el segundo'
xlabel('Current Injection (pA)')
ylabel('MFR (Hz)')

%SCN2A
figure()
hold on
errorbar(i_inj, mean_CNOonly_MFR, sem_CNOonly_MFR, 'Color', 'k', 'LineWidth', 3); %sem_CNOonly_IFR,
errorbar(i_inj, mean_Hm4di_MFR, sem_Hm4di_MFR, 'Color', 'b', 'LineWidth', 3) %sem_Hm4di_IFR, 
set(gca, 'box', 'off', 'Fontsize', 20)
legend('Control negative','DREADDs Expressing Cells')%'SCN2A (DR+) es el segundo'
xlabel('Current Injection (pA)')
ylabel('MFR (Hz)')
title('SCN2A DREADDs Control')

%WT
figure(4)
hold on
errorbar(i_inj, mean_CNOonly_MFR_Control, sem_CNOonly_MFR_Control, 'Color', 'g', 'LineWidth', 3); %sem_CNOonly_IFR,
errorbar(i_inj, mean_Hm4di_MFR_Control, sem_Hm4di_MFR_Control, 'Color', 'm', 'LineWidth', 3) %sem_Hm4di_IFR,
set(gca, 'box', 'off', 'Fontsize', 20)
legend('Control negative', 'DREADDs expressing cells')%'SCN2A (DR+) es el segundo'
xlabel('Current Injection (pA)')
ylabel('MFR (Hz)')
title('')
% plot latency 

figure()
hold on
errorbar(i_inj, mean_CNOonly_latency, sem_CNOonly_latency, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_latency, sem_Hm4di_latency, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_CNOonly_latency_Control, sem_CNOonly_latency_Control, 'Color', 'm', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_latency_Control, sem_Hm4di_latency_Control, 'Color', 'g', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)')
xlabel('Current Injection (pA)')
ylabel('Latency')

% plot AP width  (median)

figure(4)
hold on
errorbar(i_inj, mean_CNOonly_width, sem_CNOonly_width, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_width, sem_Hm4di_width, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_CNOonly_width_Control, sem_CNOonly_width_Control, 'Color', 'g', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_width_Control, sem_Hm4di_width_Control, 'Color', 'm', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)')
xlabel('Current Injection (pA)')
ylabel('AP Width (median)')

% plot V threshold (average) 

figure(5)
hold on
errorbar(i_inj, mean_CNOonly_Vth, sem_CNOonly_Vth, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_Vth, sem_Hm4di_Vth, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_CNOonly_Vth_Control, sem_CNOonly_Vth_Control, 'Color', 'g', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_Vth_Control, sem_Hm4di_Vth_Control, 'Color', 'm', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)','WT(EV+)', 'WT(DR+)')
xlabel('Current Injection (pA)')
ylabel('V threshold (average)')

% plot adaptation index 

figure(6)
hold on
errorbar(i_inj, mean_CNOonly_ADP, sem_CNOonly_ADP, 'Color', 'k', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_ADP, sem_Hm4di_ADP, 'Color', 'b', 'LineWidth', 3)
errorbar(i_inj, mean_CNOonly_ADP_Control, sem_CNOonly_ADP_Control, 'Color', 'g', 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_ADP_Control, sem_Hm4di_ADP_Control, 'Color', 'm', 'LineWidth', 3)
set(gca, 'box', 'off', 'Fontsize', 20)
legend('SCN2A (EV+)', 'SCN2A (DR+)', 'WT (EV+)', 'WT(DR+)')
xlabel('Current Injection (pA)')
ylabel('Adaptation Index')



%%
% plot AUC
% add jitter!
figure(7)
boxplot([CNOonly_IFR_AUC, Hm4di_IFR_AUC,CNOonly_IFR_AUC_Control,  Hm4di_IFR_AUC_Control ], 'BoxStyle', 'filled', 'Colors', [0.5, 0.5, 0.5])
hold on

for i = 1:length(CNOonly_IFR_AUC)
    plot(1, CNOonly_IFR_AUC(i), 'ko', 'MarkerSize', 10)
    plot(2, Hm4di_IFR_AUC(i), 'bo', 'MarkerSize', 10)
    plot(3, CNOonly_IFR_AUC_Control(i), 'go', 'MarkerSize', 10)
    plot(4, Hm4di_IFR_AUC_Control(i), 'mo', 'MarkerSize', 10)l
end

set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'CNO only', 'Hm4di'})
ylabel('IFR AUC')

% plot rheobase
% add jitter!
figure(8)
boxplot([CNOonly_rheobase, Hm4di_rheobase, CNOonly_rheobase_Control, Hm4di_rheobase_Control(i)], 'BoxStyle', 'filled', 'Colors', [0.5, 0.5, 0.5])
hold on

for i = 1:length(CNOonly_rheobase)
    plot(1, CNOonly_rheobase(i), 'ko', 'MarkerSize', 10)
    plot(2, Hm4di_rheobase(i), 'go', 'MarkerSize', 10)
    plot(3, CNOonly_rheobase_Control(i), 'go', 'MarkerSize', 10)
    plot(4, Hm4di_rheobase_Control(i), 'mo', 'MarkerSize', 10) 
end

set(gca, 'box', 'off', 'Fontsize', 20, 'XTickLabel', {'CNO only', 'Hm4di'})
ylabel('Rheobase')

