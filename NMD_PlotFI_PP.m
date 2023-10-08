% PLOT PASSIVE PROPERTIES FROM FI DATA

%FI_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/DRPthrough081022to081122.mat');
%FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/DRP_22junthrough30jun.mat');

% FI_Hm4di =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/olddata_DRP_WT_062222to063022.mat');
% FI_CNOonly =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/Mat_data_update/EVP_WT_062222to022123_alldata.mat');
% %FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/EVPthrough081022to081122.mat'); % Hm4di in V1b and BF!
% FI_BF_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023_h5_files/DRN_WT_062222to022123_alldata.mat');
% FI_BF_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023_h5_files/EVN_WT_062222to022123_alldata.mat');

FI_Hm4di=load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/WT_DRP_June22_22toMay4_23.mat');% WT DREADDs
%FI_Hm4di =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/olddata_DRP_WT_062222to063022.mat'); % WT DREADDs

FI_CNOonly =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/WT_mergeEVN_EVP_Jun22_22toMay4_23.mat');

%FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/nuevoexp2023_h5_files/mergeEVN_EVP_062222to021223.mat');
%FI_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/EVPthrough081022to081122.mat'); % Hm4di in V1b and BF!
FI_BF_Hm4di=load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/nuevoexp2023/SCN2A_DRP_Aug10_22toMay22_23');
%FI_BF_Hm4di = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/nuevoexp2023/Only_Analyzed/DRPthrou081022to030923.mat');
FI_BF_CNOonly= load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/SCN2AmergeEVN_EVP_June22_22toMay22_23.mat');
%FI_BF_CNOonly = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/SCN2AmergeEVN_EVP_062222to021223.mat');

%FI_BF_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/090822thru091322_m1_3shRNA_Hm4di_FI_DATA.mat'); % Hm4di in V1m and BF!
%FI_BF_CNOonly = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/090822thru091322_m1_3shRNA_CNOonly_FI_DATA.mat'); % Hm4di in V1m and BF!



% CCh Response cells

%%%%%% would like to adjust code to check these cells too, but currently doesn't
%%%%%% work, since data is formatted slightly differently
%FI_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/ALL_ctrl_CCh_Analysis.mat'); % ctrl
%FI_BF_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/m1_1sh_CCh_Analysis.mat'); % m1_1sh
%FI_BF_CNOonly = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/m1_3sh_CCh_Analysis.mat'); % m1_3sh



FI_Hm4di = FI_Hm4di.FI_DATA;
FI_CNOonly = FI_CNOonly.FI_DATA;
FI_BF_Hm4di = FI_BF_Hm4di.FI_DATA;
FI_BF_CNOonly = FI_BF_CNOonly.FI_DATA;

CNOonly_data_mat = NaN(length(FI_CNOonly.filename), 5);
count = 0;
for num = 1:length(FI_CNOonly.filename)
    
    if mean(FI_CNOonly.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_CNOonly.Rin{1,num}) < 90000000 || mean(FI_CNOonly.Rin{1,num}) > 350000000
        continue
%     elseif mean(FI_CNOonly.Cp{1,num}) < 0.0000000001
%         continue
    end
    
    count = count + 1;
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
    CNOonly_data_mat(num, 2) = nanmean(FI_CNOonly.Ra{1, num});
    CNOonly_data_mat(num, 3) = nanmean(FI_CNOonly.Rin{1, num});
    CNOonly_data_mat(num, 4) = nanmean(FI_CNOonly.Cp{1, num});
    CNOonly_data_mat(num, 5) = nanmean(FI_CNOonly.Vr{1, num});
end


BF_CNOonly_data_mat = NaN(length(FI_BF_CNOonly.filename), 5);
count = 0;
for num = 1:length(FI_BF_CNOonly.filename)
    
    if mean(FI_BF_CNOonly.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_BF_CNOonly.Rin{1,num}) < 90000000 || mean(FI_BF_CNOonly.Rin{1,num}) > 350000000
        continue
%     elseif mean(FI_BF_CNOonly.Cp{1,num}) < 0.0000000001
%         continue
    end
    
    count = count + 1;
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
    BF_CNOonly_data_mat(num, 2) = nanmean(FI_BF_CNOonly.Ra{1, num});
    BF_CNOonly_data_mat(num, 3) = nanmean(FI_BF_CNOonly.Rin{1, num});
    BF_CNOonly_data_mat(num, 4) = nanmean(FI_BF_CNOonly.Cp{1, num});
    BF_CNOonly_data_mat(num, 5) = nanmean(FI_BF_CNOonly.Vr{1, num});
end


Hm4di_data_mat = NaN(length(FI_Hm4di.filename), 5);
count = 0;
for num = 1:length(FI_Hm4di.filename)

    if mean(FI_Hm4di.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_Hm4di.Rin{1,num}) < 90000000 || mean(FI_Hm4di.Rin{1,num}) > 350000000
        continue
%     elseif mean(FI_Hm4di.Cp{1,num}) < 0.0000000001
%         continue
    end

    count = count + 1;
    
    Hm4di_data_mat(num, 2) = nanmean(FI_Hm4di.Ra{1, num});
    Hm4di_data_mat(num, 3) = nanmean(FI_Hm4di.Rin{1, num});
    Hm4di_data_mat(num, 4) = nanmean(FI_Hm4di.Cp{1, num});
    Hm4di_data_mat(num, 5) = nanmean(FI_Hm4di.Vr{1, num});
end

BF_Hm4di_data_mat = NaN(length(FI_BF_Hm4di.filename), 5);
count = 0;
for num = 1:length(FI_BF_Hm4di.filename)

    if mean(FI_BF_Hm4di.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_BF_Hm4di.Rin{1,num}) < 90000000 || mean(FI_BF_Hm4di.Rin{1,num}) > 350000000
        continue
%     elseif mean(FI_BF_Hm4di.Cp{1,num}) < 0.0000000001
%         continue
    end

    count = count + 1;
    
    BF_Hm4di_data_mat(num, 2) = nanmean(FI_BF_Hm4di.Ra{1, num});
    BF_Hm4di_data_mat(num, 3) = nanmean(FI_BF_Hm4di.Rin{1, num});
    BF_Hm4di_data_mat(num, 4) = nanmean(FI_BF_Hm4di.Cp{1, num});
    BF_Hm4di_data_mat(num, 5) = nanmean(FI_BF_Hm4di.Vr{1, num});
end


%%%%%
%PassProp plotting 
Colors = [0 0 0; .4660 .6740 .1880; .5 .5 .5; .30 .74 .73];

%0.4660 0.6740 0.1880
%.70 .130 .80

figure_coords = [200,800,250,400];

%%
%%
%VC Rin
Hm4di_prop_data = Hm4di_data_mat(:,3);
CNOonly_prop_data = CNOonly_data_mat(:,3);
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,3);
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,3);

Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor', Colors(2,:),'LineWidth',4,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor', Colors(1,:),'LineWidth',4,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor', Colors(4,:),'LineWidth',4,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor', Colors(3,:),'LineWidth',4,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('Input Resistance (Ohms)')
%set(gca,'XTickLabel',{'WT(DR+)', 'SCN2A(DR+)'})
%set(gca,'XTickLabel',{'WT(Control)', 'WT(DREADDs)', 'SCN2A(Control)', 'SCN2A(DREADDs)'}, 'FontSize', 15)
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'}, 'FontSize', 20)
%sigstar({[1,2],[1,3],[1,4], [2,3], [2,4], [3, 4]},[0.0091, nan, nan, nan, nan, nan])
%sigstar({[1,2],[1,4]},[0.0158, 0.0181])
%sigstar({[1,2],[1,3],[1,4], [2,3], [2,4], [3, 4]},[0.0091, nan, nan, nan, nan, nan])
%H = sigstar({[1,2], [1,4]},stats)

%sigstar(data_toplot, stats)


z = line([1 2],[300000000 300000000]);
z.Color = 'k';
z.LineWidth = 2;
q = text(1.4,300008000,'**'); %change the * according to the analysis  
q.FontSize = 30;
hold on; 

% 
% z2 = line([1 4],[33 33]);
% z2.Color = 'k';
% z2.LineWidth = 2;
% q2 = text(1.5,34,'**'); %change the * according to the analysis  
% q2.FontSize = 30;
% hold on; 

% hold on
% legend('WT', 'SCN2A +/-')
% xtickangle(45)
%legend(ax(1) , {'sin','cos'} ) ;

%%
% %Stastical analysis (Multi anova)
data_toplot_1 =padmat(CNOonly_prop_data(:), Hm4di_prop_data(:), 2);
data_toplot_2 = padmat(BF_CNOonly_prop_data(:), BF_Hm4di_prop_data(:), 2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);

%[p, tbl, stats] = kruskalwallis(data_toplot_2)
[p, tbl, stats] = kruskalwallis(data_toplot)

% [h,p_Rin,ci,stats] = ttest2(Hm4di_prop_data,CNOonly_prop_data)
% [p_anova_Rin, tbl, stats] = anova1(data_toplot) 
% [p_anova_Rin, tbl, stats] = anova1(data_toplot_1); 
Rin_stats = multcompare(stats)


%%
%%
%VC Cp
Hm4di_prop_data = Hm4di_data_mat(:,4);
CNOonly_prop_data = CNOonly_data_mat(:,4);
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,4);
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,4);

Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor',Colors(2,:),'LineWidth',4,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor',Colors(1,:),'LineWidth',4,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor',Colors(4,:),'LineWidth',4,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor',Colors(3,:),'LineWidth',4,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('Capacitance (F)')
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'}, 'FontSize', 20)
%set(gca,'XTickLabel',{'WT(EV+)', 'WT(DR+)', 'SCN2A (DR+)', 'SCN2A (EV+)'}, 'FontSize', 15)
%set(gca,'XTickLabel',{'WT(DR+)', 'SCN2A(DR+)'})
%set(gca,'XTickLabel',{'WT(DR+)', 'WT(EV+)', 'WT (DR-)', 'WT (EV-)'}, 'FontSize', 15)
% set(gca,'XTickLabel',{'CNO only', 'Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
%xtickangle(45)


[h,p_Cp,ci,stats] = ttest2(Hm4di_prop_data, CNOonly_prop_data)

%%
%VC Vr
Hm4di_prop_data = Hm4di_data_mat(:,5).*10^3;
CNOonly_prop_data = CNOonly_data_mat(:,5).*10^3;
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,5).*10^3;
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,5).*10^3;

Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor', Colors(2,:),'LineWidth',4,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor', Colors(1,:),'LineWidth',4,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor', Colors(4,:),'LineWidth',4,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor', Colors(3,:),'LineWidth',4,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('Resting Membrane Potential (mV)')
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'}, 'FontSize', 20)
%set(gca,'XTickLabel',{'WT(EV+)', 'WT(DR+)', 'SCN2A (DR+)', 'SCN2A (EV+)'}, 'FontSize', 15)
%set(gca,'XTickLabel',{'WT(DR+)', 'SCN2A(DR+)'})
%set(gca,'XTickLabel',{'WT(DR+)', 'WT(EV+)', 'WT (DR-)', 'WT (EV-)'}, 'FontSize', 15)
%set(gca,'XTickLabel',{'CNO only','Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
%xtickangle(45)


[h,p_Vr,ci,stats] = ttest2(Hm4di_prop_data,CNOonly_prop_data)

%%
%VC Ra
Hm4di_prop_data = Hm4di_data_mat(:,2).*10^3;
CNOonly_prop_data = CNOonly_data_mat(:,2).*10^3;
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,2).*10^3;
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,2).*10^3;

Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor', Colors(2,:),'LineWidth',2,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor', Colors(1,:),'LineWidth',2,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor', Colors(4,:),'LineWidth',2,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor', Colors(3,:),'LineWidth',2,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('Access Resistance (Ohms)')
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'}, 'FontSize', 20)
%set(gca,'XTickLabel',{'WT(EV+)', 'WT(DR+)', 'SCN2A (DR+)', 'SCN2A (EV+)'}, 'FontSize', 15)
%set(gca,'XTickLabel',{'WT(DR+)', 'WT(EV+)', 'WT (DR-)', 'WT (EV-)'}, 'FontSize', 15)
%set(gca,'XTickLabel',{'WT(DR+)', 'SCN2A(DR+)'})
%set(gca,'XTickLabel',{'CNO only', 'Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
%xtickangle(45)

[h,p_Ra,ci,stats] = ttest2(Hm4di_prop_data,CNOonly_prop_data)

