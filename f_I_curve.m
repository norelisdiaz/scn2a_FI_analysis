%%%% This script takes the imported f-I mat file and calculate basic stats
%%%% for each property

%%%% Under each field, data from each experimental condition (i.e. WT,
%%%% WT_CNO, DR_CNO, and, WT_saline), are saved in corresponding columns 1, 2, 3,
%%%% and 4.

%% Initialization

%location where the mat file will be saved
fp_analyzed_data = '/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/';
%'C:\Users\Andre\Documents\_Spring_2022\Patching_PC\FI_data_by_groups_updated\';

%name of the mat file
filename = {'DRP_new.mat', 'EVP_new.mat'};

%experimental conditions (stats in each column follow this order)
cond = {'SCN2A_(DR+)', 'SCN2A_(EV+)'};      %{'UT','UT_TTX','AA_UN','AA_TTX','DD_UN','DD_TTX','WT_UN','WT_TTX'}; 

%functions for calculating stats
func = {@nanmean, @nanstd, @nansem};

%field names (corresponding to properties)
fields = {'MFR','IFR','mean_IFR'};
%{'MFR', 'IFR', 'mean_IFR', 'Threshold', 'ADP_ind', 'Rheobase',...
%    'Lat', 'Rin', 'Cm', 'Udratio', 'Width_med','Width_first'};

%plot
plot_on = 1;

%save results
save_results = 1;

%% calculate stats 

%calculate stats for each current step
struct_all = cell(1,3); %results from each function are saved in corresponding cells

for sti = 1:numel(func) % per stat
    for cdi = 1:numel(cond) % per experimental condition
        curr_c = eval(cond{cdi});
        
        for fi = 1:numel(fields) % per field
            if size(curr_c.(fields{fi}),2)<2 % ignore fields that are not current-step related
                continue
            else
                for cui = 1:size(curr_c.(fields{fi}),1) % per current step
                    struct_all{sti}.(fields{fi})(cui,cdi) = func{sti}(curr_c.(fields{fi})(cui,:));
                end
            end
        end
    end
end

%current injection vector (just pick one vector, should all be the same in
%terms of plotting
current_inj = curr_c.curr_inj(:,2);

%transfer results to corresponding data structs for storage
ave = struct_all{1};
std = struct_all{2};
sem = struct_all{3};

%calculate median for each cell (med)
% and
%get values from the rheobase step for each cell (rheo)
%get values from 100 pA over rheobase for each cell (above100)

for cdii = 1:numel(cond) % per experimental condition
        curr_c = eval(cond{cdii});
        
        for fii = 1:numel(fields) % per field
            if size(curr_c.(fields{fii}),2)<2 % ignore fields that are not current-step related
                continue
            else
                for cuii = 1:size(curr_c.(fields{fii}),2) % per cell
                    med.(fields{fii})(cuii,cdii) = nanmedian(curr_c.(fields{fii})(:,cuii));
                end
                
%                 vals = get_step_values(curr_c.(fields{fii}), 2, 0);
%                 rheo.(fields{fii})(1:size(vals,1),cdii) = vals;
%                 
%                 vals_1 = get_step_values(curr_c.(fields{fii}), 2, 100);
%                 above100.(fields{fii})(1:size(vals_1,1),cdii) = vals_1;
                
            end
        end
end
 
%current injection vector (just pick one vector, should all be the same in
%terms of plotting

%% plotting

% data range for plotting
plot_range = 1:10;
plot_variable = 'IFR';
plot_stat = 'ave';
err = 'sem';
y_label = 'IFR (Hz)';

%tint factors
tint_factor = [0 0 0 0 0 0 0 0];
%tint and convert function
%x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

%color code
color1 = fun([210 43 43],tint_factor(1)); %cadmium red
color2 = fun([248 131 121],tint_factor(2)); %coral pink
color3 = fun([70 130 180],tint_factor(3)); %steel blue
color4 = fun([153 186 221],tint_factor(4)); %carolina blue
color5 = fun([138 72 141],tint_factor(5)); %purple
%color6 = fun([255 153 01], tint_factor(6)); %pumpkin orange
color6 = fun([185 185 220],tint_factor(6)); %lavender
color7 = fun([79 121 66],tint_factor(5)); %fern green
color8 = fun([147 197 114], tint_factor(8)); %pistachio green

% 
%color2 = fun([255 203 164],tint_factor(3)); %deep peach
%color4 = fun([153 186 221],tint_factor(4)); %carolina blue

%for washin
% color1 = '#EA91C2';
% color2 = '#E64AA0';
% color3 = '#94EA91';
% color4 = '#4FE649';


curr_y = eval(strcat(plot_stat,'.',plot_variable));
curr_y_err = eval(strcat(err,'.',plot_variable));

if plot_on == 1

    figure
    hold on
    %cond1
    errorbar(current_inj(plot_range),curr_y(plot_range,1),curr_y_err(plot_range,1),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'LineWidth',2,'Color',color1)
    %cond2
    errorbar(current_inj(plot_range),curr_y(plot_range,2),curr_y_err(plot_range,2),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color2, 'MarkerFaceColor',color2,'LineWidth',2,'Color',color2)    
        
    %cond3
    errorbar(current_inj(plot_range),curr_y(plot_range,3),curr_y_err(plot_range,3),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color3, 'MarkerFaceColor','none','LineWidth',2,'Color',color3)
    %cond4
    errorbar(current_inj(plot_range),curr_y(plot_range,4),curr_y_err(plot_range,4),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color4, 'MarkerFaceColor','none','LineWidth',2,'Color',color4)
    %cond5
    errorbar(current_inj(plot_range),curr_y(plot_range,5),curr_y_err(plot_range,5),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color5, 'MarkerFaceColor','none','LineWidth',2,'Color',color5)
    %cond6
    errorbar(current_inj(plot_range),curr_y(plot_range,6),curr_y_err(plot_range,6),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color6, 'MarkerFaceColor','none','LineWidth',2,'Color',color6)
    %cond7
    errorbar(current_inj(plot_range),curr_y(plot_range,7),curr_y_err(plot_range,7),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color7, 'MarkerFaceColor','none','LineWidth',2,'Color',color7)
     %cond8
    errorbar(current_inj(plot_range),curr_y(plot_range,8),curr_y_err(plot_range,8),'-o','MarkerSize',8,...
            'MarkerEdgeColor',color8, 'MarkerFaceColor','none','LineWidth',2,'Color',color8)
% 

    %WT+saline
%     errorbar(current_inj(plot_range),curr_y(plot_range,2),curr_y_err(plot_range,2),'-o','MarkerSize',8,...
%             'MarkerEdgeColor','k', 'MarkerFaceColor','k','LineWidth',2,...
%             'Color','k')

    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 2;
    ax.YLabel.String = y_label;
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = 'Injected Current (pA)';
    ax.YLabel.FontSize = 14;
    ax.YLabel.Interpreter = 'none';
    ax.XLim = [0 250];
    ax.YLim = [0 200];

    legend('UT','UT_TTX','AA_UN','AA_TTX','DD_UN','DD_TTX','WT_UN','WT_TTX','FontSize',12,'Location','northwest','Interpreter','none','Box','off') %'AA_UN','AA_TTX','DD_UN','DD_TTX','WT_UN','WT_TTX'
    %title('Bursting')

end

%% save results
if save_results == 1
    cd (fp_analyzed_data)
    save(filename{1},'ave','std','sem','med')
end
