%% Select which cell and which trace to plot

% where to save the selected trace data
save_data_path = '/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/';

%'C:\Users\Andre\Documents\_Spring_2022\Patching_PC\representative_traces\fI\';

% name of the saved data
save_data_name = 'plot_trial';

% analyzed fI results folder
fp_data = '/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/Only_Analyzed/';

%'C:\Users\Andre\Documents\_Spring_2022\Patching_PC\analyzed_FI_results';

% experimental conditions
cond = {'EV+','DR+'};

% Choose f_I analysis mat file
experiment = { 'EVPthrough081022to100522','DRPthrough081022to100522'};

% Choose a cell in each experiment for plotting
cell_num = [10 5];

% Plotting mode (0 for non-normarlized traces, 1 for normalized traces) 
plot_mode = 0;

% current step
stim = 180;

% current increment
cur_inc = 20;

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor = [0 0 0.5];

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

% color1 = '#E88BE2';
% color2 = '#78BEEB';
color1 = fun([0 0 0],tint_factor(1)); %black
% color2 = fun([150 150 150],tint_factor(2)); %gray
%color2 = fun([70 130 180],tint_factor(2)); %steel blue
color2 = fun([70 130 180],tint_factor(3)); %steel blue+50% tint

colorcode = cat(1, color1, color2);
%% data readout
% Variable initializations
cell_data = cell(1,numel(experiment));
trace_data = cell(1,numel(experiment));
trace_id = NaN(1,numel(experiment));
rheo_id = NaN(1,numel(experiment));


for fi = 1:numel(experiment)
    
%    
    curr_file = matfile(strcat(fp_data,'/',experiment{fi}, '.mat'));
    curr_file_FI_data = curr_file.FI_DATA;
    
    curr_file_data= curr_file_FI_data.aDAT;   %add this variable cause of the properties 
    curr_cell_id = curr_file_FI_data.cell_id;
    curr_rheobase_ind = curr_file_FI_data.rheobase_ind;
    
    cell_data{1,fi} = curr_file_data{1,cell_num(fi)};
    rheo_id(1,fi) = curr_rheobase_ind(fi,1);
    
    if plot_mode == 0
        trace_id(1,fi) = curr_cell_id{1,cell_num(fi)-1}(1,1)+stim/cur_inc-1;
    elseif plot_mode == 1
        trace_id(1,fi) = curr_cell_id{1,cell_num(fi)-1}(1,1)+rheo_id(1,fi)-1+stim/cur_inc;
    end
    
    trace_data{1,fi} = cell_data{1,fi}(:,trace_id(1,fi));
end

%% Plotting 

plot_range = 5000:25000;
figure('position',[4816 14802 1000 490])
hold on
for ci = 1:numel(cell_num)
    plot(trace_data{1,ci}(plot_range,1),'Color',colorcode(ci,1:3),'LineWidth',3)
end

%draw scale
plot([17500; 19500],[0; 0], '-k',[17500;17500],[0; 20], '-k', 'LineWidth',2)
text(17400, 5, '20 mV', 'HorizontalAlignment', 'right')
text(18500, -5, '0.2 s', 'HorizontalAlignment', 'center')

legend(cond)

if plot_mode == 0
    title(strcat(num2str(stim), 'pA injected'))
elseif plot_mode == 1
    title(strcat(num2str(stim),'pA over rheobase'))
end

hold off

%% Drawing current step

figure('position',[56 200 1000 490])
bl_before = zeros(1,1999);
current_step = repmat(stim/1000, 1, 5000);
bl_after = zeros(1,1999);

current_all = horzcat(bl_before,current_step,bl_after);

plot(current_all,'k', 'LineWidth',3)
ylim([-0.2 1])
    
%% save data
cd(save_data_path)

save(save_data_name, 'cond','experiment','cell_num','trace_data')
