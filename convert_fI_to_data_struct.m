% This script loops through the analyzed_fI_results folder and groups
% them by experimental conditions: 
% For example: (these indice should match 'Cat' in the cell_id excel file) 
    % DR+CNO- 1
    % CNO- 2
    % NT- 3
    % NT+saline- 4

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  
%% %% Change these accordingly based on how you want to group the data
% in each anaylzed_fI_results file, dates can be extracted from the last six digits

%save results
save_results = 1;

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'updated_WT_AA_DD';

%file name
saved_file_name = 'all_cond.mat';

%'AA_DD_WT_TTX_analysis_updated.mat';

%where to save grouped files
fp_grouped_data = ...,
    'C:\Users\Andre\Documents\_Spring_2022\Patching_PC\FI_data_by_groups_updated\';

%experimental conditions
exp_con = {'UT','UT_TTX','AA_UN','AA_TTX','DD_UN','DD_TTX','WT_UN','WT_TTX'};

%import cell_id_index table 
cd('C:\Users\Andre\Documents\_Spring_2022\Patching_PC\')
cell_id_index = readtable('NEW_AG_cell_id_index.xlsx','Sheet',exp_name);

%location of analyzed fI results
fp_analyzed_fI = ...,
    'C:\Users\Andre\Documents\_Spring_2022\Patching_PC\analyzed_FI_results_updated\';

%% extract and group data
data_struct_temp = cell(1,numel(exp_con));

cd(fp_analyzed_fI)
all_files = dir;
file_num = numel(all_files)-2;

for fi = 1:file_num
    curr_name = all_files(fi+2).name;
    if ismember('m',curr_name(4:9))
        continue
    else
        curr_date = curr_name(4:9);
    end
    
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(cell_id,2)
            if isempty(cell_id{1,ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    ti_start = cell_id{1,ci}(1,1);
                    trace_ct = size(cell_id{1,ci},1);
                    ti_end = cell_id{1,ci}(trace_ct,1);
                    
                    data_struct_temp{cond_i}.MFR(1:trace_ct,cell_ID) = MFR{1,ci}(ti_start:ti_end,1);                 
                    data_struct_temp{cond_i}.IFR(1:trace_ct,cell_ID) = IFR{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.mean_IFR(1:trace_ct,cell_ID) = IFR_ave{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Threshold(1:trace_ct,cell_ID) = V_th_1st{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.ADP_ind(1:trace_ct,cell_ID) = adp_index{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Lat(1:trace_ct,cell_ID) = lat{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Udratio(1:trace_ct,cell_ID) = udratio_median{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Width_first(1:trace_ct,cell_ID) = width_1st{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.AP_peak_1st(1:trace_ct,cell_ID) = AP_peak_1st{1,ci}(ti_start:ti_end,1);
                    data_struct_temp{cond_i}.Rheobase(cell_ID,1) = rheobase(ci,1);
                    data_struct_temp{cond_i}.Rheobase_ind(cell_ID,1) = rheobase_ind(ci,1);
                    data_struct_temp{cond_i}.Vm(cell_ID,1) = cell_stats{1,1}.Vm_stats(ci,1);
                    
                    %for Rin and Cm
                    if isnan(cell_stats{1,1}.Rin_stats(ci,5))                      
                        data_struct_temp{cond_i}.Rin(cell_ID,1) = cell_stats{1,1}.Rin_stats(ci,1);
                    else
                        data_struct_temp{cond_i}.Rin(cell_ID,1) = cell_stats{1,1}.Rin_stats(ci,5);
                    end
                    
                    if isnan(cell_stats{1,1}.Cm_stats(ci,5))                      
                        data_struct_temp{cond_i}.Cm(cell_ID,1) = cell_stats{1,1}.Cm_stats(ci,1);
                    else
                        data_struct_temp{cond_i}.Cm(cell_ID,1) = cell_stats{1,1}.Cm_stats(ci,5);
                    end                    
                        
                    %current injection vectprs should be saved separately
                    %for each cell since they might differ
                    data_struct_temp{cond_i}.curr_inj(1:trace_ct,cell_ID) = curr_inj{1,ci}(ti_start:ti_end,1);
                end     
            end
        end
    end
end
 
%%%% CHANGE FOR EACH RUN %%%%
UT = data_struct_temp{1,1};
UT_TTX = data_struct_temp{1,2};
AA_UN = data_struct_temp{1,3};
AA_TTX = data_struct_temp{1,4};
DD_UN = data_struct_temp{1,5};
DD_TTX = data_struct_temp{1,6};
WT_UN = data_struct_temp{1,7};
WT_TTX = data_struct_temp{1,8};
% WT_CNO_B = data_struct_temp{1,4};

%% save grouped data
if save_results == 1
    cd(fp_grouped_data)
    save(saved_file_name,exp_con{1},exp_con{2},exp_con{3},exp_con{4},exp_con{5},exp_con{6},exp_con{7},exp_con{8},'exp_con','cell_id_index')
end
