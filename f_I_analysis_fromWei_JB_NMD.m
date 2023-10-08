% this script is written to be run for a folder of .h5 files of evoked
% firing by current steps
% UPDATED BY JB TO RUN .MAT FILES
% Updated BY NMD to run 
%%% ****STILL NEED TO TROUBLESHOOT THIS REVISED SCRIPT!!!!***

%% these need to be changed for each new data run

%save_file_name = 'SCN2AmergeEVN_EVP_062222to052223.mat';
save_file_name = 'randomDRP_062222to021223.mat';
%save_file_name = 'EVP_WT_062222to022123_alldata.mat';
%save_file_name= 'SCM'; 

%%%%%% Name of file you would like to have data saved to

%loop through all the folders ** can only do one folder at a time - doesn't
%work for more than 1!
experiment = {'DRP'};

%filepath where you keep data folders
fp_data = '/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/';

%location of the analyzed mini data
% fp_analyzed_data = ...,
%     '/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A(fI)/';

%save results
save_results = 1;

%plot on
plot_on= 0;

sealtestend  = 1050; % in ms
sealteststart  = 0.5; % in s
sealtestlength  = 0.5; % in s
sealtest_dvolt = 0.005; % V
samplerate = 10000;

%%
close all
for jj = 1:numel(experiment)
    current_experiment = experiment{jj};
    current_folder = strcat(fp_data,current_experiment);
    cd(current_folder)
    
    all_files = dir('*.mat');
    A = {all_files.name};
    A = unique(A);
    file_num = length(A);
%     
%     
%    all_files = dir;
%    file_num = numel(all_files);


    %pre-allocate space
    %Puts names on the allocate matrix called cell
    filename = cell(1,file_num);
    cell_id = cell(1,file_num);
    cell_num = NaN(1,file_num);
    aDAT = cell(1,file_num);
    all_Ra = cell(1,file_num);
    all_Rin = cell(1,file_num);
    all_Cp = cell(1,file_num);
    all_Vr = cell(1, file_num);

    MFR = cell(1,file_num); %mean firing rate
    ISI = cell(1,file_num); %interspike interval (for each two APs)
    IFR = cell(1,file_num); %instantenous firing rate (for each ISI)
    IFR_ave = cell(1,file_num); %average of the first two IFR
    V_th = cell(1,file_num); %AP threshold (for each AP)
    V_th_ave = cell(1,file_num); %average AP threshold
    V_th_1st = cell(1,file_num); %threshold of the first AP in each trace
    V_th_4th = cell(1,file_num); %threshold of the fourth AP in each trace
    f_trough = cell(1,file_num); %AP fast trough (for each AP)
    f_trough_ave = cell(1,file_num); %average fast trough
    dV_sec = cell(1,file_num); %dV/dt for each trace
    AP_peak = cell(1,file_num); %peak amplitde of each AP
    AP_peak_ave = cell(1,file_num); %average AP peak amp
    AP_peak_1st = cell(1,file_num); %first AP amp
    Vrm = cell(1,file_num); %resting membrane potential
    Vrm_stat = NaN(file_num,4); %stats of Vrm during current injection 
    spike_count = cell(1,file_num); %number of spikes evoked
    curr_inj = cell(1,file_num); %current injected
    rheobase = NaN(file_num,1); % each row stores the rheobase of a cell
    rheobase_ind = NaN(file_num,1); %eahc row stores the index of rheobase (which stimulation)
    lat = cell(1,file_num); %latency between the start of the current step and the first AP
    adp_index = cell(1,file_num); %spike adaptation index
    udratio = cell(1,file_num); %upstroke/downstroke ratio
    udratio_slope = cell(1,file_num); %slope of the u/d ratio change across APs
    udratio_median = cell(1,file_num); %median of u/d ratio for each current step
    width = cell(1,file_num); %AP width at half height
    width_1st = cell(1,file_num); %width of the first AP for each current step
    width_slope = cell(1,file_num); %slope of the width change across APs
    width_median = cell(1,file_num); %median of width for each current step


    %%% Data readout: a cell array stores data for a certain cell, in which
    %%% each column stores a certain trace.
    curr_cell = 0;
    cellID = 0;
    FI_trace = 0;
    trace_id = 0;
    for f = 1:file_num
        num_filename = numel(filename{f});
   
            data = load([current_folder, '/', A{f}]);
             data = data.sweep0;
%          
            cellName = A{f}(1:13);
            if ~strcmp(cellName, curr_cell)
                curr_cell = cellName;
                cellID = cellID + 1 ;
                filename{cellID} = A{f}(1:13);
            end
        if length(data) == 40000
            % IS FI TRACE
            this_trace = cellName;
            disp (A{f})
            if ~strcmp(this_trace, FI_trace)
                FI_trace = this_trace;
                trace_id = 1;
            else
                trace_id = trace_id + 1;
            end
            if trace_id > 20 %cantidad de traces 
               keyboard
            end
            cell_id{1,cellID}(trace_id,1) = trace_id; %el problema es aqui maybe are the requirements
            aDAT{1,cellID}(1:length(data),trace_id) = data;
        else
            % pull passive properties!
            [Ra, Cp, Rin, Vr] = get_PassProp_VC(data.* 10^-12, sealteststart, sealtestlength, sealtest_dvolt, samplerate);
            all_Ra{1,cellID}(end+1) = Ra;
            all_Rin{1,cellID}(end+1) = Rin;
            all_Cp{1,cellID}(end+1) = Cp;
            all_Vr{1,cellID}(end+1) = Vr;
        end
        
        
        
        
        
        
        %{
        if ~strfind(filename{f}, '-')
            warning('Not a f-I trace!')
        else

            if isnan(str2double(filename{f}(6))) == 1
                cellID = str2double(filename{f}(5));
            else
                cellID = str2double(filename{f}(5:6));
            end
            %cellID = str2double(filename{f}(5));
            cell_num(f) = cellID;

            range(1) = str2double(filename{f}(num_filename - 11:num_filename - 8)); %first sweep #
            range(2) = str2double(filename{f}(num_filename - 6:num_filename - 3)); %last sweep #


            %formatting file names
            if range(1) < 10
                start = strcat('000',num2str(range(1)));
            elseif range(1) >= 10 && range(1) < 100
                start = strcat('00',num2str(range(1)));
            elseif range(1) >= 100 && range(1) < 1000
                start = strcat('0',num2str(range(1)));
            else
                start = num2str(range(1));
            end

            if range(2) < 10
                finish = strcat('000',num2str(range(2)));
            elseif range(2) >= 10 && range(2) < 100
                finish = strcat('00',num2str(range(2)));
            elseif range(2) >= 100 && range(2) < 1000
                finish = strcat('0',num2str(range(2)));
            else
                finish = num2str(range(2));
            end


            for tii = range(1):range(2)

                trace_id = tii;
                cell_id{1,cellID}(trace_id,1) = trace_id;

                if tii < 10
                    sname = strcat('000',num2str(tii));
                elseif tii >= 10 && tii < 100
                    sname = strcat('00',num2str(tii));
                elseif tii >= 100 && tii < 1000
                    sname = strcat('0',num2str(tii));
                else
                    sname = num2str(tii);
                end


                %datafile readout
                dtstruct = ws.loadDataFile(filename{f});
                prefix = strcat('dtstruct.sweep_',sname,'.analogScans');
                data = eval(prefix); 
                aDAT{1,cellID}(1:numel(data(:,1)),trace_id) = data(1:numel(data(:,1)),1);

            end
        end
        %}
    end

    %%% Replace zeroes in each cell with NaN
    for gi = 1:cellID
        cell_id{1,gi} = nonzeros(cell_id{1,gi});
        if isempty(cell_id{1,gi}) == 1
            continue
        else
            for ri = 1:cell_id{1,gi}(end,1)
                if ismember(ri,cell_id{1,gi}) == 0
                    aDAT{1,gi}(:,ri) = NaN;
                end
            end
        end
    end
%%
    %% AP detection and characterization

    for ci = 1:cellID
        %if contains(filename{ci}, 'cell_07') %deletes files 
         %   continue
        %end
        if isempty(cell_id{1,ci}) == 1
            continue
        else
            trace_start = cell_id{1,ci}(1,1);

            for ti = trace_start:cell_id{1,ci}(end,1)
                if ismember(ti,cell_id{1,ci}) == 0
                    MFR{1,ci}(ti,1) = NaN;
                    ISI{1,ci}(ti,1:20) = NaN;
                    IFR{1,ci}(ti,1:20) = NaN; 
                    V_th{1,ci}{ti} = {};
                    f_trough{1,ci}{ti} = {};
                    dV_sec{1,ci}{ti} = {};
                    AP_peak{1,ci}{ti} = NaN;
                    Vrm{1,ci}(ti,1) = NaN; 
                    spike_count{1,ci}(ti,1) = NaN;
                    curr_inj{1,ci}(ti,1) = NaN;
                    lat{1,ci}(ti,1) = NaN; 
                    adp_index{1,ci}(ti,1) = NaN;
                    udratio{1,ci}(ti,1:20) = NaN;
                    width{1,ci}(ti,1:20) = NaN;

                    continue
                else

                    data = aDAT{1,ci}(:,ti);

                    %Resting membrane potential (in mV)
                    Vrm{1,ci}(ti,1) = nanmean(data(1000:2000));

                    %current injected (in pA)
    %                 if ti-trace_start+1 > 10
    %                     trace_start = ti;
    %                 end
                    curr_inj{1,ci}(ti,1) = (ti-trace_start+1)*20;

                    %Detect AP in each trace
                    [dV_sec_temp, V_th_temp, f_trough_temp, sp_num_temp, udratio_temp] = get_vthresh(data, plot_on);


                    %number of spikes detected in each trace
                    spike_count{1,ci}(ti,1) = sp_num_temp;

                    % AP threshold
                    V_th{1,ci}{ti} = V_th_temp; %first col stores indices, second col stores potentials

                    %average AP threshold (start from the 4th AP)
                    V_th_ave{1,ci}(ti,1) = nanmean(V_th{1,ci}{1,ti}(4:end,2));
                    
                    %the first AP in each trace
                    V_th_1st{1,ci}(ti,1) = V_th{1,ci}{1,ti}(1,2);
                    
                    %the fourth AP in each trace
                    V_th_4th{1,ci}(ti,1) = V_th{1,ci}{1,ti}(4,2);


                    %AP fast trough
                    f_trough{1,ci}{ti} = f_trough_temp;%first col stores indices, second col stores potentials

                    %average fast trough (all APs)
                    f_trough_ave{1,ci}(ti,1) = nanmean(f_trough{1,ci}{1,ti}(:,2));

                    %dV/dt for each trace
                    dV_sec{1,ci}{ti} = dV_sec_temp;

                    %MFR (mean firing rate)
                    MFR{1,ci}(ti,1) = sp_num_temp; %mean FR in Hz (duration 1s)

                    %ISI (interspike interval, in ms)
                    %defined as the difference between the fast trough of the
                    %initial AP and the threshold of the next AP
                    if sp_num_temp == 0 || sp_num_temp == 1
                        ISI{1,ci}(ti,:) = NaN;
                    else
                        for ii = 1:(sp_num_temp-1)
                            ISI{1,ci}(ti,ii) = (V_th_temp(ii+1,1) - f_trough_temp(ii,1))*0.1;
                        end
                    end

                    %IFR (instantaneous firing rate, in Hz))
                    %calculated as the reciprocal of ISI
                    if sp_num_temp == 0 || sp_num_temp == 1
                        IFR{1,ci}(ti,:) = 0;
                    else
                        for fi = 1:(sp_num_temp-1)
                            IFR{1,ci}(ti,fi) = 1000/ISI{1,ci}(ti,fi);
                        end

                    end

                    %AP peak amplitude
                    %first col stores amplitudes, second col stores peak indice
                    if sp_num_temp == 0
                        AP_peak{1,ci}{ti} = NaN;
                    else
                        for pi = 1:sp_num_temp
                            [peak_amp,peak_ind] = ...
                                max(data(int32(V_th_temp(pi,1)):int32(f_trough_temp(pi,1))));
                            AP_peak{1,ci}{1,ti}(pi,1) = peak_ind +V_th_temp(pi,1); %peak index
                            AP_peak{1,ci}{1,ti}(pi,2) = peak_amp; %peak amplitude
                            %[AP_peak{1,ci}{1,ti}(pi,2),AP_peak{1,ci}{1,ti}(pi,1)]...
                            %    = max(data(int32(V_th_temp(pi,1)):int32(f_trough_temp(pi,1))));
                        end
                    end

                    %calulate average/1st AP peak amplitude for each trace
                    if isnan(AP_peak{1,ci}{1,ti}) == 0
                        AP_peak_ave{1,ci}(ti,1) = nanmean(AP_peak{1,ci}{1,ti}(:,2));
                        AP_peak_1st{1,ci}(ti,1) = AP_peak{1,ci}{1,ti}(1,2);
                    else
                        AP_peak_ave{1,ci}(ti,1) = NaN;
                        AP_peak_1st{1,ci}(ti,1) = NaN;
                    end

                    %latency (in ms)
                    if sp_num_temp == 0
                       lat{1,ci}(ti,1) = NaN;
                    else
                       lat{1,ci}(ti,1) = (V_th_temp(1,1)-10000)*0.1;
                    end

                    %spike adaptation index
                    %the rate of ISI change
                    if sp_num_temp == 0 || sp_num_temp == 1 || sp_num_temp == 2
                        adp_index{1,ci}(ti,1) = NaN;
                    else
                        adp_index{1,ci}(ti,1) = nanmean(diff(ISI{1,ci}(ti,:)))/nansum(ISI{1,ci}(ti,:));
                    end

                    %upstroke/downstroke ratio
                    %defined as the ratio of the maxmium to the minimum dV/dt
                    %for each AP
                    if sp_num_temp == 0
                        udratio{1,ci}(ti,1:sp_num_temp) = NaN;
                        udratio_median{1,ci}(ti,1) = NaN;
                    else
                        udratio{1,ci}(ti,1:sp_num_temp) = udratio_temp(1:sp_num_temp);

                        %calculate median
                        udratio_median{1,ci}(ti,1) = nanmedian(udratio{1,ci}(ti,1:sp_num_temp));

                        %calculate the slope of change in U/D ratio for each trace
                        if sp_num_temp < 2
                            udratio_slope{1,ci}(ti,1) = NaN;
                        else
                            [fm1,gof] = fit((1:sp_num_temp)',(udratio{1,ci}(ti,1:sp_num_temp))','poly1');

    %                         subplot(4,5,ti-trace_start+1);
    %                         plot((1:sp_num_temp)',(udratio{1,ci}(ti,1:sp_num_temp))','o')
    %                         hold on
    %                         plot(fm1)

                            udratio_slope{1,ci}(ti,1) = fm1.p1;
                        end
                    end


                    %AP width at half height (in ms)
                    if sp_num_temp == 0
                        width{1,ci}(ti,:) = NaN;
                        width_1st{1,ci}(ti,1) = NaN;
                        width_median{1,ci}(ti,1) = NaN;
                    else
                        for wi = 1:sp_num_temp
                            %downward halfheight
                            halfheight_d = 0.5*(AP_peak{1,ci}{1,ti}(wi,2)-f_trough{1,ci}{1,ti}(wi,2))+...
                                f_trough{1,ci}{1,ti}(wi,2);
                            halfheight_d_ind = find(data(int32(AP_peak{1,ci}{1,ti}(wi,1)):int32(f_trough{1,ci}{1,ti}(wi,1)))<=...
                                halfheight_d,1,'first')+AP_peak{1,ci}{1,ti}(wi,1);

                            %extrapolate to the upward side
                            data_extrp = data(int32(V_th{1,ci}{1,ti}(wi,1)):int32(AP_peak{1,ci}{1,ti}(wi,1)));
                            [data_extrp, index] = unique(data_extrp);
                            %data_extrp = data_extrp';
                            ind_extrp = V_th{1,ci}{1,ti}(wi,1):AP_peak{1,ci}{1,ti}(wi,1);
                            ind_extrp = ind_extrp(index)';


                            halfheight_u_ind = interp1(data_extrp, ind_extrp, halfheight_d);
    %                         %upward halfheight
    %                         halfheight_u = 0.5*AP_peak{1,ci}(ti,wi)-V_th{1,ci}{1,ti}(wi,1);
    %                         half_height_u_ind = find(data(V_th{1,ci}{1,ti}(wi,2):AP_peak{1,ci}{1,ti}(wi,2))>=...
    %                             AP_peak{1,ci}{1,ti}(wi,1)-halfheight_u,1,'first');

                            %calculate width
                            width{1,ci}(ti,wi) = (halfheight_d_ind - halfheight_u_ind)*0.1;

    %                         %troubleshooting width
    %                         figure;
    %                         hold on;
    %                         plot(data,'k')
    %                         plot(halfheight_d_ind,halfheight_d,'o','markersize',10,'markerfacecolor','r')
    %                         plot(halfheight_u_ind,halfheight_d,'o','markersize',10,'markerfacecolor','b')


                        end

                        %width of the first AP in each current step
                        width_1st{1,ci}(ti,1) = width{1,ci}(ti,1);
                        
                        %calculate median of width
                        width_median{1,ci}(ti,1) = nanmedian(width{1,ci}(ti,:));

                        %calculate width slope   
                        if sp_num_temp < 2
                            width_slope{1,ci}(ti,1) = NaN;
                        else
                            [fm2,gof] = fit((1:sp_num_temp)',(width{1,ci}(ti,1:sp_num_temp))','poly1');

    %                         subplot(4,5,ti-trace_start+1);
    %                         plot((1:sp_num_temp)',(width{1,ci}(ti,1:sp_num_temp))','o')
    %                         hold on
    %                         plot(fm2)

                            width_slope{1,ci}(ti,1) = fm2.p1;
                        end



                    end
                end
            end %per trace           
        end 

        %calculate rheobase for each cell
        for tti = 1:(trace_start+size(cell_id{1,ci},1)-1)
            if spike_count{1,ci}(tti,1) ~= 0
                rheobase(ci,1) = curr_inj{1,ci}(tti,1);
                rheobase_ind(ci,1) = tti-trace_start+1; %btw 1 and 20
                break
            end
        end

        %mean IFR (mean of the first two IFRs)
        for tiii = trace_start:ti
            if IFR{1,ci}(tiii,1) == 0
                IFR_ave{1,ci}(tiii,1) = 0;
            elseif IFR{1,ci}(tiii,1) == 0
                IFR_ave{1,ci}(tiii,1) = IFR{1,ci}(tiii,1);
            else
                IFR_ave{1,ci}(tiii,1) = nanmean(IFR{1,ci}(tiii,1:1));
            end
        end

        %mean and std of Vrm during current injections for each cell
        Vrm_stat(ci,1) = nanmean(Vrm{1,ci}(trace_start:ti,1));
        Vrm_stat(ci,2) = std(Vrm{1,ci}(trace_start:ti,1));
        Vrm_stat(ci,3) = abs(Vrm_stat(ci,1)+70)/70*100;
        %whether Vrm is acceptable
        if Vrm_stat(ci,3)>5 %error less than +-5%
            Vrm_stat(ci,4) = 0;
        else
            Vrm_stat(ci,4) = 1;
        end


        %waitforbuttonpress


    end %per cell

end %per experiment

    FI_DATA.filename = filename(1:cellID);
    FI_DATA.cell_id = cell_id(1:cellID);
    FI_DATA.aDAT = aDAT(1:cellID); 
    FI_DATA.Ra = all_Ra(1:cellID);
    FI_DATA.Rin = all_Rin(1:cellID); 
    FI_DATA.Cp = all_Cp(1:cellID);
    FI_DATA.Vr = all_Vr(1:cellID); 
    FI_DATA.MFR = MFR(1:cellID); %mean firing rate
    FI_DATA.ISI = ISI(1:cellID); %interspike interval (for each two APs)
    FI_DATA.IFR = IFR(1:cellID); %instantenous firing rate (for each ISI)
    FI_DATA.IFR_ave = IFR_ave(1:cellID); %average of the first two IFR
    FI_DATA.V_th = V_th(1:cellID); %AP threshold (for each AP)
    FI_DATA.V_th_ave = V_th_ave(1:cellID); %average AP threshold
    FI_DATA.V_th_1st = V_th_1st(1:cellID); %threshold of the first AP in each trace
    FI_DATA.V_th_4th = V_th_4th(1:cellID); %threshold of the fourth AP in each trace
    FI_DATA.f_trough = f_trough(1:cellID); %AP fast trough (for each AP)
    FI_DATA.f_trough_ave = f_trough_ave(1:cellID); %average fast trough
    FI_DATA.dV_sec = dV_sec(1:cellID); %dV/dt for each trace
    FI_DATA.AP_peak = AP_peak(1:cellID); %peak amplitde of each AP
    FI_DATA.AP_peak_ave = AP_peak_ave(1:cellID); %average AP peak amp
    FI_DATA.AP_peak_1st = AP_peak_1st(1:cellID); %first AP amp
    FI_DATA.Vrm = Vrm(1:cellID); %resting membrane potential
    FI_DATA.Vrm_stat = Vrm_stat(1:cellID,:); %stats of Vrm during current injection
    FI_DATA.spike_count = spike_count(1:cellID); %number of spikes evoked
    FI_DATA.curr_inj = curr_inj(1:cellID); %current injected
    FI_DATA.rheobase = rheobase(1:cellID); % each row stores the rheobase of a cell
    FI_DATA.rheobase_ind = rheobase_ind(1:cellID); %eahc row stores the index of rheobase (which stimulation)
    FI_DATA.lat = lat(1:cellID); %latency between the start of the current step and the first AP
    FI_DATA.adp_index = adp_index(1:cellID); %spike adaptation index
    FI_DATA.udratio = udratio(1:cellID); %upstroke/downstroke ratio
    FI_DATA.udratio_slope = udratio_slope(1:cellID); %slope of the u/d ratio change across APs
    FI_DATA.udratio_median = udratio_median(1:cellID); %median of u/d ratio for each current step
    FI_DATA.width = width(1:cellID); %AP width at half height
    FI_DATA.width_1st = width_1st(1:cellID); %width of the first AP for each current step
    FI_DATA.width_slope = width_slope(1:cellID); %slope of the width change across APs
    FI_DATA.width_median = width_median(1:cellID); %median of width for each current step

%% save to file
if save_results == 1
    cd(fp_data)
%     current_dir_name = cd;
%     total_dir_name_number = numel(current_dir_name);
%     display_name = strcat(current_dir_name(total_dir_name_number-5: end),'.mat');
    save(save_file_name,'FI_DATA')
end