
% script to plot the dv/dt v V circle plots for action potential current
% clamp recordings

%%% TODO:


% replace these strings with paths to the .mat FI data scripts for your
% conditions
cond1_path = "C:\Users\norelis\Documents\SFARISCN2AmergeEVN_EVP_062222to021223.mat";
cond2_path = "C:\Users\bcary\Downloads\SCN2AmergeEVN_EVP_062222to021223.mat";

cond1_load = load(cond1_path);
cond2_load = load(cond2_path);

% input_data_cond1 = FI_DATA;
input_data_cond1 = cond1_load.FI_DATA;
input_data_cond2 = cond2_load.FI_DATA;

%%% Plotting Params %%%
% change these color vectors if you want different colors for plotting
colors = [0.3 0.3 0.3;...
          0.2 0.9 0.4];

% replace these strings with desired condition names for plotting
cond_names = {'Cont','SCN2A+/-'};

% save_dir = '';

%%%%%%%%%%%%%%%%%%%%%%%

%%% THRESHOLD PARAMS %%%
Ra_thresh = 20e6;
Rin_low_thresh = 90e6;
Rin_high_thresh = 350e6;

Vrm_thresh = -50;

Vth_change_thresh = 0.6;
%%%%%%%%%%%%%%%%%%%%%%%%

%%% action potential analysis params %%%
first_ap_only = 1;

samp_rate = 10e3;
interp_factor = 3;
pre_vth_ms = 0.5;
post_vth_ms = 4;
ap_num_inds = ((post_vth_ms - pre_vth_ms + 1)/1000)*samp_rate*interp_factor + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main loop
input_conds = {input_data_cond1, input_data_cond2};
rawV_ap_cell_avgs_full = {};
dvdt_ap_cell_avgs_full = {};
vth_cell_avgs_full = {};
peak_dvdt_cell_avgs_full = {};
for cond_i = 1:length(input_conds)

    disp(['Working on condition: ',num2str(cond_i)])
    % cell array of cell data. each entry has one cell's data. within that is
    % cell array each entry is a trace
    dvdt_strct = input_conds{cond_i}.dV_sec;
    
    % cell of cell data. each entry one cell data. each entry is a matrix where
    % col's are traces
    rawV_strct = input_conds{cond_i}.aDAT;
    
    Vth_strct = input_conds{cond_i}.V_th;
    Ra_strct = input_conds{cond_i}.Ra;
    Rin_strct = input_conds{cond_i}.Rin;
    Vrm_strct = input_conds{cond_i}.Vrm;

    rawV_ap_cell_avgs = nan([ap_num_inds, length(dvdt_strct{1})]);
    dvdt_ap_cell_avgs = nan([ap_num_inds, length(dvdt_strct{1})]);
    vth_cell_avgs = nan([1,length(dvdt_strct{1})]);
    peak_dvdt_cell_avgs = nan([1,length(dvdt_strct{1})]);
    for cell_i = 1:length(dvdt_strct)
    %     cell_i = 1;

        disp(['Cell number: ',num2str(cell_i)])
        
        dvdt_cell = dvdt_strct{cell_i};
        rawV_cell = rawV_strct{cell_i};
        Vth_cell = Vth_strct{cell_i};
        Ra_cell = Ra_strct{cell_i};
        Rin_cell = Rin_strct{cell_i};
        Vrm_cell = Vrm_strct{cell_i};
        
        all_vths = [];
        first_vths = [];
        for i = 1:length(Vth_cell)
            first_vths(i) = Vth_cell{i}(1,2);
            all_vths = [all_vths; Vth_cell{i}(:,2)];
        end
        first_vths(isnan(first_vths)) = [];
        avg_vth_change = abs(sum(diff(first_vths)))/length(first_vths);

        if mean(Ra_cell,'omitnan') > Ra_thresh
            disp('Ra too high')
            continue
        elseif mean(Rin_cell,'omitnan') > Rin_high_thresh
            disp('Rin too high')
            continue
        elseif mean(Rin_cell,'omitnan') < Rin_low_thresh
            disp('Rin too low')
            continue
        elseif mean(Vrm_cell,'omitnan') > Vrm_thresh
            disp('Vrm too high')
            continue
        elseif avg_vth_change > Vth_change_thresh
            disp('Too much vth change')
            continue
        end

        if first_ap_only
            vth_cell_avgs(cell_i) = mean(first_vths,'omitnan');
        else
            vth_cell_avgs(cell_i) = mean(all_vths,'omitnan');
        end

        rawV_ap_trace_avgs = nan([ap_num_inds, length(dvdt_cell)]);
        dvdt_ap_trace_avgs = nan([ap_num_inds, length(dvdt_cell)]);
        peak_dvdt_trace = [];
        for trace_i = 1:length(dvdt_cell)
        %     trace_i = 7;
            
            dvdt_trace = dvdt_cell{trace_i};
            dvdt_trace = [0; dvdt_trace];
            rawV_trace = rawV_cell(:,trace_i);
            
            vth_inds_trace = Vth_cell{trace_i}(:,1);
            vth_trace = Vth_cell{trace_i}(:,2);
            
            x = 1:length(dvdt_trace);
            x_int = 1/interp_factor:1/interp_factor:length(dvdt_trace);
            
            dvdt_trace_int = interp1(x,dvdt_trace,x_int,'pchip');
            rawV_trace_int = interp1(x,rawV_trace,x_int,'pchip');
            
    %         figure;
    %         plot(x,rawV_trace)

    %         figure; hold on;
    %         % plot(rawV_trace,dvdt_trace)
    %         plot(rawV_trace_int,dvdt_trace_int)
    %         plot(mean(rawV_trace_int,'omitnan'),mean(dvdt_trace_int,'omitnan'))
    %         
            num_aps = sum(~isnan(vth_inds_trace));
            
            if first_ap_only && num_aps > 0
                num_aps = 1;
            end
    
            rawV_aps = [];
            dvdt_aps = [];
            if num_aps > 0
                for ap_i = 1:num_aps
                    vth_ind_ap = vth_inds_trace(ap_i);
                    vth_ap = vth_trace(ap_i);
            
                    ind1 = vth_ind_ap*interp_factor - (pre_vth_ms/1000)*samp_rate*interp_factor;
                    ind2 = vth_ind_ap*interp_factor + (post_vth_ms/1000)*samp_rate*interp_factor;
                    
                    rawV_aps(:,ap_i) = rawV_trace_int(ind1:ind2);
                    dvdt_aps(:,ap_i) = dvdt_trace_int(ind1:ind2);

                    peak_dvdt_trace = [peak_dvdt_trace max(dvdt_trace_int(ind1:ind2))];
                end
                
                rawV_ap_trace_avgs(:,trace_i) = mean(rawV_aps,2,'omitnan');
                dvdt_ap_trace_avgs(:,trace_i) = mean(dvdt_aps,2,'omitnan');        
            end
            
    %         figure;hold on
    %         plot(rawV_aps)
    %             
    %         figure;hold on
    %         plot(rawV_ap_trace_avgs(:,trace_i) ,dvdt_ap_trace_avgs(:,trace_i))
        
        end
        
%         figure;
%         plot(rawV_ap_trace_avgs,dvdt_ap_trace_avgs,'color',[0 0 0 0.45])
%         hold on
%         plot(mean(rawV_ap_trace_avgs,2,'omitnan'),mean(dvdt_ap_trace_avgs,2,'omitnan'),'k-','linewidth',3)
%         title(['Cell: ',num2str(cell_i)])
    % %     
        rawV_ap_cell_avgs(:,cell_i) = mean(rawV_ap_trace_avgs,2,'omitnan');
        dvdt_ap_cell_avgs(:,cell_i) = mean(dvdt_ap_trace_avgs,2,'omitnan');
        peak_dvdt_cell_avgs(cell_i) = mean(peak_dvdt_trace,'omitnan');
    
    end
    
    rawV_ap_cell_avgs_full{cond_i} = rawV_ap_cell_avgs;
    dvdt_ap_cell_avgs_full{cond_i} = dvdt_ap_cell_avgs;
    vth_cell_avgs_full{cond_i} = vth_cell_avgs;
    peak_dvdt_cell_avgs_full{cond_i} = peak_dvdt_cell_avgs;
end

spike_time_x = (1:ap_num_inds)/(samp_rate*interp_factor);
spike_time_x = spike_time_x * 1000;

ap_overlay_f = figure;
ap_ax = gca;
hold(ap_ax,'on')
box off
ylabel('Voltage (mV)')
xlabel('time (ms)')
title(['AP waveform - cell avgs'])
for cond_i = 1:length(input_conds)
    cond_data = rawV_ap_cell_avgs_full{cond_i};
    cond_color = colors(cond_i,:);

    figure;
    plot(spike_time_x,cond_data,'color',[cond_color 0.18])
    hold on
    plot(spike_time_x,mean(cond_data,2,'omitnan'),'color',cond_color,'linewidth',3)
    box off
    ylabel('Voltage (mV)')
    xlabel('time (ms)')
    title(['AP waveform - cell avgs - Condition: ',cond_names{cond_i}])

    plot(ap_ax,spike_time_x,mean(cond_data,2,'omitnan'),'color',cond_color,'linewidth',3)

end
legend(ap_ax,cond_names)


dvdt_overlay_f = figure;
dvdt_ax = gca;
hold(dvdt_ax,'on')
box off
ylabel('Voltage (mV)')
xlabel('time (ms)')
title(['AP waveform - cell avgs'])
for cond_i = 1:length(input_conds)
    cond_rawV = rawV_ap_cell_avgs_full{cond_i};
    cond_dvdt = dvdt_ap_cell_avgs_full{cond_i};
    cond_color = colors(cond_i,:);

    figure;
    plot(cond_rawV,cond_dvdt,'color',[cond_color 0.45])
    hold on
    plot(mean(cond_rawV,2,'omitnan'),mean(cond_dvdt,2,'omitnan'),...
        'color',cond_color,'linewidth',3)
    box off 
    xlabel('Voltage (mV)')
    ylabel('dV/dt (mV/ms)')
    title(['dVdt v V during AP - avgs per cell - Condition: ',cond_names{cond_i}])

    plot(dvdt_ax,mean(cond_rawV,2,'omitnan'),mean(cond_dvdt,2,'omitnan'),...
        'color',cond_color,'linewidth',3)
end
legend(dvdt_ax,cond_names)


% bee swarm plot: AP threshold, peak dvdt
data_toplot = [];
for i = 1:length(input_conds)
    data_toplot = padmat(data_toplot, vth_cell_avgs_full{i}', 2);
end

figure_coords = [301 550 (271+(100*(length(input_conds)-2))) 390];
f1 = figure('Position',figure_coords);

p1 = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',colors,...
    'PointSize',35,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',2,'MarkerEdgeColor',colors);

box off
ylabel('AP Thresh (mv)','FontSize',14)
set(gca,'XTickLabel',cond_names,'XTickLabelRotation',45,'FontSize',14);
xlim([0.5 length(input_conds)+0.5])
ylim([-45 -25])
% legend(cond_names)
title('AP Thresh - cell avgs','FontSize',12)


% peak dvdt
data_toplot = [];
for i = 1:length(input_conds)
    data_toplot = padmat(data_toplot, peak_dvdt_cell_avgs_full{i}', 2);
end

figure_coords = [301 550 (271+(100*(length(input_conds)-2))) 390];
f1 = figure('Position',figure_coords);

p1 = UnivarScatter_ATP(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',colors,...
    'PointSize',35,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',2,'MarkerEdgeColor',colors);

box off
ylabel('Peak dvdt (V/s ??)','FontSize',14)
set(gca,'XTickLabel',cond_names,'XTickLabelRotation',45,'FontSize',14);
xlim([0.5 length(input_conds)+0.5])
ylim([0 inf])
% legend(cond_names)
title('Peak dvdt - cell avgs','FontSize',12)
























% % experimental - probably need to process data differently to create sem
% % error around curve
% figure;
% edges_toplot = mean(rawV_ap_cell_avgs,2,'omitnan');
% alpha = 0.3;
% sh_color = 'k';
% plot_strcts = stdshade(dvdt_ap_cell_avgs',alpha,sh_color,edges_toplot,1,'sem');
% box off
% xlabel('Voltage (mV)')
% ylabel('dV/dt (mV/ms)')
% title('dVdt v V during AP - avgs per cell')






