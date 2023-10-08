%%%
%Takes Action potential kinectis analysis to asses files and plot the
%different action potential properties. 


%SCN2A vs WT
%Scn2a
FI_DRP_NaI = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/SCN2A_DRP_Aug10_22toMay22_23_NaI.mat');
FI_EVP_NaI = load('/Users/norelis/Documents/MATLAB/FI_DATA/SCN2A_FI_2Virus/SCN2AmergeEVN_EVP_June22_22toMay22_23_NaI.mat');

%Wildtype
FI_DRP_WT_NaI =load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/WT_DRP_June22_22toMay4_23_NaI.mat');
FI_EVP_WT_NaI=load('/Users/norelis/Documents/MATLAB/FI_DATA/WT_FI_2Virus/WT_mergeEVN_EVP_Jun22_22toMay4_23_NaI.mat');



%%
% Acess files for EVN & EVP merge for SCN2A animals

SCN2A_EVP_all_cell_ap_falls = [];
SCN2A_EVP_all_cell_ap_rise= [];
SCN2A_EVP_all_cell_ap_max_rise_rate=[];
SCN2A_EVP_all_cell_ap_max_fall_rate=[];

%AP_falls
for h = 1:size(FI_EVP_NaI.ap_fall,2)

    if mean(FI_EVP_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_NaI.Rin{1,h}) < 90000000 || mean(FI_EVP_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
    end
    
    SCN2A_EVP_ap_fall_cell = FI_EVP_NaI.ap_fall{1,h};
    
    SCN2A_EVP_all_trace_ap_falls = [];
    for trace_i = 1:length(SCN2A_EVP_ap_fall_cell)
        SCN2A_EVP_all_trace_ap_falls = [SCN2A_EVP_all_trace_ap_falls, SCN2A_EVP_ap_fall_cell{trace_i}'];
    end
    SCN2A_EVP_all_cell_ap_falls(end+1) = mean(SCN2A_EVP_all_trace_ap_falls, 'omitnan');
end 
    

%AP rise     
for h = 1:size(FI_EVP_NaI.ap_rise,2) 
   
     
       if mean(FI_EVP_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_NaI.Rin{1,h}) < 90000000 || mean(FI_EVP_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    SCN2A_EVP_ap_rise_cell=FI_EVP_NaI.ap_rise{1,h};
     
    SCN2A_EVP_all_trace_ap_rise=[];
    for trace_i = 1:length(SCN2A_EVP_ap_rise_cell)
        SCN2A_EVP_all_trace_ap_rise= [SCN2A_EVP_all_trace_ap_rise,SCN2A_EVP_ap_rise_cell{trace_i}'];
    end 
    SCN2A_EVP_all_cell_ap_rise(end+1)= mean(SCN2A_EVP_all_trace_ap_rise, 'omitnan');  
end

for h = 1:size(FI_EVP_NaI.ap_max_rise_rate,2) 
   
     
       if mean(FI_EVP_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_NaI.Rin{1,h}) < 90000000 || mean(FI_EVP_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    SCN2A_EVP_ap_max_rise_rate_cell=FI_EVP_NaI.ap_max_rise_rate{1,h};
     
    SCN2A_EVP_all_trace_ap_max_rise_rate=[];
    for trace_k = 1:length(SCN2A_EVP_ap_max_rise_rate_cell)
        SCN2A_EVP_all_trace_ap_max_rise_rate= [SCN2A_EVP_all_trace_ap_max_rise_rate,SCN2A_EVP_ap_max_rise_rate_cell{trace_k}'];
    end 
    SCN2A_EVP_all_cell_ap_max_rise_rate(end+1)= mean(SCN2A_EVP_all_trace_ap_max_rise_rate, 'omitnan');  
end

%AP max rise rate   
for h = 1:size(FI_EVP_NaI.ap_max_fall_rate,2) 
   
     
       if mean(FI_EVP_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_NaI.Rin{1,h}) < 90000000 || mean(FI_EVP_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    SCN2A_EVP_ap_max_fall_rate_cell=FI_EVP_NaI.ap_max_fall_rate{1,h};
     
    SCN2A_EVP_all_trace_ap_max_fall_rate=[];
    for trace_k = 1:length(SCN2A_EVP_ap_max_fall_rate_cell)
        SCN2A_EVP_all_trace_ap_max_fall_rate= [SCN2A_EVP_all_trace_ap_max_fall_rate,SCN2A_EVP_ap_max_fall_rate_cell{trace_k}'];
    end 
    SCN2A_EVP_all_cell_ap_max_fall_rate(end+1)= mean(SCN2A_EVP_all_trace_ap_max_fall_rate, 'omitnan');  
end



%%
%Asses files for DRP of SCN2A animals


SCN2A_DRP_all_cell_ap_falls = [];
SCN2A_DRP_all_cell_ap_rise= [];
SCN2A_DRP_all_cell_ap_max_rise_rate=[];
SCN2A_DRP_all_cell_ap_max_fall_rate=[];

%AP_falls
for k = 1:size(FI_DRP_NaI.ap_fall,2)

    if mean(FI_DRP_NaI.Ra{1,k}) > 20000000 % exclude cells witk an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_NaI.Rin{1,k}) < 90000000 || mean(FI_DRP_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,k}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
    end
    
    SCN2A_DRP_ap_fall_cell = FI_DRP_NaI.ap_fall{1,k};
    
    SCN2A_DRP_all_trace_ap_falls = [];
    for trace_i = 1:length(SCN2A_DRP_ap_fall_cell)
        SCN2A_DRP_all_trace_ap_falls = [SCN2A_DRP_all_trace_ap_falls, SCN2A_DRP_ap_fall_cell{trace_i}'];
    end
    SCN2A_DRP_all_cell_ap_falls(end+1) = mean(SCN2A_DRP_all_trace_ap_falls, 'omitnan');
end 
    

%AP rise     
for k = 1:size(FI_DRP_NaI.ap_rise,2) 
   
     
       if mean(FI_DRP_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_NaI.Rin{1,k}) < 90000000 || mean(FI_DRP_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    SCN2A_DRP_ap_rise_cell=FI_DRP_NaI.ap_rise{1,k};
     
    SCN2A_DRP_all_trace_ap_rise=[];
    for trace_i = 1:length(SCN2A_DRP_ap_rise_cell)
        SCN2A_DRP_all_trace_ap_rise= [SCN2A_DRP_all_trace_ap_rise,SCN2A_DRP_ap_rise_cell{trace_i}'];
    end 
    SCN2A_DRP_all_cell_ap_rise(end+1)= mean(SCN2A_DRP_all_trace_ap_rise, 'omitnan');  
end


for k = 1:size(FI_DRP_NaI.ap_max_rise_rate,2) 
   
     
       if mean(FI_DRP_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_NaI.Rin{1,k}) < 90000000 || mean(FI_DRP_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    SCN2A_DRP_ap_max_rise_rate_cell=FI_DRP_NaI.ap_max_rise_rate{1,k};
     
    SCN2A_DRP_all_trace_ap_max_rise_rate=[];
    for trace_k = 1:length(SCN2A_DRP_ap_max_rise_rate_cell)
        SCN2A_DRP_all_trace_ap_max_rise_rate= [SCN2A_DRP_all_trace_ap_max_rise_rate,SCN2A_DRP_ap_max_rise_rate_cell{trace_k}'];
    end 
    SCN2A_DRP_all_cell_ap_max_rise_rate(end+1)= mean(SCN2A_DRP_all_trace_ap_max_rise_rate, 'omitnan');  
end

%AP max rise rate   
for k = 1:size(FI_DRP_NaI.ap_max_fall_rate,2) 
   
     
       if mean(FI_DRP_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_NaI.Rin{1,k}) < 90000000 || mean(FI_DRP_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    SCN2A_DRP_ap_max_fall_rate_cell=FI_DRP_NaI.ap_max_fall_rate{1,k};
     
    SCN2A_DRP_all_trace_ap_max_fall_rate=[];
    for trace_k = 1:length(SCN2A_DRP_ap_max_fall_rate_cell)
        SCN2A_DRP_all_trace_ap_max_fall_rate= [SCN2A_DRP_all_trace_ap_max_fall_rate,SCN2A_DRP_ap_max_fall_rate_cell{trace_k}'];
    end 
    SCN2A_DRP_all_cell_ap_max_fall_rate(end+1)= mean(SCN2A_DRP_all_trace_ap_max_fall_rate, 'omitnan');  
end




%%
%Asses DRP files for WT animals 
WT_DRP_all_cell_ap_falls = [];
WT_DRP_all_cell_ap_rise= [];
WT_DRP_all_cell_ap_max_rise_rate= []; 
WT_DRP_all_cell_ap_max_fall_rate= []; 

%AP_falls
for h = 1:size(FI_DRP_WT_NaI.ap_fall,2)

    if mean(FI_DRP_WT_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_WT_NaI.Rin{1,h}) < 90000000 || mean(FI_DRP_WT_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
    end
    
    WT_DRP_ap_fall_cell = FI_DRP_WT_NaI.ap_fall{1,h};
    
    WT_DRP_all_trace_ap_falls = [];
    for trace_i = 1:length(WT_DRP_ap_fall_cell)
        WT_DRP_all_trace_ap_falls = [WT_DRP_all_trace_ap_falls, WT_DRP_ap_fall_cell{trace_i}'];
    end
    WT_DRP_all_cell_ap_falls(end+1) = mean(WT_DRP_all_trace_ap_falls,'omitnan');
end 
    

%AP rise     
for h = 1:size(FI_DRP_WT_NaI.ap_rise,2) 
   
     
       if mean(FI_DRP_WT_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_WT_NaI.Rin{1,h}) < 90000000 || mean(FI_DRP_WT_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    WT_DRP_ap_rise_cell=FI_DRP_WT_NaI.ap_rise{1,h};
     
    WT_DRP_all_trace_ap_rise=[];
    for trace_i = 1:length(WT_DRP_ap_rise_cell)
        WT_DRP_all_trace_ap_rise= [WT_DRP_all_trace_ap_rise,WT_DRP_ap_rise_cell{trace_i}'];
    end 
    WT_DRP_all_cell_ap_rise(end+1)= mean(WT_DRP_all_trace_ap_rise, 'omitnan');  
end


for h = 1:size(FI_DRP_WT_NaI.ap_max_rise_rate,2) 
   
     
       if mean(FI_DRP_WT_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_WT_NaI.Rin{1,h}) < 90000000 || mean(FI_DRP_WT_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    WT_DRP_ap_max_rise_rate_cell=FI_DRP_WT_NaI.ap_max_rise_rate{1,h};
     
    WT_DRP_all_trace_ap_max_rise_rate=[];
    for trace_k = 1:length(WT_DRP_ap_max_rise_rate_cell)
        WT_DRP_all_trace_ap_max_rise_rate= [WT_DRP_all_trace_ap_max_rise_rate,WT_DRP_ap_max_rise_rate_cell{trace_k}'];
    end 
    WT_DRP_all_cell_ap_max_rise_rate(end+1)= mean(WT_DRP_all_trace_ap_max_rise_rate, 'omitnan');  
end

%AP max rise rate   
for h = 1:size(FI_DRP_WT_NaI.ap_max_fall_rate,2) 
   
     
       if mean(FI_DRP_WT_NaI.Ra{1,h}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_DRP_WT_NaI.Rin{1,h}) < 90000000 || mean(FI_DRP_WT_NaI.Rin{1,h}) > 350000000 
         continue
% %    elseif mean(FI_DRP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    WT_DRP_ap_max_fall_rate_cell=FI_DRP_WT_NaI.ap_max_fall_rate{1,h};
     
    WT_DRP_all_trace_ap_max_fall_rate=[];
    for trace_k = 1:length(WT_DRP_ap_max_fall_rate_cell)
        WT_DRP_all_trace_ap_max_fall_rate= [WT_DRP_all_trace_ap_max_fall_rate,WT_DRP_ap_max_fall_rate_cell{trace_k}'];
    end 
    WT_DRP_all_cell_ap_max_fall_rate(end+1)= mean(WT_DRP_all_trace_ap_max_fall_rate, 'omitnan');  
end



%%
%Asses EVP files for WT animals
WT_EVP_all_cell_ap_falls = [];
WT_EVP_all_cell_ap_rise= [];
WT_EVP_all_cell_ap_max_rise_rate= []; 
WT_EVP_all_cell_ap_max_fall_rate= []; 

%AP_falls
for k = 1:size(FI_EVP_WT_NaI.ap_fall,2)

    if mean(FI_EVP_WT_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_WT_NaI.Rin{1,k}) < 90000000 || mean(FI_EVP_WT_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_EVP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
    end
    
    WT_EVP_ap_fall_cell = FI_EVP_WT_NaI.ap_fall{1,k};
    
    WT_EVP_all_trace_ap_falls = [];
    for trace_k = 1:length(WT_EVP_ap_fall_cell)
        WT_EVP_all_trace_ap_falls = [WT_EVP_all_trace_ap_falls, WT_EVP_ap_fall_cell{trace_k}'];
    end
    WT_EVP_all_cell_ap_falls(end+1) = mean(WT_EVP_all_trace_ap_falls, 'omitnan');
end 
    

%AP rise     
for k = 1:size(FI_EVP_WT_NaI.ap_rise,2) 
   
     
       if mean(FI_EVP_WT_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_WT_NaI.Rin{1,k}) < 90000000 || mean(FI_EVP_WT_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_EVP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    WT_EVP_ap_rise_cell=FI_EVP_WT_NaI.ap_rise{1,k};
     
    WT_EVP_all_trace_ap_rise=[];
    for trace_k = 1:length(WT_EVP_ap_rise_cell)
        WT_EVP_all_trace_ap_rise= [WT_EVP_all_trace_ap_rise,WT_EVP_ap_rise_cell{trace_k}'];
    end 
    WT_EVP_all_cell_ap_rise(end+1)= mean(WT_EVP_all_trace_ap_rise, 'omitnan');  
end

%AP max rise rate   
for k = 1:size(FI_EVP_WT_NaI.ap_max_rise_rate,2) 
   
     
       if mean(FI_EVP_WT_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_WT_NaI.Rin{1,k}) < 90000000 || mean(FI_EVP_WT_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_EVP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    WT_EVP_ap_max_rise_rate_cell=FI_EVP_WT_NaI.ap_max_rise_rate{1,k};
     
    WT_EVP_all_trace_ap_max_rise_rate=[];
    for trace_k = 1:length(WT_EVP_ap_max_rise_rate_cell)
        WT_EVP_all_trace_ap_max_rise_rate= [WT_EVP_all_trace_ap_max_rise_rate,WT_EVP_ap_max_rise_rate_cell{trace_k}'];
    end 
    WT_EVP_all_cell_ap_max_rise_rate(end+1)= mean(WT_EVP_all_trace_ap_max_rise_rate, 'omitnan');  
end

%AP max rise rate   
for k = 1:size(FI_EVP_WT_NaI.ap_max_fall_rate,2) 
   
     
       if mean(FI_EVP_WT_NaI.Ra{1,k}) > 20000000 % exclude cells with an acess resistance higher than 20Mohms 
        continue
%         %exclude cells with an input resistance lower 90Mohms and higher
%         %than 350Mohms
    elseif mean(FI_EVP_WT_NaI.Rin{1,k}) < 90000000 || mean(FI_EVP_WT_NaI.Rin{1,k}) > 350000000 
         continue
% %    elseif mean(FI_EVP.Cp{1,h}) < 0.0000000001 %capacitance is lower for
%     %younger cells 
%     %potentially memebrane voltage
%         continue
       end
       
    WT_EVP_ap_max_fall_rate_cell=FI_EVP_WT_NaI.ap_max_fall_rate{1,k};
     
    WT_EVP_all_trace_ap_max_fall_rate=[];
    for trace_k = 1:length(WT_EVP_ap_max_fall_rate_cell)
        WT_EVP_all_trace_ap_max_fall_rate= [WT_EVP_all_trace_ap_max_fall_rate,WT_EVP_ap_max_fall_rate_cell{trace_k}'];
    end 
    WT_EVP_all_cell_ap_max_fall_rate(end+1)= mean(WT_EVP_all_trace_ap_max_fall_rate, 'omitnan');  
end




%%
%Plot colors

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
%%
%AP_rise plots 

pair_WT= padmat(WT_EVP_all_cell_ap_rise', WT_DRP_all_cell_ap_rise', 2);
data_pair_het= padmat( pair_WT,SCN2A_EVP_all_cell_ap_rise', 2);
data_toplot_4= padmat(data_pair_het, SCN2A_DRP_all_cell_ap_rise',2);


UnivarScatter(data_toplot_4, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','lines','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('Action potential rise (ms)','FontSize', 20, 'LineWidth', 15)
ylim([0 1.5])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);

%%
pair_WT_apfall= padmat(WT_EVP_all_cell_ap_falls', WT_DRP_all_cell_ap_falls', 2);
data_pair_het_apfall= padmat(pair_WT_apfall,SCN2A_EVP_all_cell_ap_falls', 2);
data_toplot_5= padmat(data_pair_het_apfall, SCN2A_DRP_all_cell_ap_falls',2);


UnivarScatter(data_toplot_5, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','lines','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('Action potential fall (ms)','FontSize', 20, 'LineWidth', 15)
ylim([0 10])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);


%%
% AP Max fall rate 

pair_WT_ap_fall_rate= padmat(WT_EVP_all_cell_ap_max_fall_rate', WT_DRP_all_cell_ap_max_fall_rate', 2);
data_pair_het_ap_fall_rate= padmat(pair_WT_ap_fall_rate,SCN2A_EVP_all_cell_ap_max_fall_rate', 2);
data_toplot_6= padmat(data_pair_het_ap_fall_rate, SCN2A_DRP_all_cell_ap_max_fall_rate',2);


UnivarScatter(data_toplot_6, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','box','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('Max dV/dt during downstroke (V/s)','FontSize', 20, 'LineWidth', 15)
ylim([0 90])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);


%%
%
pair_WT_ap_rise_rate= padmat(WT_EVP_all_cell_ap_max_rise_rate', WT_DRP_all_cell_ap_max_rise_rate', 2);
data_pair_het_ap_rise_rate= padmat(pair_WT_ap_rise_rate,SCN2A_EVP_all_cell_ap_max_rise_rate', 2);
data_toplot_7= padmat(data_pair_het_ap_rise_rate, SCN2A_DRP_all_cell_ap_max_rise_rate',2);

UnivarScatter(data_toplot_7, 'BoxType','Quart',...
'Width',1,'Compression',35,'MarkerFaceColor',Colors,...
'PointSize',35,'StdColor',Colors1,'SEMColor','none',...
'Whiskers','box','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);
box off
ylabel('Max dV/dt during upstroke (V/s)','FontSize', 20, 'LineWidth', 15)
%ylim([0 170])
xticks([1 2 3 4])
set(gca,'XTickLabel',{'Control', 'DREADDs', 'Control', 'DREADDs'},'FontSize',20);
