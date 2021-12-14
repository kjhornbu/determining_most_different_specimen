function [Specimen] = function_difference_via_statsheet(data_frame,f_statistic_result)
contrast={'FA','AD','RD','MD','NQA','QA','Volume','Normalized Volume'}; 

[~,~,group_idx]=unique(data_frame.group,'stable');

    for m=1:max(group_idx)
        total_group(m)=sum(group_idx==m);
    end
        
    count_matrix_low=zeros(size(total_group,2),max(total_group));
    count_matrix_high=zeros(size(total_group,2),max(total_group));
    
%% Pull Most Changed From Scalar Study Statistics
%Load all the Data Sheets
for m=1:size(data_frame,1)
    Data{m}=readtable(char(data_frame.stat_path(m)));
end

for n=1:size(contrast,2)
    %Find significant regions from F-statistics Analysis
    f_stat_siginficant=readtable(strcat(f_statistic_result,'Scalar_',contrast{1,n},'_Reduced_DataTable_wPercentChangeTests.csv'));
    
    A=strcmp(contrast{1,n},'Volume');
    B=strcmp(contrast{1,n},'Normalized Volume');
            
    roi_significant=unique(f_stat_siginficant.roi); %Pull the ROI we will care about
    
%% how many significant ROI for a given contrast
    for o=1:size(roi_significant,1) 
        
    %Get Stat Sheets for each, with the roi signficant data
        for m=1:size(data_frame,1)
            if A==1
                Data_cells=regexpi(Data{m}.Properties.VariableNames,strcat('^',lower(contrast{1,n}),'_mm3')); %For Straight Volume values
                Data_idx=~cellfun(@isempty,Data_cells);
                Compare_Data(m,1)=sum(table2array(Data{m}([find(Data{m}.ROI==roi_significant(o)),find(Data{m}.ROI==roi_significant(o)+1000)],Data_idx))); %Get Actual Data
            elseif B==1
                Data_cells=regexpi(Data{m}.Properties.VariableNames,strcat('^',lower('Volume'),'_mm3')); %For Normalized Volume values
                Data_idx=~cellfun(@isempty,Data_cells);
                total_volume=sum(table2array(Data{m}(Data{m}.ROI>0,Data_idx))); %all regions besides exterior!
                Compare_Data(m,1)=sum(table2array(Data{m}([find(Data{m}.ROI==roi_significant(o)),find(Data{m}.ROI==roi_significant(o)+1000)],Data_idx)))/total_volume; %Get Actual Data
            else
                Data_cells=regexpi(Data{m}.Properties.VariableNames,strcat('^',contrast{1,n},'_mean')); %For non-volume values
                Data_idx=~cellfun(@isempty,Data_cells);
                Compare_Data(m,1)=mean(table2array(Data{m}([find(Data{m}.ROI==roi_significant(o)),find(Data{m}.ROI==roi_significant(o)+1000)],Data_idx))); %Get Actual Data
            end
        end
        
        %Compare Mean location for each specimen on test condition
        for m=1:max(group_idx)
            [~,Group_Order(1:total_group(m),m)]=sort(Compare_Data(group_idx==m)); %low to High order (column is group, row  is specimen)
           
            Group_mean(m)=mean(Compare_Data(group_idx==m)); %mean group value
        end
      
       [~,Group_Mean_Order]=sort(Group_mean); %low to High order
       
       %Smallest Term (there is always a 1 in the group because we load
       %from 1 to m)
       Runno_idx_low(n,o,:)=[Group_Mean_Order(1), Group_Order(1,Group_Mean_Order(1))]; %indicate first idx as the group second idx as the specimen
       count_matrix_low(Group_Mean_Order(1),Group_Order(1,Group_Mean_Order(1)))=count_matrix_low(Group_Mean_Order(1),Group_Order(1,Group_Mean_Order(1)))+1;
       
       %Largest Term
       if Group_Order(end,Group_Mean_Order(end)) == 0
           
           %find highest position with first non-zero (if number of specimen not always equal among groups)
           zero_condition=Group_Order(:,Group_Mean_Order(end))==0;
           order_data=Group_Order(~zero_condition,Group_Mean_Order(end)); %only pull from the order_data the final position
           
           Runno_idx_high(n,o,:)=[Group_Mean_Order(end) order_data(end)]; %indicate first idx as the group second idx as the specimen (r,c)
           count_matrix_high(Group_Mean_Order(end),order_data(end))=count_matrix_high(Group_Mean_Order(end),order_data(end))+1;
       else
           Runno_idx_high(n,o,:)=[Group_Mean_Order(end) Group_Order(end,Group_Mean_Order(end))]; %indicate first idx as the group second idx as the specimen (r,c)
           count_matrix_high(Group_Mean_Order(end),Group_Order(end,Group_Mean_Order(end)))=count_matrix_high(Group_Mean_Order(end),Group_Order(end,Group_Mean_Order(end)))+1;
       end
    end
end

%Finding the Lowest Data by the maximum number
[r,c]=find(count_matrix_low==max(max(count_matrix_low)));
%Picking the lowest value specimen across all conditions.
Specimen{1}=data_frame.specimen(sum(total_group(1:r-1))+c);

%Finding the Lowest Data by the maximum number
[r,c]=find(count_matrix_high==max(max(count_matrix_high)));
%Picking the highest value specimen across all conditions.
Specimen{2}=data_frame.specimen(sum(total_group(1:r-1))+c);
end

