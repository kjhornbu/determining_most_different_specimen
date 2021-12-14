clear all
close all

data_frame=readtable('/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/21_DAS_01_DataFrame.csv','Delimiter', 'comma');
data_frame=readtable('/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/21_DAS_01_DataFrame_2Group.csv','Delimiter', 'comma');

f_statistic_result='/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/F_Stat_Analysis/FeM_FeF_AirM_AirF/';
f_statistic_result='/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/F_Stat_Analysis/2Group/';

%% DO NOT EDIT BELOW THIS LINE
contrast={'FA','AD','RD','MD','NQA','QA','Volume','Normalized Volume'};

[~,~,group_idx]=unique(data_frame.group,'stable');
%Adding Gender on Top (Interested in Opposite Gender and Same Gender)
[~,~,gender_idx]=unique(data_frame.gender,'stable');

%check if gender is a subset of group by seeing if group subsets have
%gender varying within it or not... if mostly varying within it then fine.

if sum(group_idx==1)>sum(group_idx==1 & gender_idx==1) %Just trying the first one for now really want this to be a smarter checker later on...
    
else
    error('You cannot add a gender covariate when you already have it added as a grouping condition!')
end

    for m=1:max(group_idx)
        total_group(m)=sum(group_idx==m); 
    end
    for m=1:max(gender_idx)
        total_gender(m)=sum(gender_idx==m);
    end
    
    count_matrix_low=zeros(size(total_group,2),max(total_group));
    count_matrix_high=zeros(size(total_group,2),max(total_group));
    
%% Pull Most Changed From Scalar Study Statistics
%Load all the Data Sheets
for m=1:size(data_frame,1)
    Data{m}=readtable(char(data_frame.stat_path(m)));
end

N=0; %index of terms considered in dissimilar function

for n=1:size(contrast,2)
    %Find significant regions from F-statistics Analysis
    f_stat_siginficant=readtable(strcat(f_statistic_result,'Scalar_',contrast{1,n},'_Reduced_DataTable_wPercentChangeTests.csv'));
    
    A=strcmp(contrast{1,n},'Volume');
    B=strcmp(contrast{1,n},'Normalized Volume');
            
    roi_significant=unique(f_stat_siginficant.roi); %Pull the ROI we will care about
    
    N=size(roi_significant,1)+N; %number of points in a given contrast.
    
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
        
       dissimilar(n,o,:)= ((Compare_Data-mean(Compare_Data))/std(Compare_Data)).^2; %number of standard deviations above/below mean
        
       [~,Group_Mean_Order]=sort(Group_mean); %low to High order
       
       %Smallest Term Iterating (there is always a 1 in the group because we load from 1 to m)
       Runno_idx_low(n,o,:)=[Group_Mean_Order(1), Group_Order(1,Group_Mean_Order(1))]; %indicate first idx as the group second idx as the specimen
       count_matrix_low(Group_Mean_Order(1),Group_Order(1,Group_Mean_Order(1)))=count_matrix_low(Group_Mean_Order(1),Group_Order(1,Group_Mean_Order(1)))+1;
       
       %Largest Term Iterating (not always at end position because the different sized groups for various conditions)
       if Group_Order(end,Group_Mean_Order(end)) == 0
           
           %find highest position with first non-zero (if number of specimen not always
           %equal among groups)
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

RMSE=sqrt(squeeze(sum(sum(dissimilar)))/N); %finishing the RMSE calculation -> sum up difference from sum components divide by the number of possible 

%reshape RSME by group condition
for m=1:max(group_idx)
    [Group_RMSE(1:total_group(m),m),Group_RMSE_Order(1:total_group(m),m)]=sort(RMSE(group_idx==m));
    
    Group_RMSE_Order(1:total_group(m),m)=sum(total_group(1:m-1))+Group_RMSE_Order(1:total_group(m),m);
end

reshaped_gender_idx=reshape(gender_idx,[max(total_group), size(total_group,2)]); %rehaped into groups (second index) the full gender indexing.
reshaped_gender_idx=reshaped_gender_idx(Group_RMSE_Order);

for m=1:max(group_idx)
    for n=1:max(gender_idx)
        total_reshaped_gender(m,n)=sum(reshaped_gender_idx(:,m)==n); %group (r) ,gender(c)
    end
end

 Group_Gender_RMSE=zeros(max(group_idx),max(gender_idx),max(max(total_reshaped_gender)));
 Group_Gender_RMSE_Order=zeros(max(group_idx),max(gender_idx),max(max(total_reshaped_gender)));
 
%THEN separate RSME by gender condition
for m=1:max(group_idx)
    for n=1:max(gender_idx)
        Group_Gender_RMSE(m,n,1:total_reshaped_gender(m,n))=Group_RMSE(reshaped_gender_idx(:,m)==n,m);
        Group_Gender_RMSE_Order(m,n,1:total_reshaped_gender(m,n))=Group_RMSE_Order(reshaped_gender_idx(:,m)==n,m);
    end
end

%Picking the lowest value specimen across all conditions.
[r,c]=find(count_matrix_low==max(max(count_matrix_low)));
Lowest_specimen=data_frame.specimen(sum(total_group(1:r-1))+c)

%next closest for the opposite gender condition
[group,gender,~]=ind2sub([size(Group_Gender_RMSE_Order)],find(Group_Gender_RMSE_Order==sum(total_group(1:r-1))+c));

if gender==1
    Lowest_specimen{2}=data_frame.specimen(Group_Gender_RMSE_Order(Group_Gender_RMSE==max(Group_Gender_RMSE(group,2,:))));
else
    Lowest_specimen{2}=data_frame.specimen(Group_Gender_RMSE_Order(Group_Gender_RMSE==max(Group_Gender_RMSE(group,1,:))));
end

%Picking the highest value specimen across all conditions.
[r,c]=find(count_matrix_high==max(max(count_matrix_high)));
Highest_specimen=data_frame.specimen(sum(total_group(1:r-1))+c)

[group,gender,~]=ind2sub([size(Group_Gender_RMSE_Order)],find(Group_Gender_RMSE_Order==sum(total_group(1:r-1))+c));

%next closest for the opposite gender condition
if gender==1
    Highest_specimen{2}=data_frame.specimen(Group_Gender_RMSE_Order(Group_Gender_RMSE==max(Group_Gender_RMSE(group,2,:))));
else
    Highest_specimen{2}=data_frame.specimen(Group_Gender_RMSE_Order(Group_Gender_RMSE==max(Group_Gender_RMSE(group,1,:))));
end







