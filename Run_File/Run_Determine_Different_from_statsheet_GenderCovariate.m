close all
clear all

%/Users/Shared/
addpath('/Volumes/workstation/code/analysis/Determining_Most_Different_Specimen/');
addpath('/Volumes/workstation/code/shared/civm_matlab_common_utils/');

%data_frame=readtable('/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/21_DAS_01_DataFrame.csv','Delimiter', 'comma');
data_frame=readtable('/Volumes/workstation/code/analysis/Determining_Most_Different_Specimen/Run_File/21_DAS_01_DataFrame_2Group_removedlightsheet.csv','Delimiter', 'comma');

%f_statistic_result='/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/F_Stat_Analysis/FeM_FeF_AirM_AirF/';
f_statistic_result='/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/F_Stat_Analysis/2Group/';

[Specimen_low,Specimen_high] = function_difference_via_statsheet_GenderCovariate(data_frame,f_statistic_result);