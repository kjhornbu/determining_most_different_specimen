close all
clear all

addpath('/Users/Shared/workstation/code/analysis/Determining_Most_Different_Specimen/');
addpath('/Users/Shared/workstation/code/shared/civm_matlab_common_utils/');

data_frame=readtable('/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/21_DAS_01_DataFrame.csv','Delimiter', 'comma');
%data_frame=readtable('/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/21_DAS_01_DataFrame_2Group.csv','Delimiter', 'comma');

f_statistic_result='/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/F_Stat_Analysis/FeM_FeF_AirM_AirF/';
%f_statistic_result='/Volumes/dusom_civm-kjh60/All_Staff/21.DAS.01/F_Stat_Analysis/2Group/';

[Specimen] = function_difference_via_statsheet(data_frame,f_statistic_result);

