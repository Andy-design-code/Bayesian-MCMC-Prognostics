%ʹ����ʷ����ѵ��ģ�͵�������%
clear;
file_format='*.csv'; %��������.csv��ʽ�ļ�
filepath1='E:\Datasets\PHM data challenge\2010 PHM Society Conference Data Challenge-cutter\PHM2010\c4';
filepath=fullfile(filepath1,file_format);

dir_output=dir(filepath);
% file_name is a cell array including all the acceleration data .csv files.
file_name_list={dir_output.name}';

if ~isempty(file_name_list)
    number_csvfiles=numel(file_name_list);
    for i=1:number_csvfiles
        file_name=fullfile(filepath1,file_name_list{i});
        vibration=csvread(file_name,0,0);
        fx=vibration(:,1);
        fy=vibration(:,2);
        fz=vibration(:,3);
        
        vx=vibration(:,4);
        vy=vibration(:,5);
        vz=vibration(:,6);
        
        ae=vibration(:,7);
        
        [~,save_name,~]=fileparts(file_name_list{i});
        savepath='E:\Datasets\PHM data challenge\2010 PHM Society Conference Data Challenge-cutter\PHM2010\c4_mat';
        savepath_absolute=[savepath  '\' save_name '.mat'];
        save(savepath_absolute,'fx', 'fy', 'fz', 'vx', 'vy', 'vz', 'ae');
        X=sprintf('finished the %dth file', i);
        disp(X);
    end   
else
    msgbox('empty file list')
end
  