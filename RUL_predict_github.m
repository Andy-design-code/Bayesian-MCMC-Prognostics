% ==============================����ֶ�ִ�д��룬�Է������==================
clear
clc;
close all;
timeUnit='cut number';
% ���߱�ţ��Ե�4�ѵ���Ϊ��������Ҫ���������������ݣ����Ϊ"c3", "c5", "c6"
cutterNo='c4';
% �Ƚ�ԭʼ.csv��ʽ������"csv_to_mat"����תΪ.mat��ʽ�ļ�����ִ��������䣬·�����޸�Ϊ.mat�ļ�����·��
datastore=fileEnsembleDatastore(...
    fullfile(['E:\Datasets\PHM data challenge\2010 PHM Society Conference Data Challenge-cutter\PHM2010\', cutterNo '_mat']), ...
    '.mat');
% fx, fy, fz, vx, vy, vz, ae�ֱ�Ϊx,y,z�����������������񶯼��������ź�
datastore.DataVariables=["fx","fy","fz","vx","vy","vz","ae"];
datastore.IndependentVariables="Date";
datastore.SelectedVariables=["Date","fx","fy","fz","vx","vy","vz","ae"];
datastore.ReadFcn=@helperReadData;
datastore.WriteToMemberFcn=@helperWriteToHSBearing;
 Ttable=tall(datastore);
 head(Ttable);
 fs=50000;

%% ����չʾ��rawdataDisplay��������չʾԭʼ���ݣ�ע�͸��д���ɽ�Լ��Լ5����
%  rawdataDisplay(datastore, fs);
%% ������ȡ������featureExtraction��������������ȡ
% ============= fx ������ȡ================
% fx feature extraction
% ѡ��Ҫ��ȡ���ź�ͨ��������
SelectedVariables=["Date", "fx"];
features_fx=featureExtraction(datastore, SelectedVariables, fs);
% ============fy feature extraction============
% ѡ��Ҫ��ȡ���ź�ͨ��������
SelectedVariables = ["Date", "fy"];
features_fy=featureExtraction(datastore, SelectedVariables, fs);
% fz feature extraction
% ѡ��Ҫ��ȡ���ź�ͨ��������
SelectedVariables=["Date","fz"];
features_fz=featureExtraction(datastore, SelectedVariables, fs);
% vx feature extraction
SelectedVariables=["Date", "vx"];
features_vx=featureExtraction(datastore, SelectedVariables, fs);
% vy feature extraction
SelectedVariables=["Date", "vy"];
features_vy=featureExtraction(datastore, SelectedVariables, fs);
% vz feature extraction
SelectedVariables=["Date", "vz"];
features_vz=featureExtraction(datastore, SelectedVariables, fs);
% ae feature extraction
SelectedVariables=["Date", "ae"];
features_ae=featureExtraction(datastore, SelectedVariables, fs);
%% =====================���������б�==========================
for i=1:size(features_fx,2)
    features_fx.Properties.VariableNames{i}=['fx_' features_fx.Properties.VariableNames{i}];
    features_fy.Properties.VariableNames{i}=['fy_' features_fy.Properties.VariableNames{i}];
    features_fz.Properties.VariableNames{i}=['fz_' features_fz.Properties.VariableNames{i}];
    
    features_vx.Properties.VariableNames{i}=['vx_' features_vx.Properties.VariableNames{i}];
    features_vy.Properties.VariableNames{i}=['vy_' features_vy.Properties.VariableNames{i}];
    features_vz.Properties.VariableNames{i}=['vz_' features_vz.Properties.VariableNames{i}];
    
    features_ae.Properties.VariableNames{i}=['ae_' features_ae.Properties.VariableNames{i}];
end

featureTable=[features_fx, features_fy, features_fz, features_vx, features_vy, features_vz, features_ae];
% ================ add the "Date" in the first coloum of feature table=====
featureTableSmooth=varfun(@(x) movmean(x, [5,0]), featureTable);
datastore.SelectedVariables=["Date"];
Date=gather(tall(datastore));
featureTableSmooth=[Date, featureTableSmooth];
%% =======================����ɸѡ===========================
importscore=0.8;  % ��Ҫ�Է�����ֵ������˹Ƥ����ϵ��������������󣬵÷ִ���0.8��������ɸѡ����
breaktime=150;  % �״�Ԥ��ʱ��
breakpoint=find(featureTableSmooth.Date<=breaktime, 1, 'last');
trainData=featureTableSmooth(1: breakpoint, :);
trainData_withoutDate=trainData(:, 2:end);
% ���������̶�����
featureImportance=monotonicity(trainData, 'WindowSize', 0, 'method', 'rank');
helperSortedBarPlot(featureImportance, 'Monotonicity');
aa=featureImportance{:,:}>importscore;
featureTableSmooth_withoutDate=featureTableSmooth(:,2:end);
featureSelected=featureTableSmooth_withoutDate(1:end, aa);
trainDataSelected=trainData_withoutDate(:, aa);
%% =================����PCA�������ں�=========================
%���ɷַ��� (pca) ���ڽ�ά�������ں�%
meanTrain = mean(trainDataSelected{:,:});
sdTrain = std(trainDataSelected{:,:});
trainDataNormalized = (trainDataSelected{:,:} - meanTrain)./sdTrain;
coef = pca(trainDataNormalized);

%���þ�ֵ����׼ƫ��� pca ϵ�����������ݼ����д���%
PCA1 = (featureSelected{:,:} - meanTrain) ./ sdTrain * coef(:, 1);
PCA2 = (featureSelected{:,:} - meanTrain) ./ sdTrain * coef(:, 2);

%��ǰ6����Ҫ����Ŀռ��п��ӻ����ݡ�%
figure
hold on
numData = size(featureSelected, 1);

plot(1:numData,PCA1)
plot(1:numData, PCA2)

hold off
xlabel('Cut Number')
ylabel('Feature Value')
% legend('PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6')
title('PCA')

%��ǰ������Ҫ����Ŀռ��п��ӻ����ݡ�%
figure
numData = size(featureSelected, 1);
scatter(PCA1, PCA2, [], 1:numData, 'filled')
xlabel('PCA 1')
ylabel('PCA 2')
cbar = colorbar;
ylabel(cbar, ['Time (' timeUnit ')'])

%% ======================����һ��ָ��ƽ��������������===================
healthIndicator = PCA1;
%healthIndicator = movmean(healthIndicator, 10);

healthIndicator = PCA1;
m=size(featureSelected, 1);
alpha = 0.1;
s1 = zeros(m,1);
s1(1,1) = healthIndicator(1,1);
for i=2:m
s1(i) = alpha*healthIndicator(i,1)+(1-alpha)*s1(i-1);
end
healthIndicator = s1;

s2 = zeros(m,1);
s2(1,1) = healthIndicator(1,1);
for i=2:m
s2(i) = alpha*s1(i,1)+(1-alpha)*s2(i-1);
end
healthIndicator = s2;

m=size(featureSelected, 1);
% ƽ�������ķ�ΧΪ[0,1]
s3 = zeros(m,1);
s3(1,1) = healthIndicator(1,1);
for i=2:m
s3(i) = alpha*s2(i,1)+(1-alpha)*s3(i-1);
end
healthIndicator = s3;

% ���ӻ�����ָ��%
healthIndicator = healthIndicator - healthIndicator(1);
figure
numData = size(featureSelected, 1);
plot(1:numData, healthIndicator)
xlabel('cut number')
title('Health Indicator')
box
%% =======================RUL Ԥ��====================
%===============Bayesian-MCMC�����ܣ����Ĵ���==============
global y ParamName ny np
y_total=healthIndicator;

para0=[4; 0.1]; % a, b�ĳ�ֵ
weigh=[0.5; 0.5]; % ����ֲ�Ϊ���ȷֲ���������ȷ�Χ[-0.5, 0.5]
TimeUnit='Cycles';
time=[1:length(y_total)]';
threshold=healthIndicator(end);
 
 % define the name of parameters
 ParamName=['a'; 'b'];
 
 % define the start time from where to predict the RUL
 startime=breaktime;
 
 % parameters used in MCMC sampling
 ns=5e4;  % ����������Ϊ50000��
 burnIn=0.2; % �յ�ǰ20%�����ȶ�����
 
 % start the predict RUL from startime until 315
 ny=startime;
 np=size(ParamName,1);
 para_old=para0;
 RUL_cell=cell(length(time),1);
 while ny<time(end)
     y=y_total(1:ny);
     ny=length(y);
     samples(:,1)=para_old; %initial samples
     % MODEL�����м���a,b���������ϸ����ܶȺ���ֲ�
     [~, jpdf_old]=MODEL(para_old, time(1:ny)); % initial joint pdf
     t=time(1:ny);
     % MCMC��������ʵ��MCMC����
     [thetaHat_unique, para_old, jpdf_old]=MCMC(para_old, jpdf_old, weigh, ns, burnIn, t);
    
     % predict RUL at current time ny
     for i=1:size(thetaHat_unique,2)
         k=ny+1;
         a=thetaHat_unique(1,i);
         b=thetaHat_unique(2,i);
         z=0;
         flag=0;
         while z<threshold
             z=-1+a.*exp(b.*k);
             k=k+1;
             if z<0
                 flag=1;
                 break;
             end
         end
         RUL(:,i)=(k-1)-ny;
     end
     RUL_cell{ny}=RUL'; 
     RUL=[];
     ny=ny+1;
     X=sprintf('finished the %d measurement data', ny);
     disp(X);
 end

%% ==============��ӡ���=========================
RUL_statis=[];
 for i=startime:length(time)-1
     rul=RUL_cell{i};
     RUL_statis=cat(1, RUL_statis, prctile(rul, [5 50 95]));
 end
RUL_statis(:,4)=[startime:length(time)-1]';
% plot the 5% and 95% condidence interval
FaceColor=[0.1010 0.7450 0.980];
downsample_interval=5;
Y=[RUL_statis(1:downsample_interval:end-1, 1), RUL_statis(1:downsample_interval:end-1, 3)-RUL_statis(1:downsample_interval:end-1,1)];
X=[RUL_statis(1:downsample_interval:end-1,4), RUL_statis(1:downsample_interval:end-1,4)];
figure
hold on
ar=area(X,Y);
ar1=ar(1); ar2=ar(2);
ar1.FaceAlpha=0;
ar1.EdgeAlpha=0;
ar2.EdgeAlpha=0;
ar2.FaceAlpha=0.4;
ar2.FaceColor=FaceColor;

linewidth=1.5;
 % plot the real RUL
 plot([0 length(time)], [length(time), 0], 'linewidth', linewidth)

 RUL_statis(:,4)=[startime:length(time)-1];
 % plot the 50% percent estimated RUL
 plot(RUL_statis(1:end-1, 4), RUL_statis(1:end-1, 2), 'color', 'black', 'linewidth', linewidth, 'marker', '*', 'markersize',5);
 
 accuracy_zone=0.2;
 true_rul=[length(time)-startime:-1:1];
 alpha_lamda=[true_rul+accuracy_zone*true_rul; true_rul-accuracy_zone*true_rul]';
 plot(startime:time(end-1), alpha_lamda, 'linewidth', linewidth, 'color', [0 0.5 0], 'linestyle', '--');
 
 box on;
 h=get(gca,'children');
 legend([h(3), h(1), h(4), h(5)], 'RULԤ��ֵ', '\alpha=\pm20%','RUL��ʵֵ', '95%��������');
 set(gca,'fontsize',12, 'fontweight','bold')
 xlim([100,315])
 xlabel('��������')
 ylabel('RUL')

 xlim([150, time(end)])
 ylim([0,200])

 % ����rmse��mae, �ӵ�200ʱ�̿�ʼ����
loc = find(true_rul == time(end)-200);
er=[true_rul(loc : end)]'-RUL_statis(loc:end, 2);
rmse = sqrt(sum((er).^2)/length(er));
mae=sum(abs(er))/length(er);
%% MCMC ��������
function [thetaHat_unique, para_output, jpdf_output]=MCMC(para_old, jpdf_old, weigh, ns, burnIn, time)
    for i=2:ns/(1-burnIn)
        % �ӽ���ֲ������ȷֲ����в����µĲ�������
        para_new=unifrnd(para_old-weigh, para_old+weigh);
        % 
        [~, jpdf_new]=MODEL(para_new, time);
        % ����-�ܾ�׼��(acceptance/rejection
        % criterion)���������׼���������������ͬʱ��f_po(theta_new | Y{1:T})����f_po(theta_old | Y{1:T})��������һѭ��
        % ���򣬾ܾ�������������������
        if rand<min(1, jpdf_new/jpdf_old)
            para_old=para_new;
            jpdf_old=jpdf_new;
        end
        samples(:, i)=para_old;
    end  
    para_output=para_old;
    jpdf_output=jpdf_old;
    nBurn=ns/(1-burnIn)-ns;
    thetaHat=samples(:, nBurn+1:end);
    thetaHat_unique=[unique(thetaHat(1,:), 'stable'); unique(thetaHat(2,:), 'stable')];
end

%% ���㱴Ҷ˹��������ܶȺ���
function [z, poste]=MODEL(param,t)
    global y ParamName ny np
    for j=1:np
        eval([ParamName(j,:) '=param(j,:);']);
    end
    s=2.9;
    %=========�˻�ģ�Ͷ���========
    z= -1+a.*exp(b.*t);   
    if length(y)~=length(t); post=[];
    else
        % ����������������ϸ����ܶ�ֵ���μ������й�ʽ(4)
        prior=lognpdf(param(1), 1, 1e2)*normpdf(param(2), 1,1e2);
        % ������Ȼ�������μ������й�ʽ(5)
        likel=(1./(sqrt(2.*pi).*s)).^ny.*...
            exp(-0.5./s.^2.*(y-z)' * (y-z));
        % ��������ĺ������ϸ����ܶ�ֵ���μ������й�ʽ(7)
        poste=likel.*prior;  % compute the posterior value
    end
end









 