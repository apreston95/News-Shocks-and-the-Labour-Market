% KurmannSims_Part2_Figures.m
%
% Utility to plot figures 2, 4, 5 and 6 (and figures 3-5 in appendix) as they appear in the 
%
%   Kurmann, A. and E. Sims (2019). "Revisions in Utilization-Adjusted TFP and Robust Identification of News Shocks." Review of Economics and Statistics (forthcoming).
%
% Based on results generated and saved in KurmannSims_Part2.m. 
% Run after estimating appropriate VAR in KurmannSims_Part2.m.
%
% Written by Andre Kurmann (Drexel University) and Eric Sims (Notre Dame
% and NBER)
%
% Last modified: Andre Kurmann, January 7, 2020
%
% Use of code for research purposes is permitted as long as proper reference to source is given. 
%
%------------------------------------------------------------------------------------------------------------------------------------------------

close all

%% Code to plot Figures 2, 4 and 5 (and Figures 3-5 in appendix)

%load appropriate data files (after running KurmannSims_Part2.m with the
%same specification)
load('data_ms80_VARlarge_tfpu2016');   %load FEV&IRF data from VAR estimation(as saved in KurmannSims_Part2.m)
data2007=data;
load('data_ms80_VARlarge_tfpu2007');   %load FEV&IRF data from VAR estimation (as saved in KurmannSims_Part2.m)
data2016=data;

nimp=40;
R=round(nvars/2);
figure('Name','Figure 2 (or 4 or 5)'); 
for n=1:nvars;
        subplotn = subplot(R,2,n);  
        set(subplotn,'FontName','Times New Roman','FontSize',12,'Layer','top','YGrid','on');
                                    if ndims(data)==4;   %if available, plot confidence bands
                                        grpyat = [(1:nimp)' data2016(1:nimp,nvars+n,1,2); (nimp:-1:1)' data2016(nimp:-1:1,nvars+n,1,3)];
                                        patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]); 
                                    end
                                    hold on;
                                    k = plot(1:nimp,data2016(1:nimp,nvars+n,1,1),'-k','LineWidth',1.5); 
                                    r = plot(1:nimp,data2007(1:nimp,nvars+n,1,1),'-.r','LineWidth',1.5);
                                    plot(1:nimp,data2007(1:nimp,nvars+n,1,2),'--r','LineWidth',1.2);
                                    plot(1:nimp,data2007(1:nimp,nvars+n,1,3),'--r','LineWidth',1.2);
                                    plot(1:nimp,zeros(nimp),':k');
                                    grid on
                                    title(vars(n,:),'FontSize',12) 
                                    ylabel('percent','FontSize',12)
                                    xlabel('quarters','FontSize',12)
                                    axis tight
                                    hold off
  end
%   legend([k,r],'2016 TFPu vintage','2007 TFPu vintage','Location','southeast')
    legend([k,r],'2016 TFP vintage','2007 TFP vintage','Location','southeast')
%     legend([k,r],'2016 TFPu vintage','2013 TFPu unfiltered','Location','southeast')
%     legend([k,r],'Fernald labor productivity','Aggregate labor productivity','Location','southeast')   


%% Code to plot panel a) of Figure 6

%Data for innovation variables 
load('data_ms80_VARlarge_SSOICT_tfpu2007');  %load data from slope shock identification (as saved in identm.m)
data1_2007=data;
load('data_ms80_VARlarge_SSOICT_tfpu2016');   %load data from news shock identification (as saved in identm.m)
data1_2016=data;
load('data_ms80_VARannual_AlexTECH2_tfpu2007');  %load data from slope shock identification (as saved in identm.m)
data2_2007=data;
load('data_ms80_VARannual_AlexTECH2_tfpu2016');   %load data from news shock identification (as saved in identm.m)
data2_2016=data;
load('data_ms80_VARlarge_RD_tfpu2007');  %load data from slope shock identification (as saved in identm.m)
data3_2007=data;
load('data_ms80_VARlarge_RD_tfpu2016');   %load data from news shock identification (as saved in identm.m)
data3_2016=data;
load('data_ms80_VARlarge_PcPi_tfpu2007');  %load data from slope shock identification (as saved in identm.m)
data4_2007=data;
load('data_ms80_VARlarge_PcPi_tfpu2016');   %load data from news shock identification (as saved in identm.m)
data4_2016=data;
varnames = {'ICT standardization';'New technology book titles';'Real R&D expenditures';'Relative investment price (inverse)';};

varselect = 8;  %position of variable to use in 8-variable large VAR
nhorizon = size(data1_2007,1);
IRFselect_1 = zeros(nhorizon,4,1,3);
IRFselect_1(:,1,:,:) = data1_2016(:,nvars+varselect,1,:); 
IRFselect_1(1:20,2,:,:) = data2_2016(:,4+2,1,:);       %for annual 4-variable VAR, tech variable is in 2nd position 
IRFselect_1(:,3,:,:) = data3_2016(:,nvars+varselect,1,:);
IRFselect_1(:,4,:,:) = data4_2016(:,nvars+varselect,1,:); 

IRFselect_2 = zeros(nhorizon,4,1,3);
IRFselect_2(:,1,:,:) = data1_2007(:,nvars+varselect,1,:);   
IRFselect_2(1:20,2,:,:) = data2_2007(:,4+2,1,:);       %for annual 4-variable VAR, tech variable is in 2nd position    
IRFselect_2(:,3,:,:) = data3_2007(:,nvars+varselect,1,:);
IRFselect_2(:,4,:,:) = data4_2007(:,nvars+varselect,1,:); 

N=size(IRFselect_1,2);
R=round(N/2);
figure('Name','Figure 6a'); 
for n=1:N;
        subplotn = subplot(R,2,n);  
        set(subplotn,'FontName','Times New Roman','FontSize',12,'Layer','top','YGrid','on');             
                                    if n==2     
                                        nimp=10;    %for annual VAR results, choose nimp = 10 years
                                    else
                                        nimp=40;    %for quarterly VAR, choose nimp = 40 quarters as usual
                                    end
        
                                    if ndims(data)==4;   %if available, plot confidence bands
                                        grpyat = [(1:nimp)' IRFselect_1(1:nimp,n,1,2); (nimp:-1:1)' IRFselect_1(nimp:-1:1,n,1,3)];
                                        patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]); 
                                    end
                                    hold on;
                                    k = plot(1:nimp,IRFselect_1(1:nimp,n,1,1),'-k','LineWidth',1.5); 
                                    r = plot(1:nimp,IRFselect_2(1:nimp,n,1,1),'-.r','LineWidth',1.5);
                                    plot(1:nimp,IRFselect_2(1:nimp,n,1,2),'--r','LineWidth',1.2);
                                    plot(1:nimp,IRFselect_2(1:nimp,n,1,3),'--r','LineWidth',1.2);
                                    plot(1:nimp,zeros(nimp),':k');
                                    grid on;
                                    title(varnames(n,:),'FontSize',14) 
                                    ylabel('percent','FontSize',12)
                                    if n==2;
                                        xlabel('years','FontSize',12)
                                    else
                                        xlabel('quarters','FontSize',12)
                                    end
                                    axis tight
                                    hold off
 end
legend([k,r],'2016 TFP vintage','2007 TFP vintage')   


%% Code to plot panel b) of Figure 6

%Data for news variables
load('data_ms80_VARlarge_tfpu2007');  %load FEV&IRF data from VAR estimation with 2007 vintage (as saved in KurmannSims_Part2.m)
data1_2007=data;
load('data_ms80_VARlarge_tfpu2016');   %load FEV&IRF data from VAR estimation with 2006 vintage (as saved in KurmannSims_Part2.m)
data1_2016=data;
load('data_ms80_VARlarge_spread_tfpu2007');  %load FEV&IRF data from VAR estimation with 2007 vintage (as saved in KurmannSims_Part2.m)
data2_2007=data;
load('data_ms80_VARlarge_spread_tfpu2016');   %load FEV&IRF data from VAR estimation with 2016 vintage (as saved in KurmannSims_Part2.m)
data2_2016=data;
load('data_ms80_VARlarge_E5Y_tfpu2007');  %load FEV&IRF data from VAR estimation with 2007 vintage (as saved in KurmannSims_Part2.m)
data3_2007=data;
load('data_ms80_VARlarge_E5Y_tfpu2016');  %load FEV&IRF data from VAR estimation with 2016 vintage (as saved in KurmannSims_Part2.m)
data3_2016=data;
load('data_ms80_VARlarge_BOS_tfpu2007');  %load FEV&IRF data from VAR estimation with 2007 vintage (as saved in KurmannSims_Part2.m)
data4_2007=data;
load('data_ms80_VARlarge_BOS_tfpu2016');   %load FEV&IRF data from VAR estimation with 2016 vintage (as saved in KurmannSims_Part2.m)
data4_2016=data;
varnames = {'Real S&P500 index';'5year-FFR Spread';'Consumer confidence';'Business confidence'};

varselect = 8;  %position of variable to use in 8-variable large VAR
nhorizon = size(data1_2007,1);
IRFselect_1 = zeros(nhorizon,4,1,3);
IRFselect_1(:,1,:,:) = data1_2016(:,nvars+2,1,:); %S&P500 index is in 2nd position of baseline 8-variable VAR
IRFselect_1(:,2,:,:) = data2_2016(:,nvars+varselect,1,:);        
IRFselect_1(:,3,:,:) = data3_2016(:,nvars+varselect,1,:);
IRFselect_1(:,4,:,:) = data4_2016(:,nvars+varselect,1,:); 

IRFselect_2 = zeros(nhorizon,4,1,3);
IRFselect_2(:,1,:,:) = data1_2007(:,nvars+2,1,:);  %S&P500 index is in 2nd position of baseline 8-variable VAR
IRFselect_2(:,2,:,:) = data2_2007(:,nvars+varselect,1,:); 
IRFselect_2(:,3,:,:) = data3_2007(:,nvars+varselect,1,:);
IRFselect_2(:,4,:,:) = data4_2007(:,nvars+varselect,1,:); 

nimp=40;
N=size(IRFselect_1,2);
R=round(N/2);
figure('Name','Figure 6b'); 
for n=1:N;
        subplotn = subplot(R,2,n);  
        set(subplotn,'FontName','Times New Roman','FontSize',12,'Layer','top','YGrid','on');            
                                    if ndims(data)==4;   %if available, plot confidence bands
                                        grpyat = [(1:nimp)' IRFselect_1(1:nimp,n,1,2); (nimp:-1:1)' IRFselect_1(nimp:-1:1,n,1,3)];
                                        patch(grpyat(:,1),grpyat(:,2),[0.7 0.7 0.7],'edgecolor',[0.65 0.65 0.65]); 
                                    end
                                    hold on;
                                    k = plot(1:nimp,IRFselect_1(1:nimp,n,1,1),'-k','LineWidth',1.5); 
                                    r = plot(1:nimp,IRFselect_2(1:nimp,n,1,1),'-.r','LineWidth',1.5);
                                    plot(1:nimp,IRFselect_2(1:nimp,n,1,2),'--r','LineWidth',1.2);
                                    plot(1:nimp,IRFselect_2(1:nimp,n,1,3),'--r','LineWidth',1.2);
                                    plot(1:nimp,zeros(nimp),':k');
                                    grid on;
                                    title(varnames(n,:),'FontSize',14) 
                                    ylabel('percent','FontSize',12)
                                    xlabel('quarters','FontSize',12)
                                    axis tight
                                    hold off
end
legend([k,r],'2016 TFP vintage','2007 TFP vintage')   


