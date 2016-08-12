clear all, clc, close all
cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff16\am14341\Documents\RATES\R Code\data\AIS_Data\GIA_models')
importdata('GIA_REAR_EEMD.dat') %in meters per year
GIA=ans.data*1000;
importdata('GIA_REAR_EEMD_std.dat')
GIA_std = ans.data*1000;

%%Sens Analysis using low density over WA
%xlsread('\\ads.bris.ac.uk\filestore\MyFiles\Staff16\am14341\Documents\RATES\R Code\img\rock_dens\GIA_dens_lowWA.xls')
%GIA=ans;
%GIA(:,3:4)=GIA(:,3:4).*1000;


%GIA=xlsread('\\ads.bris.ac.uk\filestore\MyFiles\Staff16\am14341\Documents\RATES\R Code\img\GIA\250\GIA250.xlsx')
%Uploading GPS %in mm/yr
GPS = xlsread('C:\Users\am14341\Dropbox\WORK\GIA\GPS-velocities\GPS-velocities-MK-500.xlsx');
 
%names = GPS.textdata(2:69);
%GPS = GPS.data; 
cd('O:\Documents\RATES\R Code\matlab\polarstereo_fwd')
[x1,y1]=polarstereo_fwd(GPS(:,1), GPS(:,2),6378137.0,axes2ecc(6378137.0, 6356752.3),-71,0);
GPS_data = [GPS, [x1,y1]]; 

%DUM correction is inlcuded
 
%%%%%%%%%%%%%%%%%
%%Attributing GIA value to GPS location

F_G = TriScatteredInterp(GIA(:,1)./1000,GIA(:,2)./1000,GIA(:,3),'nearest'); %coordinates are in m, so i divide by 1000
GPS_data(:,7) = F_G(GPS_data(:,5)./1000,GPS_data(:,6)./1000);

F_E = TriScatteredInterp(GIA_std(:,1)./1000,GIA_std(:,2)./1000,GIA_std(:,3),'nearest');
GPS_data(:,8) = F_E(GPS_data(:,5)./1000,GPS_data(:,6)./1000);

%GPS_data(16,:) = []

clear F_E 
clear F_g
clear GIA
clear GIA_std
cd('O:\Documents\RATES\DATA_processing\ELASTIC_REAR\Error_Elastic')

load elastic-error2.mat
%V(61)=[]; 

% Compute WRMS
% wi = 1./(GPS_data(:,4).^2); % std gps 
error = sqrt(GPS_data(:,4).^2+V );
wi = 1./(error.^2); % std gps ,

%%%% THis is the WRMS value without spatial correlation
% D = GPS_data(:,7)-GPS_data(:,3); 
% F_E = sqrt((sum(wi.*D.^2))/sum(wi)); % gia - gps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(1,2,1)
% hist(Misf500(:,2),10)
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% %set(gca, 'FontSize', fsz);
% set(gca,'TickLength',[0.02,0.02]);
% set(gcf, 'PaperPositionMode', 'auto');
% axis([-10 10, 0 30])
% title('Misfits using l=500km')
% 
% subplot(1,2,2)
% hist(Misf250(:,2),10)
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% %set(gca, 'FontSize', fsz);
% set(gca,'TickLength',[0.02,0.02]);
% set(gcf, 'PaperPositionMode', 'auto');
% axis([-10 10, 0 30])
% title('Misfits using l=250km')

%%%%%%Compare GIA models
cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff16\am14341\Documents\RATES\R Code\data\AIS_Data\GIA_models')

  importdata('whitehouse_full.dat')
  W12=ans.data;
  F = TriScatteredInterp(W12(:,1),W12(:,2),W12(:,3),'nearest');
  GPS_W12 = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);
 
 importdata('GIA_Gunter_res.dat')
 GUN=ans.data;
 F = TriScatteredInterp(GUN(:,1),GUN(:,2),GUN(:,3),'nearest');
 GPS_GUN = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);
 
 importdata('GIA_ICE6G_res_V2.dat')
 ICE=ans.data;
 F = TriScatteredInterp(ICE(:,1),ICE(:,2),ICE(:,3),'nearest');
 GPS_ICE = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);
 
 importdata('GIA_A_res.dat');
 A=ans.data;
 F = TriScatteredInterp(A(:,1),A(:,2),A(:,3),'nearest');
 GPS_A = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);
 
 importdata('ij05_R2115.dat');
 IJ05=ans.data;
 F = TriScatteredInterp(IJ05(:,1),IJ05(:,2),IJ05(:,3),'nearest');
 GPS_IJ05 = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);
 
 importdata('ricardo_full.dat');
 RI=ans.data;
 F = TriScatteredInterp(RI(:,1),RI(:,2),RI(:,3),'nearest');
 GPS_RIVA = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);
 
 importdata('ingo_full_header.dat');
 AGE=ans.data;
 F = TriScatteredInterp(AGE(:,1),AGE(:,2),AGE(:,3),'nearest');
 GPS_AGE = F(GPS_data(:,5)/1000,GPS_data(:,6)/1000);

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%Statistics original

%%%%%%%%%%%%%%%%%%

%% Compute WRMS
%%wi = 1./(GPS_data(:,4).^2+GPS_data(:,8).^2); % std gps , std gia
%wi = 1./(GPS_data(:,4).^2); % std gps ,

%%%%%%%%%%%%%%%
%%%% THis is the WRMS value without spatial correlation
 D = GPS_data(:,7)-GPS_data(:,3); % Discrepancy in RATES
% 
% F_RATES = sqrt((sum(wi.*D.^2))/sum(wi)); % WRMS gia - gps 
 GPS_RATES=[GPS_data(:,7),D];
%  
% F_W12 = sqrt((sum(wi.*(GPS_W12(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_W12(:,2) = GPS_W12-GPS_data(:,3); 
% 
% F_GUN = sqrt((sum(wi.*(GPS_GUN(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_GUN(:,2) = GPS_GUN-GPS_data(:,3); 
% 
% F_ICE = sqrt((sum(wi.*(GPS_ICE(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_ICE(:,2) = GPS_ICE-GPS_data(:,3); 
% 
% F_A = sqrt((sum(wi.*(GPS_A(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_A(:,2) = GPS_A-GPS_data(:,3); 
% 
% F_IJ05 = sqrt((sum(wi.*(GPS_IJ05(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_IJ05(:,2) = GPS_IJ05-GPS_data(:,3); 
% 
% F_RIVA = sqrt((sum(wi.*(GPS_RIVA(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_RIVA(:,2) = GPS_RIVA-GPS_data(:,3); 
% 
% F_AGE = sqrt((sum(wi.*(GPS_AGE(:,1)-GPS_data(:,3)).^2))/sum(wi));
 GPS_AGE(:,2) = GPS_AGE-GPS_data(:,3); 
% 
% WRMS = [F_GUN, F_RIVA, F_RATES, F_A, F_AGE, F_ICE, F_IJ05, F_W12];
 summ = [GPS_IJ05,GPS_W12,GPS_AGE,GPS_A,GPS_ICE,GPS_RIVA,GPS_GUN,GPS_data(:,7),D];     

%%%%%%%%%%%THIS IS WM
WMgun = (sum(GPS_GUN(:,2).*wi))/sum(wi);
WMR = (sum(GPS_RIVA(:,2).*wi))/sum(wi);
WMRATES = (sum(GPS_RATES(:,2).*wi))/sum(wi);
WMA = (sum(GPS_A(:,2).*wi))/sum(wi);
WMAGE = (sum(GPS_AGE(:,2).*wi))/sum(wi);
WMICE = (sum(GPS_ICE(:,2).*wi))/sum(wi);
WMIJ05 = (sum(GPS_IJ05(:,2).*wi))/sum(wi);
WMW12 = (sum(GPS_W12(:,2).*wi))/sum(wi);
WMS = [WMgun,WMR, WMRATES,WMA,WMAGE,WMICE, WMIJ05,WMW12];

%% Statistics considering spatial correlation- changes weights
 cd ('O:\Documents\RATES\R Code\matlab')
% %variogram values for range -  from GSTAT in R
% %%% l=[262,398,146,587,316,348,1004,361]; 
% % pseudo ranges asumming spatial correlation up to 200 km
l=[250,250,250,250,250,250,250,250]; %range from GPS dhdt

Dij=distmat(GPS_data(:,5:6));
 
 for i=1:length(GPS_data(:,1))
     for j=1:length(l)
       ci(i,j) = sum(exp(-Dij(i,:)/l(j)));
     end
 end
for j=1:8
  %  w2i(:,j) = 1./(ci(:,j).*(GPS_data(:,4).^2 ) );
 w2i(:,j) = 1./(ci(:,j).*(sqrt(GPS_data(:,4).^2+V).^2 ) );
end
%
 WMgun = (sum(GPS_GUN(:,2).*w2i(:,1)))/sum(w2i(:,1));
 WMR = (sum(GPS_RIVA(:,2).*w2i(:,2)))/sum(w2i(:,2));
 WMRATES = (sum(GPS_RATES(:,2).*w2i(:,3)))/sum(w2i(:,3));
 WMA = (sum(GPS_A(:,2).*w2i(:,4)))/sum(w2i(:,4));
 WMAGE = (sum(GPS_AGE(:,2).*w2i(:,5)))/sum(w2i(:,5));
 WMICE = (sum(GPS_ICE(:,2).*w2i(:,6)))/sum(w2i(:,6));
 WMIJ05 = (sum(GPS_IJ05(:,2).*w2i(:,7)))/sum(w2i(:,7));
 WMW12 = (sum(GPS_W12(:,2).*w2i(:,8)))/sum(w2i(:,8));
 WMS_sp1 = [WMgun,WMR, WMRATES,WMA,WMAGE,WMICE, WMIJ05,WMW12];
 
%  %%%%%%%%%WRMS with cluster approach
WRMS_RATES = sqrt((sum(w2i(:,1).*D.^2))/sum(w2i(:,1))); % WRMS gia - gps 
WRMS_W12 = sqrt((sum(w2i(:,2).*(GPS_W12(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,2)));
WRMS_GUN = sqrt((sum(w2i(:,3).*(GPS_GUN(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,3)));
WRMS_ICE = sqrt((sum(w2i(:,4).*(GPS_ICE(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,4)));
WRMS_A = sqrt((sum(w2i(:,5).*(GPS_A(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,5)));
WRMS_IJ05 = sqrt((sum(w2i(:,6).*(GPS_IJ05(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,6)));
WRMS_RIVA = sqrt((sum(w2i(:,7).*(GPS_RIVA(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,7)));
WRMS_AGE = sqrt((sum(w2i(:,8).*(GPS_AGE(:,1)-GPS_data(:,3)).^2))/sum(w2i(:,8)));
WRMS_sp1 = [WRMS_GUN, WRMS_RIVA, WRMS_RATES,WRMS_A, WRMS_AGE,WRMS_ICE, WRMS_IJ05,WRMS_W12];

%% Statistics considering spatial correlation- changes weights using cluster
% cd ('O:\Documents\RATES\R Code\matlab')
% Dij=distmat(GPS_data(:,5:6));
% 
% for i=1:68
% ci2(i) = length(find(Dij(i,:)<250));
% end
% w3i = 1./(ci2'.*(GPS_data(:,4).^2));
% WMgun = (sum(GPS_GUN(:,2).*w3i))/sum(w3i);
% WMR = (sum(GPS_RIVA(:,2).*w3i))/sum(w3i);
% WMRATES = (sum(GPS_RATES(:,2).*w3i))/sum(w3i);
% WMA = (sum(GPS_A(:,2).*w3i))/sum(w3i);
% WMAGE = (sum(GPS_AGE(:,2).*w3i))/sum(w3i);
% WMICE = (sum(GPS_ICE(:,2).*w3i))/sum(w3i);
% WMIJ05 = (sum(GPS_IJ05(:,2).*w3i))/sum(w3i);
% WMW12 = (sum(GPS_W12(:,2).*w3i))/sum(w3i);
% 
% WMS_sp2 = [WMgun,WMR, WMRATES,WMA,WMAGE,WMICE, WMIJ05,WMW12];

%%%%%%%%%WRMS with cluster approach
% WRMS_RATES = sqrt((sum(w3i.*D.^2))/sum(w3i)); % WRMS gia - gps 
% WRMS_W12 = sqrt((sum(w3i.*(GPS_W12(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_GUN = sqrt((sum(w3i.*(GPS_GUN(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_ICE = sqrt((sum(w3i.*(GPS_ICE(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_A = sqrt((sum(w3i.*(GPS_A(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_IJ05 = sqrt((sum(w3i.*(GPS_IJ05(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_RIVA = sqrt((sum(w3i.*(GPS_RIVA(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_AGE = sqrt((sum(w3i.*(GPS_AGE(:,1)-GPS_data(:,3)).^2))/sum(w3i));
% WRMS_sp2 = [WRMS_GUN, WRMS_RIVA, WRMS_RATES,WRMS_A, WRMS_AGE,WRMS_ICE, WRMS_IJ05,WRMS_W12];
% GPS_models=[GPS_data(:,5),GPS_data(:,6),GPS_IJ05,GPS_W12, GPS_AGE, GPS_A, GPS_ICE, GPS_RIVA, GPS_GUN, GPS_RATES];



models = char( 'IJ05R2','W12', 'AGE1b', 'Geruo A', 'ICE6G', 'RIVA','GUNTER', 'RATES');
v=-20:2:20;
grey = [0.4,0.4,0.4];
orange=[1.0,0.4,0.0];
green =[0 .5 0];
yL = get(gca,'YLim');
figure
for i=1:8
    n=i*2;
subplot(2,4,i)
[nall,xall] = hist(summ(:,n),v);
bal = bar(xall,nall);
%set(bal,'facecolor','none','EdgeColor',[1.0,0.4,0.0]);
set(bal,'facecolor',grey,'EdgeColor',grey);
axis([-20 20, 0 35])
hold on
line([median(summ(:,n)),median(summ(:,n))],yL,'Color',orange)
median(summ(:,n))
end

%      
% figure
% for i=1:8
% subplot(2,4,i)
% hist(summ(:,i))
% axis([-20 20, 0 30])
% %title(models(i,:))
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% %set(gca, 'FontSize', fsz);
% set(gca,'TickLength',[0.02,0.02]);
% set(gcf, 'PaperPositionMode', 'auto');
% set(get(gca,'child'),'FaceColor','r','EdgeColor','r');
% end


%%%% Stats
mean(summ(:,16))

% %%%%%Regions defined to obtain local WRMS
% x=GPS_data(:,5);y=GPS_data(:,6);
% 
% AP = x<=-1385 & y >=287; 
% GPS_data_AP = GPS_data(AP,:);
% 
% WAIS = x<-302 & y <=288;
% 
% ROSS = 570 > x & x > 0 & y <-615;
% ROSS2 = 320 > x & x > 270 & y <-1280 & y >-1350 ;
% 
% ASE = x < -1408 & y < -9;
% 
% Sipple = x<-302 & x>-872 & y< -197 & y> -1213;
% 
% RONNE = x>-1385 & x<-600 & y<470 & y>70;
% 
% NAP = [8,17,22,31,47,48,49,54,55,64];
% SAP = [5, 23,25,30,33,62,68];
% ASE = [2, 4, 36, 60 ,61];


% EAIS = zeros(68,1); n=[1,9,14,16,39,59,63];
% EAIS(n)=1;
% figure
% plot(x(Sipple==1),y(Sipple==1),'o')
% hold on
% plot3(xant,yant,zant,'k-');
% 
% cmax = max(abs(summ(:)));cmin = min(abs(summ(:)));
% figure
% for i=1:8
%     n=i*2;
% subplot(2,4,i)
% %h=scatter(GPS_data(:,5),GPS_data(:,6), (abs(summ(:,n))-cmin)+0.00000001/(cmax-cmin)*100, sign(summ(:,n)) );
% h=scatter(GPS_data(:,5),GPS_data(:,6), abs(summ(:,n))*100, sign(summ(:,n)) );
% colormap('cool')
% caxis([-1 1])
% hold on
% plot3(xant,yant,zant,'k-');
% title(models(i,:))
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% set(gca,'TickLength',[0.02,0.02]);
% set(gcf, 'PaperPositionMode', 'auto');
% end
% 
% figure
% for i=1:8
%     n=i*2;
% subplot(2,4,i)
% h=scatter(GPS_data(:,5),GPS_data(:,6), 50, summ(:,n)*10,'fill');
% colormap('jet')
% caxis([-20 20])
% hold on
% plot3(xant,yant,zant,'k-');
% title(models(i,:))
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% %set(gca, 'FontSize', fsz);
% set(gca,'TickLength',[0.02,0.02]);
% set(gcf, 'PaperPositionMode', 'auto');
% end
% 

