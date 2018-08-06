%%
% Created by Sholto Forbes 25/7/18
% This file processes the GPOPS-2 solution output for the
% rocket-scramjet-rocket launch trajectory.
%%

function [] = Plotter(output,auxdata,mode,returnMode,namelist,M_englist,T_englist,engine_data,MList_EngineOn,AOAList_EngineOn,mRocket,mSpartan,mFuel,h0,v0,bounds)

Timestamp = datestr(now,30)
mkdir('../ArchivedResults', strcat(Timestamp, 'mode', num2str(mode)))


copyfile('CombinedProbGPOPS.m',sprintf('../ArchivedResults/%s/SecondStageProb.m',strcat(Timestamp,'mode',num2str(mode))))

save output
movefile('output.mat',sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

% % Run multiple times if output has mutltiple cells
for j = 1:length(output)


auxdata = output{j}.result.setup.auxdata

% =========================================================================
% Assign the primal variables

time1 = output{j}.result.solution.phase(1).time.';
% 
alt1 = output{j}.result.solution.phase(1).state(:,1).';
v1 = output{j}.result.solution.phase(1).state(:,2).';
m1 = output{j}.result.solution.phase(1).state(:,3).';
gamma1 = output{j}.result.solution.phase(1).state(:,4).';
alpha1 = output{j}.result.solution.phase(1).state(:,5).';
zeta1 = output{j}.result.solution.phase(1).state(:,6).';
lat1 = output{j}.result.solution.phase(1).state(:,8).';
lon1 = output{j}.result.solution.phase(1).state(:,9).';
% 

alt21 = output{j}.result.solution.phase(2).state(:,1);
alt22 = output{j}.result.solution.phase(3).state(:,1);
lon21 = output{j}.result.solution.phase(2).state(:,2);
lon22 = output{j}.result.solution.phase(3).state(:,2);
lat21 = output{j}.result.solution.phase(2).state(:,3);
lat22 = output{j}.result.solution.phase(3).state(:,3);
v21 = output{j}.result.solution.phase(2).state(:,4); 
v22 = output{j}.result.solution.phase(3).state(:,4); 
gamma21 = output{j}.result.solution.phase(2).state(:,5); 
gamma22 = output{j}.result.solution.phase(3).state(:,5); 
zeta21 = output{j}.result.solution.phase(2).state(:,6);
zeta22 = output{j}.result.solution.phase(3).state(:,6);
alpha21 = output{j}.result.solution.phase(2).state(:,7);
alpha22 = output{j}.result.solution.phase(3).state(:,7);
eta21 = output{j}.result.solution.phase(2).state(:,8);
eta22 = output{j}.result.solution.phase(3).state(:,8);
mFuel21 = output{j}.result.solution.phase(2).state(:,9); 
mFuel22 = output{j}.result.solution.phase(3).state(:,9); 

throttle22 = output{j}.result.solution.phase(3).state(:,10);

aoadot21  = output{j}.result.solution.phase(2).control(:,1); 
etadot21  = output{j}.result.solution.phase(2).control(:,2); 

aoadot22  = output{j}.result.solution.phase(3).control(:,1); 
etadot22  = output{j}.result.solution.phase(3).control(:,2); 

time21 = output{j}.result.solution.phase(2).time;
time22 = output{j}.result.solution.phase(3).time;


alt3  = output{j}.result.solution.phase(4).state(:,1);
v3    = output{j}.result.solution.phase(4).state(:,2);
gamma3  = output{j}.result.solution.phase(4).state(:,3);
m3    = output{j}.result.solution.phase(4).state(:,4);
aoa3    = output{j}.result.solution.phase(4).state(:,5);
lat3    = output{j}.result.solution.phase(4).state(:,6);
zeta3    = output{j}.result.solution.phase(4).state(:,7);
aoadot3       = output{j}.result.solution.phase(4).control(:,1);

time3 = output{j}.result.solution.phase(4).time;


figure(01)
subplot(9,1,1)
hold on
plot(time21,alt21)
plot(time22,alt22)
subplot(9,1,2)
hold on
plot(time21,v21)
plot(time22,v22)
subplot(9,1,3)
hold on
plot(time21,lon21)
plot(time22,lon22)
subplot(9,1,4)
hold on
plot(time21,lat21)
plot(time22,lat22)
subplot(9,1,5)
hold on
plot(time21,v21)
plot(time22,v22)
subplot(9,1,6)
hold on
plot(time21,gamma21)
plot(time22,gamma22)
subplot(9,1,7)
hold on
plot(time21,ones(1,length(time21)))
plot(time22,throttle22)


figure(230)
hold on
plot3(lon21,lat21,alt21)
plot3(lon22,lat22,alt22)

  
%     figure(2301)
% hold on
% 
% axesm('pcarree','Origin',[0 rad2deg(lon21(1)) 0])
% geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
% % plotm(rad2deg(lat),rad2deg(lon+lon0))
% plotm(rad2deg(lat21),rad2deg(lon21),'b')
% plotm(rad2deg(lat22),rad2deg(lon22),'r')
%     
%     cities = shaperead('worldcities', 'UseGeoCoords', true);
% lats = extractfield(cities,'Lat');
% lons = extractfield(cities,'Lon');
% geoshow(lats, lons,...
%         'DisplayType', 'point',...
%         'Marker', 'o',...
%         'MarkerEdgeColor', 'r',...
%         'MarkerFaceColor', 'r',...
%         'MarkerSize', 2)

% =========================================================================

%% Third Stage
% Optimise third stage trajectory from end point

ThirdStagePayloadMass = -output{j}.result.objective;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = length(alt21)

[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21, M21, Fd21, rho21,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21] = VehicleModelCombined(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,alpha21,eta21,1, mFuel21,mFuel21(1),mFuel21(end), 1, 0);
[~,~,~,~,~,~, q22, M22, Fd22, rho22,L22,Fueldt22,T22,Isp22,q122,flapdeflection22,heating_rate22] = VehicleModelCombined(gamma22, alt22, v22,auxdata,zeta22,lat22,lon22,alpha22,eta22,throttle22, mFuel22,0,0, 0, 0);

%% Modify throttle for return 


throttle22(M22<5.0) = throttle22(M22<5.0).*gaussmf(M22(M22<5.0),[.01,5]); % remove nonsense throttle points
throttle22(q122<20000) = throttle22(q122<20000).*gaussmf(q122(q122<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
    
Isp22(M22<5.0) = Isp22(M22<5.0).*gaussmf(M22(M22<5.0),[.01,5]); %
Isp22(q122<20000) = Isp22(q122<20000).*gaussmf(M22(q122<20000),[.1,5]); %


% Separation_LD = lift(end)/Fd(end)

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 


forward0 = [alt21(1),gamma21(1),v21(1),zeta21(1),lat21(1),lon21(1), mFuel21(1)];

[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelAscent_forward(f_t, f_y,auxdata,ControlInterp(time21,alpha21,f_t),ControlInterp(time21,eta21,f_t),1,mFuel21(1),mFuel21(end)),time21(1:end),forward0);

forward_error21 = [(alt21-f_y(:,1))/(max(alt21)-min(alt21)) (lon21-f_y(:,6))/(max(lon21)-min(lon21))...
    (lat21-f_y(:,5))/(max(lat21)-min(lat21)) (v21-f_y(:,3))/(max(v21)-min(v21))...
    (gamma21-f_y(:,2))/(max(gamma21)-min(gamma21)) (zeta21-f_y(:,4))/(max(zeta21)-min(zeta21)) (mFuel21-f_y(:,7))/(max(mFuel21)-min(mFuel21))];



figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time21,alt21);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time21,gamma21);


subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time21,v21);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time21,zeta21);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time21,mFuel21);



% Return Forward
forward0 = [alt22(1),gamma22(1),v22(1),zeta22(1),lat22(1),lon22(1), mFuel22(1)];
% 
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(1):time22(end),forward0);
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22,forward0);


[f_t221, f_y221] = ode15s(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(1:ceil(length(time22)/6)-1),forward0);

forward0 = [alt22(ceil(length(time22)/6)),gamma22(ceil(length(time22)/6)),v22(ceil(length(time22)/6)),zeta22(ceil(length(time22)/6)),lat22(ceil(length(time22)/6)),lon22(ceil(length(time22)/6)), mFuel22(ceil(length(time22)/6))];

[f_t222, f_y222] = ode15s(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(ceil(length(time22)/6):ceil(length(time22)/3)-1),forward0);


forward0 = [alt22(ceil(length(time22)/3)),gamma22(ceil(length(time22)/3)),v22(ceil(length(time22)/3)),zeta22(ceil(length(time22)/3)),lat22(ceil(length(time22)/3)),lon22(ceil(length(time22)/3)), mFuel22(ceil(length(time22)/3))];

[f_t223, f_y223] = ode15s(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time22,alpha22,f_t),ControlInterp(time22,eta22,f_t),ThrottleInterp(time22,throttle22,f_t)),time22(ceil(length(time22)/3):end),forward0);


f_t2 = [f_t221 ; f_t222 ; f_t223];
% 
f_y2 = [f_y221; f_y222 ; f_y223];

forward_error22 = [(alt22-f_y2(:,1))/(max(alt22)-min(alt22)) (lon22-f_y2(:,6))/(max(lon22)-min(lon22))...
    (lat22-f_y2(:,5))/(max(lat22)-min(lat22)) (v22-f_y2(:,3))/(max(v22)-min(v22))...
    (gamma22-f_y2(:,2))/(max(gamma22)-min(gamma22)) (zeta22-f_y2(:,4))/(max(zeta22)-min(zeta22)) (mFuel22-f_y2(:,7))/(max(mFuel22)-min(mFuel22))];


figure(213)
subplot(7,1,1)
hold on
plot(f_t2(1:end),f_y2(:,1));
plot(time22,alt22);

% gamma  = output.result.solution.phase.state(:,5);

subplot(7,1,2)
hold on
plot(f_t2(1:end),f_y2(:,2));
plot(time22,gamma22);

% latitude  = output.result.solution.phase.state(:,3);
subplot(7,1,3:5)
hold on
plot(f_y2(:,6),f_y2(:,5));
plot(lon22,lat22);

subplot(7,1,6)
hold on
plot(f_t2(1:end),f_y2(:,7));
plot(time22,mFuel22);



%% plot Ascent
addpath('addaxis')
addpath('axlabel')

 figure(211)
 fig = gcf;
set(fig,'Position',[200 0 850 1200])

% %suptitle('Second Stage Ascent')
subplot(4,1,1)
% set(gca,'xticklabels',[])
hold on
xlim([0 time21(end)-time21(1)]);
 plot(time21-time21(1),alt21/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
addaxis(time21-time21(1),M21,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Mach no.');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');


addaxis(time21-time21(1),q21/1000,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'DYnamic Pressure (kPa)');

legend(  'Altitude', 'Mach no.', 'Dynamic Pressure', 'location', 'best');


subplot(4,1,2)
% set(gca,'xticklabels',[])
hold on
xlim([0 time21(end)-time21(1)]);
 plot(time21-time21(1),rad2deg(gamma21),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Trajectory Angle (deg)');

addaxis(time21-time21(1),rad2deg(zeta21),'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Heading Angle (deg)');

addaxis(time21-time21(1),L21./Fd21,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'L/D');

legend(  'Trajectory Angle', 'Heading Angle', 'L/D', 'location', 'best');


subplot(4,1,3)
% set(gca,'xticklabels',[])
hold on
xlim([0 time21(end)-time21(1)]);
 plot(time21-time21(1),rad2deg(alpha21),'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Angle of Attack (deg)');

addaxis(time21-time21(1),rad2deg(eta21),'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Bank Angle (deg)');

addaxis(time21-time21(1),flapdeflection21,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Flap Deflection (deg)');

legend(  'Angle of Attack', 'Bank Angle', 'Flap Deflection', 'location', 'best');

subplot(4,1,4)
Isp21 = T21./Fueldt21./9.81;
IspNet21 = (T21-Fd21)./Fueldt21./9.81;
xlim([0 time21(end)-time21(1)]);
hold on

 plot(time21-time21(1),IspNet21,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('Net Isp (s)');
xlabel('Time (s)');

addaxis(time21-time21(1),T21/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Thrust (kN)');

addaxis(time21-time21(1),mFuel21,':','color','k', 'linewidth', 1.);
addaxislabel(2,'Fuel Mass (kg)');

legend(  'Net Isp', 'Thrust', 'Fuel Mass', 'location', 'best');


% figure(211)
% %suptitle('Second Stage Ascent')
% fig = gcf;
% set(fig,'Position',[200 0 850 1200])
% 
% title('Second Stage Ascent Trajectory');
% hold on
% 
% subplot(7,2,1)
% plot(time21-time21(1), alt21/1000,'Color','k')
% hold on
% title('Trajectory (km)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% dim = [.62 .7 .2 .2];
% annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Fuel Used: ' num2str(1562 - mFuel21(end)) ' kg']},'FitBoxToText','on');  
% 
% % subplot(7,2,2)
% % plot(time21-time21(1), q121/1000,'Color','k')
% % hold on
% % title('Dynamic Pressure Entering Engine (kPa)','FontSize',9)
% % set(gca,'xticklabels',[])
% % xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,3)
% plot(time21-time21(1), v21,'Color','k')
% hold on
% title('Velocity (m/s)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,4)
% plot(time21-time21(1), M21,'Color','k')
% hold on
% title('Mach no','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,5)
% plot(time21-time21(1), q21/1000,'Color','k')
% hold on
% title('Dynamic Pressure (kpa)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,6)
% plot(time21-time21(1), rad2deg(gamma21),'Color','k')
% hold on
% title('Trajectory Angle (Deg)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,7)
% plot(time21-time21(1), rad2deg(alpha21),'Color','k')
% hold on
% title('Angle of Attack (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% if returnMode == 1
% subplot(7,2,8)
% plot(time21-time21(1), rad2deg(eta21),'Color','k')
% hold on
% title('Bank Angle (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% else
%   subplot(7,2,8)
% plot(time21-time21(1), CG21,'Color','k')
% hold on
% title('Centre of Gravity (m)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);  
% end
% 
% subplot(7,2,9)
% plot(time21-time21(1), flapdeflection21,'Color','k')
% hold on
% title('Flap Deflection (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% % Isp1 = T1./Fueldt1./9.81;
% IspNet21 = (T21-Fd21)./Fueldt21./9.81;
% 
% 
% % Determine the thrust from the boat tail, and the last part of the nozzle.
% % Note: this is included in Fd in all other places for computational
% % efficiency
% T21_rear = auxdata.T_spline_Rear(M21,rad2deg(alpha21),alt21/1000);
% Cd_noengine = auxdata.Fd_spline_NoEngine(M21,rad2deg(alpha21),alt21/1000);
% Fd_noengine = 0.5*Cd_noengine.*auxdata.A.*rho21.*v21.^2*auxdata.Cdmod;
% 
% subplot(7,2,13)
% plot(time21-time21(1), T21/1000+T21_rear/1000 ,'Color','k')
% hold on
% title('Thrust (kN)','FontSize',9)
% % set(gca,'xticklabels',[])
% xlabel('Time (s)','FontSize',9);
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,12)
% plot(time21-time21(1), IspNet21,'Color','k')
% hold on
% title('Net Isp (s)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,11)
% plot(time21-time21(1), mFuel21,'Color','k')
% hold on
% title('Fuel Mass (kg)','FontSize',9)
% % xlabel('Time (s)','FontSize',9);
% xlim([0 time21(end)-time21(1)]);
% set(gca,'xticklabels',[])
% 
% subplot(7,2,10)
% plot(time21-time21(1), rad2deg(zeta21),'Color','k')
% hold on
% title('Heading Angle (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% xlim([0 time21(end)-time21(1)]);
% 
% subplot(7,2,14)
% plot(time21-time21(1), L21./Fd_noengine,'Color','k')
% hold on
% title('L/D','FontSize',9)
% xlabel('Time (s)','FontSize',9);
% xlim([0 time21(end)-time21(1)]);
% 


ThirdStageFlag = [ones(length(time21),1); zeros(length(time22),1)];

SecondStageStates = [[time21; time22] [alt21; alt22] [lon21; lon22] [lat21; lat22] [v21; v22] [gamma21; gamma22] [zeta21; zeta22] [alpha21; alpha22] [eta21; eta22] [mFuel21; mFuel22] ThirdStageFlag];
dlmwrite(strcat('SecondStageStates',namelist{j}),['time (s) ' 'altitude (m) ' 'longitude (rad) ' 'latitude (rad) ' 'velocity (m/s) ' 'trajectory angle (rad) ' 'heading angle (rad) ' 'angle of attack (rad) ' 'bank angle (rad) ' 'fuel mass (kg) ' 'Third Stage (flag)'],'');
dlmwrite(strcat('SecondStageStates',namelist{j}),SecondStageStates,'-append','delimiter',' ');
movefile(strcat('SecondStageStates',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

%% Plot Return
% figure(221)
% %suptitle('Second Stage Return')
% fig = gcf;
% set(fig,'Position',[200 0 850 1200])
% 
% title('Second Stage Return Trajectory');
% hold on
% 
% subplot(7,2,1)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), alt22/1000,'Color','k')
% title('Trajectory (km)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% dim = [.62 .7 .2 .2];
% annotation('textbox',dim,'string',{['Fuel Used: ' num2str(mFuel22(1)) ' kg']},'FitBoxToText','on');  
% 
% subplot(7,2,3)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), v22,'Color','k')
% title('Velocity (m/s)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,4)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), M22,'Color','k')
% title('Mach no','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,5)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), q22/1000,'Color','k')
% title('Dynamic Pressure (kpa)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,6)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), rad2deg(gamma22),'Color','k')
% title('Trajectory Angle (Deg)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,7)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), rad2deg(alpha22),'Color','k')
% title('Angle of Attack (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,8)
% xlim([0 time22(end)-time22(1)]);
% hold on
% ylim([rad2deg(min(eta22))-1 rad2deg(max(eta22))+1]);
% plot(time22-time22(1), rad2deg(eta22),'Color','k')
% title('Bank Angle (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,9)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), flapdeflection22,'Color','k')
% title('Flap Deflection (deg)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% T22_rear = auxdata.T_spline_Rear(M22,rad2deg(alpha22),alt22/1000);
% 
% subplot(7,2,13)
% xlim([0 time22(end)-time22(1)]);
% hold on
% ylim([0 max(T22+T22_rear)/1000]);
% plot(time22-time22(1), T22/1000+T22_rear/1000,'Color','k')
% title('Thrust (kN)','FontSize',9)
% 
% xlabel('Time (s)','FontSize',9);
% 
% subplot(7,2,12)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), Isp22,'Color','k')
% title('Potential Isp (s)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,11)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), mFuel22,'Color','k')
% title('Fuel Mass (kg)','FontSize',9)
% set(gca,'xticklabels',[])
% 
% subplot(7,2,14)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), throttle22,'Color','k')
% title('Throttle','FontSize',9)
% xlabel('Time (s)','FontSize',9);
% 
% subplot(7,2,10)
% xlim([0 time22(end)-time22(1)]);
% hold on
% plot(time22-time22(1), rad2deg(zeta22),'Color','k')
% title('Heading Angle (deg)','FontSize',9)
% set(gca,'xticklabels',[])


%% Check KKT and pontryagins minimum
% Check that the hamiltonian = 0 (for free end time)
% Necessary condition
input_test = output{j}.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);

H1 = [];
H2 = [];
H3 = [];
H4 = [];

lambda1 = output{j}.result.solution.phase(1).costate;
for i = 1:length(lambda1)-1
    H1(i) = lambda1(i+1,:)*phaseout_test(1).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda2 = output{j}.result.solution.phase(2).costate;
for i = 1:length(lambda2)-1
    H2(i) = lambda2(i+1,:)*phaseout_test(2).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda3 = output{j}.result.solution.phase(3).costate;
for i = 1:length(lambda3)-1
    H3(i) = lambda3(i+1,:)*phaseout_test(3).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda4 = output{j}.result.solution.phase(4).costate;
for i = 1:length(lambda4)-1
    H4(i) = lambda4(i+1,:)*phaseout_test(4).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

figure(2410)
hold on
plot(time1(1:end-1),H1)
plot(time21(1:end-1),H2)
plot(time22(1:end-1),H3)
plot(time3(1:end-1),H4)
ylabel('Hamiltonian')
xlabel('Time (s)')
legend('First Stage Stage Ascent','Second Stage Ascent','Second Stage Return','Third Stage Powered Ascent')

%% Check State Feasibility
% Check calculated derivatives with the numerical derivative of each
% porimal, scaled by that primal

% for i = 1:length(output{j}.result.solution.phase(2).state(1,:))
% plot(time21(1:end-1),([diff(output{j}.result.solution.phase(2).state(:,i))./diff(output{j}.result.solution.phase(2).time)] - phaseout_test(2).dynamics(1:end-1,i))./max(abs(phaseout_test(2).dynamics(1:end-1,i))),'--');
% end
% for i = 1:length(output{j}.result.solution.phase(3).state(1,:))
%     if i<= 7 % Plot different line styles when no. of colours exceeded
%     plot(time22(1:end-1),([diff(output{j}.result.solution.phase(3).state(:,i))./diff(output{j}.result.solution.phase(3).time)] - phaseout_test(3).dynamics(1:end-1,i))./max(abs(phaseout_test(3).dynamics(1:end-1,i))));
%     else
%     plot(time22(1:end-1),([diff(output{j}.result.solution.phase(3).state(:,i))./diff(output{j}.result.solution.phase(3).time)] - phaseout_test(3).dynamics(1:end-1,i))./max(abs(phaseout_test(3).dynamics(1:end-1,i))),':');
%     end
% end

% 
% for i = 1:length(output{j}.result.solution.phase(2).state(1,:))
% plot(time21(1:end),(cumtrapz(output{j}.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,i)) + output{j}.result.solution.phase(2).state(1,i)- output{j}.result.solution.phase(2).state(:,i))./(max(output{j}.result.solution.phase(2).state(:,i))),'--');
% end
% for i = 1:length(output{j}.result.solution.phase(3).state(1,:))
%     if i<= 7 % Plot different line styles when no. of colours exceeded
%     plot(time22(1:end),(cumtrapz(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i)) + output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))),'-');
%     else
%     plot(time22(1:end),(cumtrapz(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i)) + output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))),':');
%     end
% end

figure(2420)
fig = gcf;
set(fig,'Position',[200 200 700 600]);
hold on

%suptitle('Dynamics Feasibility Check')

for i = [1:6 9]
    if i == 1
        i1 = 1;
    elseif i ==2
        i1 = 9;
   elseif i ==3
        i1 = 8;
   elseif i ==4
        i1 = 2;
   elseif i ==5
        i1 = 4;
   elseif i ==6
        i1 = 6;
   elseif i ==9
        i1 = 3;
    end
   
        
    Stage1_int = cumtrapz(output{j}.result.solution.phase(1).time,phaseout_test(1).dynamics(1:end,i1));
    Stage1_error = (Stage1_int+ output{j}.result.solution.phase(1).state(1,i1)- output{j}.result.solution.phase(1).state(:,i1))./(max(output{j}.result.solution.phase(1).state(:,i1))-min(output{j}.result.solution.phase(1).state(:,i1)));
%     Stage1_int = cumtrapz([time1(1):0.1:time1(end)],pchip(output{j}.result.solution.phase(1).time,phaseout_test(1).dynamics(1:end,i),[time1(1):0.1:time1(end)]))';
%     Stage1_error = (interp1([time1(1):0.1:time1(end)],Stage1_int,time1)'+ output{j}.result.solution.phase(1).state(1,i)- output{j}.result.solution.phase(1).state(:,i))./(max(output{j}.result.solution.phase(1).state(:,i))-min(output{j}.result.solution.phase(1).state(:,i)));

%     Stage2_int = cumtrapz(output{j}.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,i));
%    Stage2_error = (Stage2_int+ output{j}.result.solution.phase(2).state(1,i)- output{j}.result.solution.phase(2).state(:,i))./(max(output{j}.result.solution.phase(2).state(:,i))-min(output{j}.result.solution.phase(2).state(:,i)));
Stage2_int = cumtrapz([time21(1):0.1:time21(end)],pchip(output{j}.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,i),[time21(1):0.1:time21(end)]))';
    Stage2_error = (interp1([time21(1):0.1:time21(end)],Stage2_int,time21)+ output{j}.result.solution.phase(2).state(1,i)- output{j}.result.solution.phase(2).state(:,i))./(max(output{j}.result.solution.phase(2).state(:,i))-min(output{j}.result.solution.phase(2).state(:,i)));
%     Return_int = cumtrapz(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i));
%     Return_error = (Return_int+ output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))-min(output{j}.result.solution.phase(3).state(:,i)));
  Return_int = cumtrapz([time22(1):0.1:time22(end)],pchip(output{j}.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,i),[time22(1):0.1:time22(end)]))';
    Return_error = (interp1([time22(1):0.1:time22(end)],Return_int,time22)+ output{j}.result.solution.phase(3).state(1,i)- output{j}.result.solution.phase(3).state(:,i))./(max(output{j}.result.solution.phase(3).state(:,i))-min(output{j}.result.solution.phase(3).state(:,i)));

    
    
    if i == 1
        i3 = 1;
    elseif i ==2
        i3 = 0;
   elseif i ==3
        i3 = 6;
   elseif i ==4
        i3 = 2;
   elseif i ==5
        i3 = 3;
   elseif i ==6
        i3 = 7;
   elseif i ==9
        i3 = 4;
    end
    
    if i3 == 0
        Stage3_error = zeros(length(time3),1);
    else
    Stage3_int = cumtrapz(output{j}.result.solution.phase(4).time,phaseout_test(4).dynamics(1:end,i3));
    Stage3_error = (Stage3_int+ output{j}.result.solution.phase(4).state(1,i3)- output{j}.result.solution.phase(4).state(:,i3))./(max(output{j}.result.solution.phase(4).state(:,i3))-min(output{j}.result.solution.phase(4).state(:,i3)));
    
%     Stage3_int = cumtrapz([time3(1):0.1:time3(end)],pchip(output{j}.result.solution.phase(4).time,phaseout_test(4).dynamics(1:end,i),[time3(1):0.1:time3(end)]))';
%     Stage3_error = (interp1([time3(1):0.1:time3(end)],Stage3_int,time3)+ output{j}.result.solution.phase(4).state(1,i)- output{j}.result.solution.phase(4).state(:,i))./(max(output{j}.result.solution.phase(4).state(:,i))-min(output{j}.result.solution.phase(4).state(:,i)));

    end

    
% plot([time1'; time21 ;time22; time3],[Stage1_error;Stage2_error; Return_error;Stage3_error]);


subplot(2,2,1)
hold on
title('First Stage')
plot(time1',Stage1_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time1(end)-time1(1)]);

subplot(2,2,2)
hold on
title('Second Stage Ascent')
plot(time21-time21(1),Stage2_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time21(end)-time21(1)]);
subplot(2,2,3)
hold on
title('Second Stage Return')
plot(time22-time22(1),Return_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time22(end)-time22(1)]);
subplot(2,2,4)
hold on
title('Third Stage')
plot(time3-time3(1),Stage3_error*100);
ylabel('Integrated Error (%)')
xlabel('Time (s)')
ylim([-1 1]);
xlim([0 time3(end)-time3(1)]);
end




% ylim([-1,1])
% legend('Alt Ascent','lon Ascent','lat Ascent','v Ascent','gamma Ascent','zeta Ascent','aoa Ascent','bank Ascent','mFuel Ascent', 'Alt Descent','lon Descent','lat Descent','v Descent','gamma Descent','zeta Descent','aoa Descent','bank Descent','mFuel Descent','throttle Descent')
legend('Alt ','lon ','lat','v ','gamma ','zeta ','Mass');


%% Plot Heating Rate

figure(701)
plot(time21,heating_rate21/1000);
hold on
plot(time22,heating_rate22/1000);
xlabel('time');
ylabel('Heat Flux (MW/m^2)');
title('Stagenation Point Heating')


%% plot engine interpolation visualiser
T0 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt21); 
T_in1 = auxdata.interp.tempgridded(M21,rad2deg(alpha21)).*T0;
M_in1 = auxdata.interp.M1gridded(M21, rad2deg(alpha21));

plotM = [min(M_englist):0.01:9];
plotT = [min(T_englist):1:550];
[gridM,gridT] =  ndgrid(plotM,plotT);
interpeq = auxdata.interp.eqGridded(gridM,gridT);
interpIsp = auxdata.interp.IspGridded(gridM,gridT);

figure(2100)
hold on
contourf(gridM,gridT,interpeq,100,'LineColor','none');
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,4),'k');
xlabel('M1')
ylabel('T1')
c=colorbar
c.Label.String = 'Equivalence Ratio';
caxis([.4 1])
plot(M_in1,T_in1,'r');

error_Isp = auxdata.interp.IspGridded(engine_data(:,1),engine_data(:,2))-engine_data(:,3);

figure(2110)
hold on
contourf(gridM,gridT,interpIsp,100,'LineColor','none');
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'k')
xlabel('M1')
ylabel('T1')
c=colorbar
c.Label.String = 'ISP';
plot(M_in1,T_in1,'r');

figure(2120)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.M1gridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'M1';

figure(2130)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.tempgridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'T1/T0';

figure(2140)
contourf(MList_EngineOn,AOAList_EngineOn,auxdata.interp.presgridded(MList_EngineOn,AOAList_EngineOn),100,'LineWidth',0)
xlabel('M')
ylabel('Angle of Attack (deg)')
c=colorbar
c.Label.String = 'P1/P0';

%%
[gridM2,gridAoA2] =  ndgrid(plotM,plotT);


%% ThirdStage


forward0 = [alt3(1),v3(1),gamma3(1),m3(1),lat3(1),zeta3(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,mode,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel3_forward(f_t, f_y,auxdata,ControlInterp(time3,aoa3,f_t),ControlInterp(time3,aoadot3,f_t)),time3(1:end),forward0);

[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDyn(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,lat3,zeta3);

lon3 = [];
lon3(1) = lon21(end);
for i = 2:length(time3)
    lon3(i) = lon3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
end


[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end),lon3(end), zeta3(end), m3(end), auxdata);


forward_error3 = [(alt3-f_y(:,1))/(max(alt3)-min(alt3)) zeros(length(alt3),0)...
    (lat3-f_y(:,5))/(max(lat3)-min(lat3)) (v3-f_y(:,2))/(max(v3)-min(v3))...
    (gamma3-f_y(:,3))/(max(gamma3)-min(gamma3)) (zeta3-f_y(:,6))/(max(zeta3)-min(zeta3)) (m3-f_y(:,4))/(max(m3)-min(m3))];



figure(312)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time3,alt3);

figure(313)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time3,v3);



figure(311)
    fig = gcf;
set(fig,'Position',[200 0 850 650])
%%suptitle('Third Stage Trajectory');
    hold on
    
    subplot(4,2,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time3; timeexo.'+time3(end)], [alt3; altexo.']/1000,'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
set(gca,'xticklabels',[])
    subplot(4,2,2)
    hold on
    title('Dynamic Pressure (kPa','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[q3;qexo.';qexo(end)]/1000,'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
set(gca,'xticklabels',[])
    subplot(4,2,3)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
set(gca,'xticklabels',[])
    subplot(4,2,4)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[v3;v3exo.'],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
set(gca,'xticklabels',[])
    subplot(4,2,5)
    hold on
    title('Mass (kg)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[ m3;mexo.';mexo(end)],'Color','k')
    xlim([time3(1) timeexo(end)+time3(end)])
set(gca,'xticklabels',[])
    subplot(4,2,6)
    hold on
    title('Thrust Vector Angle (deg)','FontSize',9);
    plot([time3; timeexo.'+time3(end)],[rad2deg(Vec_angle3);0*ones(length(timeexo),1)],'Color','k')
    xlabel('Time (s)','FontSize',9);
    xlim([time3(1) timeexo(end)+time3(end)])
set(gca,'xticklabels',[])
    subplot(4,2,7)
    hold on
    title('Thrust','FontSize',9);
    plot([time3], [T3],'Color','k')

    xlabel('Time (s)','FontSize',9);
    xlim([time3(1) timeexo(end)+time3(end)])
    subplot(4,2,8)
    hold on
    title('Drag','FontSize',9);
    plot([time3], [D3],'Color','k')

    xlabel('Time (s)','FontSize',9);
    xlim([time3(1) timeexo(end)+time3(end)])
    
    % Write data to file
    dlmwrite(strcat('ThirdStageStates',namelist{j}),['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg) ' 'dynamic pressure (Pa)' 'trajectory angle (rad) ' 'Lift (N)' 'Drag (N)' 'heading angle (rad) ' 'latitude (rad) ' 'angle of attack (rad) '],'');
    dlmwrite(strcat('ThirdStageStates',namelist{j}),[[time3; time3(end)+timeexo'], [alt3; altexo'], [v3; v3exo'], [m3; mexo'; mexo(end)],[q3; qexo'; qexo(end)] ,[gamma3; gammaexo'],[L3; Lexo'; Lexo(end)],[D3; Dexo'; Dexo(end)] ,[zeta3; zetaexo'], [lat3; latexo'], [aoa3; zeros(length(timeexo),1)]],'-append','delimiter',' ')
movefile(strcat('ThirdStageStates',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

% 
%% First Stage =========================================================

% 
FirstStageSMF = (mRocket - mFuel)/(m1(1) - mSpartan);
% 

FirstStageStates = [time1' alt1' v1' m1' gamma1' alpha1' zeta1' lat1' lon1'];

dlmwrite(strcat('FirstStageStates',namelist{j}),['time (s) ' 'altitude (m) ' 'velocity (m/s) ' 'mass (kg)' 'trajectory angle (rad) ' 'angle of attack (rad) ' 'heading angle (rad) ' 'latitude (rad)'],'');
dlmwrite(strcat('FirstStageStates',namelist{j}),FirstStageStates,'-append','delimiter',' ');
movefile(strcat('FirstStageStates',namelist{j}),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));


% Iterative Prepitch Determination ========================================
%This back determines the mass and launch altitude necessary to get to
%100m, 30m/s at the PS method determined fuel mass

interp = auxdata.interp;
Throttle = auxdata.Throttle;
Vehicle = auxdata.Vehicle;
Atmosphere = auxdata.Atmosphere;

% ntoe that launch altitude does vary, but it should only be slightly
controls = fminunc(@(controls) prepitch(controls,m1(1),interp,Throttle,Vehicle,Atmosphere),[-50,10]);

h_launch = controls(1)
t_prepitch = controls(2)
Isp1 = Vehicle.Isp.SL;
T1 = Vehicle.T.SL*Throttle;
dm1 = -T1./Isp1./9.81;
m0_prepitch = m1(1) - dm1*t_prepitch;



%% Forward Integrator
 phase = 'postpitch';
tspan = time1; 
% postpitch0_f = [y(end,1) y(end,2) y(end,3) deg2rad(89.9) phi(1) zeta(1)]; % set mass
postpitch0_f = [h0 v0 m1(1) deg2rad(89.9) lat1(1) zeta1(1)];

[t_postpitch_f, postpitch_f] = ode45(@(t_f,postpitch_f) rocketDynamicsForward(postpitch_f,ControlFunction(t_f,time1,zeta1),ControlFunction(t_f,time1,alpha1),phase,interp,Throttle,Vehicle,Atmosphere), tspan, postpitch0_f);


forward_error1 = [(alt1'-postpitch_f(:,1))/(max(alt1)-min(alt1)) zeros(length(alt1),1)...
    (lat1'-postpitch_f(:,5))/(max(lat1)-min(lat1)) (v1'-postpitch_f(:,2))/(max(v1)-min(v1))...
    (gamma1'-postpitch_f(:,4))/(max(gamma1)-min(gamma1)) (zeta1'-postpitch_f(:,6))/(max(zeta1)-min(zeta1)) (m1'-postpitch_f(:,3))/(max(m1)-min(m1))];


figure(103)
hold on
plot(postpitch_f(:,1));
plot(alt1);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = h_launch;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
gamma0_prepitch = deg2rad(90);

phase = 'prepitch';
tspan2 = [0 t_prepitch]; % time to fly before pitchover (ie. straight up)

y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, 0, 0, 0, 0];

% this performs a forward simulation before pitchover. The end results of
% this are used as initial conditions for the optimiser. 
[t_prepitch, y] = ode45(@(t,y) rocketDynamics(y,0,0,phase,interp,Throttle,Vehicle,Atmosphere), tspan2, y0);  


figure(111);
hold on
%%suptitle('First Stage Trajectory');
    fig = gcf;
set(fig,'Position',[200 0 850 600])
subplot(4,2,1)
hold on
title('Trajectory Angle (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(y(:,4).') rad2deg(gamma1)],'color','k');
subplot(4,2,2)
hold on
title('Velocity (m/s)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [y(:,2).' v1],'color','k');
subplot(4,2,3)
hold on
title('Altitude (km)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [y(:,1).'/1000 alt1/1000],'color','k');
subplot(4,2,4)
hold on
title('Angle of Attack (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [zeros(1,length(t_prepitch)) rad2deg(alpha1)],'color','k');
subplot(4,2,5)
hold on
title('Mass (kg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [y(:,3).' m1],'color','k');
subplot(4,2,6)
hold on
title('Heading Angle (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([time1+t_prepitch(end)], [rad2deg(zeta1)],'color','k');
subplot(4,2,7)
hold on
title('Latitude (deg)');
xlim([0 time1(end)+t_prepitch(end)]);
plot([t_prepitch.' time1+t_prepitch(end)], [rad2deg(lat1(1)+y(:,8).') rad2deg(lat1)],'color','k');

% plot([primal.nodes], [rad2deg(gamma)/100],'color','k','linestyle','-');
% plot([primal.nodes], [v/1000],'color','k','linestyle','--');
% plot([primal.nodes], [V/10000],'color','k','linestyle',':');
% plot([primal.nodes], [rad2deg(alpha)/10],'color','k','linestyle','-.')
xlabel('Time (s)')
xlim([0,time1(end)+t_prepitch(end)]);

saveas(figure(111),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('FirstStage',namelist{j},'.fig')]);
print(figure(111),strcat('FirstStage',namelist{j}),'-dpng');
movefile(strcat('FirstStage',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

%% Plot Forward Sim Error
figure(550)
fig = gcf;
set(fig,'Position',[200 200 700 600]);
%suptitle('Normalised Forward Simulation Error')
subplot(2,2,1)
hold on
for i = 1:length(forward_error1(1,:))
    plot(time1-time1(1),forward_error1(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time1(end)-time1(1)]);
title('First Stage')

subplot(2,2,2)
hold on
for i = 1:length(forward_error21(1,:))
    plot(time21-time21(1),forward_error21(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time21(end)-time21(1)]);
title('Second Stage Ascent')

subplot(2,2,3)
hold on
for i = 1:length(forward_error22(1,:))
    plot(time22-time22(1),forward_error22(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time22(end)-time22(1)]);
title('Second Stage Return')
legend('Alt ','lon ','lat','v ','gamma ','zeta ','Mass','Location','best');

subplot(2,2,4)
hold on
for i = 1:length(forward_error3(1,:))
    plot(time3-time3(1),forward_error3(:,i)*100)
end
ylabel('Integrated Error (%)')
xlabel('Time (s)')
xlim([0 time3(end)-time3(1)]);
title('Third Stage Powered Ascent')

% hSub = subplot(2,3,3); plot(1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');





%% Create Easy Latex Inputs

dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\PayloadToOrbitMode', num2str(mode) ,'}{ ', num2str(round(ThirdStagePayloadMass,1),'%.1f') , '}'), 'delimiter','','newline', 'pc')
dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\12SeparationAltMode', num2str(mode) ,'}{ ', num2str(round(alt21(1)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc')

dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\FirstStageSMFMode', num2str(mode) ,'}{ ', num2str(round(FirstStageSMF,3),'%.3f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\23SeparationAltMode', num2str(mode) ,'}{ ', num2str(round(alt21(end)/1000,2),'%.2f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\23SeparationvMode', num2str(mode) ,'}{ ', num2str(round(v21(end),0)) , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\23SeparationqMode', num2str(mode) ,'}{ ', num2str(round(q21(end)/1000,1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');
dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\23SeparationLDMode', num2str(mode) ,'}{ ', num2str(round(L21(end)/Fd21(end),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\2FlightTimeMode', num2str(mode) ,'}{ ', num2str(round(time21(end),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');

qlt20 = find(q3<20000);
dlmwrite(strcat('LatexInputs',namelist{j},'.txt'),strcat('\newcommand{\3qOver20Mode', num2str(mode) ,'}{ ', num2str(round(time3(qlt20(1))-time3(1),1),'%.1f') , '}'), '-append','delimiter','','newline', 'pc');



%% Bound Check
% Peform check to see if any of the states are hitting their bounds. This
% is an error if the bound is not intended to constrain the state. Fuel
% mass and throttle are not checked, as these will always hit bounds. 

for i = 1: length(output{j}.result.solution.phase(2).state(1,:))
    if any(output{j}.result.solution.phase(1).state(:,i) == bounds.phase(1).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(1).state(:,i) == bounds.phase(1).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 1 is hitting upper bound'))
    end
end

for i = 1: length(output{j}.result.solution.phase(2).state(1,:))-1
    if any(output{j}.result.solution.phase(2).state(:,i) == bounds.phase(2).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(2).state(:,i) == bounds.phase(2).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 2 is hitting upper bound'))
    end
end

for i = 1: length(output{j}.result.solution.phase(3).state(1,:))-2
    if any(output{j}.result.solution.phase(3).state(:,i) == bounds.phase(3).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(3).state(:,i) == bounds.phase(3).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 3 is hitting upper bound'))
    end
end

% Angle of attack is not checked on third stage, because angle of attack is hard constrained and should be checked manually. 
for i = [1:3 6: length(output{j}.result.solution.phase(4).state(1,:))]
    if any(output{j}.result.solution.phase(4).state(:,i) == bounds.phase(4).state.lower(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting lower bound'))
    end
    
    if any(output{j}.result.solution.phase(4).state(:,i) == bounds.phase(4).state.upper(i))
        disp(strcat('State Id: ',num2str(i),' in Phase 4 is hitting upper bound'))
    end
end

%% Plot combined trajectory




        figure(2301)
           fig = gcf;
        set(fig,'Position',[200 0 1200 1000])
hold on
% axesm('pcarree','MapLatLimit',[min(rad2deg(lat22))-5 max(rad2deg(latexo))+5],'MapLonLimit',[min(rad2deg(lon22))-5 max(rad2deg(lon21))+5])
worldmap([min(rad2deg(lat22))-4 max(rad2deg(latexo))+2],[min(rad2deg(lon22))-4 max(rad2deg(lon21))+4])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
tightmap
mlabel('south')
view([-30 35])
% plotm(rad2deg(lat),rad2deg(lon+lon0))
a = plot3m(rad2deg(lat1),rad2deg(lon1),alt1*5,'color',[0 0.8 0],'LineWidth',1.7);
b = plot3m(rad2deg(lat21),rad2deg(lon21),alt21*5,'r','LineWidth',1.7);
c = plot3m(rad2deg(lat22),rad2deg(lon22),alt22*5,'color',[1 0 1],'LineWidth',1.7);
d = plot3m(rad2deg(lat3),rad2deg(lon3)',alt3*5,'b','LineWidth',1.7);
e = plot3m(rad2deg(latexo),rad2deg(lonexo),altexo*5,'c','LineWidth',1.7);

legend([a b c d e],{'First Stage Ascent','Second Stage Ascent', 'Second Stage Return','Third Stage Powered Ascent','Third Stage Unpowered Ascent'},'Location','east');

 ht = text(-500000,-1600000,'Australia');
 set(ht,'Rotation',20)

plotm(rad2deg(lat1),rad2deg(lon1),'color',[0 0.8 0])
plotm(rad2deg(lat21),rad2deg(lon21),'r')
plotm(rad2deg(lat22),rad2deg(lon22),'color',[1 0 1])
plotm(rad2deg(lat3),rad2deg(lon3)','b')
plotm(rad2deg(latexo),rad2deg(lonexo),'c')
   
    cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)

    
for i = 0:20/(time1(end)-time1(1)):1
    time1_temp = (time1(end)-time1(1))*i + time1(1);
    alt1_temp = interp1(time1,alt1,time1_temp);
    lat1_temp = interp1(time1,lat1,time1_temp);
    lon1_temp = interp1(time1,lon1,time1_temp);

    plot3m(rad2deg([lat1_temp  ,lat1_temp ]), rad2deg([lon1_temp,  lon1_temp]),[0,alt1_temp*5],'color',[0.3 0.3 0.3]);

end
for i = 0:20/(time21(end)-time21(1)):1
    time21_temp = (time21(end)-time21(1))*i + time21(1) + rem(time21(1),20);
    alt21_temp = interp1(time21,alt21,time21_temp);
    lat21_temp = interp1(time21,lat21,time21_temp);
    lon21_temp = interp1(time21,lon21,time21_temp);

    plot3m(rad2deg([lat21_temp  ,lat21_temp ]), rad2deg([lon21_temp,  lon21_temp]),[0,alt21_temp*5],'color',[0.3 0.3 0.3]);

end
for i = 0:20/(time22(end)-time22(1)):1
    time22_temp = (time22(end)-time22(1))*i + time22(1) + rem(time22(1),20);
    alt22_temp = interp1(time22,alt22,time22_temp);
    lat22_temp = interp1(time22,lat22,time22_temp);
    lon22_temp = interp1(time22,lon22,time22_temp);

    plot3m(rad2deg([lat22_temp  ,lat22_temp ]), rad2deg([lon22_temp,  lon22_temp]),[0,alt22_temp*5],'color',[0.3 0.3 0.3]);

end
for i = 0:20/(time3(end)-time3(1)):1
    time3_temp = (time3(end)-time3(1))*i + time3(1) + rem(time3(1),20) ;
    alt3_temp = interp1(time3,alt3,time3_temp);
    lat3_temp = interp1(time3,lat3,time3_temp);
    lon3_temp = interp1(time3,lon3,time3_temp);
    plot3m(rad2deg([lat3_temp  ,lat3_temp ]), rad2deg([lon3_temp,  lon3_temp]),[0,alt3_temp*5],'color',[0.3 0.3 0.3]);

end
for i = 0:20/(timeexo(end)-timeexo(1)):1
    timeexo_temp = (timeexo(end)-timeexo(1))*i + timeexo(1) + rem(timeexo(1),20);
    altexo_temp = interp1(timeexo,altexo,timeexo_temp);
    latexo_temp = interp1(timeexo,latexo,timeexo_temp);
    lonexo_temp = interp1(timeexo,lonexo,timeexo_temp);
    plot3m(rad2deg([latexo_temp  ,latexo_temp ]), rad2deg([lonexo_temp,  lonexo_temp]),[0,altexo_temp*5],'color',[0.3 0.3 0.3]);
end

scaleruler on
setm(handlem('scaleruler1'), ...
    'XLoc',-9.5e4,'YLoc',-1.5e6,...
    'MajorTick',0:200:400,'TickDir','down','RulerStyle','patches')

% zoom(30) 
%% Plot Visualisation of Net ISP 

% MList = [5:0.1:10];
altlist = 25000:100:35000 ;
alphalist = 0:.1:5;

[alphagrid,altgrid] = ndgrid(alphalist,altlist);
figure(501)
title('Net Isp')
hold on
for temp = 1:5
    
v_temp = 1620+300*(temp-1);

for i = 1:numel(alphagrid)
        I = cell(1, ndims(alphagrid)); 
    [I{:}] = ind2sub(size(alphagrid),i);
    alpha_temp = alphagrid(I{1},I{2});
    alt_temp = altgrid(I{1},I{2});
    
    
[~,~,~,~,~,~, ~, ~, Fdgrid(I{(1)},I{(2)}), ~,Lgrid(I{(1)},I{(2)}),Fueldtgrid(I{(1)},I{(2)}),Tgrid(I{(1)},I{(2)}),Ispgrid(I{(1)},I{(2)}),~,~,~,~] = VehicleModelCombined(0, alt_temp, v_temp,auxdata,0,0,0,deg2rad(alpha_temp),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0);

end

% [~,~,~,~,~,~, ~, ~, Fd, ~,L,Fueldt,T,Isp,~,~,~,~] = VehicleModelCombined(0, 34000, 1600,auxdata,0,0,0,deg2rad(1),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0)


subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,(Tgrid-Fdgrid)./Fueldtgrid/9.81,1000,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title([num2str(v_temp) ' m/s'],'FontSize',10)
set(gca,'XTick',0:5)
caxis([0 1500]);
alt_temp = interp1(v21,alt21,v_temp);
alpha_temp = interp1(v21,alpha21,v_temp);
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')

end


%% Plot End of ascent net Isp comparison


% MList = [5:0.1:10];
altlist = 35000:100:45000 ;
alphalist = 5:.1:10;

[alphagrid,altgrid] = ndgrid(alphalist,altlist);
figure(502)
title('Net Isp')
hold on
for temp = 1:5
    
v_temp = interp1(time21,v21,time21(end)-21+5*(temp-1));

for i = 1:numel(alphagrid)
        I = cell(1, ndims(alphagrid)); 
    [I{:}] = ind2sub(size(alphagrid),i);
    alpha_temp = alphagrid(I{1},I{2});
    alt_temp = altgrid(I{1},I{2});
    
    
[~,~,~,~,~,~, ~, ~, Fdgrid(I{(1)},I{(2)}), ~,Lgrid(I{(1)},I{(2)}),Fueldtgrid(I{(1)},I{(2)}),Tgrid(I{(1)},I{(2)}),Ispgrid(I{(1)},I{(2)}),~,~,~,~] = VehicleModelCombined(0, alt_temp, v_temp,auxdata,0,0,0,deg2rad(alpha_temp),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0);

end

% [~,~,~,~,~,~, ~, ~, Fd, ~,L,Fueldt,T,Isp,~,~,~,~] = VehicleModelCombined(0, 34000, 1600,auxdata,0,0,0,deg2rad(1),0,1, mFuel21(1),mFuel21(1),mFuel21(end), 1, 0)


subplot(2,3,temp);
hold on
colormap jet
contourf(alphagrid,altgrid/1000,(Tgrid-Fdgrid)./Fueldtgrid/9.81,1000,'LineColor','none')
xlabel('Angle of Attack (deg)','FontSize',8);
ylabel('Altitude (km)','FontSize',8);
title([num2str(v_temp) ' m/s'],'FontSize',10)
set(gca,'XTick',5:10)
caxis([-1000 500]);
alt_temp = interp1(time21,alt21,time21(end)-21+5*(temp-1));
alpha_temp = interp1(time21,alpha21,time21(end)-21+5*(temp-1));
plot(rad2deg(alpha_temp),alt_temp/1000,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','m')

end

%%

c = colorbar;
set(c, 'Position', [.70 .13 .0581 .3])
ylabel(c, 'Net Isp')

% figure(503)
% contourf(alphagrid,altgrid,Fdgrid,1000,'LineWidth',0.)
% c = colorbar;
% 
% figure(504)
% contourf(alphagrid,altgrid,Tgrid,1000,'LineWidth',0.)
% title('Thrust')
% c = colorbar;
% 
% figure(505)
% contourf(alphagrid,altgrid,Fueldtgrid,1000,'LineWidth',0.)
% c = colorbar;





%% SAVE FIGS

saveas(figure(311),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('ThirdStage',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(311),strcat('ThirdStage',namelist{j}),'-dpng');
movefile(strcat('ThirdStage',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(211),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('SecondStage',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(211),strcat('SecondStage',namelist{j}),'-dpng');
movefile(strcat('SecondStage',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(221),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Return',namelist{j},'.fig')]); 
set(gcf, 'PaperPositionMode', 'auto');
print(figure(221),strcat('Return',namelist{j}),'-dpng');
movefile(strcat('Return',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(701),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('HeatFlux',namelist{j},'.fig')]); 
set(gcf, 'PaperPositionMode', 'auto');
print(figure(701),strcat('HeatFlux',namelist{j}),'-dpng');
movefile(strcat('HeatFlux',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(2410),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Hamiltonian',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2410),strcat('Hamiltonian',namelist{j}),'-dpng');
movefile(strcat('Hamiltonian',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(2420),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Verification',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2420),strcat('Verification',namelist{j}),'-dpng');
movefile(strcat('Verification',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(550),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('ForwardError',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(550),strcat('ForwardError',namelist{j}),'-dpng');
movefile(strcat('ForwardError',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(212),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Forward1',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(212),strcat('Forward1',namelist{j}),'-dpng');
movefile(strcat('Forward1',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(213),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Forward2',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(213),strcat('Forward2',namelist{j}),'-dpng');
movefile(strcat('Forward2',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
% saveas(figure(2100),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'eq.fig']);
% saveas(figure(2110),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ISP.fig']);
saveas(figure(2100),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('eq',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2100),strcat('eq',namelist{j}),'-dpng');
movefile(strcat('eq',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(2110),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('Isp',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2110),strcat('Isp',namelist{j}),'-dpng');
movefile(strcat('Isp',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(2301),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('GroundTrack',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperPositionMode', 'auto');
print(figure(2301),strcat('GroundTrack',namelist{j}),'-dpng');
movefile(strcat('GroundTrack',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(501),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('NetIsp',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(501),strcat('NetIsp',namelist{j}),'-dpng');
movefile(strcat('NetIsp',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));
saveas(figure(502),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('NetIspPullup',namelist{j},'.fig')]);
set(gcf, 'PaperPositionMode', 'auto');
print(figure(502),strcat('NetIspPullup',namelist{j}),'-dpng');
movefile(strcat('NetIspPullup',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));


close all


end

%% Plot Comparison Figures

if length(output)>1
    
   
time1nom = output{ceil(length(output)/2)}.result.solution.phase(1).time.';
% 
alt1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,1).';
v1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,2).';
m1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,3).';
gamma1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,4).';
alpha1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,5).';
zeta1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,6).';
lat1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,8).';
lon1nom = output{ceil(length(output)/2)}.result.solution.phase(1).state(:,9).';
% 

alt21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,1);
alt22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,1);
lon21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,2);
lon22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,2);
lat21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,3);
lat22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,3);
v21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,4); 
v22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,4); 
gamma21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,5); 
gamma22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,5); 
zeta21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,6);
zeta22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,6);
alpha21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,7);
alpha22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,7);
eta21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,8);
eta22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,8);
mFuel21nom = output{ceil(length(output)/2)}.result.solution.phase(2).state(:,9); 
mFuel22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,9); 

throttle22nom = output{ceil(length(output)/2)}.result.solution.phase(3).state(:,10);

aoadot21nom  = output{ceil(length(output)/2)}.result.solution.phase(2).control(:,1); 
etadot21nom  = output{ceil(length(output)/2)}.result.solution.phase(2).control(:,2); 

aoadot22nom  = output{ceil(length(output)/2)}.result.solution.phase(3).control(:,1); 
etadot22nom  = output{ceil(length(output)/2)}.result.solution.phase(3).control(:,2); 

time21nom = output{ceil(length(output)/2)}.result.solution.phase(2).time;
time22nom = output{ceil(length(output)/2)}.result.solution.phase(3).time;


alt3nom  = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,1);
v3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,2);
gamma3nom  = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,3);
m3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,4);
aoa3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,5);
lat3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,6);
zeta3nom    = output{ceil(length(output)/2)}.result.solution.phase(4).state(:,7);
aoadot3nom       = output{ceil(length(output)/2)}.result.solution.phase(4).control(:,1);

time3nom = output{ceil(length(output)/2)}.result.solution.phase(4).time;
   
    
[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3nom] = ThirdStageDyn(alt3nom,gamma3nom,v3nom,m3nom,aoa3nom,time3nom,auxdata,aoadot3nom,lat3nom,zeta3nom);

    

[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21nom, M21, Fd21, rho21,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21] = VehicleModelCombined(gamma21nom, alt21nom, v21nom,auxdata,zeta21nom,lat21nom,lon21nom,alpha21nom,eta21nom,1, mFuel21nom,mFuel21nom(1),mFuel21nom(end), 1, 0);
[~,~,~,~,~,~, q22nom, M22nom, Fd22, rho22,L22,Fueldt22,T22,Isp22,q122nom,flapdeflection22,heating_rate22] = VehicleModelCombined(gamma22nom, alt22nom, v22nom,auxdata,zeta22nom,lat22nom,lon22nom,alpha22nom,eta22nom,throttle22nom, mFuel22nom,0,0, 0, 0);
    
throttle22nom(M22nom<5.0) = throttle22nom(M22nom<5.0).*gaussmf(M22nom(M22nom<5.0),[.01,5]); % remove nonsense throttle points
throttle22nom(q122nom<20000) = throttle22nom(q122nom<20000).*gaussmf(q122nom(q122nom<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
%     

    figure(1110)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on

figure(3110)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
    
    figure(2110)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
    
    
    
   figure(2210)
    fig = gcf;
set(fig,'Position',[200 0 850 800])
    hold on
%     Colours = [flip(0.7:0.3/(floor(length(output)/2)-1):1)'  zeros(floor(length(output)/2),1) zeros(floor(length(output)/2),1); 0 0 0; zeros(floor(length(output)/2),1)  zeros(floor(length(output)/2),1) (0.7:0.3/(floor(length(output)/2)-1):1)'];
%     
    Colours = [0.8500, 0.3250, 0.0980; 1 0 0; 0 0 0; 0 0 1; 0.3010, 0.7450, 0.9330];
    
    LineStyleList = {'--','-','-','-','--'};
    
    for j = 1:length(output)
        
time1 = output{j}.result.solution.phase(1).time.';
% 
alt1 = output{j}.result.solution.phase(1).state(:,1).';
v1 = output{j}.result.solution.phase(1).state(:,2).';
m1 = output{j}.result.solution.phase(1).state(:,3).';
gamma1 = output{j}.result.solution.phase(1).state(:,4).';
alpha1 = output{j}.result.solution.phase(1).state(:,5).';
zeta1 = output{j}.result.solution.phase(1).state(:,6).';
lat1 = output{j}.result.solution.phase(1).state(:,8).';
lon1 = output{j}.result.solution.phase(1).state(:,9).';
% 

alt21 = output{j}.result.solution.phase(2).state(:,1);
alt22 = output{j}.result.solution.phase(3).state(:,1);
lon21 = output{j}.result.solution.phase(2).state(:,2);
lon22 = output{j}.result.solution.phase(3).state(:,2);
lat21 = output{j}.result.solution.phase(2).state(:,3);
lat22 = output{j}.result.solution.phase(3).state(:,3);
v21 = output{j}.result.solution.phase(2).state(:,4); 
v22 = output{j}.result.solution.phase(3).state(:,4); 
gamma21 = output{j}.result.solution.phase(2).state(:,5); 
gamma22 = output{j}.result.solution.phase(3).state(:,5); 
zeta21 = output{j}.result.solution.phase(2).state(:,6);
zeta22 = output{j}.result.solution.phase(3).state(:,6);
alpha21 = output{j}.result.solution.phase(2).state(:,7);
alpha22 = output{j}.result.solution.phase(3).state(:,7);
eta21 = output{j}.result.solution.phase(2).state(:,8);
eta22 = output{j}.result.solution.phase(3).state(:,8);
mFuel21 = output{j}.result.solution.phase(2).state(:,9); 
mFuel22 = output{j}.result.solution.phase(3).state(:,9); 

throttle22 = output{j}.result.solution.phase(3).state(:,10);



aoadot21  = output{j}.result.solution.phase(2).control(:,1); 
etadot21  = output{j}.result.solution.phase(2).control(:,2); 

aoadot22  = output{j}.result.solution.phase(3).control(:,1); 
etadot22  = output{j}.result.solution.phase(3).control(:,2); 

time21 = output{j}.result.solution.phase(2).time;
time22 = output{j}.result.solution.phase(3).time;


alt3  = output{j}.result.solution.phase(4).state(:,1);
v3    = output{j}.result.solution.phase(4).state(:,2);
gamma3  = output{j}.result.solution.phase(4).state(:,3);
m3    = output{j}.result.solution.phase(4).state(:,4);
aoa3    = output{j}.result.solution.phase(4).state(:,5);
lat3    = output{j}.result.solution.phase(4).state(:,6);
zeta3    = output{j}.result.solution.phase(4).state(:,7);
aoadot3       = output{j}.result.solution.phase(4).control(:,1);

time3 = output{j}.result.solution.phase(4).time;

        
[rdot3,xidot3,phidot3,gammadot3,vdot3,zetadot3, mdot3, Vec_angle3, AoA_max3, T3, L3, D3, q3] = ThirdStageDyn(alt3,gamma3,v3,m3,aoa3,time3,auxdata,aoadot3,lat3,zeta3);

lon3 = [];
lon3(1) = lon21(end);
for i = 2:length(time3)
    lon3(i) = lon3(i-1) + xidot3(i-1)*(time3(i)-time3(i-1));
end


[altdot21,xidot21,phidot21,gammadot21,a21,zetadot21, q21, M21, Fd21, rho21,L21,Fueldt21,T21,Isp21,q121,flapdeflection21,heating_rate21,CG21] = VehicleModelCombined(gamma21, alt21, v21,auxdata,zeta21,lat21,lon21,alpha21,eta21,1, mFuel21,mFuel21(1),mFuel21(end), 1, 0);
[~,~,~,~,~,~, q22, M22, Fd22, rho22,L22,Fueldt22,T22,Isp22,q122,flapdeflection22,heating_rate22] = VehicleModelCombined(gamma22, alt22, v22,auxdata,zeta22,lat22,lon22,alpha22,eta22,throttle22, mFuel22,0,0, 0, 0);
% 
[AltF_actual, v3F, altexo, v3exo, timeexo, mpayload, Alpha3, mexo,qexo,gammaexo,Dexo,zetaexo,latexo,incexo,Texo,CLexo,Lexo,incdiffexo,lonexo] = ThirdStageSim(alt3(end),gamma3(end),v3(end), lat3(end),lon3(end), zeta3(end), m3(end), auxdata);

throttle22(M22<5.0) = throttle22(M22<5.0).*gaussmf(M22(M22<5.0),[.01,5]); % remove nonsense throttle points
throttle22(q122<20000) = throttle22(q122<20000).*gaussmf(q122(q122<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
%     




        figure(3110)
    subplot(6,2,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)], [alt3; altexo.']/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

subplot(6,2,2)
    hold on
    title('Altitude Difference (km)','FontSize',9);

 
    plot(0:1/(length(time3)-1):1, (alt3'-interp1(0:1/(length(time3nom)-1):1,alt3nom,0:1/(length(time3)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])


    subplot(6,2,3)
    hold on
    title('Dynamic Pressure (kPa','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[q3;qexo.';qexo(end)]/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

subplot(6,2,4)
    hold on
    title('Dynamic Pressure Difference (kPa)','FontSize',9);
    plot(0:1/(length(time3)-1):1, (q3'-interp1(0:1/(length(time3nom)-1):1,q3nom,0:1/(length(time3)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])



    subplot(6,2,5)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[rad2deg(aoa3);0*ones(length(timeexo),1)],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

subplot(6,2,6)
    hold on
    title('Angle of Attack Difference (deg)','FontSize',9);
    plot(0:1/(length(time3)-1):1, rad2deg(aoa3'-interp1(0:1/(length(time3nom)-1):1,aoa3nom,0:1/(length(time3)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

%     xlim([time3(1) timeexo(end)+time3(end)])
    subplot(6,2,7)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)],[v3;v3exo.'],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

subplot(6,2,8)
    hold on
    title('Velocity Difference (m/s)','FontSize',9);
    plot(0:1/(length(time3)-1):1, (v3'-interp1(0:1/(length(time3nom)-1):1,v3nom,0:1/(length(time3)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])


%     xlabel('Time (s)','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])
    subplot(6,2,9)
    hold on
    title('Trajectory Angle (deg)','FontSize',9);
    plot([time3-time3(1); timeexo.'+time3(end)-time3(1)], [rad2deg(gamma3);rad2deg(gammaexo).'],'Color',Colours(j,:),'LineStyle',LineStyleList{j})

    xlabel('Time (s)','FontSize',9);
    
    
    subplot(6,2,10)
    hold on
    title('Flight Path Angle Difference (Deg)','FontSize',9);
    plot(0:1/(length(time3)-1):1, rad2deg(gamma3'-interp1(0:1/(length(time3nom)-1):1,gamma3nom,0:1/(length(time3)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])
xlabel('Normalised Time','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])



        figure(2110)
    subplot(6,2,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time21-time21(1)], [alt21]/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

subplot(6,2,2)
    hold on
    title('Altitude Difference (km)','FontSize',9);

 
    plot(0:1/(length(time21)-1):1, (alt21'-interp1(0:1/(length(time21nom)-1):1,alt21nom,0:1/(length(time21)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])


    subplot(6,2,3)
    hold on
    title('Dynamic Pressure (kPa','FontSize',9);
    plot([time21-time21(1)],[q21]/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time21(1) timeexo(end)+time21(end)])

subplot(6,2,4)
    hold on
    title('Dynamic Pressure Difference (kPa)','FontSize',9);
    plot(0:1/(length(time21)-1):1, (q21'-interp1(0:1/(length(time21nom)-1):1,q21nom,0:1/(length(time21)-1):1))/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])



    subplot(6,2,5)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time21-time21(1)],[rad2deg(alpha21)],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

subplot(6,2,6)
    hold on
    title('Angle of Attack Difference (deg)','FontSize',9);
    plot(0:1/(length(time21)-1):1, rad2deg(alpha21'-interp1(0:1/(length(time21nom)-1):1,alpha21nom,0:1/(length(time21)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

%     xlim([time21(1) timeexo(end)+time21(end)])
    subplot(6,2,7)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time21-time21(1)],[v21],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

subplot(6,2,8)
    hold on
    title('Velocity Difference (m/s)','FontSize',9);
    plot(0:1/(length(time21)-1):1, (v21'-interp1(0:1/(length(time21nom)-1):1,v21nom,0:1/(length(time21)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])


%     xlabel('Time (s)','FontSize',9);
%     xlim([time21(1) timeexo(end)+time21(end)])
    subplot(6,2,9)
    hold on
    title('Trajectory Angle (deg)','FontSize',9);
    plot([time21-time21(1)], [rad2deg(gamma21)],'Color',Colours(j,:),'LineStyle',LineStyleList{j})

    xlabel('Time (s)','FontSize',9);
    
    
    subplot(6,2,10)
    hold on
    title('Flight Path Angle Difference (Deg)','FontSize',9);
    plot(0:1/(length(time21)-1):1, rad2deg(gamma21'-interp1(0:1/(length(time21nom)-1):1,gamma21nom,0:1/(length(time21)-1):1))','Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])
xlabel('Normalised Time','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])


        figure(2210)
    subplot(6,2,1)
    hold on
    title('Altitude (km','FontSize',9);
    plot([time22-time22(1)], [alt22]/1000','Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])

subplot(6,2,2)
    hold on
    title('Altitude Difference (km)','FontSize',9);

 
    plot(time22-time22(1), (alt22'-interp1(time22nom,alt22nom,time22)')/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time3(1) timeexo(end)+time3(end)])


    subplot(6,2,3)
    hold on
    title('Dynamic Pressure (kPa','FontSize',9);
    plot([time22-time22(1)],[q22]/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
%     xlim([time22(1) timeexo(end)+time22(end)])

subplot(6,2,4)
    hold on
    title('Dynamic Pressure Difference (kPa)','FontSize',9);
    plot(time22-time22(1), (q22'-interp1(time22nom,q22nom,time22)')/1000,'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])



    subplot(6,2,5)
    hold on
    title('Angle of Attack (deg)','FontSize',9);
    plot([time22-time22(1)],[rad2deg(alpha22)],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

subplot(6,2,6)
    hold on
    title('Angle of Attack Difference (deg)','FontSize',9);
    plot(time22-time22(1), rad2deg(alpha22'-interp1(time22nom,alpha22nom,time22)'),'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

%     xlim([time22(1) timeexo(end)+time22(end)])
    subplot(6,2,7)
    hold on
    title('Velocity (m/s)','FontSize',9);
    plot([time22-time22(1)],[v22],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])

subplot(6,2,8)
    hold on
    title('Velocity Difference (m/s)','FontSize',9);
    plot(time22-time22(1), (v22'-interp1(time22nom,v22nom,time22)'),'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])


%     xlabel('Time (s)','FontSize',9);
%     xlim([time22(1) timeexo(end)+time22(end)])
    subplot(6,2,9)
    hold on
    title('Trajectory Angle (deg)','FontSize',9);
    plot([time22-time22(1)], [rad2deg(gamma22)],'Color',Colours(j,:),'LineStyle',LineStyleList{j})
set(gca,'xticklabels',[])
    
    
    
    subplot(6,2,10)
    hold on
    title('Flight Path Angle Difference (Deg)','FontSize',9);
    plot(time22-time22(1), rad2deg(gamma22'-interp1(time22nom,gamma22nom,time22)'),'Color',Colours(j,:),'LineStyle',LineStyleList{j})
% set(gca,'xticklabels',[])
xlabel('Time (s)','FontSize',9);
%     xlim([time3(1) timeexo(end)+time3(end)])

subplot(6,2,11)
    hold on
    title('Throttle','FontSize',9);
    plot([time22-time22(1)], [throttle22],'Color',Colours(j,:),'LineStyle',LineStyleList{j})

    xlabel('Time (s)','FontSize',9);

    end
    figure(2110)
    legend(namelist)
    figure(2210)
    legend(namelist)
    figure(3110)
    legend(namelist)
    
saveas(figure(2110),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('SecondStageComparison',namelist{j},'.fig')]);
print(figure(2110),strcat('SecondStageComparison',namelist{j}),'-dpng');
movefile(strcat('SecondStageComparison',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));

saveas(figure(2210),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('ReturnComparison',namelist{j},'.fig')]);
print(figure(2210),strcat('ReturnComparison',namelist{j}),'-dpng');
movefile(strcat('ReturnComparison',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));


saveas(figure(3110),[sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))),filesep,strcat('ThirdStageComparison',namelist{j},'.fig')]);
print(figure(3110),strcat('ThirdStageComparison',namelist{j}),'-dpng');
movefile(strcat('ThirdStageComparison',namelist{j},'.png'),sprintf('../ArchivedResults/%s',strcat(Timestamp,'mode',num2str(mode))));


end
close all
end

