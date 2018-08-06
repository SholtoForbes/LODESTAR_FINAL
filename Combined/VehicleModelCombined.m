function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp,q1,flap_deflection,heating_rate,CG] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,phi,xi,alpha,eta,throttle,mFuel,mFuelinit,mFuelend,ThirdStage,forwardflag)
%===================================================
%
% SPARTAN DYNAMICS SIMULATION
%
%===================================================

interp = auxdata.interp;
% =======================================================
% Vehicle Model
% =======================================================
A = auxdata.A; % reference area (m^2)

if ThirdStage == 1
m = auxdata.Stage2.mStruct+mFuel+auxdata.Stage3.mTot; 
else
m = auxdata.Stage2.mStruct+mFuel;
end

%% Flow =============================================================
c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
mach = v./c;
rho = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

T0 = ppval(interp.T0_spline, alt); 

P0 = ppval(interp.P0_spline, alt);
%% Thrust 


[Isp,Fueldt_max,eq,q1] = RESTint(M, rad2deg(alpha), auxdata,T0,P0); % Calculate C-REST engine properties

Isp = Isp*auxdata.Ispmod; % modify Isp by a set fraction if necessary

% Isp(q1<20000) = Isp(q1<20000).*gaussmf(q1(q1<20000),[1000,20000]); % rapidly reduce ISP to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.


 if ThirdStage == 0 && forwardflag ==0;

    % Turn off throttle at unoperable flight conditions for aerodynamic and
    % thrust purposes
    throttle(q1<20000) = throttle(q1<20000).*gaussmf(q1(q1<20000),[100,20000]); % rapidly reduce throttle to 0 after passing the lower limit of 20kPa dynamic pressure. This dynamic pressure is after the conical shock.
    throttle(M<5.0) =   throttle(M<5.0).*gaussmf(M(M<5.0),[.01,5]); % remove throttle points below operable range on return flight
end

Fueldt = Fueldt_max.*throttle; %

% T = Isp.*Fueldt*9.81.*cos(alpha).*gaussmf(throttle,[0.1,1]); % Thrust in
% direction of , modified by a gaussmf funtion to reduce thrutst rapidly
T = Isp.*Fueldt_max.*throttle*9.81.*cos(alpha); % Thrust in direction of motion

%% Aerodynamics
% Calculate aerodynamic coefficients

   
    % Interpolate between centre of gravity conditions for full, cylindrical tank empty, and empty conditions as fuel depletes
mFuel_cyltanks = 710;

CG_withFuel_noThirdStage = 14.518; % CG from CREO with no third stage but full fuel
CG_cyltanksEmpty_noThirdStage = 14.297;
CG_cyltanksEmpty_ThirdStage = (CG_cyltanksEmpty_noThirdStage*(4957.1+710)+16.63*3300)/(4957+710+3300); % CG with third stage, but no fuel
CG_noThirdStage = (CG_withFuel_noThirdStage*(4957+1562)-12.59*1562)/4957;% CG  at no fuel conditition (it is assumed that the fuel for the return is used so as to not change the CG)
EngineOn_CG_noFuel = (CG_noThirdStage*4957.1+16.63*3300)/(4957+3300); % CG with third stage, but no fuel
EngineOn_CG_Fuel = (EngineOn_CG_noFuel*8257.1 + 12.59*1562)/(8257.1+1562); % CG with third stage and full fuel
% 

if ThirdStage == 1 % Ascent with third stage

% Determine trajectory points which are using fuel from cylindrical fuel tank
index_overcyl = mFuel>(auxdata.Stage2.mFuel-mFuel_cyltanks);

% Determine the fraction of fuel in the cylindrical tank
Proportion_fulltocyl  = (mFuel(index_overcyl)-(auxdata.Stage2.mFuel-mFuel_cyltanks))./(mFuel_cyltanks);

% Calculate CG variation
CG(index_overcyl) = Proportion_fulltocyl*EngineOn_CG_Fuel + (1-Proportion_fulltocyl)*CG_cyltanksEmpty_ThirdStage;

% Calculate Cd
Cd(index_overcyl) = Proportion_fulltocyl.*auxdata.interp.Cd_spline_EngineOn.fullFuel(mach(index_overcyl),...
    rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000) + (1-Proportion_fulltocyl).*...
    auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach(index_overcyl),rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000);

% Calculate Cl
Cl(index_overcyl) = Proportion_fulltocyl.*auxdata.interp.Cl_spline_EngineOn.fullFuel(mach(index_overcyl),...
    rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000) + (1-Proportion_fulltocyl).*...
    auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach(index_overcyl),rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000);

% Calculate frap deflection
flap_deflection(index_overcyl) = Proportion_fulltocyl.*auxdata.interp.flap_spline_EngineOn.fullFuel(mach(index_overcyl),...
    rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000) + (1-Proportion_fulltocyl).*...
    auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach(index_overcyl),rad2deg(alpha(index_overcyl)),alt(index_overcyl)/1000);

% Determine trajectory points which are using fuel from the conical tank
index_undercyl = mFuel<=(auxdata.Stage2.mFuel-mFuel_cyltanks);

% Determine fraction of fuel in conical tank
Proportion_EmptytoCyl = ((auxdata.Stage2.mFuel-mFuel_cyltanks)-mFuel(index_undercyl))./(auxdata.Stage2.mFuel-mFuel_cyltanks);

% Calculate CG
CG(index_undercyl) = Proportion_EmptytoCyl*EngineOn_CG_noFuel + (1-Proportion_EmptytoCyl)*CG_cyltanksEmpty_ThirdStage;

Cd(index_undercyl) = Proportion_EmptytoCyl.*auxdata.interp.Cd_spline_EngineOn.noFuel(mach(index_undercyl),...
    rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach(index_undercyl),rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000);

Cl(index_undercyl) = Proportion_EmptytoCyl.*auxdata.interp.Cl_spline_EngineOn.noFuel(mach(index_undercyl),...
    rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach(index_undercyl),rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000);

flap_deflection(index_undercyl) = Proportion_EmptytoCyl.*auxdata.interp.flap_spline_EngineOn.noFuel(mach(index_undercyl),...
    rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000) + (1-Proportion_EmptytoCyl).*...
auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach(index_undercyl),rad2deg(alpha(index_undercyl)),alt(index_undercyl)/1000);

Cd = Cd.';
Cl = Cl.';
flap_deflection = flap_deflection.';

else  % Return, no third stage
Proportion_EmptytoCyl = ((auxdata.Stage2.mFuel-mFuel_cyltanks)-mFuel)./(auxdata.Stage2.mFuel-mFuel_cyltanks);

Cd_Engineoff = Proportion_EmptytoCyl.*auxdata.interp.Cd_spline_EngineOff.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cd_spline_EngineOff.cylTankEnd(mach,rad2deg(alpha),alt/1000);

Cl_Engineoff = Proportion_EmptytoCyl.*auxdata.interp.Cl_spline_EngineOff.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cl_spline_EngineOff.cylTankEnd(mach,rad2deg(alpha),alt/1000);

flap_deflection_Engineoff = Proportion_EmptytoCyl.*auxdata.interp.flap_spline_EngineOff.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
auxdata.interp.flap_spline_EngineOff.cylTankEnd(mach,rad2deg(alpha),alt/1000);

Cd_Engineon = Proportion_EmptytoCyl.*auxdata.interp.Cd_spline_EngineOn.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cd_spline_EngineOn.cylTankEnd(mach,rad2deg(alpha),alt/1000);

Cl_Engineon = Proportion_EmptytoCyl.*auxdata.interp.Cl_spline_EngineOn.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
    auxdata.interp.Cl_spline_EngineOn.cylTankEnd(mach,rad2deg(alpha),alt/1000);

flap_deflection_Engineon = Proportion_EmptytoCyl.*auxdata.interp.flap_spline_EngineOn.noFuel(mach,...
    rad2deg(alpha),alt/1000) + (1-Proportion_EmptytoCyl).*...
auxdata.interp.flap_spline_EngineOn.cylTankEnd(mach,rad2deg(alpha),alt/1000);

 Cd = (1-throttle).*Cd_Engineoff + throttle.*Cd_Engineon;
Cl = (1-throttle).*Cl_Engineoff + throttle.*Cl_Engineon;  
flap_deflection = (1-throttle).*flap_deflection_Engineoff + throttle.*flap_deflection_Engineon;  
   
    
    %Interpolate between engine on and engine off cases as throttle is adjusted    
% Cd = (1-throttle).*auxdata.interp.Cd_spline_EngineOff.noThirdStage(mach,rad2deg(alpha),alt/1000) + throttle.*auxdata.interp.Cd_spline_EngineOn.noThirdStage(mach,rad2deg(alpha),alt/1000);
% Cl = (1-throttle).*auxdata.interp.Cl_spline_EngineOff.noThirdStage(mach,rad2deg(alpha),alt/1000) + throttle.*auxdata.interp.Cl_spline_EngineOn.noThirdStage(mach,rad2deg(alpha),alt/1000);  
% flap_deflection = (1-throttle).*auxdata.interp.flap_spline_EngineOff.noThirdStage(mach,rad2deg(alpha),alt/1000) + throttle.*auxdata.interp.flap_spline_EngineOn.noThirdStage(mach,rad2deg(alpha),alt/1000);  

end

%%%% Compute the drag and lift:
D = 0.5*Cd.*A.*rho.*v.^2*auxdata.Cdmod;
L = 0.5*Cl.*A.*rho.*v.^2;

%Rotational Coordinates =================================================

[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoords(alt+auxdata.Re,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta,auxdata.delta);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

v_H = v.*cos(gamma);

%%heating---------------------------
% From Conceptual Shape Optimization of Entry Vehicles, Dirkx & Mooj & NASA
% lecture

%using hot wall correction

kappa = 1.7415e-4; % sutton-graves, from nasa lecture
Rn = 0.005; %effective nose radius (m) 

heating_rate = kappa*sqrt(rho./Rn).*v.^3; %W/m^2


% =========================================================================
end








