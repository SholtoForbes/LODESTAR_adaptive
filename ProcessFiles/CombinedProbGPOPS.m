%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket-Scramjet-Rocket Launch Optimiser
% By Sholto Forbes-Spyratos
% Utilises the GPOPS-2 proprietary optimisation software
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

disp('---------------------------------------------------------------------')
disp('|                                                                   |')
disp('|                 *         Initialising          *                 |')
disp('|                                                                   |')
disp('|             _     ___   ___   ___  ___  _____  _    ___           |')
disp('|            | |   / _ \ |   \ | __|/ __||_   _|/_\  | _ \          |')
disp('|            | |__| (_) || |) || _| \__ \  | | / _ \ |   /          |')
disp('|            |____|\___/ |___/ |___||___/  |_|/_/ \_\|_|_\          |')
disp('|                                                                   |')
disp('|                                                                   |')
disp('|                Authored by Sholto Forbes-Spyratos                 |')
disp('---------------------------------------------------------------------')
%% Import Vehicle and trajectory Config Data %%============================

% run Config.m

auxdata.Mission = Mission;
auxdata.Stage3 = Stage3;
auxdata.Stage2 = Stage2;



%%
disp('|                 *       Setting Run Mode        *                 |')

% Set modifiers for mode variation
auxdata.Ispmod = 1;
auxdata.Cdmod = 1;
auxdata.vCdmod = 1;
auxdata.m3mod = 1;
auxdata.Isp3mod = 1;
auxdata.Cd3mod = 1;

% addpath('.\Processing\num2words')

% Mode 1
if mode == 1

namelist{1} = 'Standard';
end

% mode 0
if mode == 0

namelist{1} = 'Alternate';
end

if mode == 90

namelist{1} = 'Constq';
end
% Mode 2
q_vars = [40000 45000 50000 55000 60000];  % Set dynamic pressures to be investigated
if mode == 2

%     for i = 1:length(q_vars)
%         namelist{i} = ['qmax' num2words(q_vars(i)/1000)];
%     end
if returnMode == 1
namelist = {'qForty' 'qFortyFive' 'qStandard' 'qFiftyFive' 'qSixty'};
else
namelist = {'qFortyNoReturn' 'qFortyFiveNoReturn' 'qStandardNoReturn' 'qFiftyFiveNoReturn' 'qSixtyNoReturn'};
end
end

% Mode 3
Isp_vars = [0.9 0.95 1 1.05 1.1] ;
if mode == 3

%     for i = 1:length(Isp_vars)
%         namelist{i} = ['Isp' num2words(Isp_vars(i)*100) '%'];
%     end
if returnMode == 1
    namelist = {'IspNinety' 'IspNinetyFive' 'IspStandard' 'IspOneHundredFive' 'IspOneHundredTen'};
    else
namelist = {'IspNinetyNoReturn' 'IspNinetyFiveNoReturn' 'IspStandardNoReturn' 'IspOneHundredFiveNoReturn' 'IspOneHundredTenNoReturn'};
end
end
% Mode 4
Cd_vars = [0.9 0.95 1 1.05 1.1] ;
if mode == 4
 
%     for i = 1:length(Cd_vars)
%         namelist{i} = ['Cd' num2words(Cd_vars(i)*100) '%'];
%     end 
if returnMode == 1
namelist = {'CdNinety' 'CdNinetyFive' 'CdStandard' 'CdOneHundredFive' 'CdOneHundredTen'};
else
namelist = {'CdNinetyNoReturn' 'CdNinetyFiveNoReturn' 'CdStandardNoReturn' 'CdOneHundredFiveNoReturn' 'CdOneHundredTenNoReturn'};
end
end

% Mode 5
vCd_vars = [0.2 0.5 1 1.07 1.15];
if mode == 5 
%     for i = 1:length(vCd_vars)
%         namelist{i} = ['vCd' num2words(vCd_vars(i)*100) '%'];
%     end  
if returnMode == 1
namelist = {'vCdTwenty' 'vCdFifty' 'vCdStandard' 'vCdOneHundredSeven' 'vCdOneHundredFifteen'};
else
namelist = {'vCdTwentyNoReturn' 'vCdFiftyNoReturn' 'vCdStandardNoReturn' 'vCdOneHundredSevenNoReturn' 'vCdOneHundredFifteenNoReturn'};
end
end
    
% Mode 6
% First Stage Structural Mass (total mass kept constant)
FirstStagem_vars = [0.90 0.95 1 1.05 1.1] ;
if mode == 6 
if returnMode == 1
namelist = {'FirstStagemNinety' 'FirstStagemNinetyFive' 'FirstStagemStandard' 'FirstStagemOneHundredFive' 'FirstStagemOneHundredTen'};
else
namelist = {'FirstStagemNinetyNoReturn' 'FirstStagemNinetyFiveNoReturn' 'FirstStagemStandardNoReturn' 'FirstStagemOneHundredFiveNoReturn' 'FirstStagemOneHundredTenNoReturn'};
end
end

% Mode 7
Cd3_vars = [0.80 0.9 1 1.1 1.2]  ;
if mode == 7
 if returnMode == 1
        namelist = {'CdThreeEighty' 'CdThreeNinety' 'CdThreeStandard' 'CdThreeOneHundredTen' 'CdThreeOneHundredTwenty'};
 else
     namelist = {'CdThreeEightyNoReturn' 'CdThreeNinetyNoReturn' 'CdThreeStandardNoReturn' 'CdThreeOneHundredTenNoReturn' 'CdThreeOneHundredTwentyNoReturn'};
 end
end   
   
% Mode 8
m3_vars = [0.9 0.95 1 1.05 1.1] ;
if mode == 8 
%     for i = 1:length(m3_vars)
%         namelist{i} = ['m3' num2words(m3_vars(i)*100) '%'];
%     end  
if returnMode == 1
namelist = {'mThreeNinety' 'mThreeNinetyFive' 'mThreeStandard' 'mThreeOneHundredFive' 'mThreeOneHundredTen'};
else
namelist = {'mThreeNinetyNoReturn' 'mThreeNinetyFiveNoReturn' 'mThreeStandardNoReturn' 'mThreeOneHundredFiveNoReturn' 'mThreeOneHundredTenNoReturn'};
end
end 
    
% Mode 9
Isp3_vars = [0.95 0.975 1 1.025 1.05];  
if mode == 9
%     for i = 1:length(Isp3_vars)
%         namelist{i} = ['T3' num2words(Isp3_vars(i)*100) '%'];
%     end 
if returnMode == 1
namelist = {'ISPThreeNinetyFive' 'ISPThreeNinetySevenFive' 'ISPThreeStandard' 'ISPThreeOneHundredTwoFive' 'ISPThreeOneHundredFive'};
else
namelist = {'ISPThreeNinetyFiveNoReturn' 'ISPThreeNinetySevenFiveNoReturn' 'ISPThreeStandardNoReturn' 'ISPThreeOneHundredTwoFiveNoReturn' 'ISPThreeOneHundredFiveNoReturn'};
end
end  
    
% mode 10
mSPARTAN_vars = [0.95 0.975 1 1.025 1.05] ;
if mode == 10 
%     for i = 1:length(mSPARTAN_vars)
%         namelist{i} = ['mSPARTAN' num2words(mSPARTAN_vars(i)*100) '%'];
%     end   
if returnMode == 1
namelist = {'mSPARTANNinetyFive' 'mSPARTANNinetySevenFive' 'mSPARTANStandard' 'mSPARTANOneHundredTwoFive' 'mSPARTANOneHundredFive'};
else
namelist = {'mSPARTANNinetyFiveNoReturn' 'mSPARTANNinetySevenFiveNoReturn' 'mSPARTANStandardNoReturn' 'mSPARTANOneHundredTwoFiveNoReturn' 'mSPARTANOneHundredFiveNoReturn'};
end
end 
 
%mode 11
mFuel_vars = [0.9 0.95 1 1.05 1.1]     ;
if mode == 11
%     for i = 1:length(mFuel_vars)
%         namelist{i} = ['mFuel' num2words(mFuel_vars(i)*100) '%'];
%     end  
if returnMode == 1
namelist = {'mFuelNinety' 'mFuelNinetyFive' 'mFuelStandard' 'mFuelOneHundredFive' 'mFuelOneHundredTen'};
else
namelist = {'mFuelNinetyNoReturn' 'mFuelNinetyFiveNoReturn' 'mFuelStandardNoReturn' 'mFuelOneHundredFiveNoReturn' 'mFuelOneHundredTenNoReturn'};
end

end

%interactionmode
if mode == 99
 namelist = {};
end

auxdata.namelist = namelist;

%% Misc Modifiers
auxdata.delta = deg2rad(0); % thrust vector angle 


%% Atmosphere Data %%======================================================
% Fetch atmospheric data and compute interpolation splines.
disp('|                 *  Importing Atmosphere Data    *                 |')


Atmosphere = dlmread('atmosphere.txt');
auxdata.Atmosphere = Atmosphere;
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
auxdata.interp.rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data
auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 
auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 



%% Import and Define Parameters
lat0 = Mission.lat0; % Import initial launch location (Also location of SPARTAN landing)
lon0 = Mission.lon0;

auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)

auxdata.A = Stage2.refA ; % SPARTAN reference area




%% FIRST STAGE
disp('|                 *  Importing First Stage Data   *                 |')



% Aerodynamics File Path
Aero = dlmread('FirstStageAero23.5262');


auxdata.Stage1 = Stage1;


%% Calculate Aerodynamic Splines
disp('|            * Interpolating First Stage Aerodynamics *             |')

interp.Lift = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,3));
interp.Drag = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));
interp.Moment = scatteredInterpolant(Aero(:,1),Aero(:,2),Aero(:,4));

M_list = unique(sort(Aero(:,1))); % create unique list of Mach numbers from engine data
M_interp = unique(sort(Aero(:,1)));

AoA_list = unique(sort(Aero(:,2))); % create unique list of angle of attack numbers from engine data 
AoA_interp = unique(sort(Aero(:,2)));

[grid.M,grid.AoA] =  ndgrid(M_interp,AoA_interp);
grid.Lift = interp.Lift(grid.M,grid.AoA);
grid.Drag = interp.Drag(grid.M,grid.AoA);
grid.Moment = interp.Moment(grid.M,grid.AoA);

auxdata.interp.LiftGridded = griddedInterpolant(grid.M,grid.AoA,grid.Lift,'spline','linear');
auxdata.interp.DragGridded = griddedInterpolant(grid.M,grid.AoA,grid.Drag,'spline','linear');
auxdata.interp.MomentGridded = griddedInterpolant(grid.M,grid.AoA,grid.Moment,'spline','linear');

mTotal = Stage2.mStruct+Stage2.mFuel+Stage3.mTot + Stage1.m;
mEmpty = Stage1.m-Stage1.mFuel;  %(kg)  %mass of the rocket (without fuel)


%% Assign Pitchover Conditions

%Define initial conditions at pitchover, these are assumed
h0 = 90;  
v0 = 15;    

gamma0 = deg2rad(89.9);    % set pitchover amount (start flight angle). This pitchover is 'free' movement, and should be kept small. 

mF = mEmpty+Stage2.mStruct+Stage2.mFuel +Stage3.mTot;  %Assume that we use all of the fuel


alpha0 = 0; %Set initial angle of attack to 0

% Altitude bounds (m)
AltMin1 = 1;   %Cannot go through the earth
AltMax1 = 40000;  % maximum first stage altitude

% Velocity bounds (m/s)
vMin1 = 0; 
vMax1 = 3000;  % Maximum first stage velocity, should be considerably higher than max possible

% Mass bounds (kg)
mMin1 = mEmpty;
mMax1 = mTotal;

% Latitude bounds (rad)  (for all stages)
latMin = lat0-1; 
latMax = lat0+1;

% Heading angle bounds (rad)  (for all stages)
zetaMin = -3*pi;
zetaMax = 3*pi;

% Angle of attack bounds (rad)
alphaMin1 = -deg2rad(5);
alphaMax1 = deg2rad(0);

% Angle of attack change rate bounds (rad)
dalphadt2Min1 = -0.1;
dalphadt2Max1 = 0.1;

% Longitude bounds (for all stages)
lonMin =  lon0-1;         lonMax = lon0+1;

% Trajectory angle limits
gammaMin1 = deg2rad(-.1);
gammaMax1 = gamma0;

% This sets the control limits, this is second derivative of AoA
uMin1 = [-.0005]; 
uMax1 = [.0005];

% time bound
tfMax1 	    = 300;     % large upper bound; do not choose Inf
		 
%%
%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%

disp('|                *  Setting-Up First Stage States  *                |')

% Define Time Bounds
bounds.phase(1).initialtime.lower = 0; % constrain initial time to 0
bounds.phase(1).initialtime.upper = 0;

bounds.phase(1).finaltime.lower = 0;
bounds.phase(1).finaltime.upper = tfMax1;

 % Designate State Names
auxdata.States_ID1 = {'Altitude','Velocity','Mass','FPA','AoA','Heading','AoAdt','Latitude','Longitude'};

% Define Initial State Bounds
bounds.phase(1).initialstate.lower = [h0, v0,  mF, gamma0, alpha0,  zetaMin, dalphadt2Min1, lat0, lon0 ];
bounds.phase(1).initialstate.upper = [h0, v0, mMax1, gamma0, alpha0, zetaMax, dalphadt2Max1, lat0, lon0];

% Define Continuous State Bounds
bounds.phase(1).state.lower = [AltMin1, vMin1, mF, gammaMin1, alphaMin1, zetaMin, dalphadt2Min1, latMin, lonMin ];
bounds.phase(1).state.upper = [ AltMax1,  vMax1, mMax1, gammaMax1, alphaMax1, zetaMax, dalphadt2Max1, latMax, lonMax];

% Define FInal State Bounds
bounds.phase(1).finalstate.lower = [AltMin1, vMin1, mF, gammaMin1, alphaMin1, zetaMin, dalphadt2Min1, latMin, lonMin ];
bounds.phase(1).finalstate.upper = [ AltMax1,  vMax1, mMax1, gammaMax1, alphaMax1, zetaMax, dalphadt2Max1, latMax, lonMax];

% Define Control Bounds
bounds.phase(1).control.lower = uMin1;
bounds.phase(1).control.upper = uMax1;

% Define Path Bounds (Path bounds values are calculated in Continuous.m)
bounds.phase(1).path.lower = 0;
bounds.phase(1).path.upper = Stage2.maxq; % set maximum dynamic pressure 

% Define Event Constraints (Events are calculated in Endpoints.m)
bounds.eventgroup(1).lower = [zeros(1,8)];
bounds.eventgroup(1).upper = [zeros(1,8)]; 

% Define Guesses
guess.phase(1).time =  [0; tfMax1];

guess.phase(1).state(:,1) = [h0; h0];
guess.phase(1).state(:,2) = [v0; 1500];
guess.phase(1).state(:,3) = [mMax1; mF];
guess.phase(1).state(:,4) = [gamma0; 0];
guess.phase(1).state(:,5) = [alpha0; 0];

if mode == 0
guess.phase(1).state(:,6) = [-pi;-pi];
else
guess.phase(1).state(:,6) = [0; 0];
end
guess.phase(1).state(:,7) = [0; 0];
guess.phase(1).state(:,8) = [lat0; lat0];
guess.phase(1).state(:,9) = [lon0; lon0];

guess.phase(1).control = [0; 0];


%% Create Interpolation Splines for the SPARTAN
% Aerodynamic and propulsion data splines are created for the SPARTAN,
% utilising SPARTANint.m within the 'Interpolators' folder

[auxdata, MaxviscAlt21, MaxviscAlt22] = SPARTANint(auxdata); 


%% Set Up Second Stage Bounds ---------------------------------------------

disp('|                * Setting-Up Second Stage States *                 |')

% Set State Bounds for Second Stage Scramjet (Those that are not not universal bounds)

aoaMin21 = 0;  aoaMax21 = 10*pi/180; % This is a physical bound, to limit the vehicle between 0-10 degrees angle of attack. AoA limits must be within the bounds of the aerodynamic data available. 

bankMin21 = -90*pi/180; % Banking is a physical constraint on the system. 
bankMax21 =   90*pi/180; 

AltMin21 = 1;
AltMax21 = MaxviscAlt21; % Altitude is bounded to the maximum altitude available from the viscous aerodynamic database, and should be within the bounds of this database to ensure extrapolation does not occur. 

vMin21 = 1000; % Velocity is constrained to encompass the reasonable velocity range of the vehicle. 
vMax21 = 3000;

mFuelMin21 = 0; % Fuel mass is constrained to the physical limits of fuel within the vehicle. 
mFuelMax21 = Stage2.mFuel;

gammaMin21 = -0.5; % The flight path angle is constrained to encompass the reeasonable flight path angle range of the vehicle. 
gammaMax21 = deg2rad(15);


% Designate State Names
auxdata.States_ID21 = {'Altitude', 'Longitude', 'Latitude','Velocity','FPA','Heading','AoA', 'Bank Angle','Fuel Mass'}; 

% Define Continuous State Bounds
bounds.phase(2).state.lower = [AltMin21, lonMin, latMin, vMin21, gammaMin21, zetaMin, aoaMin21, bankMin21, mFuelMin21];
bounds.phase(2).state.upper = [AltMax21, lonMax, latMax, vMax21, gammaMax21, zetaMax, aoaMax21, bankMax21, mFuelMax21];

% Define Initial State Bounds
bounds.phase(2).initialstate.lower = [AltMin21,lonMin, latMin, vMin21, gammaMin21, zetaMin, aoaMin21, 0, mFuelMin21] ; % initial bank angle at separation is 0
bounds.phase(2).initialstate.upper = [AltMax21,lonMax, latMax, vMax21, gammaMax21, zetaMax, aoaMax21, 0, mFuelMax21];

% Define End StateBounds
% End bounds are set slightly differently to the continuous states, even when not being constrained to a particular value, to encourage an optimal solution
bounds.phase(2).finalstate.lower = [AltMin21, lonMin, latMin, vMin21, 0, zetaMin, aoaMin21, 0, mFuelMin21]; % limit trajectory angle to minimum of 0 at SPARTAN-third stage separation
bounds.phase(2).finalstate.upper = [AltMax21, lonMax, latMax, vMax21, gammaMax21, zetaMax, aoaMax21, 0, mFuelMax21];
 

% Define Control Bounds
bounds.phase(2).control.lower = [deg2rad(-.5), deg2rad(-1)]; % Define the maximum angle of attack and bank angle change rates. These have been approximated to be physically realistic values. 
bounds.phase(2).control.upper = [deg2rad(.5), deg2rad(1)];

% Set Time Bounds
bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = 1000; % much longer time than it should fly for
bounds.phase(2).finaltime.lower = 100; % makes sure that time progresses forward
bounds.phase(2).finaltime.upper = 1000;

% Define Path Bounds (Path bounds values are calculated in Continuous.m)
bounds.phase(2).path.lower = [0];
bounds.phase(2).path.upper = Stage2.maxq;

if mode == 90 % if constant dynamic pressure flight is desired
bounds.phase(2).integral.lower = 0;
bounds.phase(2).integral.upper = 10000000;
end
%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 

% Define initial guess for start and end pt of each state. These have been
% chosen to give a very rough approximation of the beginning and end
% states of the solution. 
guess.phase(2).state(:,1)   = [24000;35000];
guess.phase(2).state(:,2)   = [lon0;lon0];
guess.phase(2).state(:,3)   = [lat0;lat0+0.2]; 
guess.phase(2).state(:,4)   = [1500; 2900];
guess.phase(2).state(:,5)   = [0.0; deg2rad(3)];
guess.phase(2).state(:,6)   = [pi;pi];
guess.phase(2).state(:,7)   = [2*pi/180; 5*pi/180];
guess.phase(2).state(:,8)   = [deg2rad(10);deg2rad(10)];
guess.phase(2).state(:,9) 	= [Stage2.mFuel; 100];

guess.phase(2).control      = [[0;0],[0;0]];
guess.phase(2).time          = [0;650];

if mode == 90
guess.phase(2).integral = 0
end

% Tie stages together
bounds.eventgroup(2).lower = [zeros(1,9)];
bounds.eventgroup(2).upper = [zeros(1,9)]; 

%% Flyback

disp('|                 *  Setting-Up Fly-Back states   *                 |')

AltMin22 = 10;  AltMax22 = MaxviscAlt21; % set maximum altitude bound to highest altitude provided in viscous database

vMin22 = 10;        vMax22 = 5000;

gammaMin22 = -80*pi/180;  gammaMax22 =  80*pi/180;

mFuelMin = 0; mFuelMax = 500;

bankMin21 = -90*pi/180; bankMax21 =   90*pi/180;  

throttleMin = 0; throttleMax = 1;

bounds.phase(3).initialtime.lower = 0;
bounds.phase(3).initialtime.upper = 3000;
bounds.phase(3).finaltime.lower = 400;
bounds.phase(3).finaltime.upper = 4000;

% Initial Bounds

 bounds.phase(3).initialstate.lower = [AltMin22, lonMin, latMin, vMin22, gammaMin22, zetaMin, aoaMin21, 0, mFuelMin, throttleMin];
bounds.phase(3).initialstate.upper = [AltMax22, lonMax, latMax, vMax22, gammaMax22, zetaMax, aoaMax21, 0, mFuelMax, throttleMax];   

% State Bounds
auxdata.States_ID22 = {'Altitude', 'Longitude', 'Latitude','Velocity','FPA','Heading','AoA', 'Bank Angle','Fuel Mass','Throttle'}; % Designate State Names

if returnMode == 0
bounds.phase(3).state.lower = [AltMin22, lonMin, latMin, vMin22, gammaMin22, zetaMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).state.upper = [AltMax22, lonMax, latMax, vMax22, gammaMax22, zetaMax, aoaMax21, bankMax21, 1, throttleMax];
else
bounds.phase(3).state.lower = [AltMin22, lonMin, latMin, vMin22, gammaMin22, zetaMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).state.upper = [AltMax22, lonMax, latMax, vMax22, gammaMax22, zetaMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
end
% End State Bounds
if returnMode == 0
bounds.phase(3).finalstate.lower = [AltMin22, lonMin, latMin, vMin22, gammaMin22, zetaMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).finalstate.upper = [AltMax22, lonMax, latMax, vMax22, gammaMax22, zetaMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
else
bounds.phase(3).finalstate.lower = [AltMin22, lon0, lat0, vMin22, gammaMin22, zetaMin, aoaMin21, bankMin21, mFuelMin, throttleMin];
bounds.phase(3).finalstate.upper = [1000, lon0, lat0, vMax22, gammaMax22, zetaMax, aoaMax21, bankMax21, mFuelMax, throttleMax];
end

% Control Bounds
if returnMode == 1 % these being different has no purpose but to keep my ascent results consistent
 bounds.phase(3).control.lower = [deg2rad(-.5), deg2rad(-1), -.1];
bounds.phase(3).control.upper = [deg2rad(.5), deg2rad(1), .1];   
else
bounds.phase(3).control.lower = [deg2rad(-.5), deg2rad(-1), -.2];
bounds.phase(3).control.upper = [deg2rad(.5), deg2rad(1), .2];
end

% Define Path Bounds (Path bounds values are calculated in Continuous.m)
bounds.phase(3).path.lower = 0;
bounds.phase(3).path.upper = 50000;

bounds.eventgroup(3).lower = [0]; 
bounds.eventgroup(3).upper = [0]; 

% Guess
tGuess              = [440; 1500];
altGuess            = [35000; 100];

if mode == 0
lonGuess            = [lon0; lon0-.1*pi/180];
latGuess            = [lat0-0.1;lat0];
else
lonGuess            = [lon0; lon0-.1*pi/180];
latGuess            = [-0.11;-0.10-0*pi/180];
end

speedGuess          = [3000; 10];
fpaGuess            = [0; 0];

if mode == 0
aziGuess            = [-pi; 0]; 
else
aziGuess            = [deg2rad(97); deg2rad(270)];
end

aoaGuess            = [8*pi/180; 5*pi/180];
bankGuess           = [89*pi/180; 89*pi/180];
mFuelGuess          = [100; mFuelMin];
if returnMode == 1
guess.phase(3).state   = [altGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[1.;1.]];
else
guess.phase(3).state   = [altGuess, lonGuess, latGuess, speedGuess, fpaGuess, aziGuess, aoaGuess, bankGuess, mFuelGuess,[0.;0.]];
end
guess.phase(3).control = [[0;0],[0;0],[0;0]];
guess.phase(3).time    = tGuess;


%% Third Stage
%% Third Stage Aerodynamic Data
disp('|            * Importing Third Stage Aerodynamic Data *             |')


% Import aerodynamic coefficients from file
auxdata.Aero3 = dlmread('ThirdStageAeroCoeffs.txt');

Aero3 = auxdata.Aero3;
auxdata.interp.Drag_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,5));
% 
auxdata.interp.Lift_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,6));

auxdata.interp.CP_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,7));

auxdata.interp.CN_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,2),Aero3(:,4));

auxdata.interp.Max_AoA_interp3 = scatteredInterpolant(Aero3(:,1),Aero3(:,4),Aero3(:,2));

%% Configure Third Stage States
disp('|                 * Setting-Up Third Stage States *                 |')

altMin3 = 30000;  altMax3 = 84000;
if mode == 0
phiMin3 = -1.5;         phiMax3 = -0.5;   
else
phiMin3 = -0.5;         phiMax3 = 0.5;
end
vMin3 = 10;        vMax3 = 8000;
gammaMin3=deg2rad(-5);  gammaMax3 =  deg2rad(30);

if mode == 0
zetaMin3 = deg2rad(-180); zetaMax3 =  deg2rad(0); 
else
zetaMin3 = deg2rad(80); zetaMax3 =  deg2rad(120);
end

aoaMin3 = 0;  aoaMax3= deg2rad(20);
aoadotMin3 = -deg2rad(1);
aoadotMax3 = deg2rad(1);

%-------------------------------------------------------------------%
%--------------- Set Up Problem Using Data Provided Above ----------%
%-------------------------------------------------------------------%
% Set Time Bounds
bounds.phase(4).initialtime.lower = 0;
bounds.phase(4).initialtime.upper = 10000;
bounds.phase(4).finaltime.lower = 1;
bounds.phase(4).finaltime.upper = 10000;

% Designate State Names
auxdata.States_ID3 = {'Altitude','Velocity','FPA', 'Mass', 'Aoa', 'Latitude','Heading'}; 

% Define Initial State Bounds
bounds.phase(4).initialstate.lower = [altMin3,vMin3, 0,  auxdata.Stage3.mTot, aoaMin3, phiMin3, zetaMin3];
bounds.phase(4).initialstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

% Define Continuous State Bounds
bounds.phase(4).state.lower = [altMin3,vMin3, gammaMin3, 0, aoaMin3, phiMin3, zetaMin3];
bounds.phase(4).state.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, aoaMax3, phiMax3, zetaMax3];

% Define Final State Bounds
bounds.phase(4).finalstate.lower = [altMin3, vMin3, 0, 0, 0, phiMin3, zetaMin3];
bounds.phase(4).finalstate.upper = [altMax3, vMax3, gammaMax3, auxdata.Stage3.mTot, 0, phiMax3, zetaMax3];

% Define Control Bounds
bounds.phase(4).control.lower = [aoadotMin3];
bounds.phase(4).control.upper = [aoadotMax3];

% Define Path Bounds (Path bounds values are calculated in Continuous.m)
bounds.phase(4).path.lower = [-deg2rad(8), -inf];
bounds.phase(4).path.upper = [deg2rad(8), 0];

% Define Event Bounds
bounds.eventgroup(4).lower = [zeros(1,7) 90000 0];
bounds.eventgroup(4).upper = [zeros(1,6) 1000 566000 0];

%% Third Stage Guess
tGuess              = [0; 150];
altGuess            = [35000; 60000];
vGuess          = [2700; 6000];
gammaGuess            = [0; deg2rad(10)];
mGuess              = [3300; 2000];
aoaGuess            = [deg2rad(20); deg2rad(20)];



if mode == 0
phiGuess = [lat0-0.1;lat0-0.1];
zetaGuess = [-deg2rad(97);-deg2rad(97)];  
else
phiGuess = [-0.11;-0.11];
zetaGuess = [deg2rad(97);deg2rad(97)];
end

guess.phase(4).state   = [altGuess, vGuess, gammaGuess, mGuess, aoaGuess, phiGuess, zetaGuess];
guess.phase(4).control = [0;0];
guess.phase(4).time    = tGuess;

%%
disp('|                 *     Configuring GPOPS-2       *                 |')
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre'; 
% mesh.method       = 'hp-LiuRao'; 
%  mesh.method       = 'hp-DarbyRao';
% mesh.maxiterations = 10;
if mode == 90
  mesh.maxiterations = 6;  
elseif returnMode == 1
 mesh.maxiterations = 4;   
else
mesh.maxiterations = 3;
end
mesh.colpointsmin = 8;
mesh.colpointsmax = 500;
mesh.tolerance    = 1e-5;

%--------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%

setup.name                           = 'SPARTAN-Combined';
setup.functions.continuous           = @CombinedContinuous;
setup.functions.endpoint             = @CombinedEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 0;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';

setup.nlp.ipoptoptions.maxiterations = 150;

setup.derivatives.supplier           = 'sparseFD';

setup.derivatives.derivativelevel    = 'first';
% setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.scales.method                  = 'automatic-guessUpdate';
setup.derivatives.dependencies      = 'full';



%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%

%% Parallel Loop
disp('|                * Setting-Up GPOPS-2 in Parallel *                 |')

num_it = 3; % Define number of parallel iterations

% Create variable setup structure
%% for interaction testing
if mode == 99
    modes_list = [2 3 4]; % define modes to generate interactions for
    runs1 = [-1 -1  zeros(1,length(modes_list))]; % defines the number of runs to perform. Need to define two -1 and two 1
    runs2 = [-1  1 zeros(1,length(modes_list))]; % defines the number of runs to perform. Need to define two -1 and two 1
    runs3 = [ 1 1 zeros(1,length(modes_list))]; % defines the number of runs to perform. Need to define two -1 and two 1
    
    % calculate all possible permutations where 2 of the given modes are
    % chosen
    k = length(modes_list);
    nk1 = nchoosek(runs1,k);
    p1=zeros(0,k);
    for i=1:size(nk1,1),
        pi1 = perms(nk1(i,:));
        p1 = unique([p1; pi1],'rows');
    end
    % calculate all possible permutations where 2 of the given modes are
    % chosen
    k = length(modes_list);
    nk2 = nchoosek(runs2,k);
    p2=zeros(0,k);
    for i=1:size(nk2,1),
        pi2 = perms(nk2(i,:));
        p2 = unique([p2; pi2],'rows');
    end
    % calculate all possible permutations where 2 of the given modes are
    % chosen
    k = length(modes_list);
    nk3 = nchoosek(runs3,k);
    p3=zeros(0,k);
    for i=1:size(nk3,1),
        pi3 = perms(nk3(i,:));
        p3 = unique([p3; pi3],'rows');
    end
    p = [p1;p2;p3];
   
     auxdata.p = p;
% Step through modes
    for j = 1:size(p,1)
        namelist{j} = num2str(j);
        setup_variations{j} = setup;
        if p(j,1) == 1
                
                setup_variations{j}.bounds.phase(1).path.upper = q_vars(end-1);
                setup_variations{j}.bounds.phase(2).path.upper = q_vars(end-1);
                setup_variations{j}.bounds.phase(3).path.upper = q_vars(end-1);
        elseif p(j,1) == -1
 
                setup_variations{j}.bounds.phase(1).path.upper = q_vars(2);
                setup_variations{j}.bounds.phase(2).path.upper = q_vars(2);
                setup_variations{j}.bounds.phase(3).path.upper = q_vars(2);
        end
        if p(j,2) == 1
          
                setup_variations{j}.auxdata.Ispmod = Isp_vars(end-1);
        elseif p(j,2) == -1
              
                setup_variations{j}.auxdata.Ispmod = Isp_vars(2);
        end
        if p(j,3) == 1
           
                setup_variations{j}.auxdata.Cdmod = Cd_vars(end-1);
        elseif p(j,3) == -1
             
                setup_variations{j}.auxdata.Cdmod = Cd_vars(2);
        end
%         if mode == 5 % not yet implemented
%                 setup_variations{j} = setup;
%                 setup_variations{j}.auxdata.vCdmod = vCd_vars(end);
%         end
%         if mode == 8 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.bounds.phase(4).initialstate.lower(4) = auxdata.Stage3.mTot*m3_vars(end);
%                 setup_variations{j}.bounds.phase(4).initialstate.upper(4) = auxdata.Stage3.mTot*m3_vars(end);
%                 setup_variations{j}.bounds.phase(4).state.upper(4) = auxdata.Stage3.mTot*m3_vars(end);
%                 setup_variations{j}.bounds.phase(4).finalstate.upper(4) = auxdata.Stage3.mTot*m3_vars(end);
%         end
%         if mode == 9 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.auxdata.Isp3mod = Isp3_vars(end);
%         end
%         if mode == 10 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.auxdata.Stage2.mStruct = auxdata.Stage2.mStruct*mSPARTAN_vars(end);
%         end
%         if mode == 11 %
%                 setup_variations{j} = setup;
%                 setup_variations{j}.bounds.phase(2).initialstate.lower(9) = Stage2.Initial.mFuel*mFuel_vars(end);
%                 setup_variations{j}.bounds.phase(2).initialstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(end);
%                 setup_variations{j}.bounds.phase(2).state.upper(9) = Stage2.Initial.mFuel*mFuel_vars(end);
%                 setup_variations{j}.bounds.phase(2).finalstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(end);
% 
%         end

    end
end

%% for independent variable testing
if mode == 1 || mode == 0
    setup_variations{1} = setup;
elseif mode == 90
    setup_variations{1} = setup;
elseif mode == 2
    for i = 1:length(q_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.bounds.phase(1).path.upper = q_vars(i);
        setup_variations{i}.bounds.phase(2).path.upper = q_vars(i);
        setup_variations{i}.bounds.phase(3).path.upper = q_vars(i);
    end
%     setup_variations{2}.guess.phase(2).state(:,4) = [1500;2800]; % modify guess for badly converging solution
elseif mode == 3
    for i = 1:length(Isp_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Ispmod = Isp_vars(i);
    end
%     setup_variations{2}.guess.phase(2).state(:,4) = [1500;2800]; % modify guess for badly converging solution
%     if returnMode == 1
%         setup_variations{1}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%         setup_variations{2}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%         setup_variations{4}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%         setup_variations{5}.guess.phase(2).state(:,4) = [1500;2700]; % modify guess for badly converging solution
%     end
elseif mode == 4
    for i = 1:length(Cd_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Cdmod = Cd_vars(i);
    end
    
elseif mode == 5 % not yet implemented
    for i = 1:length(vCd_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.vCdmod = vCd_vars(i);
    end
    elseif mode == 6 % not yet implemented
    for i = 1:length(FirstStagem_vars)  
        setup_variations{i} = setup;
        
%         setup_variations{i}.bounds.phase(1).initialstate.upper(3) = mMax1+mRocket*(FirstStagem_vars(i)-1);
%         setup_variations{i}.bounds.phase(1).state.upper(3) = mMax1+mRocket*(FirstStagem_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = mF + mEmpty*(FirstStagem_vars(i)-1);
%         setup_variations{i}.bounds.phase(1).finalstate.upper(3) = mMax1+mRocket*(FirstStagem_vars(i)-1);
    end
    elseif mode == 7 %
    for i = 1:length(Cd3_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Cd3mod = Cd3_vars(i);
    end
elseif mode == 8 %
    for i = 1:length(m3_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.bounds.phase(4).initialstate.lower(4) = auxdata.Stage3.mTot*m3_vars(i);
        setup_variations{i}.bounds.phase(4).initialstate.upper(4) = auxdata.Stage3.mTot*m3_vars(i);
        setup_variations{i}.bounds.phase(4).state.upper(4) = auxdata.Stage3.mTot*m3_vars(i);
        setup_variations{i}.bounds.phase(4).finalstate.upper(4) = auxdata.Stage3.mTot*m3_vars(i);
        
        setup_variations{i}.bounds.phase(1).initialstate.lower(3) = setup.bounds.phase(1).initialstate.lower(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).initialstate.upper(3) = setup.bounds.phase(1).initialstate.upper(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.lower(3) = setup.bounds.phase(1).state.lower(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.upper(3) = setup.bounds.phase(1).state.upper(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = setup.bounds.phase(1).finalstate.lower(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.upper(3) = setup.bounds.phase(1).finalstate.upper(3) + auxdata.Stage3.mTot*(m3_vars(i)-1);
    
        setup_variations{i}.auxdata.Stage3.mTot = auxdata.Stage3.mTot*m3_vars(i);
    end
elseif mode == 9 %
    for i = 1:length(Isp3_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Isp3mod = Isp3_vars(i);
    end
elseif mode == 10 %
    for i = 1:length(mSPARTAN_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.auxdata.Stage2.mStruct = auxdata.Stage2.mStruct*mSPARTAN_vars(i);
        setup_variations{i}.bounds.phase(1).initialstate.lower(3) = setup.bounds.phase(1).initialstate.lower(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).initialstate.upper(3) = setup.bounds.phase(1).initialstate.upper(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.lower(3) = setup.bounds.phase(1).state.lower(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.upper(3) = setup.bounds.phase(1).state.upper(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = setup.bounds.phase(1).finalstate.lower(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.upper(3) = setup.bounds.phase(1).finalstate.upper(3) + auxdata.Stage2.mStruct*(mSPARTAN_vars(i)-1);
    
    
    end
elseif mode == 11 %
    for i = 1:length(mFuel_vars)  
        setup_variations{i} = setup;
        setup_variations{i}.bounds.phase(2).initialstate.lower(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        setup_variations{i}.bounds.phase(2).initialstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        setup_variations{i}.bounds.phase(2).state.upper(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        setup_variations{i}.bounds.phase(2).finalstate.upper(9) = Stage2.Initial.mFuel*mFuel_vars(i);
        
        setup_variations{i}.bounds.phase(1).initialstate.lower(3) = setup.bounds.phase(1).initialstate.lower(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).initialstate.upper(3) = setup.bounds.phase(1).initialstate.upper(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.lower(3) = setup.bounds.phase(1).state.lower(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).state.upper(3) = setup.bounds.phase(1).state.upper(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.lower(3) = setup.bounds.phase(1).finalstate.lower(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
        setup_variations{i}.bounds.phase(1).finalstate.upper(3) = setup.bounds.phase(1).finalstate.upper(3) + Stage2.Initial.mFuel*(mFuel_vars(i)-1);
    
    end
end



for j = 1:length(setup_variations)
disp(['Starting Setup Variation ',num2str(j)])
    
for i = 1:num_it
setup_par(i) = setup_variations{j};
% setup_par(i).nlp.ipoptoptions.maxiterations = 500 + 10*i;
setup_par(i).guess.phase(2).state(:,1)   = [24000; 25000 + 1000*i]; % vary altitude guess
end

parfor i = 1:num_it
    
disp(strcat('|                * Initialising GPOPS-2 on Core_',' ', num2str(i),' *                 |'))
disp(strcat('|                Working...                                                           |'))


output_temp = gpops2(setup_par(i)); % Run GPOPS-2. Use setup for each parallel iteration.


disp(strcat('|              * Completed Run of GPOPS-2 on Core_',' ', num2str(i),' *               |'))


% Compute State Error
input_test = output_temp.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);
norm_error1 = [];
for num = [1:6 9]
Stage2_int = cumtrapz([output_temp.result.solution.phase(2).time(1):0.1:output_temp.result.solution.phase(2).time(end)],pchip(output_temp.result.solution.phase(2).time,phaseout_test(2).dynamics(1:end,num),[output_temp.result.solution.phase(2).time(1):0.1:output_temp.result.solution.phase(2).time(end)]))';
norm_error1(num,:) = (interp1([output_temp.result.solution.phase(2).time(1):0.1:output_temp.result.solution.phase(2).time(end)],Stage2_int,output_temp.result.solution.phase(2).time)+ output_temp.result.solution.phase(2).state(1,num)- output_temp.result.solution.phase(2).state(:,num))./(max(output_temp.result.solution.phase(2).state(:,num))-min(output_temp.result.solution.phase(2).state(:,num)));

end
norm_error2 = [];
if returnMode == 1
for num = [1:6 9]
Return_int = cumtrapz([output_temp.result.solution.phase(3).time(1):0.1:output_temp.result.solution.phase(3).time(end)],pchip(output_temp.result.solution.phase(3).time,phaseout_test(3).dynamics(1:end,num),[output_temp.result.solution.phase(3).time(1):0.1:output_temp.result.solution.phase(3).time(end)]))';
norm_error2(num,:) = (interp1([output_temp.result.solution.phase(3).time(1):0.1:output_temp.result.solution.phase(3).time(end)],Return_int,output_temp.result.solution.phase(3).time)+ output_temp.result.solution.phase(3).state(1,num)- output_temp.result.solution.phase(3).state(:,num))./(max(output_temp.result.solution.phase(3).state(:,num))-min(output_temp.result.solution.phase(3).state(:,num)));


end
else
 norm_error2 = 0;  
end
error(i) = max([max(abs(norm_error1)) max(abs(norm_error2))]);


PayloadMass(i) = -output_temp.result.objective;

output_store{i} = output_temp;

end

[min_error,index] = min(error); % Calculate the result which minimises the chosen error function

% [max_pl,index] = max(PayloadMass);% Calculate the result which maximises payload mass the chosen error function

disp('|               * Determining Best GPOPS-2 Solution *               |')


output{j} = output_store{index};

% Clear output store to save memory and prevent write issues
if mode ~= 1
clear output_store
end

end

%% Process Results










Plotter(output,auxdata,auxdata.mode,auxdata.returnMode,auxdata.namelist,M_englist,T_englist,engine_data,MList_EngineOn,AOAList_EngineOn,Stage1.m,Stage2.mStruct+Stage2.mFuel+Stage3.mTot ,Stage2.mFuel,h0,v0,bounds);







disp('---------------------------------------------------------------------')
disp('|                 *     LODESTAR Run Complete     *                 |')
disp('---------------------------------------------------------------------')

%% =========================================================================
% Troubleshooting Procedure
% =========================================================================
% 1: Check that you have posed your problem correctly ie. it is physically
% feasible and the bounds allow for a solution.
% 2: Check if there is any extrapolations causing bad dynamics. 
% 3: Check guess, is it reasonable?.
% 4: Increase number of iterations and collocation points.
% 5: Check for large nonlinearities, eg. atmospheric properties suddenly going to zero or thrust cutting off. These need to be eliminated or separated into phases.
% 6: Check for NaN values (check derivatives in Dynamics file while
% running).




