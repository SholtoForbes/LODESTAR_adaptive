function [auxdata, MaxviscAlt21, MaxviscAlt22] = SPARTANint(auxdata)
% This function calculates aerodynamic and propulsion interpolation splined
% for the SPARTAN


%% Conical Shock Data %%===================================================
disp('|                 * Importing Conical Shock Data  *                 |')
% Import conical shock data which defines the conditions entering the
% SPARTAN's scramjet engines. 

% Create interpolation splines for M1 (Mach
% after shock), P_rat (pressure ratio over shock) and T rat (temperature
% ratio over shock)

% Import data from file
shockdata = dlmread('ShockMat');

%change Mach and AoA data to ndgrid format
[MList_EngineOn,AOAList_EngineOn] = ndgrid(unique(shockdata(:,1)),unique(shockdata(:,2)));

% Reshape M1, P_rat and T_rat data into ndgrid form
M1_Grid = reshape(shockdata(:,3),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
Prat_Grid = reshape(shockdata(:,4),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';
Trat_Grid = reshape(shockdata(:,5),[length(unique(shockdata(:,2))),length(unique(shockdata(:,1)))]).';

% Create interpolation splines of gridded data
auxdata.interp.M1gridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,M1_Grid,'spline','linear'); 
auxdata.interp.presgridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Prat_Grid,'spline','linear');
auxdata.interp.tempgridded = griddedInterpolant(MList_EngineOn,AOAList_EngineOn,Trat_Grid,'spline','linear');

% These interpolation splines are utilised within RESTint.m, when called by
% SPARTANDynamics.m from Continuous.m

%% Equivalence Ratio %%==========================================================
disp('|                * Importing Scramjet Engine Data *                 |')


% Import engine data
auxdata.interp.engine_data = dlmread('ENGINEDATASORTED.txt');  % reads four columns; Mach no after conical shock, temp after conical shock, Isp, max equivalence ratio
engine_data = auxdata.interp.engine_data;

% Create uniform grid of Mach no. and temperature values. 
M_englist = unique(sort(engine_data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(engine_data(:,1)));

T_englist = unique(sort(engine_data(:,2))); % create unique list of angle of attack numbers from engine data
T_eng_interp = unique(sort(engine_data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);

% Set the equivalence ratio interpolation region %-------------------------
% VERY IMPORTANT

% The interpolators have trouble with equivalence ratio because its equal
% to 1 over a certain Mach no. (causes error in interpolator, as the
% interpolator will find values of equivalence ratio < 1 where they should
% not exist)

% This makes anything outside of the region where it is actually changing
% extrapolate to over 1 (which is then set to 1 by RESTM12int)

% the the maximum of this to around where equivalence ratio stops changing,
% and check the end results

eq_data = [];
j=1;
for i = 1: length(engine_data(:,1))
    if engine_data(i,1) < 5.
        eq_data(j,:) = engine_data(i,:);
        j=j+1;
    end
end

auxdata.interp.equivalence = scatteredInterpolant(eq_data(:,1),eq_data(:,2),eq_data(:,4), 'linear');
grid.eq_eng = auxdata.interp.equivalence(grid.Mgrid_eng,grid.T_eng);
auxdata.interp.eqGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,grid.eq_eng,'linear','linear');

%% Isp  %------------------------------------------------------------------

disp('|                 * Performing Isp Interpolation  *                 |')

% The engine data set is arranged to as to be difficult to interpolate
% effectively, with scattered Mach, temp and pressure data points. 
% Small imperfections in the interpolation scheme can have
% large effects on the optimal trajectory. To account for this, each 'line'
% of the C-rest data set is normalised so that the whole data set forms a
% regular 6x5 square, with even coordinate spacing. This is interpolated
% between using cubic interpolation. When this is interpolated for using a
% physical query point, the position of the data point in normalised
% coordinates is found by first finding the position of the query point
% relative to the existing data set. This position is converted to normalised 
% coordinates using the relative distance to the closest four data points. 
% See appendix of Ph.D. thesis by Sholto Forbes-Spyratos for an illustration. 


% Create coordinate system normalised to each 'line' in the engine data set
coords_noextrap = [];
temp = 1;

% create coordinates list for existing data set
for i = 2:7 % need to preallocate with extrapolation in mind
    for j = 2:6
        coords_noextrap(temp,1) = i;
        coords_noextrap(temp,2) = j;
        temp = temp+1;
        
        
    end
end

scattered_Isp_norm = scatteredInterpolant(coords_noextrap(:,1),coords_noextrap(:,2),engine_data(:,3));
scattered_M_norm = scatteredInterpolant(coords_noextrap(:,1),coords_noextrap(:,2),engine_data(:,1));
scattered_T_norm = scatteredInterpolant(coords_noextrap(:,1),coords_noextrap(:,2),engine_data(:,2));



% extrapolate in normalised coordinates and create grid
coords = [];
temp = 1;

for i = 1:8 % extrapolate by 1 'line'
    for j = 1:9
        coords(temp,1) = i;
        coords(temp,2) = j;
%         grid.Ispnorm(i,j) = engine_data(temp,3); % form an Isp grid at each coordinate
        grid.Ispnorm(i,j) = scattered_Isp_norm(i,j); % form an Isp grid at each coordinate
        M_List_data(temp) = scattered_M_norm(i,j);
        T_List_data(temp) = scattered_T_norm(i,j);
        temp = temp+1;

    end
end

M_List_data = M_List_data';
T_List_data = T_List_data';


[grid.normx,grid.normy] =  ndgrid(unique(coords(:,1)),unique(coords(:,2))); % form grids of normalised coordinates


% Create cubic spline on normalised coordinates
auxdata.interp.Ispnorm = griddedInterpolant(grid.normx,grid.normy,grid.Ispnorm,'cubic','linear'); % Cubic used as it is the best representation of the data

% Create the physical space over which to interpolate
M_eng_interp2 = min(engine_data(:,1)):0.1:max(engine_data(:,1));
% T_eng_interp2 = min(engine_data(:,2)):1:max(engine_data(:,2));
T_eng_interp2 = 220:1:max(engine_data(:,2));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp2,T_eng_interp2);

size_grid = size(grid.Mgrid_eng);

Isp_grid= [];
M_list=[];
T_list=[];
Isp_norm_list=[];
temp = 1;

% Interpolate for M and T grids, this is done external of the main vehicle
% simulation for efficiency
for i = 1:size_grid(1)
    for j = 1:size_grid(2)
M_temp = grid.Mgrid_eng(i,j);
T_temp = grid.T_eng(i,j);

in_dex = [0 0 0 0];

% Determine if the qeury point is within any four points of the data set,
% and if so, which four
for polyindex_x = 1:7
    for polyindex_y = 1:8
        % Create polygon fof four data points
        cells_index = [polyindex_y+9*(polyindex_x-1)  polyindex_y+1+9*(polyindex_x-1)  polyindex_y+1+9*(polyindex_x) polyindex_y+9*(polyindex_x)];
        % Search
        in = inpolygon(M_temp,T_temp,M_List_data(cells_index),T_List_data(cells_index));
        if in
            in_dex = cells_index; % record cells index for polygon if query point inside
        end
    end
end

if in_dex == [0 0 0 0]
Isp_grid(i,j) = 0; % if not inside data set set point to 0

else
    % if inside a polygon of data points
    % Define a polygon of the data points
pt0 = [M_temp T_temp 0]; %query pt
pt1 = [M_List_data(in_dex(1)) T_List_data(in_dex(1)) 0]; %data pt 1
pt2 = [M_List_data(in_dex(2)) T_List_data(in_dex(2))  0];
pt3 = [M_List_data(in_dex(3)) T_List_data(in_dex(3))  0];
pt4 = [M_List_data(in_dex(4)) T_List_data(in_dex(4)) 0];

% Calculate minimum distances from query pt to a line between each 'side'
% of data pts
      a1 = pt1 - pt2;
      b1 = pt0 - pt2;
      d1 = norm(cross(a1,b1)) / norm(a1);

      a2 = pt3 - pt4;
      b2 = pt0 - pt4;
      d2 = norm(cross(a2,b2)) / norm(a2);
      
      a3 = pt1 - pt4;
      b3 = pt0 - pt4;
      d3 = norm(cross(a3,b3)) / norm(a3);
      
      a4 = pt2 - pt3;
      b4 = pt0 - pt3;
      d4 = norm(cross(a4,b4)) / norm(a4);
      
      % Calculate normalised coordinate as the distance fraction from the southwest data pt to the query pt in the x
      % and y directions, added to the coordinate of the southwest data pt
x_norm(i,j) = coords(in_dex(1),1) + d1/(d1+d2);
y_norm(i,j) = coords(in_dex(1),2) + d3/(d3+d4);  

Isp_grid(i,j) = auxdata.interp.Ispnorm(x_norm(i,j),y_norm(i,j)); % interpolate for Isp and create grid corresponding to M and T

% Define lists for extrapolation purposes
M_list(temp) = M_temp;
T_list(temp) = T_temp;
Isp_norm_list(temp) = Isp_grid(i,j);
temp = temp+1;
end

    end
end

% Use the interpolated data set to generate a scattered interpolant for
% extrapolation
auxdata.interp.Ispscatterednorm = scatteredInterpolant(M_list',T_list',Isp_norm_list','linear','nearest');


% % Extrapolate on any point outside of data set
for i = 1:size_grid(1)
    for j = 1:size_grid(2)
        if Isp_grid(i,j) == 0
            
            Isp_grid(i,j) = auxdata.interp.Ispscatterednorm(grid.Mgrid_eng(i,j),grid.T_eng(i,j));
        end
    end
end
% contourf(grid.Mgrid_eng,grid.T_eng,Isp_grid,100)
% Create spline for use in vehicle model.
auxdata.interp.IspGridded = griddedInterpolant(grid.Mgrid_eng,grid.T_eng,Isp_grid,'spline','spline');


%% Aerodynamic Data
disp('|           * Importing Scramjet Stage Aerodynamic Data *           |')

% Fetch aerodynamic data and compute interpolation splines.
% Calculate the flap deflection necessary for trim.
% Each set of aero corresponds to a different CG. 

Viscousaero_EngineOn = importdata('VC3D_viscousCoefficients_ascent.dat');
Viscousaero_EngineOff = importdata('VC3D_viscousCoefficients_descent.dat');

MaxviscAlt21 = max(Viscousaero_EngineOn(:,3)); % set maximum altitude to the max provided by viscous database
MaxviscAlt22 = max(Viscousaero_EngineOff(:,3));

% These aerodynamic datasets have been created in ClicCalcCGVar.m

T_L = -1.327; % Thrust location in z, average (m), measured from CREO

% Full of fuel, with third stage

CG_z = (-0.1974*(4.9571e+03+1562) + 3300*0.547)/(4.9571e+03+1562+3300); % Centre of gravity location in z, calculated from Creo


disp('|    * Creating and Interpolating Trimmed Aedoynamic Databases *    |')


aero1.aero_EngineOff = importdata('SPARTANaero15.228');
aero1.flapaero = importdata('SPARTANaeroFlaps15.228');
aero1.aero_EngineOn = importdata('SPARTANaeroEngineOn15.228');
aero1.aero_Engine = importdata('SPARTANEngine15.228');
aero1.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero1.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.fullFuel,auxdata.interp.Cd_spline_EngineOff.fullFuel,...
    auxdata.interp.Cl_spline_EngineOn.fullFuel,auxdata.interp.Cd_spline_EngineOn.fullFuel,...
    auxdata.interp.flap_spline_EngineOff.fullFuel,auxdata.interp.flap_spline_EngineOn.fullFuel,...
    auxdata.T_spline_Rear,auxdata.Fd_spline_NoEngine,auxdata.Cd_spline_ViscousEngineOff,auxdata.Cd_spline_ViscousEngineOn,...
    auxdata.L_spline_Rear,auxdata.T_spline] = AeroInt(aero1,auxdata,T_L,CG_z);


% Cylindrical fuel tanks depleted, with third stage
CG_z = (-0.1974*(4.9571e+03+852) + 3300*0.547)/(4.9571e+03+852+3300);

aero2.aero_EngineOff = importdata('SPARTANaero15.1557');
aero2.flapaero = importdata('SPARTANaeroFlaps15.1557');
aero2.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1557');
aero2.aero_Engine = importdata('SPARTANEngine15.1557');
aero2.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero2.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.cylTankEnd,auxdata.interp.Cd_spline_EngineOff.cylTankEnd,auxdata.interp.Cl_spline_EngineOn.cylTankEnd,auxdata.interp.Cd_spline_EngineOn.cylTankEnd,auxdata.interp.flap_spline_EngineOff.cylTankEnd,auxdata.interp.flap_spline_EngineOn.cylTankEnd] = AeroInt(aero2,auxdata,T_L,CG_z);

% Empty, with third stage. 

CG_z = (-0.2134*4.9571e+03+ 3300*0.547)/(4.9571e+03+3300);

aero3.aero_EngineOff = importdata('SPARTANaero15.727');
aero3.flapaero = importdata('SPARTANaeroFlaps15.727');
aero3.aero_EngineOn = importdata('SPARTANaeroEngineOn15.727');
aero3.aero_Engine = importdata('SPARTANEngine15.727');
aero3.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero3.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noFuel,auxdata.interp.Cd_spline_EngineOff.noFuel,auxdata.interp.Cl_spline_EngineOn.noFuel,auxdata.interp.Cd_spline_EngineOn.noFuel,auxdata.interp.flap_spline_EngineOff.noFuel,auxdata.interp.flap_spline_EngineOn.noFuel] = AeroInt(aero3,auxdata,T_L,CG_z);

% Flyback, empty without third stage, empty.

CG_z = -0.2134; % calculated fom CREO

aero4.aero_EngineOff = importdata('SPARTANaero15.1255');
aero4.flapaero = importdata('SPARTANaeroFlaps15.1255');
aero4.aero_EngineOn = importdata('SPARTANaeroEngineOn15.1255');
aero4.aero_Engine = importdata('SPARTANEngine15.1255');
aero4.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero4.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStageEmpty,...
    auxdata.interp.Cd_spline_EngineOff.noThirdStageEmpty,auxdata.interp.Cl_spline_EngineOn.noThirdStageEmpty,...
    auxdata.interp.Cd_spline_EngineOn.noThirdStageEmpty,auxdata.interp.flap_spline_EngineOff.noThirdStageEmpty,...
    auxdata.interp.flap_spline_EngineOn.noThirdStageEmpty] = AeroInt(aero4,auxdata,T_L,CG_z);

% Flyback, empty without third stage, conical fuel tank full.
CG_z = -0.1974; % calculated fom CREO

aero5.aero_EngineOff = importdata('SPARTANaero14.297');
aero5.flapaero = importdata('SPARTANaeroFlaps14.297');
aero5.aero_EngineOn = importdata('SPARTANaeroEngineOn14.297');
aero5.aero_Engine = importdata('SPARTANEngine14.297');
aero5.Viscousaero_EngineOff = Viscousaero_EngineOff;
aero5.Viscousaero_EngineOn = Viscousaero_EngineOn;

[auxdata.interp.Cl_spline_EngineOff.noThirdStagecylTankEnd,...
    auxdata.interp.Cd_spline_EngineOff.noThirdStagecylTankEnd,auxdata.interp.Cl_spline_EngineOn.noThirdStagecylTankEnd,...
    auxdata.interp.Cd_spline_EngineOn.noThirdStagecylTankEnd,auxdata.interp.flap_spline_EngineOff.noThirdStagecylTankEnd,...
    auxdata.interp.flap_spline_EngineOn.noThirdStagecylTankEnd] = AeroInt(aero5,auxdata,T_L,CG_z);

clear aero1
clear aero2
clear aero3
clear aero4
clear aero5
clear Viscousaero_EngineOn
clear Viscousaero_EngineOff

end

