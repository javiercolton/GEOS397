%% Homework 5: Seafloor Subsisdence due to cooling
% GEOS 397
% Clay Roehner and Javier Colton
%
%
close all
clc

%% Part 1 Conductive heat flow
%% Step 1: A model
% Imagine an infinitely long and wide solid plate. The plate has thickness d. The temperature at the top of
% the plate is T1 and the temperature at the bottom of the plate is T2. Draw a diagram of this plate and label
% these parameters.
figure (1)
rectangle('Position',[1 3 8 3]);
axis ([0 10 0 10]);
axis off;
title('Diagram of a Solid Plate');

dThick = '\leftarrow d';
tempT1 = 'T_{1}';
tempT2 = 'T_{2}';

text(1, 4, dThick);
text(9, 6, tempT1);
text(9, 3, tempT2);

dim = [.2 .5 .3 .3];
str = {'T_{1} = Surface Temp','T_{2} = Temp at Depth (d)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% Step 2: Heat flow
%
% The rate of heat flow per unit area ($\it Q$) $\it down$ through the
% plate is
%
% $\it Q = \it -k * (T_{2}-T_{1}/\it d)$
%
% Why does heat flow $\it down$ in this equation?
%
% Heat flows down because of the negative nature of $\it k$. Q is a $\it
% rate$. At the surface a greater $\it rate$ occurs because of a greater
% temperature gradient than at depth where the gradient is less (everything
% is hotter). 
%
%% Step 3: Thermal conductivities
%
% Thermal Conductivity Values (in $W/m \bullet C$
%
% $\bullet$ Silver: 429
%
% $\bullet$ Magnesium: 156
%
% $\bullet$ Glass: 1.05
%
% $\bullet$ Rock: 2-7
%
% $\bullet$ Wood (oak): 0.17 


%% Step 4: The heat transport equation
%
% Substitute our values into heat transport equation:
%
% $Q = -k * ((T(z+\delta z)- t(z))/(Z+\delta z))$
%
% Write the derivative of the right hand side of the equation:
%
% $-k (\delta T / \delta z)$
%
% Heat flow equation using erivative
%
% $Q = -k (\delta T / \delta z)$
%
%Write the same thing using the gradient operator?
%
%
%% Step 5: The conservation equation
%
% Insert the calculated Q (above) into the conservation equation
%
% $c_{p}  \rho  (\delta T/\delta t) = A - (\delta/\delta z)(-k(\delta T/\delta z))$
%
% Assuming that internal heat generation is zero we get:
%
% $c_{p}  \rho  (\delta T/\delta t) = -k(\delta/\delta z)(\delta T/\delta z))$
%
% By dividing $c_{p}$ and $\rho$ through we get $\kappa = -K/c_{p}*\rho$ and
% finally:
%
% $(\delta T/\delta t) = k(\delta^{2} T/\delta z^{2}))$
%
% What is the value of K in this case?
%
% $\kappa = -K/c_{p}*\rho$
%% Part 2 : Oceanic lithosphere cooling
% Step 1: Setup the model domain and compute
%
%
depthOcean = linspace(0 , 100, 100); % [km]
timeInt = linspace(0 , 100, 100); %[Ma]
constantK = 31.56; %km^2/Ma
morT = 640;
%
T = zeros(length(depthOcean),length(timeInt));

tic
for it =  1 : length(timeInt);
    for iz = 1 : length(depthOcean);
       T(iz,it) =  morT * erf(depthOcean(iz)/(2*(sqrt(constantK*timeInt(it)))));
    end
end
toc
%

figure (2)

imagesc(timeInt, depthOcean, T); axis ij;
c=colorbar;
colormap(jet);
xlabel('Time [Ma]');
ylabel ('Ocean depth [km]');
c.Label.String = ('Temperature^{o}C');
title('Temperature using two "for" loops');
tic
[t1 z1] = meshgrid (timeInt, depthOcean);
meshT = morT * erf((z1)./(2*(sqrt(constantK*t1))));
toc

figure (3)

imagesc(timeInt, depthOcean, meshT)
c=colorbar;
colormap;
xlabel('Time [Ma]');
ylabel ('Ocean depth [km]');
c.Label.String = ('Temperature^{o}C');
title('Temperature using meshgrid');

% EXTRA CREDIT: 
% For loop elapsed time = 0.012884 seconds
% meshgrid elapsed time = 0.002297 seconds
%
% meshgrid is .0106 seconds faster than the for loop

%% Step 2: Analyze the model output
% Does your model make sense given the boundary conditions used to derive the solution?
% \bullet Yes the model makes sense given the boundary condition supplied.
% At T_{0} the original temp at the mid-ocean ridge is 640, as temperature
% decreases with depth, as we move further from the MOR our depth value
% increases and temperature decreases, as it does using this model.
%
% What controls the rate at which the temperature decays? (List all things you can think of.)
% \bullet The rate of temperature decay is controlled by the
% composition/heat conductivity of the rock, the temperature at the ocean
% floor, and the original heat of the material controlling the temperature
% gradient between the rock being cooled and the ocean temperature.
%
% How could we convert this model from age to distance from ridge axis?
% \bullet We could convert from age to distance from ridege axi by
% determining the relationship between age and distance (using velocity and
% time) and substitute the appropriate equation in for age related to
% distance.
% What would be a more appropriate boundary condition at T(z = 0) given what we know about the
% oceans?
% ???????????\bullet A more appropriate boundary condition T(z=0) would be a pixel by
% pixel elevation for what z actually is over space???????
% Is 640?C an appropriate value for the temperature at a mid-ocean ridge? Why or why not?
% \bullet 640^{o}C might be an appropriate temperature for certin spots along the
% mid-ocean ridge however this value can vary by up to 250^{o}C. Using one
% value for the entire length of the MOR simplifies the problem.



%% Part 3: Plate velocity and the depths of oceans
% Step 1: Load and plot sea-floor depth

load('spreadingData.mat');

fields(Bath)


subplot (2,1,1); 

plot (Bath.atlanticx,Bath.atlanticz);
ylabel ('Depth (m)');

title ('Atlantic Ocean');
grid on;
hold on;
subplot (2,1,2);
plot (Bath.pacificx,Bath.pacificz);
xlabel ('Distance from MOR (km)');
ylabel('Depth [m]');
title ('Pacific Ocean');
grid on;

%% Step 2: A half-space model


depthPredictPac = zeros(length(Bath.pacificx),1); % Create an empty vector for predicted Pacific depths
Velocity = 45; %km/Ma
for ix =  1 : length(Bath.pacificx); % Run a loop to populate depthPredictPac using the depth half-space model
       depthPredictPac(ix) = -(2.65 + .345 *(Bath.pacificx(ix)/Velocity)^(1/2)); % run Bath.pacific through the loop
end

depthPredictAtl = zeros(length(Bath.atlanticx),1); % Create an empty vector for predicted Pacific depths
VelocityAtl = 20;


for ix =  1 : length(Bath.atlanticx); % Run a loop to populate depthPredictAtl using the depth half-space model
       depthPredictAtl(ix) = -(3.12 + .345 *(Bath.atlanticx(ix)/VelocityAtl)^(1/2)); % run Bath.atlanticx through the loop
end

figure(3)


subplot (2,1,1);

plot (Bath.atlanticx,Bath.atlanticz,'b');
ylabel ('Depth (m)');

title ('Atlantic Ocean');
grid on;
hold on;
plot(Bath.atlanticx,(depthPredictAtl*1000), 'r');

subplot (2,1,2);
plot (Bath.pacificx,Bath.pacificz,'b');
hold on;
xlabel ('Distance from MOR (km)');
ylabel ('Depth (m)');
title ('Pacific Ocean');
grid on;
plot(Bath.pacificx,(depthPredictPac*1000), 'r');

legend on;
legend('Actual','Predicted' ,'Location','Best');

% What does 2.65 represent in the equation above?
% 2.65 representes the deoth of the MOR in km's. I adjusted the value to
% better fit the actual depths in the Atlantic Ocean.
%
% What are the velocity rates that best fit each ocean?
%
% Pacific Ocean = 45 km/Ma
% Atlantic Ocean = 20 km/Ma
%
% Convert to cm/year and compare to literature

velPac = Velocity*1000/1000000; % = 0.045 cm/yr

velAtl = VelocityAtl*1000/1000000; % = 0.02 cm/yr
% Measured velocities
% The velocities vary throughout both ocean but in general the PAcific is
% spreading from 40-100 mm/yr while the atlantic is spreading slower from
% 10-20 mm/yr (NOAA, 2008) . Our model shows a similar trend.



%% Part 4: Global oceanic plate ages
%
% Step 1: Load topo data and plot seafloor depths
 
load('topo.mat');

load('coastlines');

%[platelat, platelon] = importPlates('All_boundaries.txt');
%% Step 2: Kill the topography and convert to km
topo (topo>0) = 0;  % kill the topography
topo = topo/1000;   % convert depth to km



g = figure;
hold on;
g.InvertHardcopy = 'off'; %Figure background color set to black when printing
g.Color = 'k';  %Set figure window background color to black
g.Position = [100 100 1000 500]; %The location and size of figure's drawable area set to this extent
g.PaperPositionMode = 'auto'; %Printed figure matches the displayed figure size
%
ax = axesm('Mollweid', 'Frame', 'on', 'Grid', 'on'); %Defined the projection, turned frame and grid on
setm(ax,'MLabelLocation',60); %Set meridian label locations to every 60 degrees
setm(ax,'PLabelLocation',30); %Set parallel label locations to every 30 degrees
mlabel('MLabelParallel',0);  %Sets the location of the parallel labels
plabel('PLabelMeridian',-25);  %Sets the location of the meridian labels
axis('off');  %Turns off axis labeling, tickmarks and background
setm(ax,'FontColor',[0.9 0.9 0.9]); %Sets font color to nearly white
setm(ax,'GColor',[0.9 0.9 0.9]);
setm(ax,'fontweight','bold','fontsize',16);
% 
plotm(coastlat, coastlon);
LAT = topolatlim(1):topolatlim(2);
LON =  topolonlim(1):topolonlim(2);
%
idx180 = find(LON > 180, 1); % find first index of longitude array that is larger than 180 deg 
LON(idx180:end) = LON(idx180:end) - 360; % subtract 360 degrees to make lon interval [-180 180] instead of [0 360] 
%
 
[lon, lat] = meshgrid(LON,LAT); % compute the lat/lon of every grid point in topo
pcolorm(lat,lon,topo); % plot the matrix of elevations on the map
%plotm(platelat,platelon)

%
c = colorbar; 
colormap((jet(20))); % colorbar with jet scale.
c.Color = [1 1 1]; % colorbar text white

c.Label.String = ('Ocean depth [km]'); %colormap title

%% Step 3: Compute seafloor age
%
% $d = 2.65 + 0.345t^{1/2}$ Solve for t
%
% $t =((d-2.65)/0.345))^{2}$

%% Step 4: Plot the oceanic lithosphere age map
ageSeafloor= zeros(180,360);

depth = -topo;

for ix = 1:numel(depth);
    
    ageSeafloor(ix) = ((depth(ix) - 2.65)/(0.345))^2;
    
    if ageSeafloor(ix) == ((0- 2.65)/(0.345))^2;
        ageSeafloor(ix) = -10;
    end
end



l = figure;
hold on
l.InvertHardcopy = 'off'; %Figure background color set to black when printing
l.Color = 'k';  %Set figure window background color to black
l.Position = [100 100 1000 500]; %The location and size of figure's drawable area set to this extent
l.PaperPositionMode = 'auto'; %Printed figure matches the displayed figure size
%
ax = axesm('Mollweid', 'Frame', 'on', 'Grid', 'on'); %Defined the projection, turned frame and grid on
setm(ax,'MLabelLocation',60); %Set meridian label locations to every 60 degrees
setm(ax,'PLabelLocation',30); %Set parallel label locations to every 30 degrees
mlabel('MLabelParallel',0);  %Sets the location of the parallel labels
plabel('PLabelMeridian',-25);  %Sets the location of the meridian labels
axis('off');  %Turns off axis labeling, tickmarks and background
setm(ax,'FontColor',[0.9 0.9 0.9]); %Sets font color to nearly white
setm(ax,'GColor',[0.9 0.9 0.9]);
setm(ax,'fontweight','bold','fontsize',16);
% 
plotm(coastlat, coastlon);
LAT = topolatlim(1):topolatlim(2);
LON =  topolonlim(1):topolonlim(2);
%
idx180 = find(LON > 180, 1); % find first index of longitude array that is larger than 180 deg 
LON(idx180:end) = LON(idx180:end) - 360; % subtract 360 degrees to make lon interval [-180 180] instead of [0 360] 
%
 
[lon, lat] = meshgrid(LON,LAT); % compute the lat/lon of every grid point in topo
pcolorm(lat,lon,ageSeafloor); % plot the matrix of elevations on the map
%plotm(platelat,platelon,'b')



cmap = flipud( jet(20) ); % create a flipped jet colormap 
cmap = [0.5 0.5 0.5 ; cmap]; % 
cmap(end,:) = []; % (1 pt.) your comment here
colormap(cmap); % (1 pt.) your comment here
cmap = colorbar; 

cmap.Color = ([1 1 1]);
cmap.Label.String = ('Age [Ma]');

%% Step 5: Discussion

% • Does your map of ocean ages make sense given the plate boundaries?
%Yes. In general the sea floor is the youngest along spreading boundaries
%with a gradually increasing age away from the boundaries.

% • What is the oldest age in your map?
%The oldest age on our map is 195.4Ma.
%
%• Where does this oldest age occur and does this make sense geologically?
%This age occurs in the northeast Pacific. There are also some very old
%sections in the northeast Atlantic. These areas are the farthest from
%MOD's and given the relationship between time and distance this makes
%sense.
%
% • Where do the youngest ages occur? Does this conform to your knowledge of oceanic lithosphere
% generation?
%The youngest ages occur right along the midoocean ridge whihc makes sense
%because that is where new oceanic lithosphere is being created.
%
% • Are there any assumptions that have gone into this model that might not be accurate?
%We assumed that the teperature at the spreading centers is 640 degrees C
%and that there are even spreading rates everywhere. 


