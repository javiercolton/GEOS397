%JavierColton Homework 1
%% Part 1

clc
clear all

load ('topo.mat');
load('coastlines');

% setup figure properties that we want.
h = figure;
h.InvertHardcopy = 'off';
h.Color = 'k';
h.Position = [100 100 1000 500];
h.PaperPositionMode = 'auto';

%setup the map axes
ax = axesm('Mollweid', 'Frame', 'on', 'Grid','on');
setm(ax, 'MLabelLocation', 60);
setm(ax, 'PLabelLocation', 30);
mlabel('MLabelParallel', 0);
plabel('PLabelMeridian', -25);
axis('off');
setm(ax, 'FontColor', [0.9 0.9 0.9]);
setm(ax, 'GColor', [0.9 0.9 0.9]);

plotm(coastlat, coastlon);
LAT = topolatlim(1):topolatlim(2);
LON =  topolonlim(1):topolonlim(2);

idx180 = find(LON > 180, 1); % find first index of longitude array that is larger than 180 deg 
LON(idx180:end) = LON(idx180:end) - 360; % subtract 360 degrees to make lon interval [-180 180] instead of [0 360] 

[lon, lat] = meshgrid(LON,LAT); % compute the lat/lon of every grid point in topo
pcolorm(lat,lon,topo); % plot the matrix of elevations on the map
demcmap(topo); % give it a better colormap

colorbar; 
c = colorbar;
c.Color = 'w'
ylabel(c, 'Elavation [m]');
hold off

%%
distdim(topo, 'm', 'km');
kmtopo=ans; %pulls answers from function above and names it kmtopo


g = figure;
g.InvertHardcopy = 'off';
g.Color = 'k';
g.Position = [100 100 1000 500];
g.PaperPositionMode = 'auto';

%setup the map axes
ax = axesm('Mollweid', 'Frame', 'on', 'Grid','on');
setm(ax, 'MLabelLocation', 60);
setm(ax, 'PLabelLocation', 30);
mlabel('MLabelParallel', 0);
plabel('PLabelMeridian', -25);
axis('off');
setm(ax, 'FontColor', [0.9 0.9 0.9]);
setm(ax, 'GColor', [0.9 0.9 0.9]);

plotm(coastlat, coastlon);
LAT = topolatlim(1):topolatlim(2);
LON =  topolonlim(1):topolonlim(2);

idx180 = find(LON > 180, 1); % find first index of longitude array that is larger than 180 deg 
LON(idx180:end) = LON(idx180:end) - 360; % subtract 360 degrees to make lon interval [-180 180] instead of [0 360] 

[lon, lat] = meshgrid(LON,LAT); % compute the lat/lon of every grid point in topo
pcolorm(lat,lon,kmtopo); % plot the matrix of elevations on the map
demcmap(kmtopo); % give it a better colormap

colorbar; 
c = colorbar;
c.Color = 'w';

ylabel(c, 'Elavation [km]');
hold off


kmtopo(kmtopo>0)=0;
%%
f = figure;
f.InvertHardcopy = 'off';
f.Color = 'k';
f.Position = [100 100 1000 500];
f.PaperPositionMode = 'auto';

%setup the map axes
ax = axesm('Mollweid', 'Frame', 'on', 'Grid','on');
setm(ax, 'MLabelLocation', 60);
setm(ax, 'PLabelLocation', 30);
mlabel('MLabelParallel', 0);
plabel('PLabelMeridian', -25);
axis('off');
setm(ax, 'FontColor', [0.9 0.9 0.9]);
setm(ax, 'GColor', [0.9 0.9 0.9]);

plotm(coastlat, coastlon);
LAT = topolatlim(1):topolatlim(2);
LON =  topolonlim(1):topolonlim(2);

idx180 = find(LON > 180, 1); % find first index of longitude array that is larger than 180 deg 
LON(idx180:end) = LON(idx180:end) - 360; % subtract 360 degrees to make lon interval [-180 180] instead of [0 360] 

[lon, lat] = meshgrid(LON,LAT); % compute the lat/lon of every grid point in topo
pcolorm(lat,lon,kmtopo); % plot the matrix of elevations on the map
demcmap(kmtopo); % give it a better colormap
caxis([-8, 5])

colorbar; 
c = colorbar;
c.Color = 'w';

ylabel(c, 'Elavation [km]');
hold off

% To Find the volume one must also find the area of each pixel. We the have
% to multiply the area times the depth and add it up for 
%%

circ = 2*pi*6371
pixelWidth = circ/360
pixelHeight = (circ/2)/180
pixelArea = pixelWidth * pixelHeight

volumekmtopo = pixelArea * kmtopo;
negativeTotalVolume = sum(sum(volumekmtopo));
totalVolume = abs(negativeTotalVolume)

%% Part 2
%We are slightly overestimating the value of ocean volume using this
%method.

Acell = @(lat1,lat2) (pi*2*(6371^2)*abs(sin((pi/180)*lat1) - sin((pi/180)*lat2)))/360
Acell(29,30)

%%










