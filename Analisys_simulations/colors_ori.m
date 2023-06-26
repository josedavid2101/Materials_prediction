function grainColors = colors_ori(ori) 

% Using MTEX this functions generates colors in RGB space depending on
% orientations. Use "figure; plot(oM)" to plot the legend

%% Initialize MTEX library
% mtex = genpath('mtex-5.3');
%  addpath(mtex)
% startup_mtex

%% generate colors

CS = crystalSymmetry('cubic'); %Define Crystal Symmetry
oM = ipdfHSVOrientationMapping(CS); %Orientation mapping
o = orientation('quaternion',ori(:,1),ori(:,2),ori(:,3),ori(:,4),CS); %orientations in MTEX format

grainColors = oM.orientation2color(o);

% figure
%     plot(oM)
end