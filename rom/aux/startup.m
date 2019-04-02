%Show execution of startup file
disp('   Executing startup.m in ./aux')

%adding paths
addpath(genpath('../plotting'))

%Plotting defaults
set(groot,'defaultLineLineWidth',1)%default line width
screen=get(0,'ScreenSize');
%default figure position
set(groot, 'defaultFigurePosition',[10 screen(4) - 10 720 540])
clear screen
set(0,'DefaultAxesLineWidth', 1.5)%Axes line width
set(0,'DefaultAxesFontSize',18)%Axes font size
set(0,'DefaultAxesXGrid', 'on')%grid by default
set(0,'DefaultAxesYGrid', 'on')
co = distinguishable_colors(12);%distinguishable color order
set(groot,'defaultAxesColorOrder', co)
clear co
cm_inferno = inferno();%take python inferno as the default colormap
set(groot,'DefaultFigureColormap',cm_inferno)%default color map
clear cm_inferno

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
set(groot, 'defaultTextboxshapeInterpreter', 'latex');
set(groot, 'defaultTextarrowshapeInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter', 'latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesBoxStyle', 'full');
set(groot, 'defaultFigureColor', [1 1 1]); %white figure background
%default axis 'square'
set(groot,'defaultAxesPlotBoxAspectRatioMode', 'manual');
set(groot,'defaultAxesPlotBoxAspectRatio', [1 1 1]);