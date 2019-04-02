%prediction script
addpath('./mesh')
addpath('./featureFunctions/nonOverlappingPolydisperseSpheres')
addpath('./FEM')
addpath('./rom')
addpath('./comp')
addpath('./aux')

startup;        %proper plot configuration

rom = StokesROM;
testSamples = 16:32;

testData = StokesData(testSamples);
[~,predVar,~, meanSqDist, ~, mll, R] = rom.predict(testData, 'local');

save('./prediction.mat', 'meanSqDist', 'mll', 'R');

