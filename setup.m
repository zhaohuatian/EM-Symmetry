function setup
%SETUP Add production modules and runnable examples. No genpath required.
root=fileparts(mfilename('fullpath'));
names={'01_PlaneWaveMultilayer','02_DipoleMultilayer','03_DipoleSphere'};
for j=1:numel(names)
 addpath(fullfile(root,'modules',names{j}));
 addpath(fullfile(root,'modules',names{j},'examples'));
end
end
