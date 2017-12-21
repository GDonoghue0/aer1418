function setup(mode)
if (nargin < 1), mode = 1; end

dirs = {'distmesh','kernel'};
d0=fileparts([pwd,filesep]);
d0=[d0,filesep];
switch mode
 case 1 % initialization
  fprintf('   --> Initializing fem2d ...')
  for k = 1:length(dirs)
    addpath([d0,dirs{k}]);
  end
  fprintf(' Done.\n\n');
 case 0 % uninitializing
  fprintf('   --> Uninitializing fem2d ...')
  for k = 1:length(dirs)
    rmpath([d0,dirs{k}]);
  end
  fprintf(' Done.\n\n');
end

