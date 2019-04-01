function mesh = make_line_mesh(h)
% MAKE_LINE_MESH creates a square mesh
% INPUT
%   h: approximate diameter of elements
% OUTPUT
%   mesh: mesh structure
% REMARKS
%   boundary groups: 
%     1: left 
%     2: right 

% Copyright 2018 Masayuki Yano, University of Toronto

if (nargin < 1)
    h = 1;
end

ne1d = ceil(1.0/h);
mesh.coord = linspace(0,1,ne1d+1)';
mesh.tri = [1:ne1d; 2:ne1d+1]';

% left
mesh.bgrp{1} = 1;

% right
mesh.bgrp{2} = ne1d+1;

end