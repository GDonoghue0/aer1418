function grd = create_rectangular_qN(N,L,LF,porder,parityFlag);
% function grd = create_rectangular_q1(varargin);
% 
% Inputs:  (Note some or all inputs may be specified )
%    
%  N  = [Nx, Ny] :   number of elements in x and y direction
%  L  = [Lx, Ly] :   length of box in x and y direction
%  LF = [LFx, LFy] : logflag for x and y directions.
%
%  Note: LF = 0 implies linear, LF > 0 implies logarithmic 
%        increasing in size for left to right or bottom to top,
%        LF < 0 implies logarithmic in opposite direction. Magnitude
%        of LF determines amount of shifting in size of elements
%
if (nargin < 1) N = [3,3]; end
if (nargin < 2) L = [1,1]; end
if (nargin < 3) LF = [0,0]; end
if (nargin < 4) porder = 1; end
if (nargin < 5)
  % fill out parity flag
  parityFlag = zeros(N(1),N(2));
% $$$   for i = 1:N(1)
% $$$     for j = 1:(2)
% $$$       if (((i==1) && (j==N(2))) || ((i==N(1)) && (j==1)))
% $$$         parityFlag(i,j) = 0;
% $$$       else
% $$$         parityFlag(i,j) = 1;
% $$$       end
% $$$     end
% $$$   end
% symm cut
% $$$ for i = 1:N(1)
% $$$   for j = 1:N(2)
% $$$     if (i <= N(1)/2)
% $$$       parityFlag(i,j) = 1;
% $$$     else
% $$$       parityFlag(i,j) = 0;
% $$$     end
% $$$ % $$$     if ((i==1) && (j==1))
% $$$ % $$$       parityFlag(i,j) = 1;
% $$$ % $$$     end
% $$$ % $$$     if ((i==Nx) && (j==1))
% $$$ % $$$       parityFlag(i,j) = 0;
% $$$ % $$$     end
% $$$   end
% $$$ end
end

% Set number of points in x and y directions 
Nx = N(1);
Ny = N(2);

% Set length in x and y directions
Lx = L(1);
Ly = L(2);

% Get Log Flag 
LFx = LF(1);
LFy = LF(2);






% Set x Coordinates 
nxNode = porder*Nx+1;
if (LFx == 0) 
  x = linspace(0, Lx, nxNode);
elseif (LFx > 0 )
  x = Lx/(10^LFx-1)*(logspace(0, LFx, nxNode)-1);
else 
  x = Lx/(1-10^LFx)*(1-logspace(LFx, 0, nxNode));  
end

% Set y Coordinates 
nyNode = porder*Ny+1;
if (LFy == 0) 
  y = linspace(0, Ly, nyNode);
elseif (LFy > 0 )
  y = Ly/(10^LFy-1)*(logspace(0, LFy, nyNode)-1);
else 
  y = Ly/(1-10^LFy)*(1-logspace(LFy, 0, nyNode));  
end

% Create coordinates 
coordinate=zeros(nxNode*nyNode,2);
nNode = 1;
for i=1:nxNode
  for j=1:nyNode
    coordinate(nNode,:)=[x(i),y(j)];
    nNode = nNode+1;
  end
end

pp1 = porder + 1;
Node=zeros(2*Nx*Ny, (porder+1)*(porder+2)/2);
nElem = 1;
for i=1:Nx
  for j=1:Ny
    
    for k = 1:pp1
      n((k-1)*pp1+1:k*pp1,1) = (((i-1)*porder+(k-1))*nyNode + (j-1)*porder) + [1:pp1]';
    end
    
    switch porder 
     case 1
      %
      % 2-4  2-4 
      % |\|  |/|
      % 1-3  1-3
      %
      if (parityFlag(i,j) == 0)
        Node(nElem,:) = [n(1),n(3),n(2)];
        nElem = nElem+1;
        Node(nElem,:) = [n(3),n(4),n(2)];
        nElem = nElem+1;
      else
        Node(nElem,:) = [n(1),n(3),n(4)];
        nElem = nElem+1;
        Node(nElem,:) = [n(1),n(4),n(2)];
        nElem = nElem+1;
      end
      
     case 2
      %
      % 3-6-9  3-6-9
      % |\  |  |  /|
      % 2 5 8  2 5 8
      % |  \|  |/  |
      % 1-4-7  1-4-7
      %
      if (parityFlag(i,j) == 0)
        Node(nElem,:) = [n(1),n(7),n(3),n(5),n(2),n(4)];
        nElem = nElem+1;
        Node(nElem,:) = [n(7),n(9),n(3),n(6),n(5),n(8)];
        nElem = nElem+1;
      else
        Node(nElem,:) = [n(1),n(7),n(9),n(8),n(5),n(4)];
        nElem = nElem+1;
        Node(nElem,:) = [n(1),n(9),n(3),n(6),n(2),n(5)];
        nElem = nElem+1;
       end
       
     case 3
      %
      % 4- 8-12-16  4- 8-12-16
      % |\       |  |      / | 
      % 3  7 11 15  3  7 11 15 
      % |   \    |  |   /    | 
      % 2  6 10 14  2  6 10 14   
      % |      \ |  |/       |
      % 1- 5- 9-13  1- 5- 9-13
      %
      if (parityFlag(i,j) == 0)
        Node(nElem,:) = [n(1),n(13),n(4),n(10),n(7),n(3),n(2),n(5),n(9),n(6)];
        nElem = nElem+1;
        Node(nElem,:) = [n(13),n(16),n(4),n(12),n(8),n(7),n(10),n(14),n(15),n(11)];
        nElem = nElem+1;
      else
        Node(nElem,:) = [n(1),n(13),n(16),n(14),n(15),n(11),n(6),n(5),n(9),n(10)];
        nElem = nElem+1;
        Node(nElem,:) = [n(1),n(16),n(4),n(12),n(8),n(3),n(2),n(6),n(11),n(7)];
        nElem = nElem+1;
      end      
      
     case 4
      %
      % 5-10-15-20-25  5-10-15-20-25 
      % |\          |  |         / | 
      % 4  9 14 19 24  4  9 14 19 24
      % |   \       |  |      /    | 
      % 3  8 13 18 23  3  8 13 18 23 
      % |      \    |  |   /       |
      % 2  7 12 17 22  2  7 12 17 22
      % |         \ |  |/          |
      % 1- 6-11-16-21  1- 6-11-16-21 
      %
      if (parityFlag(i,j) == 0)
        Node(nElem,:) = [n(1),n(21),n(5),n(17),n(13),n(9),n(4),n(3),n(2),n(6),n(11),n(16),n(7),n(12),n(8)];
        nElem = nElem+1;
        Node(nElem,:) = [n(21),n(25),n(5),n(20),n(15),n(10),n(9),n(13),n(17),n(22),n(23),n(24),n(18),n(19),n(14)];
        nElem = nElem+1;
      else
        Node(nElem,:) = [n(1),n(21),n(25),n(22),n(23),n(24),n(19),n(13),n(7),n(6),n(11),n(16),n(12),n(17),n(18)];
        nElem = nElem+1;
        Node(nElem,:) = [n(1),n(25),n(5),n(20),n(15),n(10),n(4),n(3),n(2),n(7),n(13),n(19),n(8),n(14),n(9)];
        nElem = nElem+1;
      end        
      
     case 5
      %
      % 6-12-18-24-30-36   6-12-18-24-30-36
      % |\             |   |            / |
      % 5 11 17 23 29 35   5 11 17 23 29 35  
      % |   \          |   |         /    |
      % 4 10 16 22 28 34   4 10 16 22 28 34 
      % |      \       |   |      /       |
      % 3  9 15 21 27 33   3  9 15 21 27 33
      % |         \    |   |   /          |
      % 2  8 14 20 26 32   2  8 14 20 26 32
      % |            \ |   |/             |
      % 1- 7-13-19-25-31   1- 7-13-19-25-31
      %
      if (parityFlag(i,j) == 0)
        Node(nElem,:) = [n(1),n(31),n(6),n(26),n(21),n(16),n(11),n(5),n(4),n(3),n(2),n(7),n(13),n(19),n(25),n(8),n(14),n(20),n(9),n(15),n(10)];
        nElem = nElem+1;
        Node(nElem,:) = [n(31),n(36),n(6),n(30),n(24),n(18),n(12),n(11),n(16),n(21),n(26),n(32),n(33),n(34),n(35),n(27),n(28),n(29),n(22),n(23),n(17)];
        nElem = nElem+1;
      else
        Node(nElem,:) = [n(1),n(31),n(36),n(32),n(33),n(34),n(35),n(29),n(22),n(15),n(8),n(7),n(13),n(19),n(25),n(14),n(20),n(26),n(21),n(27),n(28)];
        nElem = nElem+1;
        Node(nElem,:) = [n(1),n(36),n(6),n(30),n(24),n(18),n(12),n(5),n(4),n(3),n(2),n(8),n(15),n(22),n(29),n(9),n(16),n(23),n(10),n(17),n(11)];
        nElem = nElem+1;
      end       
      
      
     otherwise
      error('unsuppoted order %d', porder);
    end
    
  end
end



%%%%%%%%%%%%%%%%  Boundary Groups
%        4
%    * ----- *
%    |       |
%  1 |       | 3
%    |       |
%    * ----- *
%        2
grd.nBFaceGroup = 4;
grd.nbface = [Ny,Nx,Ny,Nx];
grd.bface = zeros(max(N),2,4);

bfgrp = 1;
i=1;
for j=1:Ny
  bface = j;
  n1 = (j-1)*porder + 1;
  n2 = j*porder + 1;
  grd.bface(bface,:,bfgrp) = [n1,n2];
end

bfgrp = 2;
j=1;
for i=1:Nx
  bface = i;
  n1 = (i-1)*porder*nyNode + 1;
  n2 = i*porder*nyNode + 1;
  grd.bface(bface,:,bfgrp) = [n1,n2];
end


bfgrp = 3;
i=Nx;
for j=1:Ny
  bface = j;
  n1 = i*porder*nyNode + (j-1)*porder + 1;
  n2 = i*porder*nyNode + j*porder + 1;
  grd.bface(bface,:,bfgrp) = [n1,n2];
end

bfgrp = 4;
j=Ny;
for i=1:Nx
  bface = i;
  n1 = j*porder + (i-1)*porder*nyNode + 1;
  n2 = j*porder + i*porder*nyNode + 1;
  grd.bface(bface,:,bfgrp) = [n1,n2];
end

% create grd
grd.Dim                      = 2;
grd.nNode                    = size(coordinate,1);
grd.coordinate               = coordinate;
grd.nElement                 = size(Node,1);
grd.nElementGroup            = 1;
grd.ElementGroup{1}.nElement = size(Node,1);
grd.ElementGroup{1}.Node     = Node;
grd.ElementGroup{1}.qorder    = porder;


% write out grd
%viewmesh(grd)

