%
%	function [kern,kbu] = calckbkernel(kwidth,overgridfactor,klength)
%
%	Function calculates the appropriate Kaiser-Bessel kernel
%	for gridding, using the approach of Jackson et al.
%
%	INPUT:
%		kwidth = kernel width in grid samples.
%		overgridfactor = over-gridding factor.
%		klength = kernel look-up-table length.			
%
%	OUTPUT:
%		kern = kernel values for klength values
%			of u, uniformly spaced from 0 to kwidth/2.
%		kbu = u values.

%	B. Hargreaves



function [kern,kbu] = calckbkernel(kwidth,overgridfactor,klength)

if (nargin < 3)
	klength = 32;                   	  	% Pretty Arbitrary.
end;
if (klength < 2)
	klength = 2;	
	disp('Warning:  klength must be 2 or more - using 2.');
end;

%beta = pi*kwidth*(overgridfactor-0.5);		% From Jackson et al.

a = overgridfactor;
w = kwidth;	
beta = pi*sqrt( w^2/a^2*(a-0.5)^2-0.8 );	% From Beatty et al.



u = [0:klength-1]/(klength-1) * kwidth/2;     	% Kernel radii - grid samples. 

kern = kb(u,kwidth, beta);
kern=kern/kern(1);                 % Normalize.
kbu = u;


