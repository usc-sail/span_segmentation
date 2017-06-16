%
%	function [dat] = gridlut(ktraj,ksamps,dcf,gridsize,convhalfwidth,
%					kerneltable)
%
%	Function does gridding with a pre-calculated density function,
%		and gridding kernel from the given kerneltable (lookup table).
%
%	ktraj		K-space sample locations (complex).
%	ksamps		K-space data (complex).
%	dcf		Density compensation factors at each sample.
%	gridsize 	Size of grid onto which to grid data.
%	convhalfwidth	(optional) width of convolution kernal in grid points.
%	kerneltable	Convolution kernel lookup table.
%
%	dat		Gridded data (grid goes from -.5 to 0.5 in Kx and Ky).
%
%

function [dat] = gridlut(ktraj,ksamps,dcf,gridsize,convwidth, kerneltable)



% --------- Check Enough Inputs ------
if (nargin < 3)
	disp('Usage:  dat=gridlut(ktraj,ksamps,dcf,gridsize,cwidth,kerneltable, nkpts)'); 
	disp(' ');
	error('Not enough input arguments');
end;

% --------- Check k-space trajectory is complex-valued -------
s = size(ktraj);
if (s(2)==2)
	ktraj = ktraj(:,1)+i*ktraj(:,2);
end;
if isreal(ktraj)
	ktraj = ktraj + 2*eps*i;
	disp('Warning: k-space trajectory should be complex-valued, or Nx2');
end;


% --------- Check k-space samples are complex-valued -------
s = size(ksamps);
if (s(2)==2)
	ksamps = ksamps(:,1)+i*ksamps(:,2);
end;
if isreal(ksamps)
	ksamps = ksamps + 2*eps*i;
	disp('Warning k-space samples should be complex-valued, or Nx2');
end;


% --------- Check k-space file does not exceed 0.5 -------
if ( max(abs(real(ktraj)))>=0.5)  | (max(abs(imag(ktraj)))>=0.5 )
	disp('Warning:  k-space location radii should not exceed 0.5.');
	disp('		Some samples could be ignored.');
end;

% --------- Check DCFs do not exceed 1.0 -------
if ( max(abs(dcf)) > 1.0)
	disp('Warning:  DCF values should be between 0 and 1.0');
	disp('		(Some Recon Programs may fail)');
end;


% --------- Default Arguments -------------
if (nargin < 4)
	gridsize=256;
end;
if (nargin < 5)
	convwidth=2;
end;
if (nargin < 6)
	kerneltable=[1 0];	 
end;

% --------- Call Mex Function for this ---------

dat = gridlut_mex(ktraj,ksamps,dcf,gridsize,convwidth,kerneltable);


