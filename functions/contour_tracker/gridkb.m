%
%   [dat,kbku,kbkval]  = gridkb(ktraj,ksamps,dcf,gridsize,kwidth,overgridfactor)
%
%	OR...
%
%   [dat,kbku,kbkval]  = gridkb(kspfile,ksamps,gridsize,kwidth,overgridfactor)
%
%	Function does gridding with a pre-calculated density function,
%		and Kaiser-Bessel gridding kernel, as described in 
%		Jackson et al, 1991.  The k-space locations (ktraj) are 
%		normalized (ie  |real(ktraj)|<0.5 and |imag(k)|<0.5  ).  
%
%	ktraj		K-space sample locations (complex), normalized so
%				that |real(k)|<0.5 and |imag(k)|<0.5.  
%	ksamps		K-space data (complex).
%	dcf		Density compensation factors at each sample.
%	gridsize 	Integer size of grid onto which to grid data.
%	kwidth		Kernel width, w as described in Jackson's paper.
%	overgridfactor	Factor by which these grid parameters overgrid.
%			(This is gridsize/(FOV/resolution) )
%
%	dat		Gridded data (grid goes from -.5 to 0.5 in Kx and Ky).
%
%	kbkval		Kaiser-bessel kernel values.
%	kbku		Kaiser-bessel radius in grid points.

function [dat,kbku,kbkval] = gridkb(varargin)


% ================= First assign input arguments =====================

% 	Get k-space trajectory, DCFs and data, either from parameter 
% 		list or from file.

if (~ischar(varargin{1}) )
	ktraj = varargin{1};
	ksamps = varargin{2};
	dcf = varargin{3};
	otherargs = {varargin{4:end}};

else						% Read k-space file.
	[ktraj,dcf] = readkspace(varargin{1});
	ksamps = varargin{2};
	otherargs = {varargin{3:end}};
end;

% 	---- Check lengths ------

if (max(size(ktraj)) ~= max(size(dcf)))
	disp('Warning:  K-space trajectory and DCFs have different length.');
end;
if (max(size(ktraj(:))) ~= max(size(ksamps(:))))
	disp('Warning:  K-space trajectory and signal data have different length.');
end;



%	---- Now get other arguments, or assign defaults. ----
%

if (length(otherargs)>=1)		% -- Get gridsize. 
	gridsize = otherargs{1}; 
else
	disp('Warning:  No gridsize passed - using gridsize=256');
	gridsize = 256;
end;

if (length(otherargs)>=2)		% -- Get kernel width.
	kwidth = otherargs{2}; 
else
	disp('Warning:  No kernel width passed - using width=1.5');
	kwidth = 1.5;
end;

if (length(otherargs)>=3)		% -- Get overgridfactor
	overgridfactor = otherargs{3};
else
	disp('Warning:  No overgridfactor passed.  Using overgridfactor=1.5');
	disp('		This is an important parameter in selecting');
	disp('		the convolution kernel, so if you are actually');
	disp('		using these images, you should specifiy this.');
	overgridfactor = 1.5;
end;





% ========== Calculate Kaiser-Bessel Kernel for given parameters =========

[kerneltable,u] = calckbkernel(kwidth,overgridfactor);



% ========== Do Gridding =========

dat = gridlut(ktraj,ksamps,dcf,gridsize,kwidth/2, kerneltable);


kbkval = kerneltable;		% ---- Assign output value.
kbku = u;			% ---- Assign output value.


