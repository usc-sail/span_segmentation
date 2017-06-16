%
%	function y = kb(u,w,beta)
%
%	Computes the Kaiser-Bessel function used for gridding, namely
%
%	y = f(u,w,beta) = I0 [ beta*sqrt(1-(2u/w)^2) ]/w
%
%	where I0 is the zero-order modified Bessel function
%		of the first kind.
%
%	INPUT:
%		u = vector of k-space locations for calculation.
%		w = width parameter - see Jackson et al.
%		beta = beta parameter - see Jackson et al.
%
%	OUTPUT:
%		y = vector of Kaiser-Bessel values.
%

%	B. Hargreaves	Oct, 2003.




function y = kb(u,w,beta)

if (nargin < 3)
	error('Not enough arguments -- 3 arguments required. ');
end;


if (length(w) > 1)
	error('w should be a single scalar value.');
end;


y = 0*u;				% Allocate space.
uz = find(abs(u)< w/2);			% Indices where u<w/2.

if (length(uz) > 0)				% Calculate y at indices uz.
	x = beta*sqrt(1-(2*u(uz)/w).^2);	% Argument - see Jackson '91.
	y(uz) = besseli(0,x)./w;
end;





