function x=ntdetrend(x,order,w,basis)
%y=nt_detrend(x,order,w,basis) - remove polynomial or sinusoidal trend
% 
%  y: detrended data
%
%  x: raw data
%  order: order of polynomial
%  w: weight
%  basis: 'polynomials' [default] or 'sinusoids'

if nargin<2; error('!'); end
if nargin<3; w=[]; end
if nargin<4||isempty(basis); basis='polynomials'; end

dims=size(x);
x=x(:,:); % concatenates dims >= 2

% regressor
switch basis
    case 'polynomials'
        r=zeros(size(x,1),order);
        lin=linspace(-1,1,size(x,1));
        for k=1:order
            r(:,k)=lin.^k;
        end
    case 'sinusoids'
        r=zeros(size(x,1),order*2);
        lin=linspace(-1,1,size(x,1));
        for k=1:order
            r(:,2*k-1)=sin(2*pi*k*lin/2);
            r(:,2*k)=cos(2*pi*k*lin/2);
        end
    otherwise
        error('!');
end

x=nt_demean(x,w); % order zero, remove mean
if size(r,2)>0; 
    x=nt_tsr(x,r,0,w); % project out regressor
end

x=reshape (x,dims);




