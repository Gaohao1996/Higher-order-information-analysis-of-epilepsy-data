function I = cmi_ggg(x, y, z, biascorrect, demeaned, use_cov,unit)
% CMI_GGG Conditional Mutual information (CMI) between two Gaussian variables
%        conditioned on a third
%
%   I = cmi_ggg(x,y,z) returns the CMI between two (possibly multidimensional)
%   Gassian variables, x and y, conditioned on a third, z, with bias correction.
%   If x / y / z are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if isvector(z)
    z = z(:);
end
if ndims(x)~=2 || ndims(y)~=2 || ndims(z)~=2
    error('cmi_ggg: input arrays should be 2d')
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);
Nvarz = size(z,2);

if (size(y,1) ~= Ntrl) || (size(z,1) ~= Ntrl)
    error('cmi_ggg: number of trials do not match')
end

% default option values
if nargin<4
    biascorrect = true;
end
if nargin<5
    demeaned = false;
end
if nargin<6
    use_cov = true;
end
if nargin < 7 || isempty(unit)
    unit = 'nats';
end


% demean data if required
if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrl);
    y = bsxfun(@minus,y,sum(y,1)/Ntrl);
    z = bsxfun(@minus,z,sum(z,1)/Ntrl);
end

% joint variable
xyz = [x y z];

if use_cov
    Cxyz = (xyz'*xyz) / (Ntrl - 1);
else
    Cxyz = shrinkage_cov(xyz); %use when needed
end
% submatrices of joint covariance
Nvaryz = Nvary + Nvarz;
Nvarxyz = Nvarx + Nvaryz;
zidx = (Nvarx + Nvary + 1):Nvarxyz;
Cz = Cxyz(zidx,zidx);

idx = (Nvarx + 1):Nvarxyz;
Cyz = Cxyz(idx, idx);

Nvarxz = Nvarx + Nvarz;
Cxz = zeros(Nvarxz);
xidx = 1:Nvarx;
Cxz(xidx,xidx) = Cxyz(xidx,xidx);
zidxxz = (Nvarx+1):Nvarxz;
Cxz(xidx,zidxxz) = Cxyz(xidx,zidx);
Cxz(zidxxz,xidx) = Cxyz(zidx,xidx);
Cxz(zidxxz,zidxxz) = Cxyz(zidx,zidx);

Cz = Cz+1e-8*eye(size(Cz));
Cxz = Cxz+1e-8*eye(size(Cxz));
Cyz = Cyz+1e-8*eye(size(Cyz));
Cxyz = Cxyz+1e-8*eye(size(Cxyz));

chCz = chol(Cz);
chCxz = chol(Cxz);
chCyz = chol(Cyz);
chCxyz = chol(Cxyz);

% entropies in nats
% normalisations cancel for cmi
HZ = sum(log(diag(chCz)));% + 0.5*Nvarz*log(2*pi*exp(1));
HXZ = sum(log(diag(chCxz)));% + 0.5*(Nvarx+Nvarz)*log(2*pi*exp(1));
HYZ = sum(log(diag(chCyz)));%+ 0.5*(Nvary+Nvarz)*log(2*pi*exp(1));
HXYZ = sum(log(diag(chCxyz)));% + 0.5*(Nvarx+Nvary+Nvarz)*log(2*pi*exp(1));

ln2 = log(2);
if biascorrect
    psiterms = psi((Ntrl - (1:Nvarxyz))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    HZ = (HZ - Nvarz*dterm - sum(psiterms(1:Nvarz)));
    HXZ = (HXZ - Nvarxz*dterm - sum(psiterms(1:Nvarxz)));
    HYZ = (HYZ - Nvaryz*dterm - sum(psiterms(1:Nvaryz)));
    HXYZ = (HXYZ - Nvarxyz*dterm - sum(psiterms));
end

% convert to bits
I_nats = (HXZ + HYZ - HXYZ - HZ);
switch lower(unit)
    case 'nats'
        I = I_nats;
    case 'bits'
        I = I_nats / log(2);
    otherwise
        error('unit must be ''nats'' or ''bits''.');
end

