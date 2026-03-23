function x = icumdf(dist,y,varargin)
switch dist 
    case 'norm'
        mu = varargin{1}; sigma = varargin{2};
        x = mu  + sigma*sqrt(2)*erfinv(2*y -1); 
    case 't'
        nu = varargin{1};
        x = sqrt(nu/betaincinv(2*(1-y),nu/2,1/2) - nu);
    case 'gamma'
        alpha = varargin{1}; lambda = varargin{2};
            % different from Stats toolbox
        x = gammaincinv(y,alpha)/lambda; 
    case 'chi2'
        n = varargin{1};
        x = gammaincinv(y,n/2)*2;
    case 'F'
        m = varargin{1}; n = varargin{2};
        x = n/m/betaincinv(1-y,n/2,m/2) - n/m;
end





