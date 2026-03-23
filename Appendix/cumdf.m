function y = cumdf(dist,x,varargin)
switch dist 
    case 'norm'
        mu = varargin{1}; sigma = varargin{2};
        y = (erf(( (x - mu)/sigma )/sqrt(2)) + 1)/2;
    case 't'
        nu = varargin{1}; 
        y = 1-0.5*betainc(nu/(nu+x.^2),nu/2,1/2);
    case 'gamma'
        alpha = varargin{1}; lambda = varargin{2};
            % different from Stats toolbox
        y = gammainc(lambda*x,alpha); 
    case 'chi2'
        n = varargin{1};
        y = gammainc(x/2,n/2);
    case 'F'
        m = varargin{1}; n = varargin{2};
        y = 1 - betainc(n/(n+m*x),n/2,m/2);
end



