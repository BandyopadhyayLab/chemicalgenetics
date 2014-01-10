% modification of the matlab ttest function to incorporate a minimum
% variance measurement.
% Author: SB (modified from matlab)

function [h,p,ci,stats] = ttest_knownvariance(x,y,v_x,v_y, alpha,tail,vartype,dim,minvar)


if nargin < 2
    error('stats:ttest2:TooFewInputs','Requires at least two input arguments');
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1
    error('stats:ttest2:BadAlpha','ALPHA must be a scalar between 0 and 1.');
end

if nargin < 4 || isempty(tail)
    tail = 0;
elseif ischar(tail) && (size(tail,1)==1)
    tail = find(strncmpi(tail,{'left','both','right'},length(tail))) - 2;
end
if ~isscalar(tail) || ~isnumeric(tail)
    error('stats:ttest2:BadTail', ...
          'TAIL must be one of the strings ''both'', ''right'', or ''left''.');
end

if nargin < 5 || isempty(vartype)
    vartype = 1;
elseif ischar(vartype) && (size(vartype,1)==1)
    vartype = find(strncmpi(vartype,{'equal','unequal'},length(vartype)));
end
if ~isscalar(vartype) || ~isnumeric(vartype)
    error('stats:ttest2:BadVarType', ...
          'VARTYPE must be one of the strings ''equal'' or ''unequal''.');
end

if nargin < 6 || isempty(dim)
    % Figure out which dimension mean will work along by looking at x.  y
    % will have be compatible. If x is a scalar, look at y.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = find(size(y) ~= 1, 1); end
    if isempty(dim), dim = 1; end
    
    % If we haven't been given an explicit dimension, and we have two
    % vectors, then make y the same orientation as x.
    if isvector(x) && isvector(y)
        if dim == 2
            y = y(:)';
        else % dim == 1
            y = y(:);
        end
    end
end

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
    error('stats:ttest2:InputSizeMismatch',...
          'The data in a 2-sample t-test must be commensurate.');
end

xnans = isnan(x);
if any(xnans(:))
    nx = sum(~xnans,dim);
else
    nx = size(x,dim); % a scalar, => a scalar call to tinv
end
ynans = isnan(y);
if any(ynans(:))
    ny = sum(~ynans,dim);
else
    ny = size(y,dim); % a scalar, => a scalar call to tinv
end

%
% Incorporation of minimum variance
%

s2y = v_y;
s2x = v_x;

s2y(find(s2y<minvar)) = minvar;
s2x(find(s2x<minvar)) = minvar;

difference = nanmean(x,dim) - nanmean(y,dim);
if vartype == 1 % equal variances
    dfe = nx + ny - 2;
    sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
    se = sPooled .* sqrt(1./nx + 1./ny);
    ratio = difference ./ se;

    if (nargout>3)
        stats = struct('tstat', ratio, 'df', cast(dfe,class(ratio)), ...
                       'sd', sPooled);
        if isscalar(dfe) && ~isscalar(ratio)
            stats.df = repmat(stats.df,size(ratio));
        end
    end
elseif vartype == 2 % unequal variances
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
    se = sqrt(s2xbar + s2ybar);
    ratio = difference ./ se;

    if (nargout>3)
        stats = struct('tstat', ratio, 'df', cast(dfe,class(ratio)), ...
                       'sd', sqrt(cat(dim, s2x, s2y)));
        if isscalar(dfe) && ~isscalar(ratio)
            stats.df = repmat(stats.df,size(ratio));
        end
    end
    
    % Satterthwaite's approximation breaks down when both samples have zero
    % variance, so we may have gotten a NaN dfe.  But if the difference in
    % means is non-zero, the hypothesis test can still reasonable results,
    % that don't depend on the dfe, so give dfe a dummy value.  If difference
    % in means is zero, the hypothesis test returns NaN.  The CI can be
    % computed ok in either case.
    if se == 0, dfe = 1; end
else
    error('stats:ttest2:BadVarType',...
          'VARTYPE must be ''equal'' or ''unequal'', or 1 or 2.');
end

% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = 2 * tcdf(-abs(ratio),dfe);
    if nargout > 2
        spread = tinv(1 - alpha ./ 2, dfe) .* se;
        ci = cat(dim, difference-spread, difference+spread);
    end
elseif tail == 1 % right one-tailed test
    p = tcdf(-ratio,dfe);
    if nargout > 2
        spread = tinv(1 - alpha, dfe) .* se;
        ci = cat(dim, difference-spread, Inf(size(p)));
    end
elseif tail == -1 % left one-tailed test
    p = tcdf(ratio,dfe);
    if nargout > 2
        spread = tinv(1 - alpha, dfe) .* se;
        ci = cat(dim, -Inf(size(p)), difference+spread);
    end
else
    error('stats:ttest2:BadTail',...
          'TAIL must be ''both'', ''right'', or ''left'', or 0, 1, or -1.');
end

% Determine if the actual significance exceeds the desired significance
h = cast(p <= alpha, class(p));
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha
