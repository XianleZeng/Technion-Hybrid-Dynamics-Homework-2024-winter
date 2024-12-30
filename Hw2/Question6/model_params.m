function [m, l, I_c, d, w, b, c]=model_params(varargin)

options.damping = false;

if nargin > 0
    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
end

if options.damping
    c = 1;
else
    c = 0;
end

m = 30;
l = 0.6;
I_c = 0.15;
d = 0.15;
w = 0.25;
b = 0.2;

end

