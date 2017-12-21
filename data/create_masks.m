function [Masks,param_data] = create_masks(param_data)
% Masks : cell array of size param_data.nb_tests
%         each cell contains the mask of dimention Ni selecting
%         param_data.M frequencies 

Masks = cell(param_data.nb_tests,1) ;

%%
if ~isfield(param_data,'C'), param_data.C = 1*[1,1] ; end;  % define calibration zone

if param_data.acc==1   % if no undersampling
    mask = ones(param_data.Ni) ;
    for t = 1:param_data.nb_tests
        Masks{t} = mask ;
    end
else        % if acc>1
for t = 1:param_data.nb_tests

mask = zeros(param_data.Ni) ;
% define calibration zone if requiered
mask(param_data.Ni(1)/2 - floor(param_data.C(1)/2) + 1 : param_data.Ni(1)/2 + floor(param_data.C(1)/2) + 1, ...
     param_data.Ni(1)/2 - floor(param_data.C(1)/2) + 1 : param_data.Ni(1)/2 + floor(param_data.C(1)/2) + 1) = 1 ;
mask = mask(:) ;
% select randomly other frequencies
n = prod(param_data.C) ;
rng(t)  % initialize randomization
while n<param_data.M
ind = randperm(prod(param_data.Ni),param_data.M-n) ;
mask(ind) = 1 ;
n = sum(mask(:)) ;
end

Masks{t} = reshape(mask,param_data.Ni) ;

end
end

end