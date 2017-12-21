function [DATA, Y, DATA_clean, param_data] = create_data(param_data, im_true, sens_true, Masks)
% create param_data.nb_tests sets of MRI measurements
% DATA_clean : full Fourier measurements without noise
% DATA : cell array of size param_data.nb_tests
%        each cell contains a cube of size param_data.Ni x param_data.Ncoils
%        corresponding to the noisy uncomplete Fourier measurements
% Y : cell array of size param_data.nb_tests
%     each cell contains a vector of size M x param_data.Ncoils, with
%     M = prod(param_data.Ni)/param_data.acc

%%
if ~isfield(param_data, 'input_snr') ; param_data.input_snr = 30 ; end ;

DATA = cell(param_data.nb_tests,1) ;
Y = cell(param_data.nb_tests,1) ;
DATA_clean = zeros(param_data.Ni(1),param_data.Ni(2),param_data.Ncoils) ; 
%%
sigma_noise = zeros(param_data.Ncoils,1) ;
for coil = 1:param_data.Ncoils
DATA_clean(:,:,coil) = param_data.TF(sens_true(:,:,coil) .* im_true) ;
sigma_noise(coil) = sum(abs(reshape(DATA_clean(:,:,coil),prod(param_data.Ni),1)).^2)/ numel(DATA_clean(:,:,coil)) * 10^(-param_data.input_snr/20) ;
end


for t = 1:param_data.nb_tests
DATA{t} = zeros(param_data.Ni(1),param_data.Ni(2),param_data.Ncoils) ; 
Y{t} = zeros(param_data.M,param_data.Ncoils) ;  
rng(t) % initialize randomization
for coil = 1:param_data.Ncoils
datacoil = DATA_clean(:,:,coil) ...
           + (randn(param_data.Ni) + 1i*randn(param_data.Ni))*sigma_noise(coil)/sqrt(2) ;
DATA{t}(:,:,coil) = datacoil .* Masks{t} ;
Y{t}(:,coil) = datacoil(Masks{t}>0) ; 
end
end


end