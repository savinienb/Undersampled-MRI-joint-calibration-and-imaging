close all
clc
clear all

addpath('./data')
addpath('./data/Generate_MRI_data')
addpath('./tools')
addpath('./algorithm')

%% Choose parameters for the simulations


Ncoils = 4 ;        % number of coils
acc = 3 ;           % acceleration parameter for the undersampling mask
Ni = 256 * [1,1] ;  % dimension of the image

nb_tests = 2 ;      % number of tests (i.e. noise realizations)

disp(' ')
disp('***************************************************')
disp('Run tests for ')
disp(['Ncoils = ',num2str(Ncoils)]) 
disp(['acc = ',num2str(acc)]) 
disp(['Ni = ',num2str(Ni(1)),'x',num2str(Ni(2))]) 
disp(['nb_tests = ',num2str(nb_tests)]) 
disp('***************************************************')


%% create data

% useful parameters 
param_data = parameters_data(Ni, Ncoils, acc, nb_tests) ;

%param_data.M = N^2/2;

im_true = create_image(param_data) ;    % original image
figure, imagesc(im_true), axis image; colorbar, colormap bone, xlabel('original phantom image')

sens_true = create_sensitivity(param_data) ;    %original sensitivity maps
figure, 
subplot 211, imagesc(abs(reshape(sens_true,Ni(1), Ni(2)*Ncoils))), axis image; colorbar, colormap bone
xlabel('magnitude of sensitivity maps')
subplot 212, imagesc(angle(reshape(sens_true,Ni(1), Ni(2)*Ncoils))), axis image; colorbar, colormap bone
xlabel('phase of sensitivity maps')


% create measurements 
[Masks, param_data] = create_masks(param_data) ;    
[DATA, Y, DATA_clean, param_data] = create_data(param_data, im_true, sens_true, Masks) ;

%% Auto-calibration imaging method
%  For each experiment t = 1:nb_tests
%  Find an estimate of im_true and sens_true from DATA{t}

% -------------------------------------------------
% Define here your algorithm parameters
param_algo.NbIt = 1000 ; % number of iterations (example)
param_algo.x0 = zeros(Ni) ; % initialization (example)
param_algo.S0 = zeros(Ni(1),Ni(2),Ncoils) ; % initialization (example)
% -------------------------------------------------

im_rec = cell(nb_tests,1) ;
sens_rec = cell(nb_tests,1) ;
Time_im = zeros(nb_tests,1) ;
for t = 1:nb_tests
param_data.mask = Masks{t} ; % selection mask for test t
start_test = tic;
[im_rec{t},sens_rec{t}] = Joint_cal_im_method_prefinal(Y{t}, param_data, param_algo,im_true,sens_true) ; 
Time_im(t) = toc(start_test) ; % computation time for test t
end
