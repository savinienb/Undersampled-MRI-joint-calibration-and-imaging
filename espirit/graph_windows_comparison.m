% close all
% clc
% clear all

addpath('./data')
addpath('./data/Generate_MRI_data')
addpath(genpath('./tools'))
addpath(genpath('./Espirit'));

%% Choose parameters for the simulations


Ncoils = 8 ;        % number of coils
acc = 4 ;           % acceleration parameter for the undersampling mask
Ni = 256 * [1,1] ;  % dimension of the image

nb_tests = 1 ;      % number of tests (i.e. noise realizations)

disp(' ')
disp('***************************************************')
disp('Run tests for ')
disp(['Ncoils = ',num2str(Ncoils)])
disp(['acc = ',num2str(acc)])
disp(['Ni = ',num2str(Ni(1)),'x',num2str(Ni(2))])
disp(['nb_tests = ',num2str(nb_tests)])
disp('***************************************************')


%% create data

C=[7 7; 15 15; 25 25; 35 35];
for run=3:size(C,1)
% useful parameters
param_data = parameters_data(Ni, Ncoils, acc, nb_tests) ;
param_data.C = C(run,:);
param_data.M=Ni(1).^2/4;

im_true = create_image(param_data) ;    % original image
sens_true = create_sensitivity(param_data) ;    

% create measurements
[Masks, param_data] = create_masks(param_data) ;
[DATA, Y, DATA_clean, param_data] = create_data(param_data, im_true, sens_true, Masks) ;

%% Auto-calibration imaging method
%  For each experiment t = 1:nb_tests
%  Find an estimate of im_true and sens_true from DATA{t}

% -------------------------------------------------
% Define here your algorithm parameters
param_algo.NbIt = 140 ; % number of iterations (example)
param_algo.x0 = zeros(Ni) ; % initialization (example)
param_algo.S0 = zeros(Ni(1),Ni(2),Ncoils) ; % initialization (example)
% -------------------------------------------------

im_rec = cell(nb_tests,1) ;
sens_rec = cell(nb_tests,1) ;
im_init= cell(nb_tests,1) ;
Time_im = zeros(nb_tests,1) ;
for t = 1:nb_tests
    param_data.mask = Masks{t} ; % selection mask for test t
    start_test = tic;
    [im_rec{t},sens_rec{t},im_init{t}] = Joint_cal_im_method_prefinal(Y{t}, param_data, param_algo) ;
    Time_im(t) = toc(start_test) ; % computation time for test t
end

%% ESPIRiT method
start_test=tic;
for i=1:nb_tests
    [ESP_sens{i} temp]=ESP(DATA,Ni,Ncoils,Masks,1);
    ESP_im{i}=temp{1};
end
Time_esp = toc(start_test) ; % computation time for test t

%% image comparaison

im_snr_value=zeros(2,nb_tests);
for i=1:nb_tests
    [psnr_temp,im_snr_value(1,i)]=psnr(histeq(mat2gray(abs(im_rec{i}))),histeq(mat2gray(im_true)));
    [psnr_temp,im_snr_value(2,i)]=psnr(histeq(mat2gray(abs(ESP_im{i}))),histeq(mat2gray(im_true)));
end

temp=find(im_snr_value(1,:)==max(im_snr_value(1,:)));
max_snr(1,:)=temp(1);
temp=find(im_snr_value(2,:)==max(im_snr_value(2,:)));
max_snr(2,:)=temp(1);

mean_snr_im(:,run)=mean(im_snr_value,2);

%% Sensitivity
mask_im=im_true>0;
se = strel('disk',5);
mask_im=imclose(mask_im,se);

best_rec=-realmax;
best_ESP=-realmax;
for i=1:numel(sens_rec)
    temp_sens_rec=bsxfun(@times,sens_rec{i},mask_im);
    temp_sens_ESP=bsxfun(@times,ESP_sens{i},mask_im);
    for coils=1:param_data.Ncoils
        [psnr_temp,sens_snr_value(1,coils,i)]=psnr((temp_sens_rec(:,:,coils)),(sens_true(:,:,coils)));
        [psnr_temp,sens_snr_value(2,coils,i)]=psnr((temp_sens_ESP(:,:,coils)),(sens_true(:,:,coils)));
    end
end



mean_snr_sens(:,run)=mean(mean(sens_snr_value,2),3);
%% Residual


sens_true2=bsxfun(@times,sens_true,mask_im);
SX_true=bsxfun(@times,sens_true2,im_true);

for i=1:nb_tests

sens_rec2=bsxfun(@times,sens_rec{i},mask_im);
SX_rec=bsxfun(@times,sens_rec2,im_rec{i});
SX_rec_l2(i)=immse(SX_rec,SX_true)

ESP_sens2=bsxfun(@times,ESP_sens{i},mask_im);
SX_ESP=bsxfun(@times,ESP_sens2,ESP_im{i});
SX_ESP_l2(i)=immse(SX_ESP,SX_true)
end

mean_l2_res(1,run)=mean(SX_rec_l2);
mean_l2_res(2,run)=mean(SX_ESP_l2);



end
%% Plot

figure
hold all
plot(C,abs(mean_snr_im(1,:)))
plot(C,abs(mean_snr_im(2,:)))
plot(C,abs(mean_snr_sens(1,:)))
plot(C,abs(mean_snr_sens(2,:)))

title('MRI Image and Sensitivity maps reconstruction comparaison between our method and ESPIRiT')
legend('mean images SNR with our method','mean  images SNR with ESPIRiT','mean sensitivities SNR with our method','mean sensitivities SNR with ESPIRiT')


figure
hold all
plot(C,mean_l2_res(1,:))
plot(C,mean_l2_res(2,:))

legend('mean average residual SX with L2 norm with our method','mean average residual SX with L2 norm with ESPIRIT')