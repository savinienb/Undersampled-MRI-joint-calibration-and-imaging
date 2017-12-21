close all
clc
clear all

addpath('./tools')
addpath(genpath('./data'));
addpath(genpath('./ESPIRiT'));
%This file show and measure, for ESPIRiT and our method, the difference
%between : reconstituted MRI images, the reconstitued sensitivity map, the residual of each method and finally, their run times

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

% useful parameters
param_data = parameters_data(Ni, Ncoils, acc, nb_tests) ;
param_data.C = [35,35];
param_data.M=Ni(1).^2/4;

im_true = create_image(param_data) ;    % original image
sens_true = create_sensitivity(param_data) ;    %original sensitivity maps

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
    [im_rec{t},sens_rec{t},im_init{t}] = Joint_cal_im_method(Y{t}, param_data, param_algo) ;
    Time_rec(t) = toc(start_test) ; % computation time for test t
    t
    disp('------------------')
end

%% ESPIRiT method
for i=1:nb_tests
    start_test=tic;
    [ESP_sens{i} temp]=ESP(DATA,Ni,Ncoils,Masks,1);
    ESP_im{i}=temp{1};
    Time_esp(i) = toc(start_test) ; % computation time for test t
end

%% image comparison

%calcul the SNR between each reconstitued images (for both algorithm) and
%the original image for the n iteration
im_snr_value=zeros(2,nb_tests);
for i=1:nb_tests
    [psnr_temp,im_snr_value(1,i)]=psnr(histeq(mat2gray(abs(im_rec{i}))),histeq(mat2gray(im_true)));
    [psnr_temp,im_snr_value(2,i)]=psnr(histeq(mat2gray(abs(ESP_im{i}))),histeq(mat2gray(im_true)));
end
%Get the best image out of the 10 iterations
temp=find(im_snr_value(1,:)==max(im_snr_value(1,:)));
max_snr(1,:)=temp(1);
temp=find(im_snr_value(2,:)==max(im_snr_value(2,:)));
max_snr(2,:)=temp(1);

mean_snr_im(:)=mean(im_snr_value,2);

figure
subplot(3,2,1)
imshow((mat2gray(im_true)))
title('Groundtruth image')
subplot(3,2,2)
imshow((mat2gray(im_init{max_snr(1,1)})))
[psnr_t,snr_t]=psnr(histeq(mat2gray(im_init{max_snr(1,1)})),histeq(mat2gray(im_true)));
title(strcat('X0 initialisation image, SNR : ',num2str(snr_t)))
subplot(3,2,3)
imshow((mat2gray(abs(ESP_im{max_snr(2,1)}))))
[psnr_t,snr_t]=psnr(histeq(mat2gray(ESP_im{max_snr(2,1)})),histeq(mat2gray(im_true)));
title(strcat('ESPIRiT reconstruction, SNR : ',num2str(snr_t)))
subplot(3,2,4)
imshow((mat2gray(abs(im_rec{max_snr(1,1)}))))
[psnr_t,snr_t]=psnr(histeq(mat2gray(abs(im_rec{max_snr(1,1)}))),histeq(mat2gray(im_true)));
title(strcat('Our reconstruction, SNR : ',num2str(snr_t)))
subplot(3,2,5)
imshow(histeq(mat2gray(im_true))-histeq(mat2gray(ESP_im{max_snr(2,1)})),[])
title('Normalized histogram equalized  difference between the ESPIRiT best image result ')
subplot(3,2,6)
imshow(histeq(mat2gray(im_true))-histeq(mat2gray(abs(im_rec{1}))),[])
title('Normalized histogram equalized  difference between our method best image result ')


%% Sensitivity

%create the brain mask
mask_im=im_true>0;
se = strel('disk',5);
mask_im=imclose(mask_im,se);

best_rec=-realmax;
best_ESP=-realmax;
for i=1:nb_tests
    %apply the mask to each coil
    temp_sens_rec=bsxfun(@times,sens_rec{i},mask_im);
    temp_sens_ESP=bsxfun(@times,ESP_sens{i},mask_im);
    
    %calcul the snr for each coil for each method
    for coils=1:param_data.Ncoils
        [psnr_temp,sens_snr_value(1,coils,i)]=psnr(abs(temp_sens_rec(:,:,coils)),abs(sens_true(:,:,coils)));
        [psnr_temp,sens_snr_value(2,coils,i)]=psnr(abs(temp_sens_ESP(:,:,coils)),abs(sens_true(:,:,coils)));
    end
    
    mean_rec=mean(sens_snr_value(1,:,i));
    mean_ESP=mean(sens_snr_value(2,:,i));
    if mean_rec>best_rec
        best_rec=mean_rec;
        rec_ind=i;
    end
    if mean_ESP>best_ESP
        best_ESP=mean_ESP;
        ESP_ind=i;
    end
end

sens_true2=bsxfun(@times,sens_true,mask_im);
sens_rec2=bsxfun(@times,sens_rec{rec_ind},mask_im);
ESP_sens2=bsxfun(@times,ESP_sens{ESP_ind},mask_im);
snr_rec_disp=sens_snr_value(1,:,rec_ind);
snr_ESP_disp=sens_snr_value(2,:,ESP_ind);
figure
for i=1:param_data.Ncoils
    subplot(3,param_data.Ncoils,(i))
    imshow(abs(sens_true2(:,:,i)),[]);
    if i==4
        title('                                 True sensitivity map')
    end
    subplot(3,param_data.Ncoils,i+param_data.Ncoils)
    imshow(abs(sens_rec2(:,:,i)),[]);
    xlabel(strcat('SNR:',num2str(snr_rec_disp(1,i))))
    if i==4
        title(strcat('                                 Best sensitivity map (10 iterations) obtained with our joint calibration method, mean SNR= ',num2str(mean(snr_rec_disp(1,:)))))
    end
    subplot(3,param_data.Ncoils,i+param_data.Ncoils*2)
    imshow(abs(ESP_sens2(:,:,i)),[]);
    xlabel(strcat('SNR:',num2str(snr_ESP_disp(1,i))))
    if i==4
        title(strcat('                                 Best sensitivity map (10 iterations) obtained with ESPIRiT joint calibration method, mean SNR= ',num2str(mean(snr_ESP_disp(1,:)))));
    end
end


%% Residual



sens_true2=bsxfun(@times,sens_true,mask_im);
SX_true=bsxfun(@times,sens_true2,im_true);

sens_rec2=bsxfun(@times,sens_rec{rec_ind},mask_im);
SX_rec=bsxfun(@times,sens_rec2,im_rec{rec_ind});

ESP_sens2=bsxfun(@times,ESP_sens{ESP_ind},mask_im);
SX_ESP=bsxfun(@times,ESP_sens2,ESP_im{ESP_ind});

res_rec=sum(abs(SX_true-SX_rec));
res_ESP=sum(abs(SX_true-SX_ESP));


figure
for i=1:param_data.Ncoils
    subplot(2,param_data.Ncoils,i)
    imshow(abs(SX_rec(:,:,i)-SX_true(:,:,i)),[]);
    if i==4
        title(strcat('                                 SX residual with our method, mean residual difference= ',num2str(mean(res_rec(:)))))
    end
    subplot(2,param_data.Ncoils,i+param_data.Ncoils)
    imshow(abs(SX_ESP(:,:,i)-SX_true(:,:,i)),[]);
    if i==4
        title(strcat('                                 SX residual with ESPIRIT, mean residual difference= ',num2str(mean(res_ESP(:)))));
    end
    
end

%% MEAN VALUES

%time
mean_time_rec=mean(Time_rec)
std_time_rec=std(Time_rec)
mean_time_ESP=mean(Time_esp)
std_time_ESP=std(Time_esp)

%quality
mean_snr_im_rec=mean(mean(im_snr_value(1,:,:),2),3)
std_snr_im_rec=std(mean(im_snr_value(1,:,:),2))
mean_snr_sens_rec=mean(mean(im_snr_value(1,:,:),2),3)
std_snr_sens_rec=std(mean(sens_snr_value(1,:,:),2))


mean_snr_im_ESP=mean(mean(im_snr_value(2,:,:),2),3)
std_snr_im_ESP=std(mean(im_snr_value(2,:,:),2))
mean_snr_sens_ESP=mean(mean(sens_snr_value(2,:,:),2),3)
std_snr_ESP=std(mean(sens_snr_value(2,:,:),2))