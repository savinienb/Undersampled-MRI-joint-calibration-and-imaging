function [ ESP_sens,ESP_im] = ESP(DATA,Ni,Ncoils,Masks,nb_tests)
%ESPIRIT Summary of this function goes here
%   Detailed explanation goes here

%  In this section, you need to introduce the useful functions to estimate
%  the image and the sensitivity maps using the ESPIRiT method
%  The functions you need are available in the ESPIRiT package available at
%  https://people.eecs.berkeley.edu/~mlustig/Software.html
%  An example of the joint calibration and imaging ESPIRiT method is given
%  in the file demo_ESPIRiT_L1_recon.m of this package.
ncalib = cell(nb_tests,1) ;
DATAc = cell(nb_tests,1) ;
calib = cell(nb_tests,1) ;
resL1ESPIRiT = cell(nb_tests,1) ;

ksize = [6,6]; % ESPIRiT kernel-window-size
eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.9; % threshold of eigenvectors in image space

% parameters for L1-reconstruction with splitting
nIterCG = 5;       % number of CG iterations for the PI part
nIterSplit = 15;    % number of splitting iterations for CS part
splitWeight = 0.4;  % reasonable value
lambda = 0.0025;    % L1-Wavelet threshold

for i = 1:nb_tests
    ncalib{i} = getCalibSize(Masks{i});
    mask = repmat(Masks{i},[1,1,Ncoils]);
    DATAc{i} = DATA{i}.*mask;
    calib{i} = crop(DATAc{i},[ncalib{i},Ncoils]);
    
    % Compute eigen value maps
    [k,S] = dat2Kernel(calib{i},ksize);
    idx = max(find(S >= S(1)*eigThresh_k));
    
    % crop kernels and compute eigen-value decomposition in image space to get
    % maps
    [M,W] = kernelEig(k(:,:,:,1:idx),[Ni(1),Ni(2)]);
    
    % Compute Soft-SENSE ESPIRiT Maps
    % crop sensitivity maps according to eigenvalues==1. Note that we have to
    % use 2 sets of maps. Here we weight the 2 maps with the eigen-values
    maps = M(:,:,:,end-1:end);
    % Weight the eigenvectors with soft-senses eigen-values
    weights = W(:,:,end-1:end) ;
    weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
    weights = -cos(pi*weights)/2 + 1/2;
    
    % create and ESPIRiT operator
    ESP = ESPIRiT(maps,weights);
    
    % Reconstruction
    XOP = Wavelet('Daubechies_TI',4,6);
    FT = p2DFT(mask,[Ni(1),Ni(2),Ncoils]);
    
    disp('Performing L1-ESPIRiT reconstruction from 2 maps')
    [resL1ESPIRiT] = cgL1ESPIRiT(DATAc{i}, zeros(Ni(1),Ni(2),2), FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);
    
    %store output
    ESP_im{i}=sos(resL1ESPIRiT);
    ESP_sens=maps(:,:,:,2);
    
end


end

