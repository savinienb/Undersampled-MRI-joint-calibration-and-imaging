listing = dir('ESPIRiT');
listing(1:2)=[];
path_dir=listing.folder;
name={'/utils ' '/ESPIRiT_code /SPIRiT_code /nufft_files'}
% for i

addpath strcat(path_dir,/utils)
addpath ../ESPIRiT_code
addpath ../SPIRiT_code
addpath ../nufft_files
addpath ../data
addpath ../coilCompression_code
addpath ../SAKE_code


