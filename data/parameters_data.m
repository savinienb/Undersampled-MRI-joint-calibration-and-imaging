function param_data = parameters_data(Ni, Ncoils, acc, nb_tests) 

%%
param_data.Ni = Ni ;
param_data.Ncoils = Ncoils ;
param_data.acc = acc ;
param_data.nb_tests = nb_tests ;
param_data.M = floor(prod(param_data.Ni) / param_data.acc) ;   %number of measurements

%%
param_data.FOV = 0.28; % FOV width
param_data.pixelsize = param_data.FOV./Ni;

%%
param_data.TF =@(x) (1/sqrt(prod(Ni))) * fftshift(fft2(ifftshift( x ))) ;
param_data.TFadj =@(y) sqrt(prod(param_data.Ni)) * fftshift(ifft2(ifftshift( y ))) ;


end