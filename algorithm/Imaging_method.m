function im_rec = Imaging_method(Y, param_data, sens_true, param_algo)
% MATLAB code for your imaging method
% INPUT:
%     Y: observed noisy undersampled MRI measurements
%     param_data: useful parameters related to the data
%     sens_true: original sensitivity maps assumed to be known
%     param_algo: algorithm parameters

%% Initialization
im_rec = param_algo.x0 ;
tic

%% Masking for sampling
%Create a mask to simulate undersampling

N = param_data.Ni(1)*param_data.Ni(2);

ind = find(param_data.mask==1);
M = numel(ind);

%Mask matrix (sparse matrix in matlab)
Ma = sparse(1:M, ind, ones(M, 1), M, N);


%% Sampling operator settings
% Measurement Operator
Phit2 = @(x) Ma*reshape(param_data.TF(x),N, param_data.Ncoils);
Phi2 = @(x) param_data.TFadj(reshape(Ma'*x,[param_data.Ni(1) param_data.Ni(2) param_data.Ncoils]));

%% Sparsity operator settings
%reconstruct the image solely from the undersampled data
conj_S=conj(sens_true);
xd_coils= conj_S.*Phi2(Y);
xd=sqrt(sum(xd_coils.^2,3));

%% Definition of the TV proximal operator
%tvbp

%Estimation of the regularization parameter
[alphad1, alphad2] = gradient_op(xd);
alphad = [alphad1(:);alphad2(:)];
lambda = 1e-2*norm(alphad,Inf);

%Parameters for the TV prox
paramtv.max_iter = param_algo.NbIt;
paramtv.rel_obj = 1e-3;
paramtv.verbose = 0;
delta = 1;
rel_tol = 1e-5;

%% Iterations
disp('Iterations begin');

%Initialisation
diff=zeros(param_algo.NbIt);
tv_x_initial = 0;


fidelity_term=0.5*sum(abs( Phit2(bsxfun(@times,sens_true,im_rec))-Y ).^2,1);
f_previous=lambda*tv_x_initial+sum(fidelity_term,1);



for it = 1:param_algo.NbIt
    
    %Project the reconstructed spatial image in the measurement domain, measure the
    %difference between the reconstructed and the measured data before
    %reprojecting it back into the spatial domain to guide the tv norm
    %minimization
    gradH_coils = bsxfun(@times, conj_S, Phi2(Phit2(bsxfun(@times,sens_true,im_rec)) - Y));
    gradH=sum(gradH_coils,3);
    [x,tv_x] = prox_tv(im_rec - delta*gradH,lambda,paramtv);
    
    %for each coil calculate the L2 norm distance between the reconstructed image
    %(projected into the Fourier domain) and the measured data
    fidelity_term=0.5*sum( abs( Phit2(bsxfun(@times,sens_true,im_rec))-Y ).^2, 1);
    f_current=lambda*tv_x+sum(fidelity_term,1);
    
    
    im_rec=x;
    diff(it)=abs(f_current - f_previous)/f_current;
    f_previous = f_current;
    
    
    if diff(it) <= rel_tol
        disp('--------------------------------------------')
        disp('--------------------------------------------')
        disp(strcat('Total number of iterations :',int2str(it)));
        disp(strcat('Total time elapsed: ',int2str(toc),28,'s'));
        disp('Reconstruction over')
        break
    end
    
    
    %Plot and display
    plot_disp(im_rec,diff,it,param_algo.NbIt)
    
end
end

function plot_disp(im_rec,diff,it,NbIt)

if mod(it,100) ==0
    disp('--------------------------------------------')
    disp(strcat('Iteration number: ',int2str(it),'/',28,int2str(NbIt)));
    fprintf('Relative difference with the last iteration: %d',diff(it)*100);
    disp(' %')
    disp(strcat('Total time elapsed: ',int2str(toc),28,'s'));
    
    figure(5),plot(diff(2:it),'--b'); hold on;
    xlabel('number of iterations')
    ylabel('relative difference between two consecutives iteration')
    title('The improvement of a MRI reconstruction depending on the number of iterations');
    hold off
    
    figure(6),imshow(abs(im_rec),[])
end


end