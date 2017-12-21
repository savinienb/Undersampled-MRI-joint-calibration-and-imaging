function sens_rec = Calibration_method(Y, param_data, im_true, param_algo)
% MATLAB code for your imaging method
% INPUT:
%     Y: observed noisy undersampled MRI measurements
%     param_data: useful parameters related to the data
%     im_true: original image assumed to be known
%     param_algo: algorithm parameters

%% Initialization
sens_rec = param_algo.S0 ;
tic

%% Masking for sampling

N = param_data.Ni(1)*param_data.Ni(2);
ind = find(param_data.mask==1);
M = numel(ind);
% Mask matrix (sparse matrix in matlab)
Ma = sparse(1:M, ind, ones(M, 1), M, N);


%% Sampling operator settings
% Measurement Operator
Phit = @(x) Ma*reshape(param_data.TF(x),N, param_data.Ncoils);
Phi = @(x) param_data.TFadj(reshape(Ma'*x,[param_data.Ni(1) param_data.Ni(2) param_data.Ncoils]));

%% Sparsity operator settings
xd = bsxfun(@times, im_true, Phi(Y));

%% Definition of the TV proximal operator
%tvbp

%Estimation of the regularization parameters
lambda=zeros(1,1,param_data.Ncoils);
for c=1:param_data.Ncoils
    [alphad1, alphad2] = gradient_op(xd(:,:,c));
    alphad = [alphad1(:);alphad2(:)];
    lambda(:,:,c) = 1e-2*norm(alphad,Inf);
end

%Parameters for the TV prox
paramtv.max_iter = param_algo.NbIt;
paramtv.rel_obj = 1e-3;
paramtv.verbose = 0;
delta = 1;
rel_tol = 1e-5;

%% Iterations
disp('Iterations begin');

diff=zeros(param_algo.NbIt,param_data.Ncoils);
tv_x_initial = 0;
f_previous = lambda*tv_x_initial + 0.5*norm(Phit(bsxfun(@times,sens_rec,im_true))-Y,2).^2;

for it=1:param_algo.NbIt
    
    
    %Use the reconstructed sensitivity to recalculate the measured image and project it into the measurement domain, measure the
    %difference between the reconstructed and the measured data before
    %reprojecting it back into the spatial domain to guide the tv norm
    %minimization
    gradH = bsxfun(@times, im_true, Phi(Phit(bsxfun(@times,sens_rec,im_true)) - Y));
    
    
    %for each coil minimize the tv norm
    tv_s=zeros(1,1,param_data.Ncoils);
    for i=1:param_data.Ncoils
        [sens_rec(:,:,i), tv_s(:,:,i)] = prox_tv(sens_rec(:,:,i) - delta*gradH(:,:,i),lambda(:,:,i),paramtv);
    end
    f_current = lambda.*tv_s + 0.5*norm(Phit(bsxfun(@times, sens_rec,im_true))-Y,2)^2;
    
    
    diff(it,:)=abs(f_current - f_previous)./f_current;
    f_previous = f_current;
    
    
    if sum(diff,3) <= rel_tol
        disp('--------------------------------------------')
        disp('--------------------------------------------')
        disp(strcat('Total number of iterations :',int2str(it)));
        disp(strcat('Total time elapsed: ',int2str(toc),28,'s'));
        disp('Reconstruction over')
        break
    end
    
    
    %Plot and display
    plot_disp(sens_rec,diff,it,param_algo.NbIt)
    
end
end


function plot_disp(sens_rec,diff,it,NbIt)

if mod(it,100) ==0
    disp('--------------------------------------------')
    disp(strcat('Iteration number: ',int2str(it),'/',28,int2str(NbIt)));
    fprintf('Relative mean difference with the last iteration: %d',mean(diff(it,:),2)*100);
    disp(' %')
    disp(strcat('Total time elapsed: ',int2str(toc),28,'s'));
    
    figure(7), hold on;
    plot(diff(2:it,1),'--r');
    plot(diff(2:it,2),'--g');
    plot(diff(2:it,3),'--b');
    plot(diff(2:it,4),'--m');
    legend('coil 1','coil 2','coil 3','coil 4')
    xlabel('number of iterations')
    ylabel('relative difference between two consecutives iteration')
    title('Improvement of a MRI coils sensitivity reconstruction depending on the number of iterations');
    hold off
    
    figure(8),
    subplot 211, imagesc(abs(reshape(sens_rec,size(sens_rec,1),[]))), axis image; colorbar, colormap bone
    xlabel('magnitude of sensitivity maps')
    subplot 212, imagesc(angle(reshape(sens_rec,size(sens_rec,1),[]))), axis image; colorbar, colormap bone
    xlabel('phase of sensitivity maps')
    
    drawnow
end


end
