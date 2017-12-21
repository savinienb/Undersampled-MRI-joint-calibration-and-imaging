function [im_rec, sens_rec,im_init] = Joint_cal_im_method(Y, param_data, param_algo)
% MATLAB code for your imaging method
% INPUT: 
%     Y: observed noisy undersampled MRI measurements
%     param_data: useful parameters related to the data
%     param_algo: algorithm parameters

%% Initialization

im_rec = param_algo.x0 ;
sens_rec = param_algo.S0 ;

%% Masking for sampling

N = param_data.Ni(1)*param_data.Ni(2);

ind = find(param_data.mask==1);
M = numel(ind);
%Mask matrix (sparse matrix in matlab)
Ma = sparse(1:M, ind, ones(M, 1), M, N);

%% Sampling operator settings
% Measurement Operator
Phit = @(x) Ma*reshape(param_data.TF(x),N, param_data.Ncoils);
Phi = @(x) param_data.TFadj(reshape(Ma'*x,[param_data.Ni(1) param_data.Ni(2) param_data.Ncoils]));


%% Initial estimates of Sensitivities and MRI Image
xs = Phi(Y);
sens_rec = xs;
xd = bsxfun(@times,conj(xs),Phi(Y));
im_rec = sqrt(sum(xd.^2,3));
im_init= sqrt(sum(xd.^2,3));
imshow(im_init,[]);

%% Regularization parameter for Image - lambda
[alphad1, alphad2] = gradient_op(im_rec);
alphad = [alphad1(:);alphad2(:)];
lambda = 1e-2*norm(alphad,Inf);

%% Regularization parameter for Sensivities - beta1 , beta2 , beta3 , beta4
for c=1:param_data.Ncoils
    [alphad1, alphad2] = gradient_op(sens_rec(:,:,c));
    alphad = [alphad1(:);alphad2(:)];
    beta(:,:,c) = 1e-2*norm(alphad,Inf);
end

delta1 = 1;
delta2 = 1;
paramtv.max_iter = param_algo.NbIt;
paramtv.rel_obj = 1e-6;
paramtv.verbose = 0;
%% Iterations
disp('Iterations begin');

f_previous=lambda*tv_norm(im_rec,0)+ 0.5*sum( abs( Phit(bsxfun(@times,sens_rec,im_rec))-Y ).^2, 1);
for i=1:param_data.Ncoils
f_previous = f_previous + beta(:,:,i)*tv_norm(sens_rec(:,:,i),0);
end


for it = 1:param_algo.NbIt
    
    gradH_im = bsxfun(@times, conj(sens_rec), Phi(Phit(bsxfun(@times,sens_rec,im_rec)) - Y));
    
    sum_gradH_im = sum(gradH_im,3);
    
    [im_rec,tv_im] = prox_tv(im_rec - delta1*sum_gradH_im,lambda,paramtv);
    
    gradH_S =  bsxfun(@times, im_rec, Phi(Phit(bsxfun(@times,sens_rec,im_rec)) - Y));
    for i=1:param_data.Ncoils
        [sens_rec(:,:,i), tv_s(:,:,i)] = prox_tv(sens_rec(:,:,i) - delta2*gradH_S(:,:,i),beta(:,:,i),paramtv);
    end
    
    f_current = lambda*tv_im + sum(beta.*tv_s,3) + 0.5*sum( abs( Phit(bsxfun(@times,sens_rec,im_rec))-Y ).^2, 1);
    
    im_rec(isnan(im_rec)) ;
    diff = abs(f_current - f_previous)/f_current;
    f_previous = f_current;
    
    if mod(it,10) ==0
        fprintf('Iteration number %d\n',it);
%         figure(10),imshow(im_rec,[]);
%         figure(12),
%         subplot 211, imagesc(abs(reshape(sens_rec,size(sens_rec,1),[]))), axis image; colorbar, colormap bone
%         xlabel('magnitude of sensitivity maps')
%         subplot 212, imagesc(angle(reshape(sens_rec,size(sens_rec,1),[]))), axis image; colorbar, colormap bone
%         xlabel('phase of sensitivity maps')
        
%        drawnow
    end
    
end


end