function im = create_image(param_data)
% create phantom image of size param_data.Ni
% with field of view param_data.FOV

%%
DefineBrain;
Brain.FOV = param_data.FOV*[1, 1];
im = RasterizePhantom(Brain,param_data.Ni);




end