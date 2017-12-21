function sensitivity = create_sensitivity(param_data)
% create param_data.Ncoils sensitivity maps 
% of size param_data.Ni
% with field of view param_data.FOV

%% create sensitivity maps from Guerquin-Kern codes
sensitivity = GenerateSensitivityMap( param_data.FOV, param_data.pixelsize, param_data.Ncoils, .07, .3);



end