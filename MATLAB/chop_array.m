function chopped_array = chop_array(original_array)

% Because the profiles are different lengths, when the mean is taken, there
% are points when the mean profile will begin to decrease. We don't want
% anything after that decrease point, so this script chops off that end
% part

% Creates a difference array
d=sign([0; diff(original_array)]);

% Finds the negative 
idx=strfind(d',-1);

% Chops!
chopped_array = original_array(1:min(idx)-1);

end