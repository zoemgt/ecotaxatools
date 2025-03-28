%create a color vector for plankton groups
%17 possibles plankton groups: 17 colors associated to one plankton groups 

plankton_group_colors = zeros(numel(plankton_groups_data_group), 3);

    for k = 1:numel(plankton_groups_data_group)    
    current_value = strtrim(plankton_groups_data_group{k});

    if isequal(current_value, 'bacillariophyta')
    plankton_group_colors(k,:) = [0.10,0.61,0.46];

    elseif isequal(current_value, 'chaetognatha')
    plankton_group_colors(k,:) = [0.69,0.84,0.88];

    elseif isequal(current_value, 'ciliophora')
    plankton_group_colors(k,:) = [0.79,0.77,1.00];

    elseif isequal(current_value, 'cnidaria')
    plankton_group_colors(k,:) = [0.41,0.60,0.78];
    
    elseif isequal(current_value, 'nanoplankton')
    plankton_group_colors(k,:) = [0.40,0.56,0.23];

    elseif isequal(current_value, 'copepoda')
    plankton_group_colors(k,:) = [0.91,0.56,0.42];

    elseif isequal(current_value, 'crustacea')
    plankton_group_colors(k,:) = [0.93,0.69,0.13]; 

    elseif isequal(current_value, 'cyanobacteria')
    plankton_group_colors(k,:) = [0.44,0.59,0.65];

    elseif isequal(current_value, 'detritus')
    plankton_group_colors(k,:) = [0.95,0.95,0.95];

    elseif isequal(current_value, 'dictyochophyceae')
    plankton_group_colors(k,:) = [0.89,0.90,0.59];

    elseif isequal(current_value, 'dinoflagellata')
    plankton_group_colors(k,:) = [0.71,0.83,0.58];

    elseif isequal(current_value, 'mollusca')
    plankton_group_colors(k,:) = [0.78,0.51,0.51];

    elseif isequal(current_value, 'other')
    plankton_group_colors(k,:) = [0.98,0.98,0.73];

    elseif isequal(current_value, 'unidentified')
    plankton_group_colors(k,:) = [1.00,0.93,0.81];

    elseif isequal(current_value, 'plastics')
    plankton_group_colors(k,:) = [0.00,0.45,0.74];

    elseif isequal(current_value, 'rhizaria')
    plankton_group_colors(k,:) = [0.98,0.83,0.92];

    elseif isequal(current_value, 'tunicata')
    plankton_group_colors(k,:) = [0.23,0.39,0.54];

    elseif isequal(current_value, 'ctenophora')
    plankton_group_colors(k,:) = [0.85, 0.9, 0.95];

    clear current_value

    end
    end

% Create palet: 

plankton_group_colors_palet = [
    0.10, 0.61, 0.46; % bacillariophyta
    0.69, 0.84, 0.88; % chaetognatha
    0.79, 0.77, 1.00; % ciliophora
    0.41, 0.60, 0.78; % cnidaria
    0.40, 0.56, 0.23; % nanoplankton
    0.91, 0.56, 0.42; % copepoda
    0.93, 0.69, 0.13; % crustacea
    0.44, 0.59, 0.65; % cyanobacteria
    0.95, 0.95, 0.95; % detritus
    0.89, 0.90, 0.59; % dictyochophyceae
    0.71, 0.83, 0.58; % dinoflagellata
    0.78, 0.51, 0.51; % mollusca
    0.98, 0.98, 0.73; % other
    1.00, 0.93, 0.81; % other unidentified
    0.00, 0.45, 0.74; % plastic
    0.98, 0.83, 0.92; % rhizaria
    0.23, 0.39, 0.54; % tunicata
    0.85, 0.90, 0.95; % ctenophora
    ];

%% Active in order to test differents colors 

% figure;
% color_vector = [0.85, 0.9, 0.95];
% plot(1, 1, 's', 'MarkerSize', 50, 'MarkerFaceColor', color_vector, 'MarkerEdgeColor', 'k');
% axis off;
