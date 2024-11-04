
close all;

root_path = fileparts(mfilename('fullpath'));

% Add KFOSU folder to the root path
addpath(fullfile(root_path, 'KFOSU'));

% Add VCAandSUNSAL folder to the root path
addpath(fullfile(root_path, 'VCAandSUNSAL'));

% Add data folder to the root path
data_path = fullfile(root_path, 'data');

%-----------------------------------------------------------------------

% Load data

% Pure spectra matrix
pureSpectra_matrix_path = fullfile(data_path, 'pureSpectra_matrix.mat');
pureSpectra_matrix = load(pureSpectra_matrix_path ).w;

% Concentration matrix
concentration_path = fullfile(data_path, 'concentration_matrix_SubspaceTracking.mat');
concentration_matrix  = load(concentration_path).H;


% Set SubspaceTracking to 1, if simulation with 'concentration_matrix_SubspaceTracking.mat';
% set to 0 otherwise
SubspaceTracking=1;


    
% Number of pure spectra
  nb_pure_spectra=5;

% Signal-to-noise ratio (SNR) in dB
  SNR = 20; 
 
% Noiseless spectra 
  Y0 = concentration_matrix*pureSpectra_matrix(1: nb_pure_spectra,:)  ; 
    
% Noisy spectra
  [N,nb_spectral_channels] = size(Y0);
  variance = sum(Y0(:).^2)/10^(SNR/10)/N/nb_spectral_channels;
  Noise = sqrt(variance).*randn([N nb_spectral_channels]);
  Y = Y0 + Noise;
    
% Apply physical constraint: Y>=0
  Y = max(0,Y); 
    
   
% Randomize the order of spectra each time the script is executed
if ~SubspaceTracking
    Y = Y(randperm(N),:);
end


% KFOSU parameters, 7 out of many possible options;
% refer to KFOSU.m for details
nb_regression_spectra = 30;
period = 20;
confidence_level = "off";

params.nb_regression_spectra = nb_regression_spectra;
params.period = period; % For plot refresh
params.SubspaceTracking = "off";
params.confidence_level = confidence_level;
params.dimensionality_reduction_method = "DFT";
params.reconstruction_energy_percent = 0.85;  
params.uncertainty_level = 1;


% Ensure that N is a multiple of period when confidence_level is "on"
if  strcmp(confidence_level,"on") && mod(N,period)~=0
    error("N is not a multiple of period");
end


% Create the figure 
fig = figure(1);
text(0.5, 0.5, 'Waiting for update...', 'HorizontalAlignment', 'center', ...
     'FontSize', 12, 'Color', 'k', 'Units', 'normalized');

% 
x_tick_values = 1:nb_spectral_channels;
fill_color = [0.8, 0.2, 0.2];
face_alpha = 0.5;

exit_button = uicontrol('Style', 'pushbutton', 'String', 'Exit', ...
    'Position', [15 20 60 30], ... 
    'Callback', @(src, event) close(fig)); 


Y = Y.';
pureSpectra_matrix = pureSpectra_matrix.';
% Update pure spectra estimates (spectrum-by-spectrum basis)

for t = 1:N

    % Check if the exit button was pressed
    if ~isvalid(fig)
        break
    end
    
    % New observation: time index t
    y_t = Y(:, t);

    % Pure spectra estimates 
    %tic
    [S_t, uncertainties] = KFOSU(t, y_t, nb_pure_spectra, nb_spectral_channels, params);
    %toc  
    
    if t >= nb_regression_spectra &&  mod(t, period) == 0 
        match_indices=pure_spectra_matching(pureSpectra_matrix,S_t,nb_pure_spectra);
        if ~isempty(uncertainties)
            % 3-sigma confidence level
            sigma_matrix = sqrt(uncertainties);
            S_t_plus_confidence = S_t + 3*sigma_matrix;
            S_t_minus_confidence = S_t - 3*sigma_matrix;
        end
        
        if isvalid(fig)
            for i = 1:nb_pure_spectra

                fig_estimate = subplot(nb_pure_spectra, 2, 2*i-1);
                % Clear the current plot in the subplot
                cla;  
                plot(S_t(:, i), 'r', 'LineWidth', 1);
                
                if ~isempty(uncertainties)
                    hold on
                    fill([x_tick_values, fliplr(x_tick_values)], [S_t_plus_confidence(:,i).',...
                    fliplr(S_t_minus_confidence(:,i).')], fill_color, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');
                end

                hold on
                fig_true = subplot(nb_pure_spectra, 2, 2*i);
                cla; 
                plot(pureSpectra_matrix(:, match_indices(i)), 'b','LineWidth', 1);
                
                if i == 1
                    title(fig_estimate,['Pure spectra estimate  time index ', num2str(t), '/', num2str(N)]);
                    title(fig_true,'True pure spectra');
                end

                if i ~= nb_pure_spectra
                    set(fig_estimate, 'XTick', []);
                    set(fig_true, 'XTick', []);
                else
                    xlabel(fig_estimate,'Spectral channels');
                    xlabel(fig_true,'Spectral channels');
                    
                end
                hold off;
                xlim([1 nb_spectral_channels]);
            end
            
            drawnow;
        else
            % Exit if the figure was closed
            break;  
        end
    end
end


%% Pure spectra matching 

function match_indices = pure_spectra_matching(S_true, S_estimate, nb_pure_spectra)

    % Match each estimated pure spectrum with the true one
    % Inputs:
    %   S_true: True pure spectra matrix
    %   S_estimate: Estimated pure spectra matrix
    %   nb_pure_spectra: Number of pure spectra
    
    % Output:
    %   match_indices: indices of the best matches


    % Normalize (l2 norm) each pure spectrum in S_true and S_estimate
    S_true_normalized = S_true ./ vecnorm(S_true);  
    S_estimate_normalized = S_estimate ./ vecnorm(S_estimate);  

    % Compute the cosine similarity matrix
    cosine_similarity_matrix = S_true_normalized.' * S_estimate_normalized;

    % Find the best match for each column (pure spectrum) of S_estimate in S_true
    match_indices = zeros(1, nb_pure_spectra);  % Initialize output vector
    for i = 1:nb_pure_spectra
        % Find the maximum cosine similarity for the current column in S_estimate
        [~, max_idx] = max(cosine_similarity_matrix(:, i));
        match_indices(i) = max_idx;

        % Set the corresponding row and column in the similarity matrix to -1
        % to ensure unique matching
        cosine_similarity_matrix(max_idx, :) = -1;
    end
end
