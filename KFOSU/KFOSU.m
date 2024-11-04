
function [S_t, uncertainties] = KFOSU(t, y_t, nb_pure_spectra, nb_spectral_channels, params)


    % Update pure spectra estimates using the latest acquired spectrum
    % Inputs:
    %   t: time index
    %   y_t: new spectrum (data space)
    %   nb_pure_spectra: number of pure spectra
    %   nb_spectral_channels: number of spectral channels
    
    %   params, optional Input parameters:
    %     nb_regression_spectra: number of regressors
    %     dimensionality_reduction_method: "DFT" or "PCA"
    %     reconstruction_energy_percent: desired percentage of total energy (0 < value <= 1) used in dimensionality reduction
    %     uncertainty_level: process noise variance (> 0)
    %     SubspaceTracking: "on" or "off"
    %     probability_false_alarm: sensitivity control for archetype spectra update (0 <= value <= 1)
    %     adjust_val: sensitivity control for archetype spectra update (> 0)
    %     regression_rho: step size for solving the regression problem
    %     regression_nb_iter: number of iterations for solving the regression problem
    %     confidence_level: "on" or "off"; "on" -> return uncertainties associated with the estimates
    %     period: interval between uncertainties calculations
    %     sgolay_order: polynomial order for Savitzky-Golay filter, used for noise variance estimation
    %     sgolay_frame: frame length for Savitzky-Golay filter, used for noise variance estimation
    %     Q: window size, used for noise variance estimation
    %     concentration_positivity: ct = argmin|| y_t - S_t_minus1*c|| s.t. c >=0
    %     concentration_positivity_sum_one: ct = argmin|| y_t - S_t_minus1*c|| s.t. c >=0, c_transpose*1=1

    % Ouputs:
    %   S_t: pure spectra estimate at time index t
    %   uncertainties: uncertainties associated with the estimates 
 

    % Declare persistent variables: retained for the next iteration, t+1
    persistent S_t_minus1 S_t_reducedminus1 sigmatminus1 Y_regression covMat_reduced...
    R I_L  selected_frequencies Y_regression_reduced PCA_V mu_init S_init...
    nb_regression_spectra dimensionality_reduction_method reconstruction_energy_percent... 
    uncertainty_level SubspaceTracking probability_false_alarm adjust_val regression_rho...
    regression_nb_iter period confidence_level sgolay_order sgolay_frame Q...
    concentration_positivity concentration_positivity_sum_one dim noise_var;
 

    % Initialize persistent variables
    if t==1
        [S_t_reducedminus1, sigmatminus1, covMat_reduced,R,I_L,selected_frequencies,...
        Y_regression_reduced,PCA_V] = deal([]);
        [S_t_minus1, mu_init, S_init] = deal(zeros(nb_spectral_channels, nb_pure_spectra));
        Y_regression = zeros(nb_spectral_channels, params.nb_regression_spectra);

    
    % Initialize params with default values
    defaultParams.nb_regression_spectra = 30;
    defaultParams.dimensionality_reduction_method = "DFT";
    defaultParams.reconstruction_energy_percent = 0.95;      
    defaultParams.uncertainty_level = 1;
    defaultParams.SubspaceTracking = "off";
    defaultParams.probability_false_alarm = 0.05; % Only used if SubspaceTracking is "on"
    defaultParams.adjust_val = 15; % Only used if SubspaceTracking is "on"
    defaultParams.regression_rho = 1; 
    defaultParams.regression_nb_iter = 50; 
    defaultParams.period = 20; 
    defaultParams.confidence_level = "off";
    defaultParams.sgolay_order=3; 
    defaultParams.sgolay_frame = 5;
    defaultParams.Q = 10;
    defaultParams.concentration_positivity = "off";
    defaultParams.concentration_positivity_sum_one = "on";

    % Overwrite default values with those provided in 'params'
    if nargin > 4
        fnames = fieldnames(params);
        for i = 1:length(fnames)
            if isfield(defaultParams, fnames{i})
                defaultParams.(fnames{i}) = params.(fnames{i});
            end
        end
    end


    % Extract the parameters
    nb_regression_spectra = defaultParams.nb_regression_spectra;
    dimensionality_reduction_method = defaultParams.dimensionality_reduction_method;
    reconstruction_energy_percent = defaultParams.reconstruction_energy_percent;
    uncertainty_level = defaultParams.uncertainty_level;
    SubspaceTracking = defaultParams.SubspaceTracking;
    probability_false_alarm = defaultParams.probability_false_alarm;
    adjust_val = defaultParams.adjust_val;
    regression_rho = defaultParams.regression_rho;
    regression_nb_iter = defaultParams.regression_nb_iter;
    period = defaultParams.period;
    confidence_level = defaultParams.confidence_level;
    sgolay_order = defaultParams.sgolay_order;
    sgolay_frame = defaultParams.sgolay_frame;
    Q = defaultParams.Q;
    concentration_positivity = defaultParams.concentration_positivity;
    concentration_positivity_sum_one = defaultParams.concentration_positivity_sum_one;
    
    % Ensure concentration_positivity or concentration_positivity_sum_one is "on"
    if concentration_positivity==concentration_positivity_sum_one
        error("set concentration_positivity or concentration_positivity_sum_one to on");
    end
    [S_t, uncertainties] = deal([]);
    return;
    end

    % Initialize archetype spectra matrix and other variables
    if t<=nb_regression_spectra

        Y_regression(:,t) = y_t;
        [S_t, uncertainties] = deal([]);
    end

    if t==nb_regression_spectra
        
        % Initialize pure spectra estimates
        [S_t, ~, ~] = VCA(Y_regression, 'Endmembers',nb_pure_spectra);
        S_t_minus1 = S_t;

        % Estimate the observation noise variance 
        noise_var = estimate_noise_variance(Y_regression, sgolay_order, sgolay_frame, Q);

        % Prevent instability at very high SNR levels
        minVal = 1e-4* max(Y_regression(:));
        if noise_var < minVal 
            noise_var = minVal;
        end

        % DR operator for PCA/DFT
        if strcmp(dimensionality_reduction_method, "PCA")
            center = sum(Y_regression, 2)/nb_regression_spectra; 
            Y_regression_centered = Y_regression - center;
            [U,~,~]=svd(Y_regression_centered*Y_regression_centered.'/(nb_regression_spectra-1));
    
            % Compute reduced space dimension
            total_energy =  sum(Y_regression(:).^2);
            for reduced_dim = nb_pure_spectra:nb_spectral_channels 
                PCA_V = U(:,1:reduced_dim);
                Yproj = PCA_V*PCA_V.'*Y_regression; 
                energy = sum(Yproj(:).^2);
                if energy >=reconstruction_energy_percent*total_energy
                    dim = reduced_dim;
                    break
                end
            end
        else
            % Select frequency indices to retain (DFT method)
            selected_frequencies = select_frequencies(Y_regression, reconstruction_energy_percent,nb_spectral_channels);
            dim = 2*length(selected_frequencies);
        end


        % Extract linearly independent columns to prevent numerical instability
        % Compute reduced version of Y_regression and S_t_minus1
        
        rankY = rank(Y_regression);
        if rankY ~= nb_regression_spectra
            % Perform QR decomposition
            [~, ~, E] = qr(Y_regression, 'vector');
            Y_regression = Y_regression(:,E(1:rankY));
         end
    
         % Reduced version of Y_regression and S_t_minus1
         if strcmp(dimensionality_reduction_method, "PCA")
             Y_regression_reduced = PCA_V.'*Y_regression;  
             S_t_reduced = S_t_minus1.'*PCA_V;
         else
             Y_regression_reduced =  fft(Y_regression);
             Y_regression_reduced =  Y_regression_reduced(selected_frequencies,:);
             Y_regression_reduced = [real(Y_regression_reduced); imag(Y_regression_reduced)];
    
             S_t_reduced=  fft(S_t);
             S_t_reduced =  S_t_reduced(selected_frequencies,:);
             S_t_reduced= [real(S_t_reduced); imag(S_t_reduced)].';
          end
    
           % Vectorize S_t_reduced
           S_t_reducedminus1 = S_t_reduced(:);
    
           % Initialize covariance matrix (reduced space) and other variables
           covMat_reduced = eye(nb_pure_spectra*dim)*uncertainty_level;
           sigmatminus1 = covMat_reduced;
           R = eye(dim)*noise_var; 
           I_L = eye(nb_pure_spectra*dim);

    uncertainties = [];
    return;  
    end


    
    % Pure spectra update
    if t>nb_regression_spectra
      
        if strcmp(dimensionality_reduction_method, "DFT") 
            y_t_reduced = fft(y_t);
            y_t_reduced = y_t_reduced (selected_frequencies,:);
            y_t_reduced = [real(y_t_reduced); imag(y_t_reduced)]; 
        else
            y_t_reduced = PCA_V.'*y_t;
        end

        % Check for archetype spectra update
        if strcmp(SubspaceTracking , "on")
            [Y_regression, Y_regression_reduced] = archetype_spectra_update(Y_regression, Y_regression_reduced, y_t,y_t_reduced, nb_spectral_channels, R(1,1), adjust_val, probability_false_alarm);
        end
        
        % Estimate concentration 
        if strcmp(concentration_positivity_sum_one,"on")
            ct = sunsal(S_t_minus1, y_t).';
        else
            ct =  sunsal(S_t_minus1, y_t,'ADDONE','no').';
        end

        % Prediction
        cov_predict = sigmatminus1 +  covMat_reduced;
        Ht = zeros(dim, dim*nb_pure_spectra);
        for i=1:dim
            Ht(i, (i-1)*nb_pure_spectra + 1:i*nb_pure_spectra) = ct;
        end
       
        % Correction
        rt = y_t_reduced  - Ht*S_t_reducedminus1;
        Zt = Ht*cov_predict*(Ht.') + R;
        Kt = (cov_predict*Ht.')/Zt;

        % Compute pure spectra estimates (unconstrained, lower-dimensional/reduced space)
        S_t_reduced = S_t_reducedminus1 + Kt*rt;

        % Unvectorize S_t_reduced
        S_t_reduced_unvec =  reshape(S_t_reduced ,[nb_pure_spectra,dim]).'; 
    
        % Compute the regression matrix, S_t contains the constrained pure spectra estimates (data space) 
        [regression_matrix, S_t] = Regression(regression_rho,regression_nb_iter, mu_init, S_init,...
        Y_regression_reduced, S_t_reduced_unvec, Y_regression);
        S_t_minus1=S_t;

        % Compute the pure spectra estimates (constrained, lower-dimensional space)
        constrained_S_t_reduced_unvec = (Y_regression_reduced*regression_matrix).';

        % For the next iteration, t+1
        S_t_reducedminus1  = constrained_S_t_reduced_unvec(:);
        sigma_t = (I_L - Kt*Ht)*cov_predict;
        sigmatminus1 = sigma_t;

        % Uncertainties
        uncertainties=[];
        
        if strcmp(confidence_level , "on") && mod(t,period)==0
            % Compute uncertainties
            Z = Y_regression*pinv(Y_regression_reduced);
            C = zeros(nb_pure_spectra * nb_spectral_channels, nb_pure_spectra *dim);
            for i = 1:nb_pure_spectra
                 C((i-1)*nb_spectral_channels+1:i*nb_spectral_channels, (i-1)*dim+1:i*dim) = Z;
            end
            diag_cov = diag(C*sigma_t*C.' + noise_var*eye(nb_pure_spectra*nb_spectral_channels));
            uncertainties = reshape(diag_cov , [nb_spectral_channels, nb_pure_spectra]);
        end

    end
    
end % end of KFOSU



%% Required functions

function selected_frequencies = select_frequencies(Y_regression, reconstruction_energy_percent, nb_spectral_channels)
    
    % Select frequency indices to retain for the DFT method
    % Inputs:
    %    Y_regression: archetype spectra (data space) 
    %    reconstruction_energy_percent: desired percentage of the total energy
    %    nb_spectral_channels: number of spectral channels

    % Output:
    %    selected_frequencies: indices of the selected frequencies
    
    
    % Apply FFT to each column (spectrum) of Y_regression
    Y_fft = fft(Y_regression);

    % Compute the energy of the FFT coefficients for each spectrum
    energy_fft = abs(Y_fft).^2/nb_spectral_channels;

    % Total energy
    total_energy = sum(energy_fft(:));
    
    % Compute cumulative energy and select frequency indices
    stop_ind = floor(nb_spectral_channels/2) -1;
    cumulative_energy = sum(energy_fft(1,:)) +2*cumsum(sum(energy_fft(2:stop_ind,:),2));
    max_freq_idx = find(cumulative_energy >=  reconstruction_energy_percent*total_energy, 1);
    selected_frequencies = 1:max_freq_idx;
    
end


function noise_variance = estimate_noise_variance(Y_regression, sgolay_order, sgolay_frame, Q)

    % Estimate the observation noise variance using Savitzky-Golay filter
    % Inputs:
    %   Y_regression: archetype spectra (data space)
    %   sgolay_order: polynomial order 
    %   sgolay_frame: frame length 
    %   Q: window size

    % Ouput:
    %   noise_variance: estimated variance
    

    [nb_spectral_channels, nb_spectra] = size(Y_regression);
    
    % Apply the Savitzky-Golay filter to each column of Y_regression
    % Compute the estimated noise matrix
    Noise = Y_regression - sgolayfilt(Y_regression, sgolay_order, sgolay_frame);
    
    % Reshape the noise matrix to stack the Q-sized windows for all spectra
    T = floor(nb_spectral_channels / Q);
    reshaped_noise = reshape(Noise(1:T*Q, :), Q, T * nb_spectra);
    
    % Compute the variance of each portion
    variances = var(reshaped_noise, 0, 1); % Compute variance along columns
    
    % Compute the median 
    noise_variance = median(variances);
    
end


function [regression_matrix, S] = Regression(regression_rho, regression_nb_iter, mu_init, S_init,...
    Y_regression_reduced, S_t_matrix_reduced, Y_regression)
    
    % Solve the constrained regression problem
    % Inputs:
    %   regression_rho: step-size
    %   regression_nb_iter: number of iterations
    %   mu_init, S_init: Initialization matrices
    %   Y_regression: archetype spectra (data space)
    %   Y_regression_reduced: reduced version of Y_regression
    %   S_t_matrix_reduced: unconstrained pure spectra matrix (reduced space)
  

    % Ouputs:
    %   regression_matrix: regression matrix
    %   S: constrained pure spectra matrix (data space)
    

    term1 = (2*Y_regression_reduced.')*Y_regression_reduced + (regression_rho*Y_regression.')*Y_regression;

    % Prevent numerical instability
    term1 = term1 + 1e-4*eye(size(term1,1));

    S = S_init;
    mu = mu_init;
    Z = zeros(size(S));

    for i=1:regression_nb_iter

        % Update regression_matrix
        term2 = Y_regression.'*mu + regression_rho*Y_regression.'*S + 2*Y_regression_reduced.'*S_t_matrix_reduced;
        regression_matrix = term1\term2;

        % Update S
        S = max(Z,Y_regression*regression_matrix - (1/regression_rho)*mu);

        % Update dual variable
        mu = mu + regression_rho*(S - Y_regression*regression_matrix);

    end

end


function [Y_regression, Y_regression_reduced] = archetype_spectra_update(Y_regression, Y_regression_reduced, y_t,y_t_reduced,...
    nb_spectral_channels, noise_var, adjust_val, probability_false_alarm)
        
        % Update archetype spectra (if necessary)
        % Inputs:
        %   Y_regression: archetype spectra (data space)
        %   Y_regression_reduced: reduced version of Y_regression
        %   y_t: new spectrum (data space)
        %   y_t_reduced: reduced version of y_t
        %   nb_spectral_channels: number of spectral channels
        %   noise_var: observation noise variance
        %   adjust_val: adjustment value; > 1 when noise_var is underestimated, between 0 and 1 otherwise
        %   probability_false_alarm: probability of false alarm for the hypothesis test

        % Ouputs:
        %   Y_regression: (updated) archetype spectra 
        %   Y_regression_reduced: reduced version of Y_regression


        % Compute the projection matrix 
        P = Y_regression*pinv(Y_regression);

        % Project y_t onto orthogonal complement of span(Y_regression)
        y_proj = y_t - P* y_t;

        % Compute the critical value
        critical_val = chi2inv(1 - probability_false_alarm, nb_spectral_channels - rank(Y_regression));
        
        if norm(y_proj)^2 /(adjust_val*noise_var) > critical_val
            % Update the archetype spectra
            Y_regression = [Y_regression y_t];
            Y_regression_reduced = [Y_regression_reduced y_t_reduced];
        end

end

