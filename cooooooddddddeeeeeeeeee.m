% Constants
c = 340; % Sound velocity (m/s)
fs = 16000; % Sample frequency (samples/s)
L = [6 5 3]; % Room dimensions [x y z] (m)
n = 1600; % Number of samples
beta = 0.2; % Reverberation time (s)
microphone_pos = [2.5, 2, 1.5]; % Microphone position [x y z] (m)
P_values = [9, 36, 121, 289, 441, 676]; % Different P values
grid_resolutions = [0.2, 0.1, 0.05, 0.03, 0.025, 0.02]; % Corresponding grid resolutions

% Preallocate arrays for misalignment values
misalignments = zeros(length(P_values), 1);

% Loop over different P values
for i = 1:length(P_values)
    P = P_values(i);
    grid_resolution = grid_resolutions(i);
    
    % Generate grid for source positions
    [X, Y] = meshgrid(3.5:grid_resolution:4, 3:grid_resolution:3.5);
    Z = 1.4 * ones(size(X));
    source_positions = [X(:), Y(:), Z(:)];
    
    % Generate RAIRs for each source position
    H = zeros(n, P);
    for j = 1:P
        source_pos = source_positions(j, :);
        H(:, j) = rir_generator(c, fs, microphone_pos, source_pos, L, beta, n);
    end
    
    % Perform SVD and dimensionality reduction
    [U, S, V] = svd(H, 'econ');
%     
%     % Determine a cutoff for significant singular values
%     singular_values = diag(S);
%     plot(singular_values);
%     title(sprintf('Singular Values for P=%d', P));
%     xlabel('Index');
%     ylabel('Singular Value');
%     drawnow; % Force MATLAB to draw the plot immediately
%     
%     % Pause to inspect the singular values plot
%     disp('Press a key to continue...');
%     pause;
    
    % Selecting a rank for dimensionality reduction
    % This should be based on the singular value plot
    %     R_exp = findElbow(singular_values);
    %     disp(R_exp);
    singular_values = diag(S);
    if P<=289
        R = findRankForEnergyThreshold(singular_values, 0.90);
        disp(R);
    else 
        R = findRankForEnergyThreshold(singular_values,0.95);
        disp(R);
    end
    
    % Use the chosen R to reduce dimensions
    U_reduced = U(:, 1:R);
    
    % Generate test RAIR (h') from a random position
    test_source_pos = [3.74, 3.24, 1.4];
    h_prime = rir_generator(c, fs, microphone_pos, test_source_pos, L, beta, n);
    h_prime = h_prime(:);
    % Compute low-dimensional representation of test RAIR (h')
    h_prime_low_dim = U_reduced' * h_prime;
    
    % Reconstruct h' from the low-dimensional representation
    h_prime_reconstructed = U_reduced * h_prime_low_dim;
    
    % Calculate misalignment
    misalignment_numerator = norm(h_prime - h_prime_reconstructed, 2)^2;
    misalignment_denominator = norm(h_prime, 2)^2;
    misalignments(i) = (misalignment_numerator / misalignment_denominator);
    
    % Convert misalignment to decibels (dB)
    misalignments(i) = 10 * log10(misalignments(i));
end

% Plot normalized misalignment vs P
figure;
plot(P_values, misalignments, '-o');
xlabel('P (Number of source positions)');
ylabel('Normalized Misalignment (dB)');
title('Normalized Misalignment vs P');
grid on;
