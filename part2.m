% Constants
c = 340; % Sound velocity (m/s)
fs = 16000; % Sample frequency (samples/s)
L = [6 5 3]; % Room dimensions [x y z] (m)
n = 1600; % Length of room impulse response (RAIR)
n_obs = 4000; % Length of the observation data samples
beta = 0.2; % Reverberation time (s)
microphone_pos = [2.5, 2, 1.5]; % Microphone position [x y z] (m)
P_values = [9, 36, 121, 289, 441, 676]; % Different P values
grid_resolutions = [0.2, 0.1, 0.05, 0.03, 0.025, 0.02]; % Corresponding grid resolutions
SNR_target = 15;

% AR model for excitation signal generation
alpha = 0.9; % AR(1) coefficient
ar_coeffs = [1, -alpha];
excitation_signal = filter(1, ar_coeffs, randn(n_obs, 1)); % Colored noise with length n_obs

% Preallocate arrays for misalignment values for both methods
misalignments_RAIR_DR = zeros(length(P_values), 1);
misalignments_conventional = zeros(length(P_values), 1);

% Loop over different P values
for i = 1:length(P_values)
    P = P_values(i);
    R = P - 1; % Reduced dimensionality
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

    % Generate the test RAIR h'
    test_source_pos = [3.74, 3.24, 1.4];
    h_prime = rir_generator(c, fs, microphone_pos, test_source_pos, L, beta, n);

    % Convolve the colored excitation signal with the test RAIR
    d = conv(excitation_signal, h_prime, 'same');
    % Calculate the power of the signal
    signal_power = bandpower(d);
    
    % Calculate the required noise power for desired SNR
    noise_power = signal_power / (10^(SNR_target / 10));
    
    % Generate white Gaussian noise
    noise = sqrt(noise_power / 2) * randn(size(d));
    
    % Add noise to the signal
    d_with_noise = d + noise;
    
    % Perform SVD and dimensionality reduction
    [U, S, V] = svd(H, 'econ');
    U_reduced = U(:, 1:R);
    
    % Create the input signal vector x(k) for system identification
    % Assuming x(k) is the excitation signal and we have n_obs samples
    X_toep = toeplitz(excitation_signal(n_obs:-1:1), excitation_signal(n_obs:-1:n_obs-n_rir+1))';

    % Conventional method - full dimensionality
    % Compute cross-correlation and auto-correlation matrices
    Rdx_conventional = X_toep * d_with_noise;
    Rxx_conventional = X_toep * X_toep';
    % Wiener filter estimate
    gw_conventional = Rxx_conventional \ Rdx_conventional;
    h_estimated_conventional = X_toep' * gw_conventional;

    % Calculate error and misalignment for conventional method
    e_conventional = d_with_noise - h_estimated_conventional;
    misalignment_numerator_conventional = norm(e_conventional, 2)^2;
    misalignment_denominator_conventional = norm(d_with_noise, 2)^2;
    misalignments_conventional(i) = 10 * log10(misalignment_numerator_conventional / misalignment_denominator_conventional);

    % RAIR-DR method - reduced dimensionality
    % Use U_reduced from previous RAIR-DR calculations
    % Compute cross-correlation and auto-correlation matrices
    d_with_noise_truncated = d_with_noise(1:1600);
    Rdx_reduced = U_reduced' * d_with_noise_truncated;
    Rxx_reduced = U_reduced' * U_reduced;
    % Wiener filter estimate
    gw_reduced = Rxx_reduced \ Rdx_reduced;
    h_estimated_reduced = U_reduced * gw_reduced;

    % Calculate error and misalignment for RAIR-DR method
    e_reduced = d_with_noise_truncated - h_estimated_reduced;
    misalignment_numerator_reduced = norm(e_reduced, 2)^2;
    misalignment_denominator_reduced = norm(d_with_noise_truncated, 2)^2;
    misalignments_RAIR_DR(i) = 10 * log10(misalignment_numerator_reduced / misalignment_denominator_reduced);

    % Display current P value and corresponding misalignments
    disp(['P = ' num2str(P)]);
    disp(['Conventional Misalignment (dB): ' num2str(misalignments_conventional(i))]);
    disp(['RAIR-DR Misalignment (dB): ' num2str(misalignments_RAIR_DR(i))]);
end

% Plot normalized misalignment vs P for both methods
figure;
hold on;
plot(P_values, misalignments_RAIR_DR, '^-', 'MarkerFaceColor', 'blue', 'DisplayName', 'RAIR-DR');
plot(P_values, misalignments_conventional, 'x-', 'MarkerFaceColor', 'red', 'DisplayName', 'Conventional');
xlabel('P (Number of source positions)');
ylabel('Normalized Misalignment (dB)');
title('Normalized Misalignment vs P');
legend('show');
grid on;
hold off;
