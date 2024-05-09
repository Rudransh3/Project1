function R = findRankForEnergyThreshold(singular_values, energy_threshold)
    if nargin < 2
        energy_threshold = 0.95; % default energy threshold if not provided
    end
    
    total_energy = sum(singular_values.^2);
    energy_captured = 0;
    R = 0;
    
    for i = 1:length(singular_values)
        energy_captured = energy_captured + singular_values(i)^2;
        if energy_captured / total_energy >= energy_threshold
            R = i;
            break;
        end
    end
end
