function mu = sutherland(T)
%SUTHERLAND Uses Sutherland's law to calculate viscosity given temperature (for
%air)

    % Properties for air
    mu0 = 1.735e-5; % Ns/m^2
    T0 = 288.15; % K
    S1 = 110.4; % K

    % Calculate viscosity
    mu = mu0 * (T/T0).^(3/2) .* (T0 + S1) ./ (T + S1);

end