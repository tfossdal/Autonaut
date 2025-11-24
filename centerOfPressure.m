function x_cp = centerOfPressure(alpha)
    % Inputs:
    %   alpha: Angle of attack (rad)

    % Piecewise approximation for center of pressure
    if abs(alpha) <= 0.05 * pi
        x_cp = 0.25; % Quarter-chord
    elseif abs(alpha) <= 0.20 * pi
        x_cp = 0.2 + 0.1 / pi * abs(alpha);
    elseif abs(alpha) <= 0.80 * pi
        x_cp = 0.333 + 0.333 / pi * abs(alpha);
    elseif abs(alpha) <= 0.95 * pi
        x_cp = -0.2 + 0.1 / pi * abs(alpha);
    else
        x_cp = 0.75;
    end
end
