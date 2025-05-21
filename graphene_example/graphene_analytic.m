function rho = graphene_analytic(E, t)
    Z0 = (1.0 + abs(E / t)).^2 - ((E / t).^2 - 1.0).^2 / 4.0;
    Z1 = 4 * abs(E / t);

    if (abs(E) <= t)
        Z0 = (1.0 + abs(E / t)).^2 - ((E / t).^2 - 1.0).^2 / 4.0;
        Z1 = 4 * abs(E / t);
    elseif (abs(E) <= 3 * t)
        Z1 = (1.0 + abs(E / t)).^2 - ((E / t).^2 - 1.0).^2 / 4.0;
        Z0 = 4 * abs(E / t);
    else
        rho = 0;
        return;
    end

    rho = 4/pi^2 * abs(E)/t^2 / sqrt(Z0) * ellipke(Z1./Z0);
    
    rho = rho/4; % Normalize so that integral of rho is 1
end