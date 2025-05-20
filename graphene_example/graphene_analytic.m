% Analytic formula for nearest-neighbor tight-binding graphene density of states
% function rho = graphene_analytic(E, t)
%     e = abs(E/t);
%     Fe = (e+1)^3*(3-e)/16;
%     if (e <= 1)
%         rho = e/(pi*t^2)/sqrt(Fe)*ellipke(e/Fe);
%     elseif e <= 3
%         rho = e/(pi*t^2)/sqrt(e)*ellipke(Fe/e);
%     else
%         rho = 0;
%     end
% end

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