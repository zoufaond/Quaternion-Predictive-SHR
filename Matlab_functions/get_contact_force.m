dist = linspace(-0.1,0.1,100);
forces = get_contact_forcee(dist);
forces_lig = get_conoid_force(dist,0);
figure
plot(dist,forces)

figure
plot(dist,forces_lig)

function F = get_contact_forcee(distance)
    k = 20000;
    eps = 0.005;
    F = 1/2.*k.*(distance-sqrt(distance.^2+eps.^2));
end