function F = get_conoid_force(length,model)
% lopt = model.conoid_length;
lopt = 0.017;
% k = model.conoid_stiffness;
k = 1e3;
% eps = model.conoid_eps;
eps = 1e-3;
d = length-lopt;
F = k ./ 2 .* (d + sqrt(d.^2 + eps));
end