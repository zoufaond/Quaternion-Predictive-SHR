function dim = optimize_thorax_dim(trajectory_simulation,OS_model)
fun = @(x) objective(x,trajectory_simulation,OS_model);
x0 = [OS_model.thorax_radii];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.1;0.1;0.05];
ub = [0.2;0.25;0.15];

x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
dim = x;
end

function res = objective(x,trajectory_simulation,OS_model)
    numdata = size(trajectory_simulation,1);
    Epos = OS_model.thorax_center;
    for i= 1:numdata
        [iTS_pos,iAI_pos] = contact_points_position(trajectory_simulation(i,:),OS_model);
        eq_TS(i) = elips_eq(iTS_pos,Epos,x);
        eq_AI(i) = elips_eq(iAI_pos,Epos,x);
    end

    res = sum((eq_TS.^2 + eq_AI.^2).^(1/2));
end