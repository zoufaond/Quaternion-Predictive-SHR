function res = scale_dof(dof,scale)

    max_val = max(dof);
    res = scale * (dof - max_val) + dof;

end