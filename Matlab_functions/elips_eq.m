function res = elips_eq(IP,Epos,Edim)
    res = ((IP(1)-Epos(1))/Edim(1))^2+((IP(2)-Epos(2))/Edim(2))^2+((IP(3)-Epos(3))/Edim(3))^2-1;
end