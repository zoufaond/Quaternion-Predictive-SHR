clearvars
N = 10000;
vnorm = linspace(-1,1,N);
f_gauss = 0.25;
kpe = 5;
epsm0 = 0.6;
% f = switcher(lnorm,100).*(exp(kpe*(lnorm - 1)/epsm0)-1)/(exp(kpe)-1);
% f = switcher(lnorm,100);
for i=1:N
    % f_ch(i) = fvce_chad(vnorm(i));
    f_th(i) = fvce_T03(vnorm(i));
    f_gr(i) = fvce_groot(vnorm(i));
end
plot(vnorm,f_th, vnorm, f_gr)
legend('Thelen','groot')

% pp = polyfit(vnorm, f_ch,2);
% fitted = polyval(pp,vnorm);
% 
% 
% % figure
% plot(vnorm,f_ch); hold on
% plot(vnorm,fitted)
% legend('ch','fitted')

function res = switcher(val,thresh,scale)
    res = (tanh(scale*(val-thresh))+1)/2;

end

function fvce = fvce_groot(vce)
    d1 = -0.318;
    d2 = -8.149;
    d3 = -0.374;
    d4 = 0.886;
    fvce = d1 * log(d2 * vce + d3 + sqrt((d2 * vce + d3)^2 + 1)) + d4;
end

function fvce = fvce_chad(vce)
    vmax = 10;
    A = 0.25;
    gmax = 1.5;
    if vce <= 0
        fvce = (vmax+vce)/(vmax-vce/A);
    else
        c3 = vmax*A*(gmax-1)/(A+1);
        fvce = (gmax*vce+c3)/(vce+c3);
    end
    % fvce_cc = (vmax+vce)/(vmax-vce/A);
    % c3 = vmax*A*(gmax-1)/(A+1);
    % fvce_exc = (gmax*vce+c3)/(vce+c3);

    % fvce =  fvce_exc; %switcher(-vce,0,5)*fvce_cc


end

function fvce=fvce_T03(vm)
flce=1;
lmopt=1;
a=1;
vmmax = 1;
fmlen=2;
af=2;

if vm <= 0
    fvce = af*a*flce*(4*vm + vmmax*(3*a + 1))/(-4*vm + vmmax*af*(3*a + 1));
else
    fvce = a*flce*(af*vmmax*(3*a*fmlen - 3*a + fmlen - 1) + 8*vm*fmlen*(af + 1))/(af*vmmax*(3*a*fmlen - 3*a + fmlen - 1) + 8*vm*(af + 1));
end

end