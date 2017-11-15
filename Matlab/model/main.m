
f_min = ;
f_max = ;

I_modeled = @(f) temporal_psd(f);
I_measured = @(f) f;

aux = @(f) (log10(I_modeled(f)) - log10(I_measured(f)))^2;
