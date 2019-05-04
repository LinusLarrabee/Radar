function sample_carrierwave = RF_Modulation()
fs = 5e11;
Ts = 1/fs;
f0 = 9e10;
Tend = 1e-7;
t = 0:Ts:Tend-Ts;
sample_carrierwave = cos(2*pi*f0*t);
plot(real(sample_carrierwave));
