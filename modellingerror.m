%%
error = V-ModeledVoltage_Linear;
id_discharge=find(I<0);
id_charge=find(I>0);
ModErrorDischarge=error(id_discharge);
ModErrorCharge=error(id_charge);
s_d=SOC(id_discharge);
%plot(s_d,ModErrorDischarge);
s_c=SOC(id_charge);
%plot(s_c,ModErrorCharge);

ss=0:0.01:1;
error_C=zeros(1,length(ss));
error_d=zeros(1,length(ss));
for i=2:length(ss)
    idd=find(s_d<ss(i) & s_d>ss(i-1));
    idc=find(s_c<ss(i) & s_c>ss(i-1));
    error_d(i)=mean(ModErrorDischarge(idd));
    error_c(i)=mean(ModErrorCharge(idc));
end
error_both=0.5*(error_d.^2+error_c.^2);
plot(ss,error_both);
xlabel('SOC');
ylabel('Error');
title("Modelling Error - Linear Model")