data=xlsread('data.xls');
T=data(:,1);% Time in hours
I=data(:,2);% Current in Ampere
V=data(:,3);% Voltage in Volts
SOC=data(:,4);% SOC in [0 1]
s=(1-2*0.175).*SOC+0.175;%scaling to avoid computation erros

%% Linear Model
p_linear=[ones([length(SOC),1]),SOC,I];
k_linear=(p_linear'*p_linear)\(p_linear'*V);
OCV_Linear=p_linear(:,1:2)*k_linear(1:2,:);     
ModeledVoltage_Linear=p_linear*k_linear; 
ko_linear=[k_linear(1) k_linear(2)];
fprintf("\n The Parameters of Linear Model are [");
fprintf("%g, ", ko_linear(1:end-1));  
fprintf("%g]\n", ko_linear(end));
fprintf("\n The value of R_0,h is %f ohms\n",k_linear(3));

% Estimation of correlation coefficient R^2            

Sr=sum((V-p_linear*k_linear).^2); % sum of errors
meanY=(sum(V))*(1/length(V)); %mean
St=sum((V-meanY).^2); %sum of errors around mean
r=sqrt((St-Sr)./St);
r2=r.^2;
fprintf("\n The value of Correlation Coefficient r^2 of Linear Model is %f  \n",r2);

%OCV-SOC CURVE

figure; hold on; box on; grid on
plot(SOC, V, SOC, OCV_Linear, SOC, ModeledVoltage_Linear, '-', 'linewidth', 2); grid on;
legend('Measured Voltage', 'OCV', 'Modeled Voltage')
xlabel('SOC')
ylabel('OCV');
title('OCV-SOC Curve - Linear Model. ECM Parameter R0,h = ',num2str(k_linear(3)))
set(gca, 'fontsize', 14)  

%% Polynomial Model
p_polynomial=[ones([length(s),1]),s,s.^2,s.^3,s.^4,s.^5,1./s,1./s.^2,1./s.^3,1./s.^4,1./s.^5,I];
k_polynomial=(p_polynomial'*p_polynomial)\(p_polynomial'*V);
OCV_Polynomial=p_polynomial(:,1:11)*k_polynomial(1:11,:);
ModeledVoltage_Poly=p_polynomial*k_polynomial;  
ko_polynomial=[k_polynomial(1) k_polynomial(2) k_polynomial(3) k_polynomial(4)...
    k_polynomial(5) k_polynomial(6) k_polynomial(7) k_polynomial(8)...
    k_polynomial(9) k_polynomial(10) k_polynomial(11)];
fprintf("\n The Parameters of Polynomial Model are [");
fprintf('%g, ', ko_polynomial(1:end-1));  
fprintf("%g]\n", ko_polynomial(end));
fprintf("\n The value of R_0,h is %f ohms \n",k_polynomial(12));

% Estimation of correlation coefficient R^2            

Sr=sum((V-p_polynomial*k_polynomial).^2); % sum of errors
meanY=(sum(V))*(1/length(V)); %mean
St=sum((V-meanY).^2); %sum of errors around mean
r=sqrt((St-Sr)./St);
r2=r^2;
fprintf("\n The value of Correlation Coefficient r^2 of Polynomial Model is %f  \n",r2);

%OCV-SOC

figure; hold on; box on; grid on
plot(SOC, V, SOC, OCV_Polynomial, SOC, ModeledVoltage_Poly, '-', 'linewidth', 2); grid on
legend('Measured Voltage', 'OCV', 'Modeled Voltage')
xlabel('SOC')
ylabel('OCV');
title('OCV-SOC Curve - Polynomial Model. ECM Parameter R0,h = ',num2str(k_polynomial(12)))
set(gca, 'fontsize', 14)  

%% Combined Model
p_combined=[ones([length(s),1]),1./s,s,log(s),log(1-s),I];
k_combined=(p_combined'*p_combined)\(p_combined'*V);
OCV_Combined=p_combined(:,1:5)*k_combined(1:5,:);
ModeledVoltage_Combined=p_combined*k_combined;
ko_combined=[k_combined(1) k_combined(2) k_combined(3)...
    k_combined(4) k_combined(5)];
fprintf("\n The Parameters of Combined Model are [");
fprintf('%g, ', ko_combined(1:end-1));  
fprintf("%g]\n", ko_combined(end));
fprintf("\n The value of R_0,h is %f ohms \n",k_combined(6));

% Estimation of correlation coefficient R^2            

Sr=sum((V-p_combined*k_combined).^2); % sum of errors
meanY=(sum(V))*(1/length(V)); %mean
St=sum((V-meanY).^2); %sum of errors around mean
r=sqrt((St-Sr)./St);
r2=r^2;
fprintf("\n The value of correlation coefficient r^2 of Combined Model is %f \n",r2);

%OCV-SOC

figure; hold on; box on; grid on
plot(SOC, V, SOC, OCV_Combined, SOC, ModeledVoltage_Combined, '-', 'linewidth', 2); grid on
legend('Measured Voltage', 'OCV', 'Modeled Voltage')
xlabel('SOC')
ylabel('OCV');
title('OCV-SOC Curve - Combined Model. ECM Parameter R0,h = ',num2str(k_combined(6)))
set(gca, 'fontsize', 14)   

%% Combined+3 Model
p_combined3=[ones([length(s),1]),1./s,1./(s.^2),1./(s.^3),1./(s.^4),s,log(s),log(1-s),I];
k_combined3=(p_combined3'*p_combined3)\(p_combined3'*V);%least square estimate to find parameters
OCV_Combined3=p_combined3(:,1:8)*k_combined3(1:8,1);%ocv model
ModeledVoltage_Combined3= p_combined3*k_combined3;
ko_combined3=[k_combined3(1) k_combined3(2) k_combined3(3) k_combined3(4)...
    k_combined3(5) k_combined3(6) k_combined3(7) k_combined3(8)];
fprintf("\n The Parameters of Combined+3 Model are [");
fprintf('%g, ',ko_combined3(1:end-1));  
fprintf("%g]\n", ko_combined3(end));
fprintf("\n The value of R_0,h is %f ohms \n",k_combined3(9));

% Estimation of correlation coefficient R^2            

Sr=sum((V-p_combined3*k_combined3).^2); % sum of errors
meanY=(sum(V))*(1/length(V)); %mean
St=sum((V-meanY).^2); %sum of errors around mean
r=sqrt((St-Sr)./St);
r2=r^2; 
fprintf("\n The value of correlation coefficient r^2 of Combined Model+3 is %f \n",r2);
  
%OCV-SOC

figure; hold on; box on; grid on
plot(SOC, V, SOC, OCV_Combined3, SOC, ModeledVoltage_Combined3, '-', 'linewidth', 2); grid on
legend('Measured Voltage', 'OCV', 'Modeled Voltage')
ylabel('OCV');
xlabel('SOC')
title('OCV-SOC Curve - Combined+3 Model. ECM Parameter R0,h = ',num2str(k_combined3(9)))
set(gca, 'fontsize', 14)    

