%% EXERCISE 6: Multiple Linear Regression Model Estimation
% Exercises 3-4 required costly experiments. To avoid this, Van Cauter et al.
% proposed a new approach leveraging the "limited" variability of kinetics
% across different subjects, given certain anthropometric parameters: 
% gender, health status, BMI, age, etc.
% Study on over 200 subjects, covariates collected and followed by C-peptide
% analysis (previous i-o). Covariate model (regressors) vs kinetic 
% parameters (dependent variables).

clear all
close all
clc

run dati_vc.m

nomi_regressori={'Health Status' 'Gender' 'Age' 'Height' 'Weight' 'BMI' 'BSA'}; 
nomi_parametri={'Volume' 'Short Half-life' 'Long Half-life' 'Fraction'};

regressori=dati(:,1:7); % Covariates
parametri=dati(:,8:11);

%% UNIVARIATE ANALYSIS
% Iteratively compare each parameter with every covariate.
% One figure for each parameter and a subplot for each regressor.
% X-axis: regressor; Y-axis: i-th parameters. 
% Goal: Estimate a subject's parameters based on their characteristics.

[m, n]=size(regressori);
[o, p]=size(parametri);

for k=1:p
    figure(k)
    sgtitle(nomi_parametri{k})
    for i=1:n
        subplot(3,3,i)
        if i<3
            % Use boxplots for categorical data
            boxplot(dati(:,k+7),dati(:,i))
            title(nomi_regressori{i})
        else
            % Use scatter plots for interval data
            scatter(dati(:,i),dati(:,k+7))
            title(nomi_regressori{i})
        end
    end
end

% From the plots, the dispersion between parameters and their respective 
% regressor is nearly constant along the x-axis with no specific pattern 
% (higher data values do not exhibit higher variance). 
% Assuming a constant SD error model, we apply the Least Squares (LS) method.

%% STEPWISE FORWARD SELECTION:
% Start with an intercept-only model (b0) and iteratively add regressors.
% At each iteration, choose the regressor with the highest adjusted R-squared,
% then compare the adjusted R-squared of the new model against the old one
% to verify if adding the regressor increases the explained variance.

for i=2:7
    % Cell array containing values of individual regressors
    regressoricel{i}=regressori(:,i);
end

% Encode health status using a 2-level dummy variable:
% DIABETIC = [1 0], OBESE = [0 1], NORMAL = [0 0]
stato_salute=zeros(length(dati(:,1)),2); 

indici_normali=find(regressori(:,1)==0);
stato_salute(indici_normali,1)=0; stato_salute(indici_normali,2)=0;

indici_obesi=find(regressori(:,1)==1);
stato_salute(indici_obesi,1)=0; stato_salute(indici_obesi,2)=1;

indici_diabetici=find(regressori(:,1)==2);
stato_salute(indici_diabetici,1)=1; stato_salute(indici_diabetici,2)=0;

regressoricel{1}=stato_salute;

% Add interaction terms with gender to regressoricel:
for i=3:length(regressoricel) % Multiply age onwards by gender
    interazioni_sesso=cell2mat(regressoricel(:,i)).*cell2mat(regressoricel(:,2));
    regressoricel=[regressoricel interazioni_sesso];
end

% Create cubic interaction variables with gender.
% This allows evaluating whether the effect of age, height, weight, BMI, BSA
% depends on the subject's gender.
for i=3:7
    continue_cubo=cell2mat(regressoricel(:,i)).^3; % Cube values and add regressors
    regressoricel=[regressoricel continue_cubo];
end

% Add non-linear terms to model complex relationships between variables.
% If the link between regressor and covariate is non-linear, higher powers
% improve model fitting.
nomi_regressori={'Health Status' 'Gender' 'Age' 'Height' 'Weight' 'BMI' 'BSA' 'Gender*Age' 'Gender*Height' 'Gender*Weight' 'Gender*BMI' 'Gender*BSA' 'Age^3' 'Height^3' 'Weight^3' 'BMI^3' 'BSA^3'};
% Define a cell array of strings containing original and newly created 
% regressor names to identify them once selected.

% Perform the 4 regression problems (one for each parameter)
% using the stepwise.m function.

% VOLUME Regression Model (assuming constant standard deviation error model)
[r2_Volume,bhat_selezionato_Volume,Hselezionato_Volume,indici_selezionati_Volume] = stepwise(regressoricel,dati,parametri(:,1));
yhat_Volume = [ones(207,1) Hselezionato_Volume]*bhat_selezionato_Volume; % y = H * p
res_Volume = parametri(:,1) - yhat_Volume; % res = y - yhat
RSS_Volume = sum(res_Volume.^2);
Hselezionato_Volume=[ones(207,1),Hselezionato_Volume]; % Reconstruct complete regressor matrix with intercept
Alfa_V=RSS_Volume/(207-size(bhat_selezionato_Volume,1)-1); % s2 -> WRSS
SIGMAV_Volume=Alfa_V.*eye(207);
s2_V=sqrt(Alfa_V); % s

% Short Half-life (T12corto) Regression Model
[r2_T12corto,bhat_selezionato_T12corto,Hselezionato_T12corto,indici_selezionati_T12corto] = stepwise(regressoricel,dati,parametri(:,2));
yhat_T12corto = [ones(207,1) Hselezionato_T12corto]*bhat_selezionato_T12corto;
res_T12corto = parametri(:,2) - yhat_T12corto;
RSS_T12corto = sum(res_T12corto.^2);
Hselezionato_T12corto=[ones(207,1),Hselezionato_T12corto];
Alfa_T12corto=RSS_T12corto/(207-size(bhat_selezionato_T12corto,1)-1); 
SIGMAV_T12corto=Alfa_T12corto.*eye(207);
s2_T12corto=sqrt(Alfa_T12corto);

% Long Half-life (T12lungo) Regression Model
[r2_T12lungo,bhat_selezionato_T12lungo,Hselezionato_T12lungo,indici_selezionati_T12lungo] = stepwise(regressoricel,dati,parametri(:,3));
yhat_T12lungo = [ones(207,1) Hselezionato_T12lungo]*bhat_selezionato_T12lungo;
res_T12lungo = parametri(:,3) - yhat_T12lungo;
RSS_T12lungo = sum(res_T12lungo.^2);
Hselezionato_T12lungo=[ones(207,1),Hselezionato_T12lungo];
Alfa_T12lungo=RSS_T12lungo/(207-size(bhat_selezionato_T12lungo,1)-1); 
SIGMAV_T12lungo=Alfa_T12lungo.*eye(207);
s2_T12lungo=sqrt(Alfa_T12lungo);

% Fraction Regression Model
[r2_Fraction,bhat_selezionato_Fraction,Hselezionato_Fraction,indici_selezionati_Fraction] = stepwise(regressoricel,dati,parametri(:,4));
yhat_Fraction = [ones(207,1) Hselezionato_Fraction]*bhat_selezionato_Fraction;
res_Fraction = parametri(:,4) - yhat_Fraction;
RSS_Fraction = sum(res_Fraction.^2); % SSE
Hselezionato_Fraction=[ones(207,1),Hselezionato_Fraction];
Alfa_Fraction=RSS_Fraction/(207-size(bhat_selezionato_Fraction,1)-1); 
SIGMAV_Fraction=Alfa_Fraction.*eye(207);
s2_Fraction=sqrt(Alfa_Fraction);

% Display results of the 4 regression problems (selected optimal regressors):
disp(['Regressors for ' nomi_parametri{1} ' are: ' strjoin(nomi_regressori(indici_selezionati_Volume), ', ')])
disp(['Regressors for ' nomi_parametri{2} ' are: ' strjoin(nomi_regressori(indici_selezionati_T12corto), ', ')])
disp(['Regressors for ' nomi_parametri{3} ' are: ' strjoin(nomi_regressori(indici_selezionati_T12lungo), ', ')])
disp(['Regressors for ' nomi_parametri{4} ' are: ' strjoin(nomi_regressori(indici_selezionati_Fraction), ', ')])

% SIGMA_P = inv(H' * inv(SIGMA_V) * H)
SIGMAP_Volume=diag(inv(Hselezionato_Volume'*inv(SIGMAV_Volume)*Hselezionato_Volume));
SIGMAP_T12corto=diag(inv(Hselezionato_T12corto'*inv(SIGMAV_T12corto)*Hselezionato_T12corto));
SIGMAP_T12lungo=diag(inv(Hselezionato_T12lungo'*inv(SIGMAV_T12lungo)*Hselezionato_T12lungo));
SIGMAP_Fraction=diag(inv(Hselezionato_Fraction'*inv(SIGMAV_Fraction)*Hselezionato_Fraction));

sigmap_Volume=sqrt(SIGMAP_Volume);
sigmap_T12corto=sqrt(SIGMAP_T12corto);
sigmap_T12lungo=sqrt(SIGMAP_T12lungo);
sigmap_Fraction=sqrt(SIGMAP_Fraction);

% CALCULATE COEFFICIENTS OF VARIATION (CV) FOR THE 4 PARAMETER MODELS (%CV = sigma / mu):
CV_Volume=abs(sigmap_Volume./bhat_selezionato_Volume);
CV_T12corto=abs(sigmap_T12corto./bhat_selezionato_T12corto);
CV_T12lungo=abs(sigmap_T12lungo./bhat_selezionato_T12lungo);
CV_Fraction=abs(sigmap_Fraction./bhat_selezionato_Fraction);

% TWO-SIDED 95% CONFIDENCE INTERVALS FOR THE ESTIMATED BETA COEFFICIENTS:
% Since data are normal, beta coefficients are linear combinations of normal 
% variables with unknown variance, hence distributed as Student's t.
% Formula: IC = bhat +/- t(n - k - 1, alfa / 2) * std(bhat)
alfa=0.05;

% Confidence Interval for Volume:
IC_Volume=[bhat_selezionato_Volume-tinv(alfa/2,207-size(bhat_selezionato_Volume,1)-1)*sigmap_Volume, bhat_selezionato_Volume+tinv(alfa/2,207-size(bhat_selezionato_Volume,1)-1)*sigmap_Volume];
% Each row represents a CI for each bhat found for the Volume parameter

% Confidence Interval for Short Half-life:
IC_tcorto=[bhat_selezionato_T12corto-tinv(alfa/2,207-size(bhat_selezionato_T12corto,1)-1)*sigmap_T12corto, bhat_selezionato_T12corto+tinv(alfa/2,207-size(bhat_selezionato_T12corto,1)-1)*sigmap_T12corto];

% Confidence Interval for Long Half-life:
IC_tlungo=[bhat_selezionato_T12lungo-tinv(alfa/2,207-size(bhat_selezionato_T12lungo,1)-1)*sigmap_T12lungo, bhat_selezionato_T12lungo+tinv(alfa/2,207-size(bhat_selezionato_T12lungo,1)-1)*sigmap_T12lungo];

% Confidence Interval for Fraction:
IC_Fraction=[bhat_selezionato_Fraction-tinv(alfa/2,207-size(bhat_selezionato_Fraction,1)-1)*sigmap_Fraction, bhat_selezionato_Fraction+tinv(alfa/2,207-size(bhat_selezionato_Fraction,1)-1)*sigmap_Fraction];

%% PREDICTIONS VS RESIDUALS PLOT
% Good fitting plot: random pattern expected
figure('Name','Residual Analysis')
subplot(2,2,1)
scatter(yhat_Volume,res_Volume)
title('Volume')
xlabel('Estimated Values')
ylabel('Residuals')

subplot(2,2,2)
scatter(yhat_T12corto,res_T12corto)
title('Short Half-life')
xlabel('Estimated Values')
ylabel('Residuals')

subplot(2,2,3)
scatter(yhat_T12lungo,res_T12lungo)
title('Long Half-life')
xlabel('Estimated Values')
ylabel('Residuals')

subplot(2,2,4)
scatter(yhat_Fraction,res_Fraction)
title('Fraction')
xlabel('Estimated Values')
ylabel('Residuals')

%% RESIDUALS VS REGRESSORS PLOT FOR EACH MODEL
% Check for any systematic trends relative to specific regressors.
% We expect error assumptions to hold, hence no systematic patterns should appear.
figure('Name','Volume Analysis')
subplot(1,2,1)
scatter(regressoricel{:,17},res_Volume)
title('BSA^3')
xlabel('BSA^3 Data')
ylabel('Volume Residuals')
subplot(1,2,2)
scatter(regressoricel{:,13},res_Volume)
title('Age^3')
xlabel('Age^3 Data')
ylabel('Volume Residuals')

figure('Name','Short Half-life Analysis')
subplot(1,2,1)
scatter(regressoricel{:,6},res_T12corto)
title('BMI')
xlabel('BMI Data')
ylabel('Short Half-life Residuals')
subplot(1,2,2)
scatter(regressoricel{:,16},res_T12corto)
title('BMI^3')
xlabel('BMI^3 Data')
ylabel('Short Half-life Residuals')

figure('Name','Long Half-life Analysis')
subplot(2,4,1)
scatter(regressoricel{:,13},res_T12lungo)
title('Age^3')
xlabel('Age^3 Data')
ylabel('Long Half-life Residuals')
subplot(2,4,2)
scatter(regressoricel{:,14},res_T12lungo)
title('Height^3')
xlabel('Height^3 Data')
ylabel('Long Half-life Residuals')
subplot(2,4,3)
scatter(regressoricel{:,3},res_T12lungo)
title('Age')
xlabel('Age Data')
ylabel('Long Half-life Residuals')
subplot(2,4,4)
scatter(regressoricel{:,15},res_T12lungo)
title('Weight^3')
xlabel('Weight^3 Data')
ylabel('Long Half-life Residuals')
subplot(2,4,5)
scatter(regressoricel{:,5},res_T12lungo)
title('Weight')
xlabel('Weight Data')
ylabel('Long Half-life Residuals')
subplot(2,4,6)
scatter(regressoricel{:,1},res_T12lungo)
title('Health Status')
xlabel('Health Status Data')
ylabel('Long Half-life Residuals')
subplot(2,4,7)
scatter(regressoricel{:,16},res_T12lungo)
title('BMI^3')
xlabel('BMI^3 Data')
ylabel('Long Half-life Residuals')
subplot(2,4,8)
scatter(regressoricel{:,6},res_T12lungo)
title('BMI')
xlabel('BMI Data')
ylabel('Long Half-life Residuals')

figure('Name','Fraction Analysis')
scatter(regressoricel{:,6},res_Fraction)
title('BMI')
xlabel('BMI Data')
ylabel('Fraction Residuals')

%% POINT PREDICTION USING PERSONAL PARAMETERS:
disp('Personal Prediction');
stato=input('Health status: 1=Diabetic / 2=Obese / Other=Normal: ');
if stato==1
    stato1=1; stato2=0;
elseif stato==2
    stato1=0; stato2=1;
else
    stato1=0; stato2=0;
end

eta=input('Age: ');
altezza=input('Height: ');
sesso=input('Gender M=0/F=1: ');
peso=input('Weight: ');

BMI=peso/altezza^2;
BSA=0.20247*(altezza^0.725)*(peso^0.425);

% PARAMETERIZATION WITH MULTIPLE REGRESSION MODEL:
Volume=[1 BSA^3 eta^3]*bhat_selezionato_Volume; % y = H * p
Emivita_corto=[1 BMI BMI^3]*bhat_selezionato_T12corto;
Emivita_lungo=[1 eta^3 altezza^3 eta peso^3 peso stato1 stato2 BMI^3 BMI]*bhat_selezionato_T12lungo;
Fraction=[1 BMI]*bhat_selezionato_Fraction;

disp(['V=',num2str(Volume),'; Short Half-life=',num2str(Emivita_corto),'; Long Half-life=',num2str(Emivita_lungo),'; Fraction=',num2str(Fraction)]);

my_H_emivita_lungo=[1 eta^3 altezza^3 eta peso^3 peso stato1 stato2 BMI^3 BMI];
s_pred_t12lungo=sqrt(Alfa_T12lungo+my_H_emivita_lungo*inv(Hselezionato_T12lungo'*inv(SIGMAV_T12lungo)*Hselezionato_T12lungo)*my_H_emivita_lungo');

my_H_emivita_corto=[1 BMI BMI^3];
s_pred_t12corto=sqrt(Alfa_T12corto+my_H_emivita_corto*inv(Hselezionato_T12corto'*inv(SIGMAV_T12corto)*Hselezionato_T12corto)*my_H_emivita_corto');

my_H_Fraction=[1 BMI];
s_pred_fraction=sqrt(Alfa_Fraction+my_H_Fraction*inv(Hselezionato_Fraction'*inv(SIGMAV_Fraction)*Hselezionato_Fraction)*my_H_Fraction');

my_H_Volume=[1 BSA eta];
s_pred_volume=sqrt(Alfa_V+my_H_Volume*inv(Hselezionato_Volume'*inv(SIGMAV_Volume)*Hselezionato_Volume)*my_H_Volume');

% CONFIDENCE INTERVALS FOR POINT REGRESSION MODELS:
% Formula: IC = personal_betahat +/- tinv(alfa / 2, n - p - 1) * personal_s_pred
% YHAT follows Student's t distribution because it is a combination of normal variables with unknown variance.
IC_my_Volume=[Volume-tinv(0.05/2,207-size(Hselezionato_Volume,2)-1)*s_pred_volume; Volume+tinv(0.05/2,207-size(Hselezionato_Volume,2)-1)*s_pred_volume];
IC_my_T12corto=[Emivita_corto-tinv(0.05/2,207-size(Hselezionato_T12corto,2)-1)*s_pred_t12corto; Emivita_corto+tinv(0.05/2,207-size(Hselezionato_T12corto,2)-1)*s_pred_t12corto];
IC_my_T12lungo=[Emivita_lungo-tinv(0.05/2,207-size(Hselezionato_T12lungo,2)-1)*s_pred_t12lungo; Emivita_lungo+tinv(0.05/2,207-size(Hselezionato_T12lungo,2)-1)*s_pred_t12lungo];
IC_my_Fraction=[Fraction-tinv(0.05/2,207-size(Hselezionato_Fraction,2)-1)*s_pred_fraction; Fraction+tinv(0.05/2,207-size(Hselezionato_Fraction,2)-1)*s_pred_fraction];

%% PARAMETERIZATION WITH BIEXPONENTIAL MODEL (INPUT-OUTPUT MODEL):
D=1; % Unit impulse response, or D = AUC * Beta * Volume
c0=D/Volume;

% Distribution half-life = log(2) / alpha, assumed to be t-distributed
Alfa=log(2)/Emivita_corto;

% Elimination half-life = log(2) / beta, assumed to be t-distributed
Beta=log(2)/Emivita_lungo;

A=Fraction/Volume; % A is t-distributed: A = Fraction * c0
B=(1-Fraction)/Volume; % B is t-distributed: B = c0 - A

t=0:0.1:100;
figure
plot(t,A*exp(-Alfa*t)+B*exp(-Beta*t));
xlabel('Time [sec]')
ylabel('Concentration [kg/L]')
disp(['My Biexponential Model: A=',num2str(A),'; B=',num2str(B),'; Alpha=',num2str(Alfa),'; Beta=',num2str(Beta)]);

% Assuming independent parameters when calculating confidence intervals.
% CI for A at 95%: A depends on Fraction and Volume non-linearly. 
% Generate Fraction and Volume via simulation and calculate A for each pair to find the empirical CI.
for i=1:10000
    Fra(i)=Fraction+trnd(207-size(bhat_selezionato_Fraction,1),1)*s_pred_fraction;
    Vol(i)=Volume+trnd(207-size(bhat_selezionato_Volume,1),1)*s_pred_volume;
    A_camp(i)=Fra(i)./Vol(i);
end
IC_A=prctile(A_camp,[0.05/2,1-(0.05/2)]);

% CI for B at 95%: B depends on Fraction and Volume non-linearly.
for i=1:10000
    Fra(i)=Fraction+trnd(207-size(bhat_selezionato_Fraction,1),1)*s_pred_fraction;
    Vol(i)=Volume+trnd(207-size(bhat_selezionato_Volume,1),1)*s_pred_volume;
    B_camp(i)=(1-Fra(i))./Vol(i);
end
IC_B=prctile(B_camp,[0.05/2,1-(0.05/2)]);

% CI for Alpha at 95%: Alpha depends non-linearly on Short Half-life.
for i=1:10000
    T12corto(i)=Emivita_corto+trnd(207-size(bhat_selezionato_T12corto,1),1)*s_pred_t12corto;
    Alfa_camp(i)=log(2)./T12corto(i);
end
IC_Alfa=prctile(Alfa_camp,[0.05/2,1-(0.05/2)]);

% CI for Beta at 95%: Beta depends non-linearly on Long Half-life.
for i=1:10000
    T12lungo(i)=Emivita_lungo+trnd(207-size(bhat_selezionato_T12lungo,1),1)*s_pred_t12lungo;
    Beta_camp(i)=log(2)./T12lungo(i);
end
IC_Beta=prctile(Beta_camp,[0.05/2,1-(0.05/2)]);

%% TWO-COMPARTMENT MODEL PARAMETERIZATION:
k12=(Beta*A+Alfa*B)/(A+Beta);
k01=(Alfa*Beta)/k12;
v1=1/(A+B);
k21=Alfa+Beta-k01-k12;

disp(['My Two-Compartment Model: k12=',num2str(k12),'; k01=',num2str(k01),'; k21=',num2str(k21),'; V1=',num2str(v1)]);

% CI for k12: k12 depends non-linearly on A, B, Alpha, Beta. Simulate them to build an empirical distribution for k12.
for i=1:10000
    Fra(i)=Fraction+trnd(207-size(bhat_selezionato_Fraction,1),1)*s_pred_fraction;
    Vol(i)=Volume+trnd(207-size(bhat_selezionato_Volume,1),1)*s_pred_volume;
    T12corto(i)=Emivita_corto+trnd(207-size(bhat_selezionato_T12corto,1),1)*s_pred_t12corto;
    T12lungo(i)=Emivita_lungo+trnd(207-size(bhat_selezionato_T12lungo,1),1)*s_pred_t12lungo;
    A_camp(i)=Fra(i)./Vol(i);
    B_camp(i)=(1-Fra(i))./Vol(i);
    Alfa_camp(i)=log(2)./T12corto(i);
    Beta_camp(i)=log(2)./T12lungo(i);
    k12_camp(i)=(Beta_camp(i).*A_camp(i)+Alfa_camp(i).*B_camp(i))./(A_camp(i)+Beta_camp(i));
end
IC_k12=prctile(k12_camp,[0.05/2,1-(0.05/2)]);

% CI for k01: k01 depends non-linearly on A, B, Alpha, Beta, k12.
for i=1:10000
    Fra(i)=Fraction+trnd(207-size(bhat_selezionato_Fraction,1),1)*s_pred_fraction;
    Vol(i)=Volume+trnd(207-size(bhat_selezionato_Volume,1),1)*s_pred_volume;
    T12corto(i)=Emivita_corto+trnd(207-size(bhat_selezionato_T12corto,1),1)*s_pred_t12corto;
    T12lungo(i)=Emivita_lungo+trnd(207-size(bhat_selezionato_T12lungo,1),1)*s_pred_t12lungo;
    A_camp(i)=Fra(i)./Vol(i);
    B_camp(i)=(1-Fra(i))./Vol(i);
    Alfa_camp(i)=log(2)./T12corto(i);
    Beta_camp(i)=log(2)./T12lungo(i);
    k12_camp(i)=(Beta_camp(i).*A_camp(i)+Alfa_camp(i).*B_camp(i))./(A_camp(i)+Beta_camp(i));
    k01_camp(i)=(Alfa_camp(i).*Beta_camp(i))./k12_camp(i);
end
IC_k01=prctile(k01_camp,[0.05/2,1-(0.05/2)]);

% CI for k21: k21 depends non-linearly on A, B, Alpha, Beta, k12, k01.
for i=1:10000
    Fra(i)=Fraction+trnd(207-size(bhat_selezionato_Fraction,1),1)*s_pred_fraction;
    Vol(i)=Volume+trnd(207-size(bhat_selezionato_Volume,1),1)*s_pred_volume;
    T12corto(i)=Emivita_corto+trnd(207-size(bhat_selezionato_T12corto,1),1)*s_pred_t12corto;
    T12lungo(i)=Emivita_lungo+trnd(207-size(bhat_selezionato_T12lungo,1),1)*s_pred_t12lungo;
    A_camp(i)=Fra(i)./Vol(i);
    B_camp(i)=(1-Fra(i))./Vol(i);
    Alfa_camp(i)=log(2)./T12corto(i);
    Beta_camp(i)=log(2)./T12lungo(i);
    k12_camp(i)=(Beta_camp(i).*A_camp(i)+Alfa_camp(i).*B_camp(i))./(A_camp(i)+Beta_camp(i));
    k01_camp(i)=(Alfa_camp(i).*Beta_camp(i))./k12_camp(i);
    k21_camp(i)=k12_camp(i).*k01_camp(i)./Beta_camp(i)+Beta_camp(i)-k01_camp(i)-k12_camp(i);
end
IC_k21=prctile(k21_camp,[0.05/2,1-(0.05/2)]);

% CI for v1: v1 depends non-linearly on A, Alpha, Beta, k12.
for i=1:10000
    Fra(i)=Fraction+trnd(207-size(bhat_selezionato_Fraction,1),1)*s_pred_fraction;
    Vol(i)=Volume+trnd(207-size(bhat_selezionato_Volume,1),1)*s_pred_volume;
    T12corto(i)=Emivita_corto+trnd(207-size(bhat_selezionato_T12corto,1),1)*s_pred_t12corto;
    T12lungo(i)=Emivita_lungo+trnd(207-size(bhat_selezionato_T12lungo,1),1)*s_pred_t12lungo;
    A_camp(i)=Fra(i)./Vol(i);
    B_camp(i)=(1-Fra(i))./Vol(i);
    Alfa_camp(i)=log(2)./T12corto(i);
    Beta_camp(i)=log(2)./T12lungo(i);
    k12_camp(i)=(Beta_camp(i).*A_camp(i)+Alfa_camp(i).*B_camp(i))./(A_camp(i)+Beta_camp(i));
    v1_camp(i)=(1./A_camp(i)).*(k12_camp(i)-Alfa_camp(i))./(Beta_camp(i)-Alfa_camp(i));
end
IC_v1=prctile(v1_camp,[0.05/2,1-(0.05/2)]);