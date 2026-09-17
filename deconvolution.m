%% EXERCISE 7: Deconvolution
% Up & down glycemic response test in a subject. Plasma C-peptide 
% concentration is measured over time.

clear all
close all
clc

run dati_es7.m

% Leverage the multiple regression models found for the 4 parameters 
% of the Van Cauter model
bhat_Volume=[3.306637282268909;0.150737862067233;-1.743683226744852e-06];
bhat_T12corto=[7.358661701650782;-0.114385477646278;2.377542309886689e-05];
bhat_T12lungo=[60.750707655580330;5.282465913821193e-05;-2.146151466108115;-0.113079330604868;-1.254884869209628e-07;0.090397542469936;3.370565944307208;4.070281266592681;1.837471820484636e-04;-1.041882067862100];
bhat_Fraction=[0.740055359072289;0.001126646841145];

% yb is the basal concentration in pmol/ml
t=DATA(:,1);
y=DATA(:,2)-yb; % Known output (plasma concentration in pmol/ml) minus basal value

% List the subject's covariates (regressors):
stato_salute=[0,0]; % Dummy variable, coding for normal status
altezza=1.64;
peso=59.6;
eta=32;
sesso=0; % Binary variable: 0=F, 1=M 
BMI=peso/altezza^2;
BSA=0.20247*altezza^0.725*peso^0.425;

%% ESTIMATION OF VAN CAUTER PARAMETERS FOR THE SUBJECT:
% Find the subject's 4 Van Cauter parameters using the regression models 
% identified in the previous exercise:
Volume=[1 BSA^3 eta^3]*bhat_Volume*1000; % y = H * p (yhat = X * bhat), V = D / (A + B)
Emivita_corto=[1 BMI BMI^3]*bhat_T12corto;
Emivita_lungo=[1 eta^3 altezza^3 eta peso^3 peso stato_salute BMI^3 BMI]*bhat_T12lungo;
Fraction=[1 BMI]*bhat_Fraction; % F = A / (A + B)

AUC=trapz(t,y);
A=Fraction/Volume; % A = F * D / V
B=(1-Fraction)/Volume; % B = D * (1 - F) / V
Alfa=log(2)/Emivita_corto;
Beta=log(2)/Emivita_lungo;
D=1; % Impulsive dose (bolus) in mg

% Calculate the impulse response g:
g= @(t) A*exp(-Alfa*t)+B*exp(-Beta*t);

% Assuming uniform sampling with step T, calculate gtilde by definition:
% integral(gtilde_m-th, m*T, (m+1)*T)
% Solving the integral yields:
% gtilde_m-th = A/alpha * (-exp(-alpha*(m+1)*T) + B/beta * (-exp(-beta*(m+1)*T) + exp(-beta*m*t)) where m = 0, ..., n-1
for m=1:length(t)-1
     gtilde(m+1)=integral(g,t(m),t(m+1));
end
gtilde(1)=integral(g,0,t(1));
gtilde=gtilde';

% Knowing that G is a lower triangular matrix, we derive it as:
G=tril(toeplitz(gtilde));

%% INSULIN SECRETION ESTIMATION USING RAW DECONVOLUTION:
% Ignores measurement errors: z = y + v = y => res = z - y = 0 and uhat = G^-1 * y
uhat_raw=(G^-1)*y;
y_raw=G*uhat_raw;
residui_raw=y_raw-y; % Notice they are null, exactly as expected

% Calculate the condition number k: expresses the maximum amplification 
% that the measurement error introduced on the data can undergo 
% when translating into an error on the uhat_raw estimate
k=norm(G)*norm(inv(G)); % By definition of condition number

figure('Name','Raw-Deconvolution')
subplot(3,1,1)
stairs(t,uhat_raw)
title('Input Signal Reconstruction')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')

subplot(3,1,2)
plot(t,G*uhat_raw)
hold on 
plot(t, y, '*')
legend('G*u', 'data'), title("Reconvolved signal")
% The reconvolved y fits the data too well, even though they are noisy.

subplot(3,1,3)
plot(t, residui_raw, '*');
ylim([-2 2])
title("Residuals")
% Null residuals => overfitting

%% INSULIN SECRETION ESTIMATION USING REGULARIZATION:
% It is no longer required that the virtual grid is as dense 
% as the sampling grid (z = y + v = G*u + v)
% Make the virtual grid denser (ideally we want a very dense grid 
% to properly reconstruct the unknown input u)
griglia_virtuale=[0:1:max(t)]; % Virtual grid from 0 to max(t) with step 1 
for n=1:length(griglia_virtuale)-1
    gtilde_reg(n+1)=integral(g,griglia_virtuale(n),griglia_virtuale(n+1));
end
gtilde_reg(1)=integral(g,0,griglia_virtuale(1));
gtilde_reg=gtilde_reg';
G_r=tril(toeplitz(gtilde_reg));

% Find the instants of vector t that coincide with the virtual grid:
for i=1:length(t)
    indici(i)=find(t(i)==griglia_virtuale);
end

% Find the G corresponding to the created virtual grid:
% (G_regularization = G)
G_regolarizzazione=G_r(indici-1,:);

% CV and SIGMA_V are copied from exercise 3
CV=0.04;
SIGMA_V=diag((CV.*(y+yb)).^2); % Constant CV error model
RSS=trace(SIGMA_V); % By definition

d=[1;-1;zeros(length(griglia_virtuale)-2,1)]; % Difference between consecutive values 
% m=1 -> consider only the first derivative
P=tril(toeplitz(d)); % Penalty matrix (contains input derivatives, 
% i.e., jumps from one value to the next)
sum_res2=0;

% REGULARIZATION AS GAMMA VARIES
gamma=[0 10 300 500 50000];
for i=1:length(gamma)
    ur(i)={pinv(G_regolarizzazione'*pinv(SIGMA_V)*G_regolarizzazione+gamma(i)*P'*P)*G_regolarizzazione'*pinv(SIGMA_V)*y};
end

tdati=t(2:end);
figure ('Name','Regularization as gamma varies')
stairs(t,[ur{1,1}(indici),ur{1,2}(indici),ur{1,3}(indici),ur{1,4}(indici),ur{1,5}(indici)])
legend({'0', '10', '300', '500', '50000'})
title('Regularization as gamma varies')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')
% Small gamma: model closely fits data, higher noise sensitivity 
% (overfitting) gamma=0, raw-deconvolution. Large gamma: smoother and 
% more stable solution, lower data fit, high penalization. 

% Search for optimal gamma 
gamma0=1000; % Choose an arbitrary gamma value

% TWOMEY'S DISCREPANCY CRITERION (Protected against zero/negative gamma)
gammaottimo_tw=fminsearch(@(gamma)discrepanza(G_regolarizzazione,SIGMA_V,P,y,max(1e-5,abs(gamma))),gamma0);

% Find uhat referring to the optimal gamma:
[uhat_tw,RSShat_tw,residuihat_tw]=discrepanza_post(G_regolarizzazione,SIGMA_V,P,y,gammaottimo_tw); 
discr=RSShat_tw-trace(SIGMA_V); % Turns out to be approximately zero -> correct
residuitw_standard=residuihat_tw./sqrt(diag(SIGMA_V));

figure('Name','Regularization with Twomey')
subplot(3,1,1)
stairs(t,uhat_tw(indici))
title('Input Signal Reconstruction')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')

subplot(3,1,2)
plot(t,G_regolarizzazione*uhat_tw)
hold on 
plot(t, y, '*')
xlabel('Time [min]')
legend('G*u', 'data')

subplot(3,1,3)
plot(t,residuitw_standard,"*")
title("Residuals")
hold on
plot(t,ones(1,length(t))*0)
ylim([-2 2])
title('Residuals')
xlabel('Time [min]')
% Residuals improve compared to before; data are followed less closely 
% and the signal appears smoother 

% GCV CRITERION (Protected against zero/negative gamma)
gammaottimo_gcv=fminsearch(@(gamma)crossvalidation(G_regolarizzazione,SIGMA_V,P,y,max(1e-5,abs(gamma))),gamma0);
% Find uhat referring to the optimal gamma:
[uhat_gcv,RSShat_gcv,residuihat_gcv]=discrepanza_post(G_regolarizzazione,SIGMA_V,P,y,gammaottimo_gcv);
residuigcv_standard=residuihat_gcv./sqrt(diag(SIGMA_V));

figure('Name','Regularization with GCV')
subplot(3,1,1)
stairs(t,uhat_gcv(indici))
title('Input Signal Reconstruction')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')

subplot(3,1,2)
plot(t,G_regolarizzazione*uhat_gcv)
hold on 
plot(t, y, '*')
legend('G*u', 'data')

subplot(3,1,3)
plot(t,residuigcv_standard,"*")
title("Residuals")
hold on
plot(t,ones(1,length(t))*0)
ylim([-2 2])
title('Residuals')
xlabel('Time [min]')

% MAXIMUM LIKELIHOOD CRITERION (Protected against zero/negative gamma)
gammaottimo_ml=fminsearch(@(gamma)mlikelihood(G_regolarizzazione,SIGMA_V,P,y,max(1e-5,abs(gamma))),gamma0);
% Find uhat referring to the optimal gamma:
[uhat_ml,RSShat_ml,residuihat_ml]=discrepanza_post(G_regolarizzazione,SIGMA_V,P,y,gammaottimo_ml);
residuiml_standard=residuihat_ml./sqrt(diag(SIGMA_V));

figure('Name','Regularization with Maximum Likelihood')
subplot(3,1,1)
stairs(t,uhat_ml(indici))
title('Input Signal Reconstruction')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')

subplot(3,1,2)
plot(t,G_regolarizzazione*uhat_ml)
hold on 
plot(t, y, '*')
legend('G*u', 'data')

subplot(3,1,3)
plot(t,residuiml_standard,"*")
title("Residuals")
hold on
plot(t,ones(1,length(t))*0)
ylim([-2 2])
title('Residuals')
xlabel('Time [min]')

%% BAYESIAN ESTIMATION
a=100; 
b=100e-6; % 0
% a and b chosen such that CV is approximately 4%
m=length(gtilde);
sinv=inv(SIGMA_V);
P=zeros(length(gtilde));
for i=1:length(gtilde)-1
    P(i,i)=1;
    P(i+1,i)=-1;
end
P(end,end)=1;

for i=1:2000
    lambda_2=gamrnd(a,b,1);
    Bw=inv(inv(P)'*G'*sinv*G*inv(P)+eye(m)*lambda_2);
    Aw=Bw*inv(P)'*G'*sinv*y;
    w(:,i)=(mvnrnd(Aw,(Bw+Bw')/2))'; % Posterior PDF of U
    u_bayes_t(:,i)=P\w(:,i); % U = W * P^-1
    a=a+m/2;
    b=(1/b+w(:,i)'*w(:,i)/2)^(-1);
end

u_bayes=u_bayes_t(:,i/2+1:end); % Take the second half of samples to avoid burn-in 

figure
subplot(1,2,1)
stairs(t,u_bayes)
title('Bayesian estimate of u')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')

percentili=prctile(u_bayes',[5 50 95]);
subplot(1,2,2)
stairs(t,percentili')
title('Bayesian estimate of u: 5th, 50th, 95th percentile')
xlabel('Time [min]')
ylabel('Concentration [pmol/ml]')
legend('5th perc','50th perc','95th perc')