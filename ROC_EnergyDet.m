% Basic example on the creation of a set of ROC curves for the Energy detector
% based on analytical expressions and MonteCarlo simulation
% Furher reading: 
%   - A Survey of Spectrum Sensing Algorithms for Cognitive Radio Applications, Yucek and Arslan, 
%   - Energy Detection for Spectrum Sensing in Cognitive Radio, S. Atapattu et al
%
% Germán Augusto Ramírez Arroyave
% Universidad Nacional de Colombia 
% CMUN 2017

close all; clear all; clc;
%% Creation of a set of plots based on the analytical expressions
N=31;                   % number of time samples
W=-60;                  % noise power level in dB
S=-65:2:-55;            % set of signal power levels in dB
minEnTh=-50;            % minimum energy threshold
maxEnTh=-30;            % maximum energy threshold

sigma_w_2=10.^(W/10);   % noise variance
sigma_s_2=10.^(S/10);   % signal variance
lambda=10.^(linspace(minEnTh,maxEnTh,101)/10); % power threshold for detection (linear scale)

x=0.5*lambda./sigma_w_2;            
xx=zeros(length(S),length(lambda));
for cont=1:length(sigma_s_2)
    xx(cont,:)=0.5*lambda./(sigma_w_2+sigma_s_2(cont));
end

Pf=gammainc(x,N,'upper');   %this is due to the way in which Matlab calculates the incomplete gamma function. Equivalent to: Pf=gammainc(x,N,'upper')/gamma(N) 
Pfa=qfunc((lambda-N*2*sigma_w_2)/(sqrt(N)*2*sigma_w_2));    %alternative (approximate) expression for Pf by means of central limit theorem
Pd=gammainc(xx,N,'upper');  %Pd=1-gammainc(xx,N,'lower');	%equivalent, use in case of numerical stability issues

semilogx(lambda,Pd,'linewidth',2);
grid on; title('Detection probability in terms of threshold'); 
xlabel('\lambda: Detection threshold value'); ylabel('Pd: detection probability');
legendStr=cellstr(num2str((S-W)', 'SNR= %-d dB')); legend(legendStr,'Location','SouthWest');

figure, semilogx(lambda,Pf,lambda,Pfa,'r','linewidth',2);
grid on; title(['False alarm probability in terms of threshold for SNR=', num2str(W), 'dB']); 
xlabel('\lambda: Detection threshold value'); ylabel('Pf: false alarm probability');
legend('Accurate','CLT approximation');

figure, semilogx(Pf,Pd,'-*','linewidth',2);
grid on; title('Analytical ROC curve for Energy detector'); axis([1e-6 1 0 1]);
xlabel('Pf: false alarm probability'); ylabel('Pd: detection probability');
legendStr=cellstr(num2str((S-W)', 'SNR= %-d dB')); legend(legendStr,'Location','SouthEast');

%% Monte Carlo simulation, approach 1: make some signal realizations varying the target PF and create an empirical ROC curve
L=10000;            %number of MonteCarlo simulations
h=ones(1,N);        %dummy transfer function (later iterations should consider channel model)
K=101;              %number of samples for the target Pf in the simulation

ind=1;              % initialize variables to be used within the loop
simPd=zeros(length(S),K);
thr=zeros(length(S),K);
En=zeros(length(S),K);

targPf=logspace(-6,1,K)/10;
for simPf=logspace(-6,1,K)/10
    thr(:,ind)=gammaincinv(simPf,N,'upper')*2*sigma_w_2; %Detection threshold according to the analytical expression for Pf
%    thr(:,ind)=(qfuncinv(simPf)+sqrt(N))*sqrt(N)*2*sigma_w_2; %Detection threshold according to the approximate expression for Pf (enable to see the differences)
    contPd=zeros(length(S),1);
    for cont=1:L %Do a series of Monte-Carlo (MC) runs
        n=sqrt(sigma_w_2)*randn(length(S),N)+1i*sqrt(sigma_w_2)*randn(length(S),N); 
        s=sqrt(sigma_s_2).'*randn(1,N)+1i*sqrt(sigma_s_2).'*randn(1,N);     %external product to create a set of signals
        y=bsxfun(@times,h,s)+n;                                             %fastest way to do this in Matlab
        En(:,ind)=sum(abs(y).^2,2);                                         %calculate energy for each realization
        contPd=contPd+(En(:,ind)>thr(:,ind));                               %count the number of true detections
    end
    simPd(:,ind)=contPd/L;                                                  %get empirical probability
    ind=ind+1;
end

figure; semilogx(targPf,En,targPf,thr,'k','linewidth',2);
grid on; title('Simulated Energy and ''optimum'' threshold value for a set of SNR');
xlabel('Pf: false alarm probability'); ylabel('Energy and threshold');
legendStr=cellstr(num2str((S-W)', 'Energy (SNR= %-d)'));legend(legendStr,'Threshold');

figure;semilogx(targPf,simPd,'o', Pf,Pd,'-*', 'linewidth',2); %generate the ROC curve
grid on; title('Simulated and Analytical ROC curve for Energy detector'); axis([1e-6 1 0.5 1]);
xlabel('Pf: false alarm probability'); ylabel('Pd: detection probability');
legendStr=cellstr(num2str((S-W)', 'SNR= %-d dB S')); legendStr=[legendStr; cellstr(num2str((S-W)', 'SNR= %-d dB A'))];
legend(legendStr,'Location','SouthEast');

%% Monte Carlo simulation, approach 2: make some signal realizations varying threshold and create empirical ROC curve
L=10000;            %number of MonteCarlo simulations
h=ones(1,N);        %dummy transfer function, later iterations can consider channel impairments
K=101;              %number of samples in the threshold sweep

ind=1;
simPd=zeros(K,length(S)); %simulated probabilities are arranged row-wise because of the format accepted by the plot function
simPf=zeros(K,length(S));
En=zeros(length(S),K);

for thr=10.^(linspace(minEnTh,maxEnTh,K)/10) % power threshold for detection (linear scale)
    contPd=zeros(length(S),1);
    contPf=zeros(length(S),1);
    contTrue=0;
    for cont=1:L %Do a series of Monte-Carlo (MC) runs
        n=sqrt(sigma_w_2)*randn(length(S),N)+1i*sqrt(sigma_w_2)*randn(length(S),N); 
        theta=randi(2)-1;       %(on average) half samples will be '1' and half will be '0'
        s=theta*sqrt(sigma_s_2).'*randn(1,N)+theta*1i*sqrt(sigma_s_2).'*randn(1,N);
        y=bsxfun(@times,h,s)+n;             
        En(:,ind)=sum(abs(y).^2,2);
        if theta
            contTrue=contTrue+1;
            contPd=contPd+(En(:,ind)>thr);	%count the number of true detections
        else
            contPf=contPf+(En(:,ind)>thr);  %count the number of false detections
        end
    end
%     simPd(ind,:)=2*contPd/L;          %get empirical probability (2 factor: statistically L/2 true and false examples)
%     simPf(ind,:)=2*contPf/L;          %NOTE: it seems that L=10000 is not large enough for this assumption to hold
	simPd(ind,:)=contPd/contTrue;       %get empirical probability takinng int account the actual number of true examples
    simPf(ind,:)=contPf/(L-contTrue);	
    ind=ind+1;    
end

figure;semilogx(simPf,simPd,'o', Pf,Pd,'-*', 'linewidth',2); %generate the ROC curve
grid on; title('Simulated and Analytical ROC curve for Energy detector'); axis([10e-3 1 0.5 1]);
xlabel('Pf: false alarm probability'); ylabel('Pd: detection probability');
legendStr=cellstr(num2str((S-W)', 'SNR= %-d dB S')); legendStr=[legendStr; cellstr(num2str((S-W)', 'SNR= %-d dB A'))]; 
legend(legendStr,'Location','SouthEast');