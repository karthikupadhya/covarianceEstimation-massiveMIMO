% Simulation code for "Covariance Matrix Estimation for Massive MIMO"
% by Karthik Upadhya and Sergiy Vorobyov
% Published in IEEE Signal Processing Letters, vol. 25, no. 4, pp. 546-550, April 2018.
%
% This code generates Figures 1-3.
% ref14 is the method in [14] - E. Björnson, L. Sanguinetti and M. Debbah,
% "Massive MIMO with imperfect channel covariance information," 50th
% Asilomar Conf. Signals, Syst. Computers, Pacific Grove, CA, 2016, pp. 974-978.


clc
clear all
close all

numCell                     = 7;
numUser                     = 10;
%%%%%%%% - Plot Figs 1-3 with covariance matrix estimated from the method in [14] - %%%%%%%%%
covarianceEstimationMethod  = 'ref14';
channelEstimationMethodRange= {'leastSquares','exactcovariancematrix','estimatedcovariancematrix'};
nRRange = [10,50,100,200,250]; % Range of values for the parameter N_r in [14]
for nn = 1:length(channelEstimationMethodRange);
    channelEstimationMethod = channelEstimationMethodRange{nn};
    for ii = 1:length(nRRange)
        nR = nRRange(ii);
        nQ = 10 * nR;
        [sumAchRate(nn,ii),mseVal(nn,ii),matrixErrorNorm(nn,ii)] = code123(covarianceEstimationMethod,channelEstimationMethod,(nQ + numCell * nR));        
    end
end
numPilotRange = nRRange * (numCell + 10) * numUser;

figure(1); plot(numPilotRange,repmat(10*log10(mseVal(1,1)),1,length(nRRange)),'cd-'); hold on;
figure(1); plot(numPilotRange,repmat(10*log10(mseVal(2,1)),1,length(nRRange)),'r-');
figure(1); plot(numPilotRange,10*log10(mseVal(3,:)),'ko-');

figure(2); plot(numPilotRange,repmat(sumAchRate(1,1),1,length(nRRange)),'cs-'); hold on;
figure(2); plot(numPilotRange,repmat(sumAchRate(2,1),1,length(nRRange)),'r-');
figure(2); plot(numPilotRange,sumAchRate(3,:),'ko-');

figure(3); plot(numPilotRange,10*log10(matrixErrorNorm(3,:)),'ko-'); hold on; drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% - Plot Figs 1-3 with covariance matrix estimated using the proposed method with regular pilots - %%%%%%%%%
covarianceEstimationMethod  = 'proposed-regularpilot';
channelEstimationMethodRange= {'estimatedcovariancematrix'};
nRRange = [10,50,100,200,250]; % Range of values for the parameter N_r in [14]
for nn = 1:length(channelEstimationMethodRange);
    channelEstimationMethod = channelEstimationMethodRange{nn};
    for ii = 1:length(nRRange)
        nR = nRRange(ii);
        nQ = 10 * nR;
        N  = ( nQ + numCell * nR ) / 2;
        [sumAchRate(nn,ii),mseVal(nn,ii),matrixErrorNorm(nn,ii)] = code123(covarianceEstimationMethod,channelEstimationMethod,N);        
    end
end

figure(1); plot(numPilotRange,10*log10(mseVal(1,:)),'ms-');

figure(2); plot(numPilotRange,sumAchRate(1,:),'ms-');

figure(3); plot(numPilotRange,10*log10(matrixErrorNorm(1,:)),'ms-'); drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% - Plot Figs 1-3 with covariance matrix estimated using the proposed method with staggered pilots - %%%%%%%%%
covarianceEstimationMethod  = 'proposed-staggeredpilot';
channelEstimationMethodRange= {'estimatedcovariancematrix'};
nRRange = [10,50,100,200,250]; % Range of values for the parameter N_r in [14]
for nn = 1:length(channelEstimationMethodRange);
    channelEstimationMethod = channelEstimationMethodRange{nn};
    for ii = 1:length(nRRange)
        nR = nRRange(ii);
        nQ = 10 * nR;
        N  = ( nQ + numCell * nR ) / 2;
        [sumAchRate(nn,ii),mseVal(nn,ii),matrixErrorNorm(nn,ii)] = code123(covarianceEstimationMethod,channelEstimationMethod,N);        
    end
end

figure(1); plot(numPilotRange,10*log10(mseVal(1,:)),'b*-');

figure(2); plot(numPilotRange,sumAchRate(1,:),'b*-');

figure(3); plot(numPilotRange,10*log10(matrixErrorNorm(1,:)),'b*-');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); 
ylabel('MSE of the Channel Estimate (dB)'); 
xlabel('Number of UL pilots'); 
legend('LS','LMMSE','Method in [14]','Proposed Method - L_t = 1','Proposed Method - Staggered Pilot');

figure(2); 
ylabel('Sum Rate'); 
xlabel('Number of UL pilots'); 
legend('LS','LMMSE','Method in [14]','Proposed Method - L_t = 1','Proposed Method - Staggered Pilot');

figure(3); 
ylabel('MSE of the Covariance Matrix Estimate (dB)'); 
xlabel('Number of UL pilots'); 
legend('Method in [14]','Proposed Method - L_t = 1','Proposed Method - Staggered Pilot');
