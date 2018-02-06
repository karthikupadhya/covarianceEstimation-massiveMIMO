% Simulation code for "Covariance Matrix Estimation for Massive MIMO"
% by Karthik Upadhya and Sergiy Vorobyov
% Submitted to IEEE SPL.
%
% This code generates Figures 1-3.

clc
clear all
close all

numCell                     = 7;
numUser                     = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covarianceEstimationMethod  = 'ref14';
channelEstimationMethodRange= {'leastSquares','exactcovariancematrix','estimatedcovariancematrix'};
nRRange = [10,50,100,200,250];
for nn = 1:length(channelEstimationMethodRange);
    channelEstimationMethod = channelEstimationMethodRange{nn};
    for ii = 1:length(nRRange)
        nR = nRRange(ii);
        nQ = 10 * nR;
        [sumAchRate(nn,ii)] = code4(covarianceEstimationMethod,channelEstimationMethod,(nQ + numCell * nR));        
    end
end
numPilotRange = nRRange * (numCell + 10) * numUser;

figure(4); plot(numPilotRange,repmat(sumAchRate(1,1),1,length(nRRange)),'cs-'); hold on;
figure(4); plot(numPilotRange,repmat(sumAchRate(2,1),1,length(nRRange)),'r-');
figure(4); plot(numPilotRange,sumAchRate(3,:),'ko-');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covarianceEstimationMethod  = 'proposed-staggeredpilot';
channelEstimationMethodRange= {'exactcovariancematrix','estimatedcovariancematrix'};
nRRange = [10,50,100,200,250];
for nn = 1:length(channelEstimationMethodRange);
    channelEstimationMethod = channelEstimationMethodRange{nn};
    for ii = 1:length(nRRange)
        nR = nRRange(ii);
        nQ = 10 * nR;
        N  = ( nQ + numCell * nR ) / 2;
        [sumAchRate(nn,ii)] = code4(covarianceEstimationMethod,channelEstimationMethod,N);        
    end
end

figure(4); plot(numPilotRange,repmat(sumAchRate(1,1),1,length(nRRange)),'b--');
figure(4); plot(numPilotRange,sumAchRate(2,:),'b*-');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); 
ylabel('Sum Rate'); 
xlabel('Number of UL pilots'); 
legend('LS','LMMSE-Regular Pilot','Method in [14]','LMMSE-Staggered Pilot','Proposed Method - Staggered Pilot');