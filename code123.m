% Simulation code for "Covariance Matrix Estimation for Massive MIMO"
% by Karthik Upadhya and Sergiy Vorobyov
% Published in IEEE Signal Processing Letters, vol. 25, no. 4, pp. 546-550, April 2018.
%
% This code generates Figures 1-3.

%%%% - List of Arguments - %%%%
% Argument 1  : covarianceEstimationMethod
% Description : Select method for estimating the covariance matrices
% Valid choices for Argument 1 : 
% (1) ref14 - Use method in [14] 
% (2) proposed-regularpilot - Use the proposed method with regular pilots 
% (3) proposed-staggeredpilot - Use the proposed method with staggered pilots

% Argument 2 : channelEstimationMethod
% Description: Select method for estimating the channel
% Valid choices for Argument 2  : 
% (1) leastSquares - compute the LS channel estimate
% (2) exactCovarianceMatrix - compute the exact LMMSE channel estimate
% (3) estimatedCovarianceMatrix - use estimated covariance matrices to compute the LMMSE channel estimate

% Argument 3 : numCoherenceBlock 
% Description: number of coherence blocks used for estimating the covariance matrix (denoted as N in the paper).

%%%% - List of Outputs - %%%%
% Output 1 : sumAchRate
% Description : Sum achievable rate

% Output 2    : channelNmseOut
% Description : Normalized MSE of channel estimate

% Output 3    : matrixErrorNormOut
% Description : MSE of covariance matrix estimate

%%%% - Important note - %%%%
% ref14 is the method in E. Björnson, L. Sanguinetti and M. Debbah, "Massive MIMO with imperfect channel covariance information," 
% 50th Asilomar Conf. Signals, Syst. Computers, Pacific Grove, CA, 2016, pp. 974-978.
%
% To get the same results in the paper. Uncomment the rng function in line 79

function [sumAchRate,channelNmseOut,matrixErrorNormOut] = code123(covarianceEstimationMethod,channelEstimationMethod,numCoherenceBlock)

% - Simulation Parameters - %
numTrial                    = 1000; % Number of monte carlo trials
% - System Parameters - %
numBsAntenna                = 100; % Number of antennas at the base station (notation in paper - M)
numCell                     = 7; % Number of cells in the network (notation in paper - L)
numUser                     = 10; % - Number of users per cell (notation in paper - K)
pP                          = 1; % Pilot power (denoted as \rhoP^2 in the paper)
pD                          = 1; % Data power (denoted as \rhoD^2 in the paper)
Nr                          = numCoherenceBlock / (10 + numCell); % The parameter N_r in [14]
Nq                          = 10 * Nr; % The parameter N_q in [14]

% - Channel Parameters - %
userRadius                  = 120; % User radius in m (Users located on a ring around the BS)
cellRadius                  = 150; % Distance between BSs in m
normAntennaSpacing          = 1/2; % Antenna spacing normalized by wavelength
ulSlotLength                = 200; % Number of symbols in the uplink time-slot (notation in paper - C_u)
numStationaryBlock          = 25000;% Number of coherence blocks over which the channel is stationary (notation in paper - \tau_s)

% - Performance metrics - %
channelMse                 = zeros(numTrial,numUser);

% - Generate user locations and parameters for the covariance matrices - %
center  = zeros(1,numCell);
for ii = 1:numCell
    if ii == 1
        center(ii) = 0; % Location of reference BS
    else
        center(ii) = ceil((ii-1)/6) * 2 * sqrt(3)/2 * cellRadius * exp(1i * (pi/3 * mod(ii-1,6) ) ); % Location of the remaining BSs
    end
    userLocations  = userRadius * exp(1i * 2 * pi / numUser * (0:numUser-1)) + center(ii); % User locations in cell ii
    rxPowerInDb    = 78.7 - 37.6 * log10(abs(userLocations)); % Received power in dB corresponding to each user in cell ii
    rxPower(ii,:)  = 10.^(rxPowerInDb/10); % Received power corresponding to each user in cell ii
    meanAngle(ii,:)= angle(userLocations) * 180 / pi; % Mean angle of channel cluster
end

noiseVar        = 1; % Noise variance (notation in paper - \sigma^2)
angleSpread     = 20 * ones(numCell,numUser); % Angular spread of the channel cluster in Degrees ( sys.numCell x sys.numUser )


% - Generate Covariance Matrix - %
% Covariance matrix is computed as E[\alpha * a(\theta) * a(\theta)^H] where E[.] is
% the expectation operator. The expectation is computed by integrating over
% the variables \theta and \alpha.
covarianceMatrix = cell(numCell,numUser); 
[covarianceMatrix{:}] = deal(zeros(numBsAntenna)); % Initialize covariance matrix to all zeros
sqrtCovarianceMatrix = covarianceMatrix;
for jj = 1:numCell
    for mm = 1:numUser
        % - Compute covariance matrix of user (jj,mm) at the reference BS - % 
        userMeanAngle   = meanAngle(jj,mm);
        userAngleSpread = angleSpread(jj,mm);
        
        integralSpacing = 0.01;
        thetaRange      = -180:integralSpacing:360; % Account for wrap around
        pTheta          = 1/userAngleSpread * ((thetaRange >= userMeanAngle - userAngleSpread/2) & (thetaRange <= userMeanAngle + userAngleSpread/2));
        
        aVector         = exp(1i * 2 * pi * normAntennaSpacing * (0:numBsAntenna-1).' * cos(pi/180 * thetaRange(pTheta>0)));
        covarianceMatrix{jj,mm}         = rxPower(jj,mm) * (aVector * diag(pTheta(pTheta>0)) * aVector') * integralSpacing;
        sqrtCovarianceMatrix{jj,mm}     = covarianceMatrix{jj,mm}^(1/2);
        % covarianceMatrix{jj,mm} contains the covariance matrix of user (jj,mm) at
        % the reference BS.
    end
end

% - Simulation Begins - %
for ii = 1:numTrial
    % rng(ii); % Uncomment to get the results in the paper.
    
    % - Generate Channel Vectors- %
    channelMatrix       = cell(numCell,numCoherenceBlock);
    [channelMatrix{:}]  = deal(zeros(numBsAntenna,numUser));
    for jj = 1:numCell
        for mm = 1:numUser
            % - Generate N realizations of the channel vector for user (jj,mm) - %
            xii                 = complex(normrnd(0,1/sqrt(2),numBsAntenna,numCoherenceBlock),normrnd(0,1/sqrt(2),numBsAntenna,numCoherenceBlock));
            Hii                 = sqrtCovarianceMatrix{jj,mm} * xii;
            for nn = 1:numCoherenceBlock
                channelMatrix{jj,nn}(:,mm) = Hii(:,nn); % Store the realizations of the channel vector into the variable channelMatrix
            end
        end
    end
    
    % - Obtain covariance matrices for channel estimation - %
    % estIndivCovarianceMatrix - Variable containing the estimates of the individual user covariance matrices (denoted as R_{jcu} in the paper)
    % estSumCovarianceMatrix - Variable containing the estimates of the sum covariance matrices (denoted as Q_{jk}) in the paper_
    switch lower(channelEstimationMethod)
        case 'leastsquares'
        % - estIndivCovarianceMatrix and estIndivCovarianceMatrix is the identity matrix for the LS method - %
        estIndivCovarianceMatrix = cell(numCell,numUser); 
        [estIndivCovarianceMatrix{:}] = deal(speye(numBsAntenna)); 
        estSumCovarianceMatrix = cell(1,numUser);
        [estSumCovarianceMatrix{:}] = deal(speye(numBsAntenna));
        case 'exactcovariancematrix'
            % - estIndivCovarianceMatrix and estIndivCovarianceMatrix contain the exact covariance matrices - %
            estIndivCovarianceMatrix    = covarianceMatrix;
            estSumCovarianceMatrix      = cell(1,numUser);
            
            %%% - Estimation of Sum Covariance Matrix - %%%
            for mm = 1:numUser
                Rsum = zeros(numBsAntenna);
                for ll = 1:numCell
                    Rsum = Rsum + estIndivCovarianceMatrix{ll,mm};
                end
                Rsum = Rsum + noiseVar/numUser * eye(numBsAntenna);
                estSumCovarianceMatrix{1,mm} = Rsum;
            end
        case 'estimatedcovariancematrix'
        % - estIndivCovarianceMatrix and estIndivCovarianceMatrix contain the estimated the covariance matrices - %
        if strcmpi(covarianceEstimationMethod,'ref14')
            % - Compute covariance matrix estimates as in ref14 - %
            exactCovarianceMatrix       = covarianceMatrix;
            estIndivCovarianceMatrix    = cell(numCell,numUser);
            estSumCovarianceMatrix      = cell(1,numUser);
            
            %%% - Computing \widehat{Q}_{jm} - %%%
            for nn = 1:Nq
                % - Compute channel estimate \widehat{h}_{jjm} - %
                hHat1{nn} = zeros(numBsAntenna,numUser);
                for ll = 1:numCell
                    hHat1{nn} = hHat1{nn} + channelMatrix{ll,nn};
                end
                hHat1{nn} = hHat1{nn} + complex(normrnd(0,sqrt(noiseVar/(2*numUser)),numBsAntenna,numUser),normrnd(0,sqrt(noiseVar/(2*numUser)),numBsAntenna,numUser));
            end
            
            for mm = 1:numUser
                % - Compute \widehat{Q}_{jm} = 1/N * \sum_{n=0}^{N-1} \widehat{h}_{jjk} \widehat{h}_{jjm}^H - %
                tempH1 = zeros(numBsAntenna,Nq);
                for nn = 1:Nq
                    tempH1(:,nn)         = hHat1{nn}(:,mm);
                end
                estSumCovarianceMatrix{1,mm} = tempH1 * tempH1' / Nq ; 
            end
            
            %%% - Computing \widehat{R}_{jcu} - %%%
            jj = 1; % Evaluate \widehat{R}_{jcu} only for the users in the reference cell.
            cnt = 0;
            for nn = Nq+1:numCoherenceBlock % Use the Nr = N - Nq coherence blocks for computing \widehat{R}_{jjm}
                if jj == mod(nn-Nq-1,numCell)+1
                    % - Compute channel estimate of users in the reference cell - %
                    cnt = cnt + 1;
                    hHat2{jj,cnt} = zeros(numBsAntenna,numUser);
                    for ll = 1:numCell
                        if ll~=jj
                            hHat2{jj,cnt} = hHat2{jj,cnt} + channelMatrix{ll,nn};
                        end
                        hHat2{jj,cnt} = hHat2{jj,cnt} + complex(normrnd(0,sqrt(noiseVar/(2*numUser)),numBsAntenna,numUser),normrnd(0,sqrt(noiseVar/(2*numUser)),numBsAntenna,numUser));
                    end
                end
            end
            
            % - Compute sample covariance matrix and regularize according to [14] - %
            for mm = 1:numUser
                % - Compute Qexact = Q_{jm} for calculating the regularization coefficient - %
                Qexact = zeros(numBsAntenna);
                for ll = 1:numCell
                    Qexact = Qexact + exactCovarianceMatrix{ll,mm};
                end
                Qexact = Qexact + noiseVar/numUser * eye(numBsAntenna);
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                tempH2 = zeros(numBsAntenna,Nr);
                cnt    = 0;
                for nn = Nq+1:numCoherenceBlock
                    if jj == mod(nn-Nq-1,numCell)+1
                        cnt         = cnt + 1;
                        tempH2(:,cnt)= hHat2{jj,cnt}(:,mm); % - Load the channel estimate into a dummy variable
                    end
                end
                estIndivCovarianceMatrix{jj,mm} = estSumCovarianceMatrix{1,mm} - tempH2 * tempH2'/ cnt ; % Compute \widehat{R}_{jjm}
                muRange = 0:0.05:1; % \mu is the regularization coefficient for R_{jjm}. Evaluate MSE for \mu in the range [0:1] in steps of 0.05
                etaRange = 0:0.05:1; % \eta is the regularization coefficient for Q_{jm}. Evaluate MSE for \eta in the range [0:1] in steps of 0.05
                for nn = 1:length(muRange);
                    muVal = muRange(nn);
                    for pp = 1:length(etaRange);
                        % Compute the MSE for different values of \eta and mu in the range [0,1] in steps of 0.05 - %
                        etaVal           = etaRange(pp);
                        tempQMatrix      = etaVal * estSumCovarianceMatrix{1,mm} + (1 - etaVal) * diag(diag(estSumCovarianceMatrix{1,mm})); % Regularize \widehat{Q}_{jm}
                        tempRMatrix      = muVal * estIndivCovarianceMatrix{jj,mm} + (1 - muVal) * diag(diag(estIndivCovarianceMatrix{jj,mm})); % Regularize \widehat{R}_{jjm}
                        tempW            = tempRMatrix/tempQMatrix; % Compute W = \widehat{R} \widehat{Q}^{-1}
                        fVal(nn,pp)      = real(trace((eye(numBsAntenna) - tempW - tempW') * exactCovarianceMatrix{jj,mm} ) + trace( tempW *  Qexact * tempW' )); % Compute MSE with W
                    end
                end
                [~,idx]                 = min(fVal(:));
                [muMinIdx,etaMinIdx]    = ind2sub(size(fVal),idx); % Select pair of \mu and \eta that results in the least MSE
                etaOpt                  = etaRange(etaMinIdx); 
                muOpt                   = muRange(muMinIdx);
                if jj == 1
                    estSumCovarianceMatrix{1,mm}    = etaOpt * estSumCovarianceMatrix{1,mm} + (1 - etaOpt) * diag(diag(estSumCovarianceMatrix{1,mm})); % Regularize \widehat{Q}_{jm} with optimal value of \eta
                end
                estIndivCovarianceMatrix{jj,mm}     = muOpt * estIndivCovarianceMatrix{jj,mm} + (1 - muOpt) * diag(diag(estIndivCovarianceMatrix{jj,mm})); % Regularize \widehat{R}_{jjm} with optimal value of \mu
            end
        else
            % - Proposed method for estimating the covariance matrices - %
            randPhase                   = exp(1i * 2 * pi * rand(numCell,numCoherenceBlock)); % Generate N realizations of the random variable e^{j\Theta}
            switch lower(covarianceEstimationMethod)
                case 'proposed-staggeredpilot'
                    % - Use the proposed method with staggered pilots - %
                    hHat1 = cell(numCell,numCoherenceBlock); % hHat1 is the placeholder for channel estimates from the first pilot
                    [hHat1{:}] = deal(zeros(numBsAntenna,numUser));
                    hHat2 = hHat1; % hHat2 is the placeholder for channel estimates from the second pilot
                    for ll = 1:numCell
                        txData1    = complex(normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock),normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock)); % txData1 contains the data transmitted by users in the cells not in L_t when the first pilot is transmitted
                        txData2    = complex(normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock),normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock)); % txData1 contains the data transmitted by users in the cells not in L_t when the second pilot is transmitted
                        noiseVec1  = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock)); % noiseVec1 contains the AWGN in the observations corresponding to the first pilot transmission
                        noiseVec2  = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock)); % noiseVec2 contains the AWGN in the observations corresponding to the second pilot transmission
                        for nn = 1:numCoherenceBlock
                            txDataIdx    = (nn-1) * numUser + 1 : nn * numUser;
                            hInterf1     = sqrt(pD) / (sqrt(pP*numUser)) * cell2mat(channelMatrix([1:ll-1 ll+1:end],nn).') * txData1(:,txDataIdx) + 1/(sqrt(pP*numUser))* noiseVec1(:,txDataIdx); % hInterf1 contains the component of the channel estimate corresponding to interfering pilot and data transmissions when the first pilot sequence is transmitted
                            hHat1{ll,nn} = channelMatrix{ll,nn} + hInterf1; % hHat1{ll,nn} contains the channel estimates of users in cell ll in coherence block nn obtained from the first pilot sequence
                            hInterf2     = conj(randPhase(ll,nn)) * sqrt(pD) / (sqrt(pP*numUser)) * cell2mat(channelMatrix([1:ll-1 ll+1:end],nn).') * txData2(:,txDataIdx) + conj(randPhase(ll,nn)) / sqrt(pP*numUser) * noiseVec2(:,txDataIdx); % hInterf2 contains the component of the channel estimate corresponding to interfering pilot and data transmissions when the second pilot sequence is transmitted
                            hHat2{ll,nn} = channelMatrix{ll,nn} + hInterf2; % hHat2{ll,nn} contains the channel estimates of users in cell ll in coherence block nn obtained from the second pilot sequence (after compensating for the random phase shift)
                        end
                    end
                case 'proposed-regularpilot'
                    % - Use the proposed method with regular pilots - %
                    hHat1 = cell(numCell,numCoherenceBlock); % hHat1 is the placeholder for channel estimates from the first pilot
                    [hHat1{:}] = deal(zeros(numBsAntenna,numUser));
                    hHat2 = hHat1; % hHat2 is the placeholder for channel estimates from the second pilot
                    for jj = 1:numCell
                        noiseVec1   = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock)); % noiseVec1 contains the AWGN in the observations corresponding to the first pilot transmission
                        noiseVec2   = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock)); % noiseVec2 contains the AWGN in the observations corresponding to the second pilot transmission
                        for nn = 1:numCoherenceBlock
                            txDataIdx    = (nn-1) * numUser + 1 : nn * numUser;
                            for ll = 1:numCell
                                hHat1{jj,nn} = hHat1{jj,nn} + channelMatrix{ll,nn}; % hHat1{jj,nn} contains the channel estimates of users in cell jj in coherence block nn obtained from the first pilot sequence. This channel estimate is a sum of the channel vectors of all users transmitting pilots alongside users in cell jj. 
                                hHat2{jj,nn} = hHat2{jj,nn} + conj(randPhase(jj,nn)) * randPhase(ll,nn) * channelMatrix{ll,nn}; % hHat2{jj,nn} contains the channel estimates of users in cell jj in coherence block nn obtained from the second pilot sequence. This channel estimate is a sum of the channel vectors of all users transmitting pilots alongside users in cell jj. 
                            end
                            hHat1{jj,nn} = hHat1{jj,nn} + 1/sqrt(numUser) * noiseVec1(:,txDataIdx); % Add noise to the channel estimate
                            hHat2{jj,nn} = hHat2{jj,nn} + conj(randPhase(jj,nn))/sqrt(numUser) * noiseVec2(:,txDataIdx); % Add noise to the channel estimate
                        end
                    end
            end
            
            %%% - Compute \widehat{R}_{jcu} using the proposed method - %%%
            estIndivCovarianceMatrix       = cell(numCell,numUser);
            for jj = 1:numCell
                for mm = 1:numUser
                    tempH1 = zeros(numBsAntenna,numCoherenceBlock);
                    tempH2 = tempH1;
                    for nn = 1:numCoherenceBlock
                        tempH1(:,nn) = hHat1{jj,nn}(:,mm); % Load channel estimates into a dummy variable
                        tempH2(:,nn) = hHat2{jj,nn}(:,mm); % Load channel estimates into a dummy variable
                    end
                    tempCovMat  = tempH1 * tempH2' / numCoherenceBlock; % Obtain \widehat{R}_{jcu} through the cross correlation of \widehat{h}_{jcu}^{(1)} and \widehat{h}_{jcu}^{(2)}
                    % - Regularize the estimated covariance matrix by replacing it with the closest positive semidefinite matrix in Frobenius norm - %
                    tempCovMat  = (tempCovMat + tempCovMat')/2;
                    [U,D]       = eig(tempCovMat);
                    d           = real(diag(D));
                    estIndivCovarianceMatrix{jj,mm} = U(:,d>0) * sparse(diag(d(d>0))) * U(:,d>0)'; % estIndivCovarianceMatrix{jj,mm} contains the regularized covariance matrix
                end
            end
            
            %%% - Compute \widehat{Q}_{jm} - %%%
            % For Figs. 1-3 and with the proposed method, covariance
            % matrices are estimated using either regular or staggered
            % pilots. However, coherence blocks which are not used for
            % covariance matrix estimation are assumed to contain regular pilots.
            % Hence the estimation of \widehat{Q}_{jm} takes the following
            % form.
            estSumCovarianceMatrix         = cell(1,numUser);
            for mm = 1:numUser
                Rsum            = zeros(numBsAntenna);
                for jj = 1:numCell
                    Rsum        = Rsum + estIndivCovarianceMatrix{jj,mm};
                end
                Rsum = Rsum + noiseVar/numUser * eye(numBsAntenna);
                estSumCovarianceMatrix{1,mm} = Rsum;
            end
        end
    end
    
    % - Generate channel vectors for obtaining the channel estimate - %
    channelVec = cell(numCell,1);
    for jj = 1:numCell
        for mm = 1:numUser;
            channelVec{jj,1}(:,mm) = sqrtCovarianceMatrix{jj,mm} * complex(normrnd(0,1/sqrt(2),numBsAntenna,1),normrnd(0,1/sqrt(2),numBsAntenna,1)); % channelVec{jj,1} contains the channel vectors of users in cell jj
        end
    end
    
    % - Generate channel estimate for regular pilot - %
    hHat = zeros(numBsAntenna,numUser);
    for mm = 1:numUser;
        % User mm is a user in the reference cell (cell 1 is the reference cell)
        for jj = 1:numCell
            hHat(:,mm) = hHat(:,mm) + channelVec{jj,1}(:,mm); % The channel estimate of user mm contains channel vectors from users in other cells (cell jj in this case) that reuse the same pilots
        end
    end
    hHat = hHat + complex(normrnd(0,sqrt(noiseVar / (2 * numUser)),numBsAntenna,numUser),normrnd(0,sqrt(noiseVar / (2 * numUser)),numBsAntenna,numUser)); % Add noise to the channel estimate
    
    % - Compute LMMSE channel estimate - %
    channelEstimate     = cell(1,1);
    [channelEstimate{:}]= deal(zeros(numBsAntenna,numUser));
    for mm = 1:numUser
        Rsum            = estSumCovarianceMatrix{1,mm};
        R0              = estIndivCovarianceMatrix{1,mm};
        channelEstimate{1}(:,mm) = R0 * (Rsum \ hHat(:,mm)); % Compute LMMSE channel estimate of user mm in the reference cell (cell 1 is the reference cell)
    end
    
    % - Compute RZF combining vector - %
    term1 = channelEstimate{1} * channelEstimate{1}';
    term2 = noiseVar * eye(numBsAntenna);
    combiningVec        = ((term1 + term2) \ channelEstimate{1}); % combiningVec contains the RZF combining vector
    
    % - Evaluate terms for SINR calculation - %
    for mm = 1:numUser
        interfTerm1(ii,mm) = 0; % interfTerm1 is the first term in (21) in the paper. Note that \rhoD^2 = \rhoP^2 = 1 with regular pilots
        for ll = 1:numCell
            for pp = 1:numUser
                interfTerm1(ii,mm) = interfTerm1(ii,mm) + abs( combiningVec(:,mm)' * channelVec{ll,1}(:,pp) ).^2; 
            end
        end
        interfTerm2(ii,mm) = combiningVec(:,mm)' * channelVec{1,1}(:,mm); % interfTerm2 is the second term in (21) in the paper. 
        interfTerm3(ii,mm) = norm(combiningVec(:,mm))^2 * noiseVar; % interfTerm3 is the third term in (21) in the paper.
    end
    
    for mm = 1:numUser
        for ll = 1:numCell
            for pp = 1:numUser
                interfTermDebug1{mm,ll,pp}(ii,1) = abs( combiningVec(:,mm)' * channelVec{ll,1}(:,pp) )^2; % The elements of InterfTermDebug1 correspond to the inner terms in the summation in the first term in (21) in the paper.
            end
        end
        interfTermDebug2{mm}(ii,1) = combiningVec(:,mm)' * channelVec{1,1}(:,mm); % InterfTermDebug2 corresponds to the second term in (21) in the paper.
        interfTermDebug3{mm}(ii,1) = norm(combiningVec(:,mm))^2; % InterfTermDebug3 corresponds to the third term in (21) in the paper
    end
    
    % - Compute NMSE of channel estimate - %
    for mm = 1:numUser
        channelMse(ii,mm) = norm(channelEstimate{1}(:,mm) - channelVec{1,1}(:,mm))^2 / real(trace(covarianceMatrix{1,mm}));
    end
    
    % - Compute square root of the MSE of covariance matrix estimate - %
    for mm = 1:numUser
        matrixErrorNorm(ii,mm) = norm(estIndivCovarianceMatrix{1,mm} - covarianceMatrix{1,mm},'fro');
    end
    
end

% - Compute prelog factor - %
switch lower(covarianceEstimationMethod)
    case 'ref14'
        prelogFactor = (1 - numUser/ulSlotLength - Nr * numUser * numCell / (numStationaryBlock*ulSlotLength) ); % Prelog factor for the method in [14].
    case 'proposed-staggeredpilot'
        prelogFactor1 = 2 * numCoherenceBlock * numUser * ( numCell - 1) / (ulSlotLength * numStationaryBlock); % Prelog factor in the N coherence blocks containing staggered pilots
        prelogFactor2 = (1 - numUser / ulSlotLength + numCoherenceBlock * numUser / (ulSlotLength * numStationaryBlock) - 2 * numCell * numCoherenceBlock * numUser / (ulSlotLength * numStationaryBlock) ); % Prelog factor in the \tau_s - N coherence blocks containing regular pilots
    case 'proposed-regularpilot'
        prelogFactor = (1 - numUser/ulSlotLength - numCoherenceBlock * numUser / (numStationaryBlock*ulSlotLength) ); % Prelog factor when all the pilots transmitted are regular pilots
end

% - Compute SINR and achievable rate - %
switch lower(covarianceEstimationMethod)
    case {'ref14','proposed-regularpilot'}
        % - The following code evaluates Equation (15) in Reference 14 - %
        interferencePower = mean(interfTerm1,1) - abs(mean(interfTerm2,1)).^2 + mean(interfTerm3,1);
        sigPower          = abs(mean(interfTerm2,1)).^2;
        sinr              = sigPower./interferencePower;
        achRate           = prelogFactor * log2(1 + sinr);
    case 'proposed-staggeredpilot'
        % - The following code evaluates Equation (21) in published manuscript - %
        for mm = 1:numUser
            for ll = 1:numCell
                for pp = 1:numUser
                    interfVal1{mm}(ll,pp) = mean(interfTermDebug1{mm,ll,pp}(:,1));
                end
            end
            
            t1(mm) = sum(reshape(interfVal1{mm},1,[])); % Term1 in (21) in paper.
            
            t2(mm) = abs( mean(interfTermDebug2{mm}(:,1)) ).^2; % Term2 in (21) in paper.
            
            t3(mm) = noiseVar / pD * mean(interfTermDebug3{mm}(:,1)); % Term3 in (21) in paper.
            
            for ll = 1:numCell-1
                t4(ll,mm) = (pP - pD) / pD * sum(interfVal1{mm}(ll+1,:)); % Term4 in (21) in paper.
            end
        end
        t4     = mean(t4);
        
        sinrStag      = t2 ./ (t1 - t2 + t3 + t4);
        sinrReg       = t2 ./ (t1 - t2 + t3);
        
        achRate = prelogFactor1 * log2(1 + sinrStag) + prelogFactor2 * log2(1 + sinrReg);
end

sumAchRate = sum(achRate); % Compute and output the sum achievable rate
channelNmseOut= mean(mean(channelMse)); % Compute and output the NMSE of the channel estimate
matrixErrorNormOut = mean(mean(matrixErrorNorm.^2)); % Compute and output the MSE of the covariance matrix estimate
