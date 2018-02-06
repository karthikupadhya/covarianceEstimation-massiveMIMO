% Code for generating figure 4
% Takes the covariance estimation method (options - ref14,
% proposed-staggeredpilot and proposed-regularpilot),
% channelEstimationMethod (options - leastSquares,exactCovarianceMatrix,
% estimatedCovarianceMatrix), and number of coherence blocks used for
% covariance matrix estimation as input. Returns the sum achievable rate, and mse of
% the channel and covariance matrix as output.
%
% ref14 corresponds to E. Björnson, L. Sanguinetti and M. Debbah, "Massive MIMO with imperfect channel covariance information," 
% 50th Asilomar Conf. Signals, Syst. Computers, Pacific Grove, CA, 2016, pp. 974-978.
% To get the same results in the paper. Uncomment the rng function in line
% 43 and 80.

function [achRateOut] = code4(covarianceEstimationMethod,channelEstimationMethod,numCoherenceBlock)

% - Simulation Parameters - %
numUserLocations            = 100;
numTrial                    = 100;
% - System Parameters - %
numBsAntenna                = 50;
numCell                     = 7;
numUser                     = 5; % - Number of users per cell
% covarianceEstimationMethod  = 'proposed-staggeredpilot';
% channelEstimationMethod     = 'estimatedcovariancematrix';
pP                          = 1; % Pilot power
pD                          = 1; % Data power
% numCoherenceBlock           = 170;
Nr                          = numCoherenceBlock / (10 + numCell);
Nq                          = 10 * Nr;

% - Channel Parameters - %
cellRadius                  = 150; % Cell radius in m
userRadius                  = 120; % User radius in m
normAntennaSpacing          = 1/2; % Antenna spacing normalized by wavelength;
ulSlotLength                = 100;
numStationaryBlock          = 25000;

% - Performance Metrics - %
errPower                    = cell(numUserLocations,1);

% - Simulation Begins - %
for kk = 1:numUserLocations
    % rng(kk); % Uncomment to get the results in the paper.
    center  = zeros(1,numCell);
    for ii = 1:numCell
        if ii == 1
            center(ii)      = 0;
        else
            center(ii)      = ceil((ii-1)/6) * 2 * sqrt(3)/2 * cellRadius * exp(1i * (pi/3 * mod(ii-1,6) ) );
        end
        msLocations         = genHexSample(numUser,0.001) * cellRadius + center(ii);
        rxPowerInDb         = 78.7 - 37.6 * log10(abs(msLocations));
        rxPower(ii,:)       = 10.^(rxPowerInDb/10);
        angleMean(ii,:)     = angle(msLocations) * 180 / pi;
    end
    
    noiseVar                = 1;
    angleSpread             = 20 * ones(numCell,numUser); % in Degrees ( sys.numCell x sys.numUser )
    
    % - Generate Covariance Matrix - %
    covarianceMatrix        = cell(numCell,numUser);
    [covarianceMatrix{:}]   = deal(zeros(numBsAntenna));
    sqrtCovarianceMatrix    = covarianceMatrix;
    for jj = 1:numCell
        for mm = 1:numUser
            userAngleMean   = angleMean(jj,mm);
            userAngleSpread = angleSpread(jj,mm);
            
            integralSpacing = 0.01;
            thetaRange      = -180:integralSpacing:360; % Account for wrap around
            pTheta          = 1/userAngleSpread * ((thetaRange >= userAngleMean - userAngleSpread/2) & (thetaRange <= userAngleMean + userAngleSpread/2));
            
            aVector         = exp(1i * 2 * pi * normAntennaSpacing * (0:numBsAntenna-1).' * cos(pi/180 * thetaRange(pTheta>0)));
            covarianceMatrix{jj,mm}         = rxPower(jj,mm) * (aVector * diag(pTheta(pTheta>0)) * aVector') * integralSpacing;
            sqrtCovarianceMatrix{jj,mm}     = covarianceMatrix{jj,mm}^(1/2);
        end
    end
    
    for ii = 1:numTrial
        % rng(ii); % Uncomment to get the results in the paper.
        % - Generate Channel Vectors- %
        channelMatrix               = cell(numCell,numCoherenceBlock);
        [channelMatrix{:}]          = deal(zeros(numBsAntenna,numUser));
        for jj = 1:numCell
            for mm = 1:numUser
                xii                 = complex(normrnd(0,1/sqrt(2),numBsAntenna,numCoherenceBlock),normrnd(0,1/sqrt(2),numBsAntenna,numCoherenceBlock));
                Hii                 = sqrtCovarianceMatrix{jj,mm} * xii;
                for nn = 1:numCoherenceBlock
                    channelMatrix{jj,nn}(:,mm) = Hii(:,nn);
                end
            end
        end
        
        switch lower(channelEstimationMethod)
            case 'leastsquares'
                estIndivCovarianceMatrix = cell(numCell,numUser);
                [estIndivCovarianceMatrix{:}] = deal(speye(numBsAntenna));
                estSumCovarianceMatrix = cell(1,numUser);
                [estSumCovarianceMatrix{:}] = deal(speye(numBsAntenna));
            case 'exactcovariancematrix'
                % - Use exact covariance matrices - %
                switch lower(covarianceEstimationMethod)
                    case 'proposed-staggeredpilot'
                        estIndivCovarianceMatrix = covarianceMatrix;
                        estSumCovarianceMatrix   = cell(1,numUser);                        
                        for mm = 1:numUser
                            Rsum     = estIndivCovarianceMatrix{1,mm} + noiseVar/(pP * numUser) * eye(numBsAntenna);
                            for ll = 2:numCell
                                for uu = 1:numUser
                                    Rsum = Rsum + pD/(pP*numUser) * estIndivCovarianceMatrix{ll,uu};
                                end
                            end
                            estSumCovarianceMatrix{1,mm} = Rsum;
                        end
                    case 'ref14'
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
                end
            case 'estimatedcovariancematrix'
                % - Estimate the covariance matrices - %
                if strcmpi(covarianceEstimationMethod,'ref14')
                    exactCovarianceMatrix       = covarianceMatrix;
                    estIndivCovarianceMatrix    = cell(numCell,numUser);
                    estSumCovarianceMatrix      = cell(1,numUser);
                    
                    %%% - Estimation of Q - %%%
                    for nn = 1:Nq
                        hHat1{nn} = zeros(numBsAntenna,numUser);
                        for ll = 1:numCell
                            hHat1{nn} = hHat1{nn} + channelMatrix{ll,nn};
                        end
                        hHat1{nn} = hHat1{nn} + complex(normrnd(0,sqrt(noiseVar/(2*numUser)),numBsAntenna,numUser),normrnd(0,sqrt(noiseVar/(2*numUser)),numBsAntenna,numUser));
                    end
                    
                    for mm = 1:numUser
                        tempH1 = zeros(numBsAntenna,Nq);
                        for nn = 1:Nq
                            tempH1(:,nn)         = hHat1{nn}(:,mm);
                        end
                        estSumCovarianceMatrix{1,mm} = tempH1 * tempH1' / Nq ; % Desired + interfering
                    end
                    
                    %%% - Estimation of R - %%%
                    jj = 1; % Evaluate R only for the center cell.
                    cnt = 0;
                    for nn = Nq+1:numCoherenceBlock
                        if jj == mod(nn-Nq-1,numCell)+1
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
                    
                    for mm = 1:numUser
                        % - Evaluate exact Q - %
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
                                tempH2(:,cnt)= hHat2{jj,cnt}(:,mm);
                            end
                        end
                        estIndivCovarianceMatrix{jj,mm} = estSumCovarianceMatrix{1,mm} - tempH2 * tempH2'/ cnt ; % Desired
                        muRange = 0:0.05:1;
                        etaRange = 0:0.05:1;
                        for nn = 1:length(muRange);
                            muVal = muRange(nn);
                            for pp = 1:length(etaRange);
                                etaVal           = etaRange(pp);
                                tempQMatrix      = etaVal * estSumCovarianceMatrix{1,mm} + (1 - etaVal) * diag(diag(estSumCovarianceMatrix{1,mm}));
                                tempRMatrix      = muVal * estIndivCovarianceMatrix{jj,mm} + (1 - muVal) * diag(diag(estIndivCovarianceMatrix{jj,mm}));
                                tempW            = tempRMatrix/tempQMatrix;
                                fVal(nn,pp)      = real(trace((eye(numBsAntenna) - tempW - tempW') * exactCovarianceMatrix{jj,mm} ) + trace( tempW *  Qexact * tempW' ));
                            end
                        end
                        [~,idx]                 = min(fVal(:));
                        [muMinIdx,etaMinIdx]    = ind2sub(size(fVal),idx);
                        etaOpt                  = etaRange(etaMinIdx);
                        muOpt                   = muRange(muMinIdx);
                        if jj == 1
                            estSumCovarianceMatrix{1,mm}   = etaOpt * estSumCovarianceMatrix{1,mm} + (1 - etaOpt) * diag(diag(estSumCovarianceMatrix{1,mm}));
                        end
                        estIndivCovarianceMatrix{jj,mm} = muOpt * estIndivCovarianceMatrix{jj,mm} + (1 - muOpt) * diag(diag(estIndivCovarianceMatrix{jj,mm}));
                    end
                else
                    randPhase                   = exp(1i * 2 * pi * rand(numCell,numCoherenceBlock));
                    hHat1 = cell(numCell,numCoherenceBlock);
                    [hHat1{:}] = deal(zeros(numBsAntenna,numUser));
                    hHat2 = hHat1;
                    for ll = 1:numCell
                        txData1    = complex(normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock),normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock));
                        txData2    = complex(normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock),normrnd(0,1/sqrt(2),numUser*(numCell-1),numUser*numCoherenceBlock));
                        noiseVec1  = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock));
                        noiseVec2  = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser*numCoherenceBlock));
                        for nn = 1:numCoherenceBlock
                            txDataIdx    = (nn-1) * numUser + 1 : nn * numUser;
                            hInterf1     = sqrt(pD) / (sqrt(pP*numUser)) * cell2mat(channelMatrix([1:ll-1 ll+1:end],nn).') * txData1(:,txDataIdx) + 1/(sqrt(pP*numUser))* noiseVec1(:,txDataIdx);
                            hHat1{ll,nn} = channelMatrix{ll,nn} + hInterf1;
                            hInterf2     = conj(randPhase(ll,nn)) * sqrt(pD) / (sqrt(pP*numUser)) * cell2mat(channelMatrix([1:ll-1 ll+1:end],nn).') * txData2(:,txDataIdx) + conj(randPhase(ll,nn)) / sqrt(pP*numUser) * noiseVec2(:,txDataIdx);
                            hHat2{ll,nn} = channelMatrix{ll,nn} + hInterf2;
                        end
                    end                    
                    
                    %%% - Estimate R - %%%
                    estIndivCovarianceMatrix       = cell(numCell,numUser);
                    for jj = 1:numCell
                        for mm = 1:numUser
                            tempH1 = zeros(numBsAntenna,numCoherenceBlock);
                            tempH2 = tempH1;
                            for nn = 1:numCoherenceBlock
                                tempH1(:,nn) = hHat1{jj,nn}(:,mm);
                                tempH2(:,nn) = hHat2{jj,nn}(:,mm);
                            end
                            tempCovMat  = tempH1 * tempH2' / numCoherenceBlock;
                            % - Replace with closest psd matrix - %
                            tempCovMat  = (tempCovMat + tempCovMat')/2;
                            [U,D]       = eig(tempCovMat);
                            d           = real(diag(D));
                            estIndivCovarianceMatrix{jj,mm} = U(:,d>0) * sparse(diag(d(d>0))) * U(:,d>0)';
                        end
                    end
                    
                    %%% - Estimate Q - %%%
                    estSumCovarianceMatrix         = cell(1,numUser);
                    for mm = 1:numUser
                        Rsum     = estIndivCovarianceMatrix{1,mm} + noiseVar/(pP * numUser) * eye(numBsAntenna);
                        for ll = 2:numCell
                            for uu = 1:numUser
                                Rsum = Rsum + pD/(pP*numUser) * estIndivCovarianceMatrix{ll,uu};
                            end
                        end
                        estSumCovarianceMatrix{1,mm} = Rsum;
                    end
                end
        end
        
        % - Generate channel vectors for obtaining the channel estimate - %
        channelVec = cell(numCell,1);
        for jj = 1:numCell
            for mm = 1:numUser;
                channelVec{jj,1}(:,mm) = sqrtCovarianceMatrix{jj,mm} * complex(normrnd(0,1/sqrt(2),numBsAntenna,1),normrnd(0,1/sqrt(2),numBsAntenna,1));
            end
        end
        
        switch lower(covarianceEstimationMethod)
            case 'ref14'
                % - Generate channel estimate for regular pilot - %
                hHat = zeros(numBsAntenna,numUser);
                for mm = 1:numUser;
                    for jj = 1:numCell
                        hHat(:,mm) = hHat(:,mm) + channelVec{jj,1}(:,mm);
                    end
                end
                hHat = hHat + complex(normrnd(0,sqrt(noiseVar / (2 * numUser)),numBsAntenna,numUser),normrnd(0,sqrt(noiseVar / (2 * numUser)),numBsAntenna,numUser));
            case 'proposed-staggeredpilot'                
                % - Generate channel estimate for staggered pilot - %
                noiseVec    = complex(normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser),normrnd(0,sqrt(noiseVar/2),numBsAntenna,numUser));
                txData      = complex(normrnd(0,1/sqrt(2),(numCell-1)*numUser,numUser),normrnd(0,1/sqrt(2),(numCell-1)*numUser,numUser));
                hInterf     = sqrt(pD) / (sqrt(pP*numUser)) * cell2mat(channelVec(2:end).') * txData + 1/(sqrt(pP*numUser))* noiseVec;
                hHat        = channelVec{1} + hInterf;            
        end
        
        % - Compute LMMSE channel estimate - %
        channelEstimate     = cell(1,1);
        [channelEstimate{:}]= deal(zeros(numBsAntenna,numUser));
        for mm = 1:numUser
            Rsum            = estSumCovarianceMatrix{1,mm};
            R0              = estIndivCovarianceMatrix{1,mm};
            channelEstimate{1}(:,mm) = R0 * (Rsum \ hHat(:,mm));
        end
        
        % - Compute RZF combining vector - %
        term1 = channelEstimate{1} * channelEstimate{1}';
        term2 = noiseVar * eye(numBsAntenna);
        combiningVec        = ((term1 + term2) \ channelEstimate{1});        
        
        % - Evaluate terms for SINR calculation - %
        for mm = 1:numUser
            interfTerm1(ii,mm) = 0;
            for ll = 1:numCell
                for pp = 1:numUser
                    interfTerm1(ii,mm) = interfTerm1(ii,mm) + abs( combiningVec(:,mm)' * channelVec{ll,1}(:,pp) ).^2;
                end
            end
            interfTerm2(ii,mm) = combiningVec(:,mm)' * channelVec{1,1}(:,mm);
            interfTerm3(ii,mm) = norm(combiningVec(:,mm))^2 * noiseVar;
        end
        
        for mm = 1:numUser
            for ll = 1:numCell
                for pp = 1:numUser
                    interfTermDebug1{mm,ll,pp}(ii,kk) = abs( combiningVec(:,mm)' * channelVec{ll,1}(:,pp) )^2;
                end
            end
            interfTermDebug2{mm}(ii,kk) = combiningVec(:,mm)' * channelVec{1,1}(:,mm);
            interfTermDebug3{mm}(ii,kk) = norm(combiningVec(:,mm))^2;
        end        
    end
    
    % - Compute prelog factor - %
    switch lower(covarianceEstimationMethod)
        case 'ref14'
            prelogFactor = (1 - numUser/ulSlotLength - Nr * numUser * numCell / (numStationaryBlock*ulSlotLength) );
        case 'proposed-staggeredpilot'
            prelogFactor1 = (numCell-1) / (ulSlotLength * numStationaryBlock) * numUser * (numCoherenceBlock  + numStationaryBlock);
            prelogFactor2 = 1 - (numCell * numUser) / ulSlotLength - (numCoherenceBlock * numCell * numUser) / (numStationaryBlock * ulSlotLength);
    end
    
    % - Compute SINR and achievable rate - %
    switch lower(covarianceEstimationMethod)
        case 'ref14'
            % - Evaluate Equation (15) in Reference 14 - %
            interferencePower{kk} = mean(interfTerm1,1) - abs(mean(interfTerm2,1)).^2 + mean(interfTerm3,1);
            sigPower{kk}          = abs(mean(interfTerm2,1)).^2;
            sinr{kk}              = sigPower{kk}./interferencePower{kk};
            achRate(kk,:)          = prelogFactor * log2(1 + sinr{kk});
        case 'proposed-staggeredpilot' 
            % - Evaluate Equation (20) in submitted manuscript - %
            for mm = 1:numUser
                for ll = 1:numCell
                    for pp = 1:numUser
                        interfVal1{mm}(ll,pp) = mean(interfTermDebug1{mm,ll,pp}(:,kk));
                    end
                end
                
                t1(mm) = sum(reshape(interfVal1{mm},1,[]));
                
                t2(mm) = abs( mean(interfTermDebug2{mm}(:,kk)) ).^2;
                
                t3(mm) = noiseVar / pD * mean(interfTermDebug3{mm}(:,kk));
                
                for ll = 1:numCell-1
                    t4(ll,mm) = (pP - pD) / pD * sum(interfVal1{mm}(ll+1,:));
                end
            end
            t4     = mean(t4);
            
            sinrStag{kk}        = t2 ./ (t1 - t2 + t3 + t4);
            sinrReg{kk}         = t2 ./ (t1 - t2 + t3);
            
            achRate(kk,:)       = prelogFactor1 * log2(1 + sinrStag{kk}) + prelogFactor2 * log2(1 + sinrReg{kk});
    end
end

% - Average over all user locations - %
achRateOut = sum(mean(achRate,1));
