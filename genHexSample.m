% Generates random user locations within a hexagonal cell. Users have a minimum distance of 10*d0 from the BS. 
% Input - numUser, d0.

function genSample = genHexSample(numUser,d0)

numSample = 0;
while numSample < numUser
    testSample  = 2 * ( rand(1,numUser - numSample) + 1i * rand(1,numUser - numSample) - 0.5 - 0.5i );
    idx1        = imag(testSample) < sqrt(3)/2;
    idx2        = imag(testSample) > -sqrt(3)/2;
    idx3        = imag(testSample) < ( -sqrt(3) * real(testSample) + sqrt(3) );
    idx4        = imag(testSample) > ( -sqrt(3) * real(testSample) - sqrt(3) );
    idx5        = imag(testSample) > ( sqrt(3) * real(testSample)  - sqrt(3) );
    idx6        = imag(testSample) < ( sqrt(3) * real(testSample)  + sqrt(3) );
    idx7        = abs(testSample)  > 10 * d0; %Minimum distance
    testSample  = testSample(idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx7);
    genSample(numSample+1:numSample+length(testSample))  = testSample;
    numSample   = length(genSample);
end

genSample       = genSample(1:numUser);

end
