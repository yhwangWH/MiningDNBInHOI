function [maxEmerge,maxEmergeTime, DNBChorm] = Fitness(dataset, DNBChorm,timePoints, numOfGenes, numOfTimePoints)

CI = zeros(numOfTimePoints,1);
for j = 1:numOfTimePoints
    datasetOfTime = dataset(:,timePoints{j});
    DNBSize = size(DNBChorm,2)/2;
    
    sumSD = 0;
    for i = 1:2:2*DNBSize
        sumSD = sumSD + std(datasetOfTime(DNBChorm(i),:));
    end
    SD = sumSD/(DNBSize);
    
    rho= corrcoef(datasetOfTime'); 
    hyperEdge = cell(DNBSize,2); 
    
    for i = 1:2:2*DNBSize
        absRho = abs(rho(:,DNBChorm(i)));
        [maxPCC, maxPCCIndex] = sort(absRho,'descend');
        hyperEdgeNodes = maxPCCIndex(2: 2+DNBChorm(i+1)-1); 
        hyperEdge{ (i+1)/2 ,1} = hyperEdgeNodes;
        hyperEdgeNodesPCC = sum(maxPCC(2: 2+DNBChorm(1+1)-1)); 
        hyperEdge{ (i+1)/2 ,2} = hyperEdgeNodesPCC;
    end
    
	% GroupIn
    GroupIn = 0;
    for i = 1:DNBSize
        GroupIn = GroupIn + hyperEdge{i,2};
    end
    numOfPccIn = 0;
    for i = 2:2:2*DNBSize+1
        numOfPccIn = numOfPccIn + DNBChorm(i);
    end
    GroupIn = GroupIn/numOfPccIn; 
    
	% GroupOut
    hyperEdgeOut = cell(DNBSize,1);
    for i = 1:2:2*DNBSize
        absRho = abs(rho(:,DNBChorm(i)));
        OutHyperEdge = setdiff(1: size(absRho,1), [hyperEdge{(i+1)/2,1}' DNBChorm(i)]); % 删除他自己和超边里面的东西
        absRhoOutHyperEdge = absRho(OutHyperEdge,:);
        hyperEdgeOut{(i+1)/2,1} = sum(absRhoOutHyperEdge);
    end
    GroupOut = 0;
    for i = 1:DNBSize
        GroupOut = GroupOut + hyperEdgeOut{i,1};
    end
    
    numOfPccOut = 0;
    for i = 2:2:2*DNBSize+1
        numOfPccOut = numOfPccOut + (numOfGenes - DNBChorm(i) -1);
    end
    
    GroupOut = GroupOut/numOfPccOut;
    % CI
    CI(j) = (SD*GroupIn)/GroupOut;
end


% detaCI
detaCI = zeros(numOfTimePoints-1,1);
for i = 1: numOfTimePoints-1
    detaCI(i) = CI(i+1)-CI(i);
end

[maxEmerge, maxEmergeTime] = sort(detaCI, 'descend');
maxEmerge = maxEmerge(1);
maxEmergeTime = maxEmergeTime(1);

end

