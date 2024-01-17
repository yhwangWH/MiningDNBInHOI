function [maxEmerge,maxEmergeTime, DNBChorm] = FitnessForDyadic(dataset, DNBChorm,timePoints, numOfGenes, numOfTimePoints)


CI = zeros(numOfTimePoints,1);

for t = 1:numOfTimePoints
    
    datasetOfTime = dataset(:,timePoints{t});
    DNBSize = size(DNBChorm,2);

    sumSD = 0;
    for i = 1:DNBSize
        sumSD = sumSD + std(datasetOfTime(DNBChorm(i),:));
    end
    SD = sumSD/(DNBSize);
    
    % GroupIn
    rho= corrcoef(datasetOfTime'); 

    allGenesCombInDNB = nchoosek(DNBChorm,2);
    GroupIn = 0;
    for i = 1:size(allGenesCombInDNB,1)
        GroupIn =  GroupIn + rho( allGenesCombInDNB(i,1), allGenesCombInDNB(i,2) );
    end
    GroupIn = GroupIn/(DNBSize*(DNBSize-1));
    
    % GroupOut
    allGenes = [1:numOfGenes];
    genesNotInDNB = setdiff(allGenes, DNBChorm);
    GroupOut = 0; 
    for i = 1:DNBSize
        for j = 1:size(genesNotInDNB,2)
            GroupOut = GroupOut + rho( DNBChorm(i) , genesNotInDNB(j) );
        end
    end
    GroupOut = GroupOut/(DNBSize*(numOfGenes - DNBSize));
    
    % CI
    CI(t) = (SD*GroupIn)/GroupOut;
end

detaCI = zeros(numOfTimePoints-1,1);
for i = 1: numOfTimePoints-1
    detaCI(i) = CI(i+1)-CI(i);
end


[maxEmerge, maxEmergeTime] = sort(detaCI, 'descend');
maxEmerge = maxEmerge(1);
maxEmergeTime = maxEmergeTime(1);



end

