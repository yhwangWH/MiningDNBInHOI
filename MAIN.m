% clc; clear all;
% Parameters setting
popSize = 200; 
maxIteration = 100;
pCross = 0.5;
pMuteGene = 0.1;
pMuteRadius = 0.1;
maxRadius = 25;
maxDNB = 50;

dataset = table2array(NormedDEGCaseSamples1);
numOfGenes = size(dataset,1);
numOfSamples = size(dataset,2);
T05 = [1:6];
T1 = [7:12];
T4 = [13:18];
T8 = [19:24];
T12 = [25:30];
numOfTimePoints = 5;
timePoints = cell(numOfTimePoints,1);
timePoints{1} = T05;
timePoints{2} = T1;
timePoints{3} = T4;
timePoints{4} = T8;
timePoints{5} = T12;

%% Initiation
GAPop = cell(popSize,1);
for i = 1: popSize
    randLength = randi([2, maxDNB]); 
    DNB = randperm(numOfGenes ,randLength);
    
    DNBChorm = zeros(1, 2*randLength);
    for j = 1: randLength 
        DNBChorm(1, 2*j-1)= DNB(j);
    end
    for j = 1: randLength 
        DNBChorm (1, 2*j) = randi(randLength);
    end
    GAPop{i} = DNBChorm;
end

globalFitnessVal = cell(maxIteration,3);
m = 1;

while m < maxIteration+1
    tic;
    %% Fitness evaluation
    FitnessVal = cell(popSize,3);
    for i = 1:popSize
        DNBChorm = GAPop{i};
        [maxEmerge,maxEmergeTime, DNBChorm] = Fitness(dataset, DNBChorm,timePoints, numOfGenes, numOfTimePoints);
        FitnessVal{i,1} = maxEmerge;
        FitnessVal{i,2} = maxEmergeTime;
        FitnessVal{i,3} = DNBChorm;
    end
    orderTemp = [];
    for i = 1:popSize
        temp = FitnessVal{i,1};
        orderTemp = [orderTemp temp];
    end
    [max1, maxIndex] = max(orderTemp);
    globalFitnessVal{m,1} = FitnessVal{maxIndex,1};
    globalFitnessVal{m,2} = FitnessVal{maxIndex,2};
    globalFitnessVal{m,3} = FitnessVal{maxIndex,3};
    
    %% binary tournament selection strategy
    % Each time, two chromosomes are randomly selected from the population for competition
    % and the fittest one wins and is selected for variation.
    % We repeat this process until the population size is reached.
    FitnessValSelect = cell(popSize,3);
    for i = 1:popSize
        chorm1Index = randi(popSize);
        chorm2Index = randi(popSize);
        
        if FitnessVal{chorm1Index,1} >= FitnessVal{chorm2Index,1}
            FitnessValSelect{i,1} = FitnessVal{chorm1Index,1};
            FitnessValSelect{i,2} = FitnessVal{chorm1Index,2};
            FitnessValSelect{i,3} = FitnessVal{chorm1Index,3};
        else
            FitnessValSelect{i,1} = FitnessVal{chorm2Index,1};
            FitnessValSelect{i,2} = FitnessVal{chorm2Index,2};
            FitnessValSelect{i,3} = FitnessVal{chorm2Index,3};
        end
    end
    
    %% Crossover
    GAPopCross = cell(popSize,1);
    for i = 1:2:popSize
        p = rand();
        if p<=pCross 
            lengthOfChorm1 = size(FitnessValSelect{i,3} ,2)/2;
            lengthOfChorm2 = size(FitnessValSelect{i+1,3} ,2)/2;
            seq = randi( min(lengthOfChorm1, lengthOfChorm2) );
            
            if seq == 1
                newChorm1 = [FitnessValSelect{i+1,3}(seq:seq+1) FitnessValSelect{i,3}(2*seq+1:end)];
                newChorm2 = [FitnessValSelect{i,3}(seq:seq+1) FitnessValSelect{i+1,3}(2*seq+1:end)];
            else
                newChorm1 = [FitnessValSelect{i,3}(1:2*seq-2) FitnessValSelect{i+1,3}(2*seq-1:2*seq) FitnessValSelect{i,3}(2*seq+1:end) ];
                newChorm2 = [FitnessValSelect{i+1,3}(1:2*seq-2) FitnessValSelect{i,3}(2*seq-1:2*seq) FitnessValSelect{i+1,3}(2*seq+1:end) ];
            end
            
            GAPopCross{i}= newChorm1;
            GAPopCross{i+1} = newChorm2;
        else
            GAPopCross{i}= FitnessValSelect{i,3};
            GAPopCross{i+1} = FitnessValSelect{i+1,3};
        end
        
    end
    
    %% Mutation
    GAPopMut = GAPopCross;
    for j = 1:popSize
        for i = 1: 2: size(GAPopMut{j},2)
            p1 = rand();
            if p1 <= pMuteGene
                lambda1 = randi(numOfGenes);
            else
                lambda1 = GAPopCross{j}(i);
            end
            GAPopMut{j}(i) = lambda1;
        end
    end
    
    GAPopMutRadius = GAPopMut;
    for j=1:popSize
        for i = 2: 2: size(GAPopMut{j},2)
            p2 = rand();
            if p2 <= pMuteRadius
                lambda2 = randi(maxRadius);
            else
                lambda2 = GAPopMut{j}(i);
            end
            GAPopMutRadius{j}(i) = lambda2;
        end
    end
    
    GAPop = GAPopMutRadius;
    t = toc;
    fprintf('%d iteration is done for %f sec\n', m, t)
    m = m+1;
end


