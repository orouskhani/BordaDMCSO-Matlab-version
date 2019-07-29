function [ particle_copies, GlobalBest ] = Seeking( particle , GlobalBest , VelMax , VelMin , VarMax , VarMin , VarSize , nVar , it)
SRD = 0.3;
SMP = 4;
numOfVars = sum(nVar);
CDC = floor(0.8 * numOfVars);

%Fitness_Values = zeros(2 , SMP + 1);
particle_copies = repmat(particle , 1, SMP+1);
%Fitness_Values(1) = particle_copies(1).Cost;

for i = 2 : length(particle_copies)
    permIndices = randperm(numOfVars , CDC);
    
    for index = 1 : length(permIndices)
        booleanFlag = mod(randi(2),2);
        if(booleanFlag == 0)
            particle_copies(i).Position(permIndices(index)) = particle_copies(i).Position(permIndices(index)) * (1 + SRD) ;
        else
            particle_copies(i).Position(permIndices(index)) = particle_copies(i).Position(permIndices(index)) * (1 - SRD) ;
        end
    end
    % Apply Position Limits
    particle_copies(i).Position(1:nVar(1)) = max(particle_copies(i).Position(1:nVar(1)),VarMin(1));
    particle_copies(i).Position(1:nVar(1)) = min(particle_copies(i).Position(1:nVar(1)),VarMax(1));
    
    particle_copies(i).Position(nVar(1)+1:nVar(2)+1) = max(particle_copies(i).Position(nVar(1)+1:nVar(2)+1),VarMin(2));
    particle_copies(i).Position(nVar(1)+1:nVar(2)+1) = min(particle_copies(i).Position(nVar(1)+1:nVar(2)+1),VarMax(2));
    
    particle_copies(i).Cost = Benchmark(particle_copies(i).Position , it);
    
    
    if Dominates(particle_copies(i),GlobalBest)
        GlobalBest.Position=particle_copies(i).Position;
        GlobalBest.Cost=particle_copies(i).Cost;
        
    elseif Dominates(GlobalBest,particle_copies(i))
        % Do Nothing
    else
        if rand<0.5
            GlobalBest.Position = particle_copies(i).Position;
            GlobalBest.Cost = particle_copies(i).Cost;
        end
    end
    
    %Fitness_Values(i) = particle_copies(i).Cost;
end


% maxF = max(Fitness_Values);
% minF = min(Fitness_Values);
%
% particle_probability = abs(Fitness_Values - maxF) / abs(maxF - minF);
% randIndex = rand();
% IsGreater = particle_probability >= randIndex;
% GreaterIndex = find(IsGreater == 1);
%
% particle = particle_copies(GreaterIndex(1));
