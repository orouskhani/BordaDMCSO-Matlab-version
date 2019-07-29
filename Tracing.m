function [ particle , GlobalBest ] = Tracing( particle , GlobalBest , VelMax , VelMin , VarMax , VarMin , VarSize , nVar , w , c1 , it)
% Update Velocity
particle.Velocity(1:nVar(1)) = w*particle.Velocity(1:nVar(1)) +c1*rand(VarSize(1,:)).*(GlobalBest.Position(1:nVar(1))-particle.Position(1:nVar(1)));
particle.Velocity(nVar(1)+1 : nVar(2) + 1) = w*particle.Velocity(nVar(1)+1:nVar(2)+1) +c1*rand(VarSize(2,:)).*(GlobalBest.Position(nVar(1)+1:nVar(2)+1)-particle.Position(nVar(1)+1:nVar(2)+1));

% Apply Velocity Limits
particle.Velocity(1:nVar(1)) = max(particle.Velocity(1:nVar(1)),VelMin(1));
particle.Velocity(1:nVar(1)) = min(particle.Velocity(1:nVar(1)),VelMax(1));

particle.Velocity(nVar(1)+1:nVar(2)+1) = max(particle.Velocity(nVar(1)+1:nVar(2)+1),VelMin(2));
particle.Velocity(nVar(1)+1:nVar(2)+1) = min(particle.Velocity(nVar(1)+1:nVar(2)+1),VelMax(2));

% Update Position
particle.Position = particle.Position + particle.Velocity;

% Velocity Mirror Effect
IsOutside1=(particle.Position(1:nVar(1))<VarMin(1) | particle.Position(1:nVar(1))>VarMax(1));
IsOutside2=(particle.Position(nVar(1)+1:nVar(2)+1)<VarMin(2) | particle.Position(nVar(1)+1:nVar(2)+1)>VarMax(2));
IsOutside = [IsOutside1 , IsOutside2];

particle.Velocity(IsOutside)=-particle.Velocity(IsOutside);

% Apply Position Limits
particle.Position(1:nVar(1)) = max(particle.Position(1:nVar(1)),VarMin(1));
particle.Position(1:nVar(1)) = min(particle.Position(1:nVar(1)),VarMax(1));

particle.Position(nVar(1)+1:nVar(2)+1) = max(particle.Position(nVar(1)+1:nVar(2)+1),VarMin(2));
particle.Position(nVar(1)+1:nVar(2)+1) = min(particle.Position(nVar(1)+1:nVar(2)+1),VarMax(2));

% Evaluation
particle.Cost = Benchmark(particle.Position , it);

if Dominates(particle,GlobalBest)
    GlobalBest.Position=particle.Position;
    GlobalBest.Cost=particle.Cost;
    
elseif Dominates(GlobalBest,particle)
    % Do Nothing
else
    if rand<0.5
        GlobalBest.Position = particle.Position;
        GlobalBest.Cost = particle.Cost;
    end
end
end

