clc;
clear all;
close all;

%% Problem Definition

nVar=[3 0];             % Number of Decision Variables

VarSize=[1 , nVar(1) ; 1 , nVar(2)];   % Size of Decision Variables Matrix

VarMin=[-5 , 0];         % Lower Bound of Variables
VarMax=[5 , 0];         % Upper Bound of Variables


%% CSO Parameters

MaxIt=100;      % Maximum Number of Iterations

nPop=50;        % Population Size (Swarm Size)
nTest = 0.1 * nPop;

w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=2;           % Global Learning Coefficient
MR = 0.9;         % 1 = all cats in seeking and 0 = all cats in tracing

nGrid=5;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate
beta=2;             % Leader Selection Pressure

reInitPer = 0.2;

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.flag = false;
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.Rank=[];
empty_particle.DominationSet=[];
empty_particle.DominatedCount=[];
empty_particle.CrowdingDistance=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];
empty_particle.index = 0;

particles=repmat(empty_particle,nPop,1);
TestParticles = repmat(empty_particle , nTest , 1); % Test Swarm

%% Initialize Population
for i=1:nPop
    
    % Initialize Position
    particles(i).Position(1:nVar(1))=unifrnd(VarMin(1),VarMax(1),VarSize(1,:));
    particles(i).Position(nVar(1)+1:nVar(2)+1)=unifrnd(VarMin(2),VarMax(2),VarSize(2,:));
    
    % Initialize Velocity
    particles(i).Velocity(1:nVar(1))=zeros(VarSize(1,:));
    particles(i).Velocity(nVar(1)+1:nVar(2)+1)=zeros(VarSize(2,:));
    % Evaluation
    particles(i).Cost=Benchmark(particles(i).Position , 0);
    particles(i).index = i;
    
    % Update Personal Best
    particles(i).Best.Position=particles(i).Position;
    particles(i).Best.Cost=particles(i).Cost;
    
end

%% Initialize Test Population
for i=1:nTest
    % Initialize Position
    TestParticles(i).Position(1:nVar(1))=unifrnd(VarMin(1),VarMax(1),VarSize(1,:));
    TestParticles(i).Position(nVar(1)+1:nVar(2)+1)=unifrnd(VarMin(2),VarMax(2),VarSize(2,:));
    
    % Evaluation
    TestParticles(i).Cost=Benchmark(TestParticles(i).Position , 0);
end
%%

BestCost=zeros(MaxIt,1);

% Non-Dominated Sorting
[particles , F]=NonDominatedSorting(particles);

% Calculate Crowding Distance
particles=CalcCrowdingDistance(particles,F);

% Sort Population
[particles , F]=SortPopulation(particles);
for j = 1 : length(particles)
    particles(j).index = j;
end

% Determine Domination
F1 = particles(F{1});
Grid=CreateGrid(F1,nGrid,alpha);

for i=1:numel(F1)
    F1(i)=FindGridIndex(F1(i),Grid);
end


%% CSO Movement
for it=1:MaxIt
    TestBackupCosts = [TestParticles.Cost];
    
    TestNewCosts = zeros(2 , nTest);
    for i=1:nTest
        TestParticles(i).Cost = Benchmark(TestParticles(i).Position , it);
        TestNewCosts(: , i) = TestParticles(i).Cost;
    end
    tf = isequal(TestBackupCosts , TestNewCosts);
    if(tf == 0)
        OldIndex = [particles.index];
        OldIndex = OldIndex';
        for i=1:nPop
            particles(i).Cost = Benchmark(particles(i).Position , it);
        end
        [particles , F]=NonDominatedSorting(particles);
        particles = CalcCrowdingDistance(particles,F);
        [particles , F]=SortPopulation(particles);
        
        NewIndex = [particles.index];
        NewIndex = NewIndex';
        
        OptIndex = BordaCountAlgorithm([OldIndex , NewIndex]);
        
        NumOfNoChange = floor((1 - reInitPer) * length(OptIndex));
        
        for i = NumOfNoChange : length(OptIndex)
            particles(OptIndex(i)).Position(1:nVar(1))=unifrnd(VarMin(1),VarMax(1),VarSize(1,:));
            particles(OptIndex(i)).Position(nVar(1)+1:nVar(2)+1)=unifrnd(VarMin(2),VarMax(2),VarSize(2,:));
            
            % Evaluation
            particles(OptIndex(i)).Cost=Benchmark(particles((OptIndex(i))).Position , it);
            
        end
        
        [particles , F]=NonDominatedSorting(particles);
        particles = CalcCrowdingDistance(particles,F);
        [particles , F]=SortPopulation(particles);
        
        for j = 1 : length(particles)
            particles(j).index = j;
        end
        
        
        F1 = particles(F{1});
        Grid=CreateGrid(F1,nGrid,alpha);
        
        for i=1:numel(F1)
            F1(i)=FindGridIndex(F1(i),Grid);
        end
    end
    
    indices = randperm(nPop , floor(nPop*MR));
    SeekingPop = [];
    TracingPop = [];
    
    leader=SelectLeader(F1,beta);
    
    for i=1:nPop
        ind = find(indices == i);
        
        if(~isempty(ind))
            %% Seeking Mode
            [SeekingParticles, leader] = Seeking( particles(i) , leader , VelMax , VelMin , VarMax , VarMin , VarSize , nVar , it);
            SeekingPop = [SeekingPop , SeekingParticles];
        else
            %% Tracing Mode
            [ particles(i) , leader ] = Tracing( particles(i) , leader , VelMax , VelMin , VarMax , VarMin , VarSize , nVar , w , c1 , it);
            TracingPop = [TracingPop , particles(i)];
        end
    end
    particles = [SeekingPop , TracingPop];
    
    % Non-Dominated Sorting
    [particles , F]=NonDominatedSorting(particles);
    
    % Calculate Crowding Distance
    particles=CalcCrowdingDistance(particles,F);
    
    % Sort Population
    [particles , F]=SortPopulation(particles);
    
    % Truncate
    particles=particles(1:nPop);
    
    % Non-Dominated Sorting
    [particles , F]=NonDominatedSorting(particles);
    
    % Calculate Crowding Distance
    particles=CalcCrowdingDistance(particles,F);
    
    % Sort Population
    [particles , F]=SortPopulation(particles);
    for j = 1 : length(particles)
        particles(j).index = j;
    end
    
    % Store F1
    F1=particles(F{1});
    
    Grid=CreateGrid(F1,nGrid,alpha);
    for i=1:numel(F1)
        F1(i)=FindGridIndex(F1(i),Grid);
    end
    
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1);
    
    w=w*wdamp;
end



