function FootruleRankingInteger = OptimalSpearmanFootRuleRanking(GroupsRanking)

[C , ~] = size(GroupsRanking);
W = zeros(C);

for i = 1 : C
    for j = 1 : C
        W(i , j) = FindSpearmanCoefficient(i , j , GroupsRanking);
    end
end

FreeStatus = zeros(C);

SumOfCharacterStatuses = sum(sum(FreeStatus , 2));

ProposedPosition = zeros(C);

while SumOfCharacterStatuses ~= C
    CharactersStatuses = sum(FreeStatus , 2);
    ind = find(CharactersStatuses == 0);
    ind = ind(1);
    Positions = W(ind , :);
    
    [ ~ , RankedPosition ] = sort( Positions , 'ascend');
    for i = 1 : length(RankedPosition)
        PositionIsFree = sum(FreeStatus(: , RankedPosition(i)));
        if PositionIsFree == 0
            FreeStatus(ind , RankedPosition(i) ) = 1;
            ProposedPosition(ind , RankedPosition(i)) = 1;
            break;
        else
            AnotherCharacter = find(FreeStatus(: , RankedPosition(i)));
            if W(ind , RankedPosition(i) ) <= W( AnotherCharacter , RankedPosition(i) )
                if ProposedPosition(ind , RankedPosition(i) ) == 0
                    ProposedPosition(ind , RankedPosition(i)) = 1;
                    FreeStatus(ind , RankedPosition(i) ) = 1;
                    FreeStatus(AnotherCharacter , RankedPosition(i) ) = 0;
                    break;
                end
            end
        end
        
    end
    
    SumOfCharacterStatuses = sum(sum(FreeStatus , 2));
end

FootruleRankingInteger = zeros(C , 1);
for i = 1 : C
    ind = find(FreeStatus(i , :));
    FootruleRankingInteger(i) = ind;
end

end


function SpearmanCoeff = FindSpearmanCoefficient(harf , position , GroupsRanking)

SpearmanCoeff = 0;
[~ , CountsOfRankings] = size(GroupsRanking);

for i = 1 : CountsOfRankings
    RankingI = GroupsRanking(: , i);
    
    ind = find(RankingI == harf);
    Distance = abs(position - ind);
    
    SpearmanCoeff = SpearmanCoeff + Distance;
end
end