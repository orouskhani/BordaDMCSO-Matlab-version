function BordaCountRankingInteger = BordaCountAlgorithm(GroupsRanking)
[Row Col] = size(GroupsRanking);
Scores = zeros(Row , 1);
for j = 1 : Col
    SizeOfGroupJ = length(GroupsRanking(: , j));
    for i = 1 : Row
        Scores(GroupsRanking(i,j)) = Scores(GroupsRanking(i,j)) + ( Row - i ) * SizeOfGroupJ;
    end
end

[~ , BordaCountRankingInteger] = sort(Scores , 'descend');    
end