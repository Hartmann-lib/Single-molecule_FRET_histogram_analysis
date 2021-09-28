function randArray=randi_hist(values,occurence,N1,N2)

    cumDistr=cumsum(occurence./sum(occurence));

    randNum=rand(N1,N2);

    p=@(randNum) find(randNum<cumDistr,1,'first');

    index=arrayfun(p,randNum);
    randArray=values(index);
end

