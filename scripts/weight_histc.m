function counts=weight_histc(x,weights,edges)

    x=x(:);
    weights=weights(:);
    edges=edges(:);

    % define histogram array
    counts=zeros(length(edges),1);
    
    % extrapolate to the left for full bins
    bins=[2*edges(1)+edges(2);edges];
    
    % count weights
    for iter=1:length(edges)
        
        counts(iter)=sum(weights(bins(iter)<x&x<=bins(iter+1)));
    end
end