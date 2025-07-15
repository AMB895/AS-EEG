function [max_perm_tfce] = calc_perm_tfce_2d(tFvals,nPerm)
size1 = size(tFvals,1);
size2 = size(tFvals,2);
for n=1:nPerm
    % getting total size of t/F matrix
    totalSize = size1*size2;
    
    % reshape t/F matrix
    reshape_tF = reshape(tFvals,[1,totalSize]);
    
    % permute t/F matrix
    tF = reshape_tF(randperm(totalSize)); 
    
    % reshape back to original size
    tF = reshape(tF,[size1,size2]);
    
    % display permutation number
%     disp(n)
    
    % run limo_tfce on permuted t/F values
    perm_tfce_scores = limo_tfce(2,tF,[],0);
    
    % save max permutation tfce score per iteration
    max_perm_tfce(n) = max(perm_tfce_scores,[],[1 2]);
end

end