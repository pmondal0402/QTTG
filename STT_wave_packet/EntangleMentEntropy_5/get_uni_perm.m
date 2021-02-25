function [result] = get_uni_perm(n, k)
% Returns unique permutations of k numbers chosen from n.
% Note k <=n 
% n (total no. of elements) is total number of elements in set
% k (subset of elelments) is subset of n

% Note : If the user is interested in all possible combination of numbers
% of a given set, choose k = n

A = npermutek( n, k ) ;
result = unique(A, 'rows') ;

