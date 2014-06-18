function [D, Idx] = knn(X, Y, k)
% function [D, Idx] = knn(X, Y, k)
%   This function is a matlab verion of the knn algorithm to be used as
%   reference
% inputs: 
%   X: A [d x Q] martix containing the query points
%   Y: A [d x N] matrix containing the corpus
%   k: The number of required neighbors 
%
% output: 
%   D: A [k x Q] matrix containing the distances of the knn
%   Idx: A [k x Q] matrix containing the indices of the knn
%
%
% author: Nikos Sismanis
% date: October 2013

if(size(Y, 2)==1)
	Idx = ones(1, size(X, 2));
	D = cellfun(@(x) norm(Y - x), mat2cell(X, size(X, 1), ones(1, size(X, 2))));
	return;
end

 [D, Idx] = sort(bsxfun(@plus,(-2)*(Y'*X),dot(Y,Y,1)'));
 % if there are no enough neighbours take the minimun number
 k=min(k,size(Idx,1));       
 Idx = Idx(1:k,:);
 D = D(1:k,:);
 D = bsxfun(@plus,D,dot(X,X,1));


