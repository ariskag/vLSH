function [rec] = compute_recall(index_base, index)
% Function that computes the recall of knn search
% input:
%   index_base: The knn index of the ground truth stored in a [k x N]
%   matrix
%   index: The knn index of the approximate result strored in a [k x N]
%   matrix
%
% output:
%   rec: The recall value 
%
% author: Nikos Sismanis
% date: October 2013

% The size of the input matrix
k = size(index_base, 1);
n = size(index_base, 2);

found = 0;

% process each line separately
for i=1:k    
   f =  any(bsxfun(@eq,index_base, index(i,:)), 1); % fund the true matches
   found = found + sum(double(f)); % count how many neighbors were successfully retrieved
end

rec = found / prod([n k]); % normalize over the total number of neighbors retrieved

end
