function [ ids cand_sizes binsFound_ratios times ] = lshSearch( queries, X, lshStruct, K, T )
% function [ ids cand_sizes binsFound_ratios times ] = lshSearch( queries, X, lshStruct, K, T )
%
% Search for queries' K nearest neighbors from dataset X using LSH
%
% INPUT:
%	queries:   [D x Q] matrix containing the queries
%	X:         [D x N]	 matrix containing dataset
%	lshStruct: [1 x L] array of structs produced by lshConstruct
%	K:         # of nearest neighbors to search
%	T:         {optional} # of additional probing buckets for
%	             multiprobe searching (default=0 for simple LSH)
%
% OUTPUT:
%	ids:       [K x Q] columns of K ids of nearest neighbors for each query
%	cand_size: [1 x Q] total number of distinct candidates checked for each query
%	binsFound_ratios: ratio of existing probing bins found for queries
%	times: struct containing the following fields:
%		retrieval: [1 x L] time to retrieve K candidates
%		                   for all queries from each table
%		knn: time to compute final queries' K nearest neighbors
%		             among their L*K candidates retrieved earlier
%		cand_size: time to compute queries' candidate sizes
%		total: total time elapsed to search for nearest neighbors
%
% Projection of queries is done this way:
%	Y = floor( X'*A - b) / w);
%	Afterwards, the minimum values at each dimension(stored at lshStruct.low) are subtracted
%
% The function calls lshFindQueryBin, which makes use of 2nd level hashing
% by using Gregory Shakhnarovich's lshhash.m function, included in his LSH
% toolkit (http://ttic.uchicago.edu/~gregory/download.html)
% lshhash.m produces a scalar value for each M-digit code.
% lshStruct.bhash(key) contains indices to buckets for which lshhash
% produces hash value <key>.

%{
	Copyright (C) 2014. All Rights Reserved.
	Aris Kagias <ariskagias@gmail.com>
	Ioannis Markopoulos <markioan@gmail.com>
	Nikolaos Sismanis <nsismani@auth.gr>
	
	This file is part of vLSH.
	vLSH is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with cuLSH.  If not, see <http://www.gnu.org/licenses/>.
%}

% If T is not given, perform simple LSH
if(~exist('T', 'var'))
	T = 0;
end

tstart = tic;

[D Q] = size(queries);
L = length(lshStruct);
	% D: # of dimensions at original space
	% Q: # of queries
	% L: # of tables (LSH projections)

ids_all = zeros(L*K, Q);
candidates_all = cell(L, Q);
binsFound = zeros(1, Q);
	% ids_all: [L*K x Q] ids of queries' nearest neighbors retrieved from all the tables
	% candidates_all: [L x Q] cell matrix, cell (l, q) contains candidates ids for query q at table l
	% binsFound: [1 x Q] # of total buckets matching query's projection codes

queries_cellArray = mat2cell(queries, D, ones(1, Q));
	% queries_cellArray: [1 x Q] cell vector, each cell contains [D x 1] query data

% Find for each table the queries' corresponding buckets and search for their nearest neighbors in these buckets
for table=1:L
	tic;
%{	
	if(T==0)
		queryCodes = bsxfun(@minus, floor( bsxfun(@minus, queries'*lshStruct(table).A, lshStruct(table).b) / lshStruct(table).w), lshStruct(table).low);
			% [Q x M]
		keys = lshhash(queryCodes)';
			% [1 x Q]
		possiblyMatchingBins = lshStruct(table).bhash(keys);
			% [1 x Q] cell array
		fh_getMatchingIndexFromPossibleBins = @(queryCode, possibleBins) possibleBins( all( bsxfun(@eq, lshStruct(table).codes(possibleBins, :), queryCode), 2 ) );
		binIndices = cellfun( fh_getMatchingIndexFromPossibleBins, mat2cell(queryCodes, ones(1, Q), size(queryCodes, 2))', possiblyMatchingBins );
			% [1 x Q]
		candidates_all(table, :) = lshStruct(table).pointId(binIndices);
	else
%}
		fprintf('Retrieving candidates from table %d\n', table);
		fh_findBinIndex = @(query) lshFindQueryBins(query, lshStruct(table), T);
		binIndices = cellfun(fh_findBinIndex, queries_cellArray, 'UniformOutput', false);
			% binIndices: [1 x Q] cell vector, each cell contains bucket ids matching queries' codes
		binsFound = binsFound + cellfun(@length, binIndices);

		% Retrieve K nearest neighbors for queries from this table
		[ids_all( ((table-1)*K+1):(table*K), :), candidates_all(table, :)] = lshSearchTable(queries, X, K, lshStruct(table).pointId, binIndices);

		times.retrieval(table) = toc;
		fprintf('\t-----> %f seconds\n', toc);
%	end
end

% For each query there are (T+1) bucket codes produced for each table (original code + T probing codes)
binsFound_ratios = binsFound / (L * (T+1));	% [1 x Q] ratios of existing buckets found for each query

% Retrieve for each query its K nearest neighbors among its L*K candidates
tic;

tic;
fprintf('Removing duplicate candidates...\n');
% Find absolute nearest neighbors of queries by checking the neighbors retrieved from all the tables
% Determine unique nearest neighbor candidates for each query out of the ids retrieved from all tables
queryCandidateIds = mat2cell(ids_all, L*K, ones(1, Q));	% [1 x Q] cell vector (might contain zero values)
clear ids_all
	% If for a query KK < L*K candidates were retrieved, its last L*K - KK indices are equal to zero.
	% Remove indices to zero, since they will cause an error at knn function called later
queryCandidateIds = cellfun( @nonzeros, queryCandidateIds, 'UniformOutput', false );
	% Remove duplicate entries
queryCandidateIds = cellfun( @unique, queryCandidateIds, 'UniformOutput', false );
	% queryCandidateIds: [1 x Q] cell vector, each cell contains the union
	% of query's nearest neighbor candidates retrieved from all the tables\
metrics.times.duplicateRemoval = toc;
fprintf('\t-----> %f seconds\n', toc);

queryCandidateSizes = cellfun(@length, queryCandidateIds);
%queryCandidateSizes_lessThanK = arrayfun(@lt, queryCandidateSizes, K*ones(1, Q));	% [1 x Q]
KK = K*ones(1, Q);	% [1 x Q]
queryCandidateSizes_lessThanK = queryCandidateSizes < KK;
KK(queryCandidateSizes_lessThanK) = queryCandidateSizes(queryCandidateSizes_lessThanK);	% [1 x Q]
queryDatasets_cellArray = cellfun( @(queryCandidateIds) X(:, queryCandidateIds), queryCandidateIds, 'UniformOutput', false );
	% queryCandidateSizes: [1 x Q] # of candidates for each query
	% queryCandidateSizes_lessThanK: [1 x Q] indices to queries with less than K candidates (an error will be produced at knn if it tries to search K nearest neighbors for these queries)
	% KK: [1 x Q] # of nearest neighbors to be retrieved by knn (equals to K, unless candidate size is smaller than K)
	% queryDatasets_cellArray: [1 x Q] cell vector, each cell contains a [D x QCS] matrix containing the columns of X corresponding to query's candidates, where QCS is query's candidate size

tic;
fprintf('Performing knn search...\n');
%# queryIds = cellfun( @(queries, data, K) knnsearch(data', queries', 'K', K)', queries_cellArray, queryDatasets_cellArray, num2cell(KK), 'UniformOutput', false );
[~, queryNearestNeighborIds] = cellfun( @knn, queries_cellArray, queryDatasets_cellArray, num2cell(KK), 'UniformOutput', false );	% [1 x Q] cell vector (indices referring to queryIdsNNuniq)
	% queryNearestNeighborIds: [1 x Q] cell vector, each cell is a [KK x 1] column containing query's nearest neighbors retrieved from candidates
	% Indices refer to queries' candidate sets. Transform them to refer to initial dataset X
times.knn = toc;
fprintf('\t-----> %f seconds\n', toc);

queryNearestNeighborIds = cellfun( @(vector, index) vector(index), queryCandidateIds, queryNearestNeighborIds, 'UniformOutput', false );
	% In case KK < K neighbors were retrieved for a query, pad K - KK zeros to its results
zerosNeeded = K*ones(1, Q) - KK;
queryNearestNeighborIds = cellfun(@(x, y) cat(1, x, zeros(y, 1)), queryNearestNeighborIds, num2cell(zerosNeeded), 'UniformOutput', false);	% [1 x Q] cell vector, all cells have the same size of [K x 1]
	% Each cell of queryNearestNeighborIds now is a [K x 1] column
ids = cell2mat(queryNearestNeighborIds);	% [K x Q]
	% ids: [K x Q] ids of queries' nearest neighbors retrieved from all tables

tic;
fprintf('Calculating candidate sizes...\n');
cand_sizes = arrayfun( @(query) length(unique([candidates_all{:, query}])), 1:Q );
	% cand_sizes: [1 x Q] # of unique candidates for each query
times.cand_size = toc;
fprintf('\t-----> %f seconds\n', toc);

times.total = toc(tstart);

end
