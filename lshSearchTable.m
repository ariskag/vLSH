function [ ids candidates ] = lshSearchTable( queries, X, K, bins, binIndices )
% function [ ids ] = lshSearchTable( queries, X, K, bins, binIndices )
%
% Search for queries' nearest neighbors among specific buckets
%
% INPUT:
%	queries:    [D x Q] query data
%	X:          [D x N] dataset
%	K:          # of nearest neighbors to retrieve
%	bins:       [1 x B] cell array, cell (b) is a row containing
%	              point indices belonging to bucket <b>
%	binIndices: [1 x Q] cell array, cell (q) is a row containing
%	             indices to the buckets to which query <q> belongs
%
% OUTPUT:
%	ids:        [K x Q] indices to nearest neighbors of queries

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

[D Q] = size(queries);
	% D: # of dimensions
	% Q: # of queries

searchBins = unique( cell2mat(binIndices) );
SB = length(searchBins);
	% searchBins: [1 x SB] unique search bucket ids contained in binIndices
	% SB: # of unique bucket ids

queriesBelonging2SearchBins = cellfun(@(queryBinIndices) ismember(searchBins, queryBinIndices), binIndices, 'UniformOutput', false);
searchBinsQueriesMatrix = cell2mat(queriesBelonging2SearchBins');
searchBinsQueryIndices = cellfun(@find, mat2cell(searchBinsQueriesMatrix, Q, ones(1, SB)), 'UniformOutput', false);
queriesSearchBinsIndices = cellfun(@find, mat2cell(searchBinsQueriesMatrix, ones(1, Q), SB), 'UniformOutput', false)';
	% queriesBelonging2SearchBins: [1 x Q] cell vector, each cell is [1 x SB] vector, which has 1 in position <i> if query belongs to search bin <i>
	% searchBinsQueriesMatrix: [Q x SB], if (i, j)==1 then bin <j> is a probing bin of query <i>
	% searchBinsQueryIndices: [1 x SB] cell vector, each cell contains [QQ x 1] query indices belonging to corresponding bucket
	% queriesSearchBinsIndices: [1 x Q] cell vector, each cell contains search bin ids to which query belongs

candidates = cellfun( @(binIndices) [bins{binIndices}], binIndices, 'UniformOutput', false );
	% candidates: [1 x Q] cell vector, containing the candidate ids for all queries

binSizes = cellfun(@length, bins);
searchingBinSizes = binSizes(searchBins);
KK = K * ones(1, SB);	% [1 x SB]
searchingBinSizes_lessThanK = searchingBinSizes < KK;
KK(searchingBinSizes_lessThanK) = searchingBinSizes(searchingBinSizes_lessThanK);
	% binSizes: [1 x B] # of points belonging to each bucket
	% searchingBinSizes: [1 x SB] # of points belonging to each search bucket
	% KK: [1 x SB] # of nearest neighbors to be retrieved by knn (equals to K, unless search bucket size is smaller than K)

queryGroups_cellArray = cellfun( @(queryIndices) queries(:, queryIndices), searchBinsQueryIndices, 'UniformOutput', false );
bucketDataset_cellArray = cellfun( @(bucketPoints) X(:, bucketPoints), bins(searchBins), 'UniformOutput', false);
[~, queriesKnnIds_searchBins] = cellfun(@knn, queryGroups_cellArray, bucketDataset_cellArray, num2cell(KK), 'UniformOutput', false);
	% queryGroups_cellArray: [1 x SB] cell vector, each cell contains [D x QQ] query data
	% bucketDataset_cellArray: [1 x SB] cell vector, each cell contains [D x NN] data of dataset X
	% queriesKnnIds_searchBins: [1 x SB] cell vector, each cell contains [KK x QQ] indices, where QQ is # of points belonging to bucket
	% Indices refer to the points belonging to search bucket

% Transform queriesKnnIds_searchBins to indices in the initial dataset X
%	If ids is [K x QQ] containing ids of nearest neighbors of QQ queries
%		belonging to the same bucket B, bins{B}(ids) is also [K x QQ],
%		whether bins{B} is a row or a column.
%	However, if only one query point belongs to a bucket, the K-length vector
%		returned by bins{B}(ids) (ids [K x 1]) is a vector of the same shape
%		with bins{B}.
%	Therefore, transform row vectors bins{B} to column vectors beforehand.
binsRowIndices = cellfun(@isrow, bins(searchBins));
bins(searchBins(binsRowIndices)) = cellfun(@transpose, bins(searchBins(binsRowIndices)), 'UniformOutput', false);
queriesKnnIds_bins = cellfun( @(vector, index) vector(index), bins(searchBins), queriesKnnIds_searchBins, 'UniformOutput', false );
	% queriesKnnIds_bins: [1 x SB] cell vector, each cell contains [KK x QQ] indices
	% Indices now refer to initial dataset

% In case KK < K neighbors were retrieved for a search bucket, pad K - KK zeros to its results
zerosNeeded = K*ones(1, SB) - KK;
searchingBinsQueriesContained = cellfun(@length, searchBinsQueryIndices);
queriesKnnIds_bins = cellfun(@(x, dimension1, dimension2) cat(1, x, zeros(dimension1, dimension2)), queriesKnnIds_bins, num2cell(zerosNeeded), num2cell(searchingBinsQueriesContained), 'UniformOutput', false);
	% queriesKnnIds_bins: [1 x SB] cell vector, each cell contains [K x QQ] indices

queriesKnnIds_bins = cell2mat(queriesKnnIds_bins);
correspondingQueries = cell2mat( searchBinsQueryIndices' )';
queriesIndicesToKnnIds = arrayfun( @(query) correspondingQueries==query, 1:Q, 'UniformOutput', false );
	% queriesKnnIds_bins: [K x QQQ]
	% correspondingQueries: [1 x QQQ], queriesKnnIds_bins(:, q) corresponds to query correspondingQueries(q)
	% queriesIndicesToKnnIds: [1 x Q] cell vector, cell (q) contains indices to columns of queriesKnnIds_bins belonging to query <q>

queriesKnnIds = cellfun( @(queryIndices) queriesKnnIds_bins(:, queryIndices), queriesIndicesToKnnIds, 'UniformOutput', false );
queriesKnnIds = cellfun( @(x) x(:), queriesKnnIds, 'UniformOutput', false );
queriesKnnIds = cellfun( @nonzeros, queriesKnnIds, 'UniformOutput', false );
	% queriesKnnIds: [1 x Q] cell vector, each cell is a column with query's knn indices

% Retrieve for each query its K nearest neighbors among the candidates retrieved from its matching buckets
KK = K * ones(1, Q);
queriesKnnIdsSizes = cellfun( @length, queriesKnnIds );
queriesKnnIdsSizes_lessThanK = queriesKnnIdsSizes < KK;
KK(queriesKnnIdsSizes_lessThanK) = queriesKnnIdsSizes(queriesKnnIdsSizes_lessThanK);
	
%queries_cellArray = arrayfun( @(queryIndex) queries(:, queryIndex), 1:Q, 'UniformOutput', false );
queries_cellArray = mat2cell(queries, D, ones(1, Q));
queryDatasets_cellArray = cellfun(@(indices) X(:, indices), queriesKnnIds, 'UniformOutput', false);
[~, knnIds] = cellfun(@knn, queries_cellArray, queryDatasets_cellArray, num2cell(KK), 'UniformOutput', false);
	% knnIds: [1 x Q] cell vector, each cell contains [KK x 1] indices
	% Indices refer to queries' candidates

% Transform knnIds to indices in the initial dataset X
knnIds = cellfun(@(vector, index) vector(index), queriesKnnIds, knnIds, 'UniformOutput', false);
% In case KK < K neighbors were retrieved for a query, pad K - KK zeros to its results
zerosNeeded = K*ones(1, Q) - KK;
knnIds = cellfun(@(x, dimension1) cat(1, x, zeros(dimension1, 1)), knnIds, num2cell(zerosNeeded), 'UniformOutput', false);
ids = cell2mat(knnIds);
%	ids: [K x Q] nearest neighbor indices

end

