function [ binIndices ] = lshFindQueryBins( query, lshStructTable, T)
% function [ binIndices ] = lshFindQueryBins( query, lshStructTable, T)
%
% Given a query, return the probing bins to which its nearest neighbors belong
%
% INPUT:
%	query:          [D x 1] query data
%	lshStructTable: struct belonging to lshStruct, which is produced by lshConstruct
%	T:              {optional} # of additional probing buckets to return for
%	                    multiprobe searching (default=0, returns only one bin)
%
% lshStruct is a [1 x L] struct array
% lshStructTable is lshStruct(l), where l is the id of the
%   table for which we want to get the query's matching buckets
%
% OUTPUT:
%	binIndices: [1 x B] indices to buckets matching query
%
% The function makes use of 2nd level hashing by using Gregory
% Shakhnarovich's lshhash.m function, included in his LSH toolkit
% (http://ttic.uchicago.edu/~gregory/download.html)
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

% If T is not given, return only query's original bucket id
if(~exist('T', 'var'))
	T = 0;
end

queryCodeNonFloored = bsxfun(@minus, query'*lshStructTable.A, lshStructTable.b) / lshStructTable.w - lshStructTable.low;
	% queryCodeNonFloored: [1 x M] query's projection
queryCode = floor(queryCodeNonFloored);
	% queryCode: [1 x M] query's code (floored projection)
delta = queryCodeNonFloored - queryCode;
	% delta: [1 x M] difference between query's projection and code

M = length(queryCode);
	% M: # of dimensions at projection space

positions = [1:M 1:M];
deltas = [delta 1-delta];
borders = [-1*ones(1, M) ones(1, M)];
[~, ii] = sort(deltas);

% Generate T pertubation vectors with smallest overall delta scores
pb_sets = lshGenPertubationVectors(M, T, deltas(ii));

% Generate codes corresponding to original code and its pertubated versions
codes = repmat(queryCode, T+1, 1);
	% codes: [(T+1) x M], first row is original code
	% Generate pertubated versions of the original code
for i=1:T
	codes(i+1, positions(ii(pb_sets{i}))) = codes(i+1, positions(ii(pb_sets{i}))) + borders(ii(pb_sets{i}));
end

%{
if(size(codes, 1)~=size(unique(codes, 'rows'), 1))
	fprintf('!!!!!!!!!!!!!!!!!\n\tSIZE(CODES,1) = %d\n\tSIZE(UNIQUE(CODES, ''ROWS''), 1) = %d\N!!!!!!!!!!!!!!!!!\n', size(codes, 1), size(unique(codes, 'rows'), 1));
	filename = sprintf('log-%d', randi(10000));
	save(filename, 'codes', 'probes', 'M', 'T');
	fprintf('SAVED UNDER %s\n\n\n', filename);
	codes = unique(codes, 'rows');
end
%}

% Generate 2nd level hash values of generated codes
keys = lshhash(codes)';
	% keys: [1 x B] hash values of generated codes
	% Take into account the probability of keys<=0 or keys>length(lshStruct.bhash)
keys(keys<=0) = 1;
lshStructTable.bhash{end+1} = [];
keys(keys>length(lshStructTable.bhash)) = length(lshStructTable.bhash);

%codes(4, 3) = 100;

possiblyMatchingBins = lshStructTable.bhash(keys);
	% possiblyMatchingBins: [1 x (T+1)] cell vector, cell at position i contains indices to buckets possibly matching code codes(i, :)
possiblyMatchingBinCodes = cellfun( @(bucketIndices) lshStructTable.codes(bucketIndices, :), possiblyMatchingBins, 'UniformOutput', false );
	% possiblyMatchingBinCodes: [1 x (T+1)] cell vector, cell at position i contains codes of buckets indexed at possiblyMatchingBins(i)
fh_getMatchingIndexFromPossibleBins = @(queryCode, possibleBins, possibleBinCodes) possibleBins( all( bsxfun(@eq, possibleBinCodes, queryCode), 2 ) );
binIndices = cellfun( fh_getMatchingIndexFromPossibleBins, mat2cell(codes, ones(1, T+1), M)', possiblyMatchingBins, possiblyMatchingBinCodes, 'UniformOutput', false );
	% binIndices: [1 x (T+1)] cell array, each cell is a unique bucket index, or an empty array in case no matching bucket was found
	% Remove empty elements of binIndices
binIndices_emptyCells = cellfun(@isempty, binIndices);
binIndices = cell2mat(binIndices(~binIndices_emptyCells));
	% binIndices: [1 x B] bucket indices, where B is the # of matching buckets found

end

