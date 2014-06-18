function [ lshStruct ] = lshConstruct( X, L, M, W, A, b )
% function [ lshStruct ] = lshConstruct( X, L, M, W, A, b )
%
% Construct LSH index structure for dataset matrix X
%
% INPUT:
%	X: [D x N] matrix containing N columns of D-dimensional data
%	L: # of projections to perform for indexing
%	M: # of dimensions at projection spaces
%
%	W: {optional} bucket width
%	A: {optional} [1 x L] cell array, each cell is a [D x M] matrix
%	b: {optional} [1 x L] cell array, each cell is a [1 x M] vector with values in range [0, W)
%
% If A, b are not given, they are generated randomly.
% If W is not given, it is computed like
%	Gregory Shakhnarovich does in his LSH toolkit.
%
% OUTPUT:
%	lshStruct: [1 x L] array of structs(one for each projection) with fields:
%		A:       [D x M] multiplication matrix
%		w:       bucket width (=W)
%		b:       [1 x M] bucket quantization values
%		low:     [1 x M] minimum values of data for each dimension at projection space
%		codes:   [B x M] rows of codes corresponding to the buckets
%		pointId: [1 x B] cell array with point indices for each bucket
%		sizes:   [1 x B] sizes of buckets
%		bhash:   cell array used for 2nd level hashing
%
% Projection is done this way:
%	Y = floor( X'*A - b) / w);
%	Afterwards, the minimum values at each dimension are subtracted, so that Y>=0
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

D = size(X, 1);

% A must be [1 x L] cell vector, each cell must be a [D x M] matrix
if(exist('A', 'var'))
	if(~iscell(A) || length(A)~=L)
		fprintf(2, 'A must be [1 x L] cell vector.\n');
		return;
	end
	for i=1:L
		[dim1 dim2] = size(A{i});
		if(dim1~=D || dim2~=M)
			fprintf(2, 'Each cell of A must be a [D x M] matrix.\n');
			return;
		end
	end
end

% b must be [1 x L] cell vector, each cell must be a [1 x M] vector
if(exist('b', 'var'))
	if(~iscell(b) || length(b)~=L)
		fprintf(2, 'b must be [1 x L] cell vector.\n');
		return;
	end
	for i=1:L
		[dim1 dim2] = size(b{i});
		if(dim1~=1 || dim2~=M)
			fprintf(2, 'Each cell of b must be a [1 x M] matrix.\n');
			return;
		end
		if(min(b{i})<0 || max(b{i})>=W)
			fprintf(2, 'Each cell of b must have values in range [0, W)\n');
			return;
		end
	end
end

% If W was not given, compute it in Gregory's way
if(~exist('W', 'var'))
	range = [min(X, [], 2)'; max(X, [], 2)'];
	limits = max(abs(range), [], 1);
%	limits = max(abs(range(1,:)),abs(range(2,:)));
    rangeAct=mean(diff([-limits; limits]*2*sqrt(D)));
    W = rangeAct/16;
	fprintf('Computed W: %f\n', W);
end

% Create L tables
tstart = tic;
for l=1:L
	fprintf('Creating table %d/%d\n', l, L);
	% Define projection parameter A ( [D x M] )
	if(exist('A', 'var'))
		lshStruct(l).A = A{l};
	else
		lshStruct(l).A = randn(D, M);
	end
	% Define projection parameter w ( scalar )
	lshStruct(l).w = W;
	% Define projection parameter b ( [1 x M] )
	if(exist('b', 'var'))
		lshStruct(l).b = b{l};
	else
		lshStruct(l).b = unifrnd(0, lshStruct(l).w, 1, M);
	end
	
	% Project data
	fprintf('\tProjecting data...\n');
	Y = floor( bsxfun(@minus, X'*lshStruct(l).A, lshStruct(l).b) / lshStruct(l).w);	% [N x M]
	
		lshStruct(l).low = min(Y, [], 1);	% low: [1 x D]
	
	Y = bsxfun(@minus, Y, lshStruct(l).low);	% [N x M] >= 0
	% Now Y>=0 => keys=lshhash(Y)>0)
	
	% Identify unique bucket codes and their belonging points
	[Ycodes, ~, Ybin_groups] = unique(Y, 'rows');
		% Ycodes: [B x M] rows of unique codes
		% Ybin_groups: [N x 1], Y = Ycodes(Ybin_groups, :)
	[~, Ybi] = sort(Ybin_groups);	% [N x 1]
		% Ybi: [N x 1] bucket ids in order of ascending id of unique bucket codes
	Ysizes = histc(Ybin_groups, unique(Ybin_groups));
		% Ysizes: [B x 1] sizes of buckets
	YpointId = mat2cell(Ybi', 1, Ysizes);
		% YpointId: [1 x B] cell vector of belonging point indices for each bucket

	fprintf('\t%d buckets\n', length(Ysizes));
	fprintf('\tmin / max size: %d / %d\n', min(Ysizes), max(Ysizes));
	fprintf('\tmean / median: %d / %d\n', mean(Ysizes), median(Ysizes));
	
		lshStruct(l).codes = Ycodes;		% [B x M] rows of codes
		lshStruct(l).pointId = YpointId;	% [1 x B] cell vector of point ids per bucket
		lshStruct(l).sizes = Ysizes';		% [1 x B] sizes of buckets
	
	% Perform 2nd level hashing of the codes
	keys = lshhash(lshStruct(l).codes);
		% keys: [B x 1] scalar 2nd level hash keys
	keysUnique = unique(keys);
		% keysUnique: [UK x 1] unique keys
	keysSizes = histc(keys, keysUnique);
		% keysSizes: [UK x 1] # of occurences of each unique key
	[~, keysIds] = sort(keys);
		% keyIds: [B x 1] key ids in order of ascending value of unique key
	keysGroups = mat2cell(keysIds', 1, keysSizes);
		% keyGroups: [1 x UK] cell vector of indices to buckets with same key
	
	% lshStruct(l).bhash{value} contains ids to buckets with 2nd level hash key = <value>
		lshStruct(l).bhash = cell(1, max(keys));
		lshStruct(l).bhash(keysUnique) = keysGroups;
end

fprintf(2, 'Total time: %f\n', toc(tstart));

end

