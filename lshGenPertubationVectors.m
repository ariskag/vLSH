function [ pertubation_sets pertubation_scores ] = lshGenPertubationVectors( M, T, deltas )
% function [ pertubation_sets pertubation_scores ] = lshGenPertubationVectors( M, T, probes )
%
% Generate T pertubation vectors for M-valued code
%
% INPUT:
%	M:      length of code
%	T:      # of pertubation sets to be produced
%	probes: [1 x 2*M] struct array with fields:
%		position: scalar value in range [1, M]
%		value:    query's code value at specified position
%		delta:    difference between query's projection and border's code at specified position
%		border:   -1 when referring to the left (i.e. query's code value at
%		             specified position minus 1),
%		          +1 when referring to the right (i.e. query's code value
%		             at specified position plus 1)
% Probes must already be sorted in order of increasing DELTA values!!!
%
% Each element of query's M-valued code can be either incremented or subtracted by 1
% Therefore probes has 2*M elements:
%	-> M with deltas corresponding to the incremented version of each value,
%	-> M with deltas corresponding to the subtracted version of each value
%
% OUTPUT:
%	pertubation_sets:   [1 x T] cell array of pertubation vectors, each cell
%	                     contains indices to probes belonging to pertubation vector
%	pertubation_scores: [1 x T] scores, i.e. sums of deltas of
%	                     probes belonging to each pertubation vector
%
% The pertubation vectors are computed according to paper
%  MultiProbe LSH: Efficient Indexing for High-Dimensional Similarity Search,
%  by Qin Lv, William Josephson, Zhe Wang, Moses Charikar, Kai Li.

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

A0 = {1};
A0_scores = deltas(1);

pertubation_sets = cell(1, T);
pertubation_scores = zeros(1, T);

% Generate T pertubation vectors
for i=1:T
	min_score_pos = find(A0_scores==min(A0_scores), 1);
		% min_score_pos: index to element in A0 with smallest score
	
	% Element with smallest score is the next pertubation vector
	pertubation_sets(i) = A0(min_score_pos);
	pertubation_scores(i) = A0_scores(min_score_pos);
	
	% Remove element from A0, replace it with its shifted and its expanded versions
	A0(min_score_pos) = [];
	A0_scores(min_score_pos) = [];
	
	vector2add = pertubation_sets{i};
	vector2add_maxValue = max(vector2add);
	if(vector2add_maxValue<2*M)
		vector2add(vector2add==vector2add_maxValue) = vector2add_maxValue + 1;
		score2add = pertubation_scores(i) - deltas(vector2add_maxValue) + deltas(vector2add_maxValue + 1);
		A0{end+1} = vector2add;
		A0_scores(end+1) = score2add;
	
		vector2add = pertubation_sets{i};
		vector2add(end+1) = vector2add_maxValue + 1;
		score2add = pertubation_scores(i) + deltas(vector2add_maxValue+1);
		A0{end+1} = vector2add;
		A0_scores(end+1) = score2add;
	end
end

end

