vLSH
====

Description
-----------

This directory contains a simple implementation of a Vectorized  Multiprobing Locality-Sensitive Hashing (LSH) algorithm based on [Greg Shakhnarovich](http://ttic.uchicago.edu/~gregory/download.html)'s algorithm. Its functions are:

* _lshConstruct.m_: Constructs the LSH index structure for dataset matrix X.
* _lshSearch.m_: Searches for queries' K nearest neighbors from dataset X using LSH.
* _lshFindQueryBins.m_: Given a query, it returns the codes of the probing bins to search for its nearest neighbors. It is called by lshSearch.
* _lshGenPertubationVectors.m_: Generates T pertubation vectors for a query, in order to compute its probing bins. It is called by lshFindQueryBins.
* _lshSearchTable.m_: Searches for queries' nearest neighbors among specific buckets. It is called by lshSearch for each one of the tables of an LSH structure, after their probing bins have been found via lshFindQueryBins.
* _lshhash.m_: Simple one-level hashing function to speed up bucket search in LSH. It is copied by [Gregory Shakhnarovich](http://ttic.uchicago.edu/~gregory/download.html)'s LSH toolkit.
* _knn.m_: Perform brute force search of queries' nearest neighbors.
* _compute\_recall.m_: Given two nearest neighbor sets, of which one contains the true nearest neighbor ids and the other LSH's returned nearest neighbor ids, it returns the recall of the LSH search (i.e. the size of the intersection of the two sets divided by K, where K is the number of nearest neighbors computed).

If for a query only KK < K neighbors were found (because it had less than K candidates), then 0 is returned as neighbor id for the (K - KK) entries in the id matrix.
For more details refer to each function's source code and comments.

Usage example
-------------

Suppose X is the dataset on which we want to perform Locality Sensitive Hashing.


	Q = 1000; % # of queries
	K = 32;   % # of nearest neighbors to search for each query
	L = 2;    % # of tables
	M = 10;   % # of dimensions at projection space
	W = 1000; % bucket width
	T = 10;   % # of additional probing bins

	% Construct index tables
	lshStruct = lshConstruct(X, L, M, W);

	% Perform simple LSH search on the first Q elements of X
	[idsSIMPLE cand_sizeSIMPLE binsFound_ratiosSIMPLE] = lshSearch(X(:, 1:Q), X, lshStruct, K);
	% is the same with: [idsSIMPLE cand_sizeSIMPLE binsFound_ratiosSIMPLE] = lshSearch(X(:, 1:Q), X, lshStruct, K, 0);

	% Perform multiprobe LSH search on the first Q elements of X
	[idsMULTIPROBE cand_sizeMULTIPROBE binsFound_ratiosMULTIPROBE] = lshSearch(X(:, 1:Q), X, lshStruct, K, T);

	% Perform brute force KNN on the first Q elements of X
	[~, idsTRUE] = knn(X(:, 1:Q), X, K);

	% Compute recalls
	recallSIMPLE = compute_recall(idsTRUE, idsSIMPLE);
	recallMULTIPROBE = compute_recall(idsTRUE, idsMULTIPROBE);

	% Compute selectivities (i.e. mean size of distinct candidates for each query)
	selectivitySIMPLE = sum(cand_sizeSIMPLE) / (Q * size(X, 2));
	selectivityMULTIPROBE = sum(cand_sizeMULTIPROBE) / (Q * size(X, 2));

