//Number of population samples (demes)
5
//Population effective sizes (number of genes)
NE1
NE2
NE3
NE4
NE5
//Sample sizes
88
92
78
98
78
//Growth rates	: negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//migration rates matrix 0:
0.0 M12 0.0 0.0 0.0
M21 0.0 M23 0.0 0.0
0.0 M32 0.0 0.0 0.0
0.0 0.0 0.0 0.0 M54
0.0 0.0 0.0 M45 0.0
//migration rates matrix 1:
0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
4 historical event
TDIV 1 0 1 1 0 1
TDIV 2 0 1 1 0 1
TDIV 3 0 1 1 0 1
TDIV 4 0 1 1 0 1
//Number of independent loci [chromosome] 
10 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
MICROSAT 1 0 MUTR 0 0
