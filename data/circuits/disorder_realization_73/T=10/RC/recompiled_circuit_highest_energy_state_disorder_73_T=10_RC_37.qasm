OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4740648) q[0];
sx q[0];
rz(-0.83905554) q[0];
sx q[0];
rz(0.43247217) q[0];
rz(-2.4617545) q[1];
sx q[1];
rz(-1.6814657) q[1];
sx q[1];
rz(-2.3753994) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941152) q[0];
sx q[0];
rz(-1.9355398) q[0];
sx q[0];
rz(0.095679521) q[0];
x q[1];
rz(1.4031448) q[2];
sx q[2];
rz(-2.9580287) q[2];
sx q[2];
rz(-2.6324384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7442465) q[1];
sx q[1];
rz(-2.098171) q[1];
sx q[1];
rz(0.1583734) q[1];
rz(-pi) q[2];
rz(-0.12285439) q[3];
sx q[3];
rz(-1.4335504) q[3];
sx q[3];
rz(0.48396971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.55521518) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(-2.4797454) q[2];
rz(-1.9072748) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(2.9890649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9041651) q[0];
sx q[0];
rz(-2.4894042) q[0];
sx q[0];
rz(-2.7320614) q[0];
rz(-0.25853363) q[1];
sx q[1];
rz(-1.9347609) q[1];
sx q[1];
rz(-0.55005598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15274238) q[0];
sx q[0];
rz(-2.0439195) q[0];
sx q[0];
rz(0.097448601) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0070877) q[2];
sx q[2];
rz(-1.9795784) q[2];
sx q[2];
rz(-0.43540149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0747968) q[1];
sx q[1];
rz(-2.6535712) q[1];
sx q[1];
rz(-1.1244494) q[1];
rz(-0.34724473) q[3];
sx q[3];
rz(-0.85111516) q[3];
sx q[3];
rz(-1.5517476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6253769) q[2];
sx q[2];
rz(-1.7552525) q[2];
sx q[2];
rz(-2.5499482) q[2];
rz(2.6243465) q[3];
sx q[3];
rz(-2.0375242) q[3];
sx q[3];
rz(-0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2194694) q[0];
sx q[0];
rz(-1.2930433) q[0];
sx q[0];
rz(-1.0009694) q[0];
rz(-1.4283904) q[1];
sx q[1];
rz(-0.75894558) q[1];
sx q[1];
rz(-3.0975814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.509217) q[0];
sx q[0];
rz(-0.80501834) q[0];
sx q[0];
rz(-0.53765996) q[0];
rz(1.6268544) q[2];
sx q[2];
rz(-0.73768697) q[2];
sx q[2];
rz(2.6083168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7928501) q[1];
sx q[1];
rz(-0.89952911) q[1];
sx q[1];
rz(2.491374) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3408362) q[3];
sx q[3];
rz(-0.17647568) q[3];
sx q[3];
rz(-0.0214487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22237805) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(-0.10981336) q[2];
rz(-0.63012704) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(-0.28756791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6269864) q[0];
sx q[0];
rz(-1.0543178) q[0];
sx q[0];
rz(-1.0940254) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(2.4213562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0648616) q[0];
sx q[0];
rz(-2.0429049) q[0];
sx q[0];
rz(1.2272627) q[0];
rz(-1.5237775) q[2];
sx q[2];
rz(-1.2628619) q[2];
sx q[2];
rz(1.5297254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1853237) q[1];
sx q[1];
rz(-0.8517304) q[1];
sx q[1];
rz(-1.1683502) q[1];
x q[2];
rz(-0.31836001) q[3];
sx q[3];
rz(-2.3009319) q[3];
sx q[3];
rz(2.1206253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56095162) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(0.15216039) q[2];
rz(-1.2946607) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(-1.367182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0691567) q[0];
sx q[0];
rz(-2.5788588) q[0];
sx q[0];
rz(-0.38134545) q[0];
rz(2.5533679) q[1];
sx q[1];
rz(-1.7344079) q[1];
sx q[1];
rz(0.064195976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8255257) q[0];
sx q[0];
rz(-0.41393201) q[0];
sx q[0];
rz(2.6878306) q[0];
x q[1];
rz(0.61582404) q[2];
sx q[2];
rz(-1.9069789) q[2];
sx q[2];
rz(1.6662963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5178495) q[1];
sx q[1];
rz(-2.2961943) q[1];
sx q[1];
rz(-2.7571667) q[1];
x q[2];
rz(-2.1023687) q[3];
sx q[3];
rz(-2.9039249) q[3];
sx q[3];
rz(2.8157734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3885145) q[2];
sx q[2];
rz(-1.3103176) q[2];
sx q[2];
rz(0.7927967) q[2];
rz(0.13475284) q[3];
sx q[3];
rz(-2.5760791) q[3];
sx q[3];
rz(1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(-2.1886574) q[0];
rz(-1.2512402) q[1];
sx q[1];
rz(-1.6798991) q[1];
sx q[1];
rz(-2.9073263) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4933245) q[0];
sx q[0];
rz(-0.52487322) q[0];
sx q[0];
rz(0.29294818) q[0];
rz(-pi) q[1];
rz(1.1646005) q[2];
sx q[2];
rz(-1.9941386) q[2];
sx q[2];
rz(-0.51423954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3418048) q[1];
sx q[1];
rz(-1.4569286) q[1];
sx q[1];
rz(-2.1569685) q[1];
rz(-pi) q[2];
rz(-1.4552507) q[3];
sx q[3];
rz(-2.503814) q[3];
sx q[3];
rz(-1.4664283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.970152) q[2];
sx q[2];
rz(-2.5415387) q[2];
sx q[2];
rz(-2.1947412) q[2];
rz(1.4417449) q[3];
sx q[3];
rz(-1.3274679) q[3];
sx q[3];
rz(0.12137752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22987188) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(-2.7666996) q[0];
rz(-0.74294576) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(-0.65623647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5798168) q[0];
sx q[0];
rz(-0.72951376) q[0];
sx q[0];
rz(-0.39493409) q[0];
rz(0.63323449) q[2];
sx q[2];
rz(-1.0960964) q[2];
sx q[2];
rz(-2.1556319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23327574) q[1];
sx q[1];
rz(-1.9413949) q[1];
sx q[1];
rz(-0.45062391) q[1];
x q[2];
rz(-2.5498059) q[3];
sx q[3];
rz(-1.7370708) q[3];
sx q[3];
rz(0.70032728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6132505) q[2];
sx q[2];
rz(-2.2684147) q[2];
sx q[2];
rz(0.64995107) q[2];
rz(2.2925099) q[3];
sx q[3];
rz(-2.504118) q[3];
sx q[3];
rz(-1.1879296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7717302) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(-0.318203) q[0];
rz(-2.8089306) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(2.7645848) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.862826) q[0];
sx q[0];
rz(-1.5246848) q[0];
sx q[0];
rz(-1.9079303) q[0];
rz(-pi) q[1];
rz(-3.105806) q[2];
sx q[2];
rz(-1.1038968) q[2];
sx q[2];
rz(-2.4344079) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1313427) q[1];
sx q[1];
rz(-0.90908644) q[1];
sx q[1];
rz(2.0493839) q[1];
rz(-pi) q[2];
x q[2];
rz(1.351736) q[3];
sx q[3];
rz(-2.6483834) q[3];
sx q[3];
rz(-1.0194609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0074924) q[2];
sx q[2];
rz(-1.8635805) q[2];
sx q[2];
rz(-2.210614) q[2];
rz(1.603294) q[3];
sx q[3];
rz(-1.337991) q[3];
sx q[3];
rz(-1.3163756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304921) q[0];
sx q[0];
rz(-0.40754023) q[0];
sx q[0];
rz(2.7000309) q[0];
rz(0.26511296) q[1];
sx q[1];
rz(-1.6997063) q[1];
sx q[1];
rz(-1.5107059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.961542) q[0];
sx q[0];
rz(-1.3494216) q[0];
sx q[0];
rz(0.44261296) q[0];
rz(-2.052102) q[2];
sx q[2];
rz(-0.64844744) q[2];
sx q[2];
rz(2.5221053) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35263044) q[1];
sx q[1];
rz(-1.6704541) q[1];
sx q[1];
rz(-1.3369192) q[1];
rz(-pi) q[2];
rz(-0.36863736) q[3];
sx q[3];
rz(-2.1876699) q[3];
sx q[3];
rz(-1.6807238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.536671) q[2];
sx q[2];
rz(-1.8646381) q[2];
sx q[2];
rz(-0.46860179) q[2];
rz(-2.129715) q[3];
sx q[3];
rz(-1.7185271) q[3];
sx q[3];
rz(-2.0603777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6877947) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(1.7048365) q[0];
rz(2.1935678) q[1];
sx q[1];
rz(-1.8016305) q[1];
sx q[1];
rz(0.2737793) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38588148) q[0];
sx q[0];
rz(-1.4997002) q[0];
sx q[0];
rz(-3.0174251) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4953734) q[2];
sx q[2];
rz(-1.5685967) q[2];
sx q[2];
rz(-0.96680582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.791509) q[1];
sx q[1];
rz(-2.4285762) q[1];
sx q[1];
rz(-2.3886695) q[1];
rz(-pi) q[2];
rz(2.2446452) q[3];
sx q[3];
rz(-1.9292783) q[3];
sx q[3];
rz(1.4523539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56129366) q[2];
sx q[2];
rz(-1.0729921) q[2];
sx q[2];
rz(-0.44798526) q[2];
rz(2.5042846) q[3];
sx q[3];
rz(-2.1154805) q[3];
sx q[3];
rz(-0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7359903) q[0];
sx q[0];
rz(-1.525282) q[0];
sx q[0];
rz(-1.4646144) q[0];
rz(-2.6276656) q[1];
sx q[1];
rz(-0.75427873) q[1];
sx q[1];
rz(-2.8819081) q[1];
rz(2.9630255) q[2];
sx q[2];
rz(-2.6926276) q[2];
sx q[2];
rz(-1.2603029) q[2];
rz(-0.00058978422) q[3];
sx q[3];
rz(-1.9851709) q[3];
sx q[3];
rz(-2.6125534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
