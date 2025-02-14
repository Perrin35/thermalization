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
rz(0.67983812) q[1];
sx q[1];
rz(-1.4601269) q[1];
sx q[1];
rz(2.3753994) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1575398) q[0];
sx q[0];
rz(-1.6601642) q[0];
sx q[0];
rz(-1.2045226) q[0];
rz(0.030969663) q[2];
sx q[2];
rz(-1.3898347) q[2];
sx q[2];
rz(2.8029) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70464221) q[1];
sx q[1];
rz(-0.54848236) q[1];
sx q[1];
rz(1.30634) q[1];
rz(-pi) q[2];
rz(2.2965374) q[3];
sx q[3];
rz(-0.18394205) q[3];
sx q[3];
rz(-1.9234274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.55521518) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(-0.66184723) q[2];
rz(1.2343179) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(2.9890649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2374275) q[0];
sx q[0];
rz(-2.4894042) q[0];
sx q[0];
rz(-2.7320614) q[0];
rz(-0.25853363) q[1];
sx q[1];
rz(-1.2068318) q[1];
sx q[1];
rz(-2.5915367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058592168) q[0];
sx q[0];
rz(-2.6592836) q[0];
sx q[0];
rz(1.7586209) q[0];
rz(-pi) q[1];
rz(-1.0070877) q[2];
sx q[2];
rz(-1.9795784) q[2];
sx q[2];
rz(0.43540149) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2376783) q[1];
sx q[1];
rz(-1.774607) q[1];
sx q[1];
rz(-2.0173748) q[1];
rz(-pi) q[2];
rz(-1.9411381) q[3];
sx q[3];
rz(-2.3562288) q[3];
sx q[3];
rz(2.0919652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5162158) q[2];
sx q[2];
rz(-1.3863401) q[2];
sx q[2];
rz(0.59164444) q[2];
rz(-0.51724616) q[3];
sx q[3];
rz(-1.1040684) q[3];
sx q[3];
rz(0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.2194694) q[0];
sx q[0];
rz(-1.2930433) q[0];
sx q[0];
rz(-2.1406232) q[0];
rz(1.7132022) q[1];
sx q[1];
rz(-0.75894558) q[1];
sx q[1];
rz(0.044011291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3428033) q[0];
sx q[0];
rz(-2.2384423) q[0];
sx q[0];
rz(-1.081388) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5147382) q[2];
sx q[2];
rz(-0.73768697) q[2];
sx q[2];
rz(-0.53327582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7928501) q[1];
sx q[1];
rz(-2.2420635) q[1];
sx q[1];
rz(2.491374) q[1];
rz(-pi) q[2];
rz(0.040626133) q[3];
sx q[3];
rz(-1.7425795) q[3];
sx q[3];
rz(-0.21200997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9192146) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(3.0317793) q[2];
rz(-2.5114656) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(0.28756791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6269864) q[0];
sx q[0];
rz(-2.0872748) q[0];
sx q[0];
rz(-2.0475673) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(-0.72023645) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65533017) q[0];
sx q[0];
rz(-1.2661326) q[0];
sx q[0];
rz(0.49691864) q[0];
rz(0.14671756) q[2];
sx q[2];
rz(-2.8302022) q[2];
sx q[2];
rz(1.3757094) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1853237) q[1];
sx q[1];
rz(-0.8517304) q[1];
sx q[1];
rz(1.9732425) q[1];
rz(-2.8232326) q[3];
sx q[3];
rz(-0.84066072) q[3];
sx q[3];
rz(2.1206253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56095162) q[2];
sx q[2];
rz(-1.6683942) q[2];
sx q[2];
rz(-2.9894323) q[2];
rz(1.8469319) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(1.7744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072435943) q[0];
sx q[0];
rz(-2.5788588) q[0];
sx q[0];
rz(0.38134545) q[0];
rz(-2.5533679) q[1];
sx q[1];
rz(-1.4071848) q[1];
sx q[1];
rz(0.064195976) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4669111) q[0];
sx q[0];
rz(-1.3935601) q[0];
sx q[0];
rz(-2.7655274) q[0];
rz(1.1663127) q[2];
sx q[2];
rz(-2.1475386) q[2];
sx q[2];
rz(-2.8167644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5178495) q[1];
sx q[1];
rz(-0.84539834) q[1];
sx q[1];
rz(-2.7571667) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0392239) q[3];
sx q[3];
rz(-0.23766773) q[3];
sx q[3];
rz(0.32581926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75307816) q[2];
sx q[2];
rz(-1.831275) q[2];
sx q[2];
rz(0.7927967) q[2];
rz(-0.13475284) q[3];
sx q[3];
rz(-0.56551352) q[3];
sx q[3];
rz(-1.6511065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(-2.1886574) q[0];
rz(-1.8903525) q[1];
sx q[1];
rz(-1.6798991) q[1];
sx q[1];
rz(2.9073263) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15794663) q[0];
sx q[0];
rz(-1.0704213) q[0];
sx q[0];
rz(1.7364794) q[0];
rz(-pi) q[1];
rz(1.9769922) q[2];
sx q[2];
rz(-1.9941386) q[2];
sx q[2];
rz(-2.6273531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7997879) q[1];
sx q[1];
rz(-1.4569286) q[1];
sx q[1];
rz(0.98462414) q[1];
rz(-3.0563596) q[3];
sx q[3];
rz(-0.93795034) q[3];
sx q[3];
rz(-1.3229602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17144063) q[2];
sx q[2];
rz(-2.5415387) q[2];
sx q[2];
rz(-0.94685143) q[2];
rz(1.4417449) q[3];
sx q[3];
rz(-1.3274679) q[3];
sx q[3];
rz(-3.0202151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117208) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(0.3748931) q[0];
rz(-0.74294576) q[1];
sx q[1];
rz(-0.82287794) q[1];
sx q[1];
rz(-2.4853562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0716055) q[0];
sx q[0];
rz(-0.90815571) q[0];
sx q[0];
rz(1.2394942) q[0];
x q[1];
rz(1.0032651) q[2];
sx q[2];
rz(-2.1250688) q[2];
sx q[2];
rz(-2.2330333) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9083169) q[1];
sx q[1];
rz(-1.9413949) q[1];
sx q[1];
rz(0.45062391) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8493622) q[3];
sx q[3];
rz(-0.61200313) q[3];
sx q[3];
rz(-2.0296869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6132505) q[2];
sx q[2];
rz(-0.87317792) q[2];
sx q[2];
rz(2.4916416) q[2];
rz(-2.2925099) q[3];
sx q[3];
rz(-0.63747469) q[3];
sx q[3];
rz(1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717302) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(0.318203) q[0];
rz(-2.8089306) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(-0.37700787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.862826) q[0];
sx q[0];
rz(-1.5246848) q[0];
sx q[0];
rz(1.9079303) q[0];
rz(0.035786619) q[2];
sx q[2];
rz(-1.1038968) q[2];
sx q[2];
rz(0.70718471) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.711447) q[1];
sx q[1];
rz(-2.3466169) q[1];
sx q[1];
rz(-2.6076016) q[1];
rz(-0.11628233) q[3];
sx q[3];
rz(-2.0512037) q[3];
sx q[3];
rz(-1.8745599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13410021) q[2];
sx q[2];
rz(-1.8635805) q[2];
sx q[2];
rz(0.93097869) q[2];
rz(-1.603294) q[3];
sx q[3];
rz(-1.8036017) q[3];
sx q[3];
rz(-1.3163756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285434) q[0];
sx q[0];
rz(-2.7340524) q[0];
sx q[0];
rz(-2.7000309) q[0];
rz(2.8764797) q[1];
sx q[1];
rz(-1.6997063) q[1];
sx q[1];
rz(1.5107059) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1800507) q[0];
sx q[0];
rz(-1.7921711) q[0];
sx q[0];
rz(-0.44261296) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8042107) q[2];
sx q[2];
rz(-2.1357015) q[2];
sx q[2];
rz(3.1021822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9471252) q[1];
sx q[1];
rz(-1.8034916) q[1];
sx q[1];
rz(-0.10242771) q[1];
rz(-pi) q[2];
rz(-1.1006894) q[3];
sx q[3];
rz(-0.70610922) q[3];
sx q[3];
rz(2.2694893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6049217) q[2];
sx q[2];
rz(-1.8646381) q[2];
sx q[2];
rz(0.46860179) q[2];
rz(2.129715) q[3];
sx q[3];
rz(-1.7185271) q[3];
sx q[3];
rz(2.0603777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6877947) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(1.4367562) q[0];
rz(-0.94802481) q[1];
sx q[1];
rz(-1.3399622) q[1];
sx q[1];
rz(-0.2737793) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.439246) q[0];
sx q[0];
rz(-2.9986023) q[0];
sx q[0];
rz(2.6197394) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5735515) q[2];
sx q[2];
rz(-2.2170137) q[2];
sx q[2];
rz(0.60233145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35008365) q[1];
sx q[1];
rz(-2.4285762) q[1];
sx q[1];
rz(-0.75292315) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1115361) q[3];
sx q[3];
rz(-0.74990898) q[3];
sx q[3];
rz(-0.53242079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56129366) q[2];
sx q[2];
rz(-2.0686006) q[2];
sx q[2];
rz(2.6936074) q[2];
rz(2.5042846) q[3];
sx q[3];
rz(-1.0261122) q[3];
sx q[3];
rz(0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4056024) q[0];
sx q[0];
rz(-1.6163106) q[0];
sx q[0];
rz(1.6769782) q[0];
rz(-0.51392705) q[1];
sx q[1];
rz(-2.3873139) q[1];
sx q[1];
rz(0.25968459) q[1];
rz(2.6988637) q[2];
sx q[2];
rz(-1.6479658) q[2];
sx q[2];
rz(-2.992291) q[2];
rz(0.00058978422) q[3];
sx q[3];
rz(-1.1564217) q[3];
sx q[3];
rz(0.52903928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
