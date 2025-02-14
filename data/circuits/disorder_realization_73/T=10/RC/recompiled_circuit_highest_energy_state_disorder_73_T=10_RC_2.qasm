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
rz(4.8230584) q[1];
sx q[1];
rz(11.800177) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941152) q[0];
sx q[0];
rz(-1.2060529) q[0];
sx q[0];
rz(0.095679521) q[0];
rz(-1.3897498) q[2];
sx q[2];
rz(-1.5403325) q[2];
sx q[2];
rz(-1.2265282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.39734617) q[1];
sx q[1];
rz(-2.098171) q[1];
sx q[1];
rz(-2.9832193) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12285439) q[3];
sx q[3];
rz(-1.4335504) q[3];
sx q[3];
rz(0.48396971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5863775) q[2];
sx q[2];
rz(-2.5368498) q[2];
sx q[2];
rz(0.66184723) q[2];
rz(-1.9072748) q[3];
sx q[3];
rz(-1.581278) q[3];
sx q[3];
rz(0.1525277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9041651) q[0];
sx q[0];
rz(-0.65218848) q[0];
sx q[0];
rz(2.7320614) q[0];
rz(-2.883059) q[1];
sx q[1];
rz(-1.2068318) q[1];
sx q[1];
rz(2.5915367) q[1];
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
rz(-1.3829718) q[0];
rz(-0.47357719) q[2];
sx q[2];
rz(-1.0583631) q[2];
sx q[2];
rz(0.88922516) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0747968) q[1];
sx q[1];
rz(-2.6535712) q[1];
sx q[1];
rz(-2.0171432) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34724473) q[3];
sx q[3];
rz(-2.2904775) q[3];
sx q[3];
rz(1.5517476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6253769) q[2];
sx q[2];
rz(-1.7552525) q[2];
sx q[2];
rz(0.59164444) q[2];
rz(0.51724616) q[3];
sx q[3];
rz(-2.0375242) q[3];
sx q[3];
rz(0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9221233) q[0];
sx q[0];
rz(-1.2930433) q[0];
sx q[0];
rz(2.1406232) q[0];
rz(-1.4283904) q[1];
sx q[1];
rz(-2.3826471) q[1];
sx q[1];
rz(-0.044011291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.509217) q[0];
sx q[0];
rz(-0.80501834) q[0];
sx q[0];
rz(-2.6039327) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3077009) q[2];
sx q[2];
rz(-1.5331037) q[2];
sx q[2];
rz(2.1455763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90755573) q[1];
sx q[1];
rz(-2.2436972) q[1];
sx q[1];
rz(-0.91958811) q[1];
rz(-pi) q[2];
rz(3.1009665) q[3];
sx q[3];
rz(-1.3990132) q[3];
sx q[3];
rz(2.9295827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9192146) q[2];
sx q[2];
rz(-1.5488449) q[2];
sx q[2];
rz(3.0317793) q[2];
rz(-2.5114656) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(-2.8540247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5146062) q[0];
sx q[0];
rz(-1.0543178) q[0];
sx q[0];
rz(2.0475673) q[0];
rz(-0.47359723) q[1];
sx q[1];
rz(-2.4351599) q[1];
sx q[1];
rz(2.4213562) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65533017) q[0];
sx q[0];
rz(-1.8754601) q[0];
sx q[0];
rz(0.49691864) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14671756) q[2];
sx q[2];
rz(-0.31139049) q[2];
sx q[2];
rz(-1.7658833) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2536972) q[1];
sx q[1];
rz(-1.2716537) q[1];
sx q[1];
rz(0.76051416) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81497808) q[3];
sx q[3];
rz(-1.3354104) q[3];
sx q[3];
rz(-2.8081364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56095162) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(-2.9894323) q[2];
rz(-1.2946607) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(-1.367182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0691567) q[0];
sx q[0];
rz(-0.56273383) q[0];
sx q[0];
rz(0.38134545) q[0];
rz(-2.5533679) q[1];
sx q[1];
rz(-1.4071848) q[1];
sx q[1];
rz(0.064195976) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4669111) q[0];
sx q[0];
rz(-1.7480326) q[0];
sx q[0];
rz(2.7655274) q[0];
x q[1];
rz(0.54406329) q[2];
sx q[2];
rz(-0.69101101) q[2];
sx q[2];
rz(0.34073433) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0654352) q[1];
sx q[1];
rz(-0.80428752) q[1];
sx q[1];
rz(1.1706513) q[1];
rz(3.0194123) q[3];
sx q[3];
rz(-1.3664277) q[3];
sx q[3];
rz(2.271657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3885145) q[2];
sx q[2];
rz(-1.831275) q[2];
sx q[2];
rz(2.348796) q[2];
rz(-3.0068398) q[3];
sx q[3];
rz(-2.5760791) q[3];
sx q[3];
rz(-1.6511065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7029) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(-0.95293522) q[0];
rz(1.2512402) q[1];
sx q[1];
rz(-1.4616936) q[1];
sx q[1];
rz(-2.9073263) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3327985) q[0];
sx q[0];
rz(-1.7160128) q[0];
sx q[0];
rz(0.50620075) q[0];
x q[1];
rz(2.6855747) q[2];
sx q[2];
rz(-1.9393183) q[2];
sx q[2];
rz(-0.8816661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7997879) q[1];
sx q[1];
rz(-1.684664) q[1];
sx q[1];
rz(2.1569685) q[1];
x q[2];
rz(-0.085233099) q[3];
sx q[3];
rz(-2.2036423) q[3];
sx q[3];
rz(1.8186325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.970152) q[2];
sx q[2];
rz(-0.60005391) q[2];
sx q[2];
rz(2.1947412) q[2];
rz(1.4417449) q[3];
sx q[3];
rz(-1.8141247) q[3];
sx q[3];
rz(3.0202151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.22987188) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(-0.3748931) q[0];
rz(-2.3986469) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(0.65623647) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0716055) q[0];
sx q[0];
rz(-0.90815571) q[0];
sx q[0];
rz(-1.2394942) q[0];
x q[1];
rz(2.1383275) q[2];
sx q[2];
rz(-1.0165239) q[2];
sx q[2];
rz(-2.2330333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23327574) q[1];
sx q[1];
rz(-1.9413949) q[1];
sx q[1];
rz(-2.6909687) q[1];
rz(2.8493622) q[3];
sx q[3];
rz(-2.5295895) q[3];
sx q[3];
rz(-2.0296869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6132505) q[2];
sx q[2];
rz(-0.87317792) q[2];
sx q[2];
rz(-2.4916416) q[2];
rz(-2.2925099) q[3];
sx q[3];
rz(-0.63747469) q[3];
sx q[3];
rz(-1.1879296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36986247) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(-0.318203) q[0];
rz(0.33266208) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(2.7645848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2787667) q[0];
sx q[0];
rz(-1.5246848) q[0];
sx q[0];
rz(1.2336624) q[0];
rz(-pi) q[1];
rz(-2.0379534) q[2];
sx q[2];
rz(-1.5388414) q[2];
sx q[2];
rz(-0.87972537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.01025) q[1];
sx q[1];
rz(-0.90908644) q[1];
sx q[1];
rz(2.0493839) q[1];
x q[2];
rz(-2.0539862) q[3];
sx q[3];
rz(-1.6738664) q[3];
sx q[3];
rz(-2.7838992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0074924) q[2];
sx q[2];
rz(-1.2780122) q[2];
sx q[2];
rz(2.210614) q[2];
rz(1.603294) q[3];
sx q[3];
rz(-1.8036017) q[3];
sx q[3];
rz(-1.8252171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285434) q[0];
sx q[0];
rz(-0.40754023) q[0];
sx q[0];
rz(0.44156179) q[0];
rz(2.8764797) q[1];
sx q[1];
rz(-1.4418863) q[1];
sx q[1];
rz(1.6308867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.961542) q[0];
sx q[0];
rz(-1.7921711) q[0];
sx q[0];
rz(0.44261296) q[0];
x q[1];
rz(0.33738193) q[2];
sx q[2];
rz(-1.0058912) q[2];
sx q[2];
rz(-0.039410465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7889622) q[1];
sx q[1];
rz(-1.6704541) q[1];
sx q[1];
rz(-1.8046735) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2208414) q[3];
sx q[3];
rz(-1.8691321) q[3];
sx q[3];
rz(-0.32978299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.536671) q[2];
sx q[2];
rz(-1.2769545) q[2];
sx q[2];
rz(2.6729909) q[2];
rz(-2.129715) q[3];
sx q[3];
rz(-1.4230655) q[3];
sx q[3];
rz(-1.0812149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6877947) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(1.4367562) q[0];
rz(2.1935678) q[1];
sx q[1];
rz(-1.3399622) q[1];
sx q[1];
rz(-0.2737793) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7557112) q[0];
sx q[0];
rz(-1.6418925) q[0];
sx q[0];
rz(3.0174251) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5735515) q[2];
sx q[2];
rz(-0.92457891) q[2];
sx q[2];
rz(2.5392612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35008365) q[1];
sx q[1];
rz(-0.71301642) q[1];
sx q[1];
rz(-0.75292315) q[1];
rz(-2.1115361) q[3];
sx q[3];
rz(-0.74990898) q[3];
sx q[3];
rz(2.6091719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.580299) q[2];
sx q[2];
rz(-1.0729921) q[2];
sx q[2];
rz(0.44798526) q[2];
rz(2.5042846) q[3];
sx q[3];
rz(-2.1154805) q[3];
sx q[3];
rz(-0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7359903) q[0];
sx q[0];
rz(-1.6163106) q[0];
sx q[0];
rz(1.6769782) q[0];
rz(0.51392705) q[1];
sx q[1];
rz(-0.75427873) q[1];
sx q[1];
rz(-2.8819081) q[1];
rz(0.17856715) q[2];
sx q[2];
rz(-0.44896508) q[2];
sx q[2];
rz(1.8812897) q[2];
rz(-3.1410029) q[3];
sx q[3];
rz(-1.1564217) q[3];
sx q[3];
rz(0.52903928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
