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
rz(-0.66752783) q[0];
sx q[0];
rz(-2.3025371) q[0];
sx q[0];
rz(2.7091205) q[0];
rz(0.67983812) q[1];
sx q[1];
rz(-1.4601269) q[1];
sx q[1];
rz(2.3753994) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941152) q[0];
sx q[0];
rz(-1.2060529) q[0];
sx q[0];
rz(3.0459131) q[0];
rz(-1.7384479) q[2];
sx q[2];
rz(-2.9580287) q[2];
sx q[2];
rz(0.50915424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4369504) q[1];
sx q[1];
rz(-2.5931103) q[1];
sx q[1];
rz(1.30634) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0187383) q[3];
sx q[3];
rz(-1.4335504) q[3];
sx q[3];
rz(-0.48396971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.55521518) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(-2.4797454) q[2];
rz(1.9072748) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(0.1525277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9041651) q[0];
sx q[0];
rz(-2.4894042) q[0];
sx q[0];
rz(2.7320614) q[0];
rz(0.25853363) q[1];
sx q[1];
rz(-1.2068318) q[1];
sx q[1];
rz(-0.55005598) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058592168) q[0];
sx q[0];
rz(-0.48230902) q[0];
sx q[0];
rz(1.3829718) q[0];
rz(-pi) q[1];
rz(-2.2520355) q[2];
sx q[2];
rz(-2.4585138) q[2];
sx q[2];
rz(1.6967333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0747968) q[1];
sx q[1];
rz(-0.48802146) q[1];
sx q[1];
rz(-1.1244494) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2004545) q[3];
sx q[3];
rz(-2.3562288) q[3];
sx q[3];
rz(1.0496275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6253769) q[2];
sx q[2];
rz(-1.3863401) q[2];
sx q[2];
rz(-0.59164444) q[2];
rz(0.51724616) q[3];
sx q[3];
rz(-2.0375242) q[3];
sx q[3];
rz(0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2194694) q[0];
sx q[0];
rz(-1.8485494) q[0];
sx q[0];
rz(-2.1406232) q[0];
rz(1.4283904) q[1];
sx q[1];
rz(-2.3826471) q[1];
sx q[1];
rz(-3.0975814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6881268) q[0];
sx q[0];
rz(-1.9489053) q[0];
sx q[0];
rz(-0.72909683) q[0];
rz(-0.83389177) q[2];
sx q[2];
rz(-1.6084889) q[2];
sx q[2];
rz(2.1455763) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2340369) q[1];
sx q[1];
rz(-2.2436972) q[1];
sx q[1];
rz(2.2220045) q[1];
rz(-pi) q[2];
rz(0.040626133) q[3];
sx q[3];
rz(-1.3990132) q[3];
sx q[3];
rz(0.21200997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9192146) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(3.0317793) q[2];
rz(-0.63012704) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(-0.28756791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(1.6269864) q[0];
sx q[0];
rz(-2.0872748) q[0];
sx q[0];
rz(-1.0940254) q[0];
rz(0.47359723) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(2.4213562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0648616) q[0];
sx q[0];
rz(-1.0986878) q[0];
sx q[0];
rz(-1.2272627) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8333387) q[2];
sx q[2];
rz(-1.5259907) q[2];
sx q[2];
rz(-0.055331478) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7590649) q[1];
sx q[1];
rz(-0.80611496) q[1];
sx q[1];
rz(-2.7208946) q[1];
rz(-pi) q[2];
rz(1.9071753) q[3];
sx q[3];
rz(-2.3569538) q[3];
sx q[3];
rz(1.6616846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56095162) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(-0.15216039) q[2];
rz(-1.2946607) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(1.7744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691567) q[0];
sx q[0];
rz(-0.56273383) q[0];
sx q[0];
rz(2.7602472) q[0];
rz(2.5533679) q[1];
sx q[1];
rz(-1.7344079) q[1];
sx q[1];
rz(-3.0773967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17339029) q[0];
sx q[0];
rz(-1.2009092) q[0];
sx q[0];
rz(-1.3805519) q[0];
rz(-pi) q[1];
rz(1.97528) q[2];
sx q[2];
rz(-0.99405406) q[2];
sx q[2];
rz(-2.8167644) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.209268) q[1];
sx q[1];
rz(-1.8552244) q[1];
sx q[1];
rz(-0.80764215) q[1];
x q[2];
rz(3.0194123) q[3];
sx q[3];
rz(-1.7751649) q[3];
sx q[3];
rz(0.86993566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3885145) q[2];
sx q[2];
rz(-1.3103176) q[2];
sx q[2];
rz(-0.7927967) q[2];
rz(3.0068398) q[3];
sx q[3];
rz(-2.5760791) q[3];
sx q[3];
rz(-1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7029) q[0];
sx q[0];
rz(-1.5543992) q[0];
sx q[0];
rz(0.95293522) q[0];
rz(1.8903525) q[1];
sx q[1];
rz(-1.4616936) q[1];
sx q[1];
rz(2.9073263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.983646) q[0];
sx q[0];
rz(-2.0711714) q[0];
sx q[0];
rz(-1.7364794) q[0];
x q[1];
rz(0.71990196) q[2];
sx q[2];
rz(-2.5635519) q[2];
sx q[2];
rz(-1.3225916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30429669) q[1];
sx q[1];
rz(-2.1526744) q[1];
sx q[1];
rz(3.0051662) q[1];
rz(-pi) q[2];
rz(-0.085233099) q[3];
sx q[3];
rz(-2.2036423) q[3];
sx q[3];
rz(-1.3229602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.970152) q[2];
sx q[2];
rz(-2.5415387) q[2];
sx q[2];
rz(-2.1947412) q[2];
rz(1.4417449) q[3];
sx q[3];
rz(-1.8141247) q[3];
sx q[3];
rz(-0.12137752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117208) q[0];
sx q[0];
rz(-1.0026824) q[0];
sx q[0];
rz(-0.3748931) q[0];
rz(-2.3986469) q[1];
sx q[1];
rz(-0.82287794) q[1];
sx q[1];
rz(-0.65623647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5617758) q[0];
sx q[0];
rz(-2.4120789) q[0];
sx q[0];
rz(-0.39493409) q[0];
x q[1];
rz(-2.5083582) q[2];
sx q[2];
rz(-1.0960964) q[2];
sx q[2];
rz(0.98596078) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6306047) q[1];
sx q[1];
rz(-1.9888249) q[1];
sx q[1];
rz(1.1633148) q[1];
rz(-pi) q[2];
rz(-1.3712758) q[3];
sx q[3];
rz(-2.1533416) q[3];
sx q[3];
rz(0.75967805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6132505) q[2];
sx q[2];
rz(-2.2684147) q[2];
sx q[2];
rz(-2.4916416) q[2];
rz(-2.2925099) q[3];
sx q[3];
rz(-0.63747469) q[3];
sx q[3];
rz(1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36986247) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(2.8233897) q[0];
rz(-0.33266208) q[1];
sx q[1];
rz(-2.5883364) q[1];
sx q[1];
rz(-0.37700787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57721389) q[0];
sx q[0];
rz(-2.8014392) q[0];
sx q[0];
rz(1.7094014) q[0];
x q[1];
rz(1.6416574) q[2];
sx q[2];
rz(-0.46816816) q[2];
sx q[2];
rz(2.5137794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.01025) q[1];
sx q[1];
rz(-0.90908644) q[1];
sx q[1];
rz(-1.0922088) q[1];
rz(-2.0539862) q[3];
sx q[3];
rz(-1.4677262) q[3];
sx q[3];
rz(2.7838992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0074924) q[2];
sx q[2];
rz(-1.2780122) q[2];
sx q[2];
rz(-2.210614) q[2];
rz(1.5382986) q[3];
sx q[3];
rz(-1.337991) q[3];
sx q[3];
rz(1.3163756) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304921) q[0];
sx q[0];
rz(-2.7340524) q[0];
sx q[0];
rz(-2.7000309) q[0];
rz(0.26511296) q[1];
sx q[1];
rz(-1.4418863) q[1];
sx q[1];
rz(-1.6308867) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8245659) q[0];
sx q[0];
rz(-2.6500033) q[0];
sx q[0];
rz(-0.48382171) q[0];
rz(-2.8042107) q[2];
sx q[2];
rz(-1.0058912) q[2];
sx q[2];
rz(3.1021822) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1944674) q[1];
sx q[1];
rz(-1.3381011) q[1];
sx q[1];
rz(-0.10242771) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2208414) q[3];
sx q[3];
rz(-1.2724605) q[3];
sx q[3];
rz(0.32978299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.536671) q[2];
sx q[2];
rz(-1.8646381) q[2];
sx q[2];
rz(-0.46860179) q[2];
rz(-1.0118777) q[3];
sx q[3];
rz(-1.4230655) q[3];
sx q[3];
rz(-2.0603777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6877947) q[0];
sx q[0];
rz(-0.84306222) q[0];
sx q[0];
rz(-1.7048365) q[0];
rz(-2.1935678) q[1];
sx q[1];
rz(-1.8016305) q[1];
sx q[1];
rz(-0.2737793) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1760491) q[0];
sx q[0];
rz(-1.6946486) q[0];
sx q[0];
rz(1.6424422) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5735515) q[2];
sx q[2];
rz(-0.92457891) q[2];
sx q[2];
rz(-0.60233145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5999881) q[1];
sx q[1];
rz(-2.0683824) q[1];
sx q[1];
rz(1.0367839) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.580299) q[2];
sx q[2];
rz(-1.0729921) q[2];
sx q[2];
rz(2.6936074) q[2];
rz(-0.63730803) q[3];
sx q[3];
rz(-1.0261122) q[3];
sx q[3];
rz(0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4056024) q[0];
sx q[0];
rz(-1.6163106) q[0];
sx q[0];
rz(1.6769782) q[0];
rz(2.6276656) q[1];
sx q[1];
rz(-2.3873139) q[1];
sx q[1];
rz(0.25968459) q[1];
rz(2.6988637) q[2];
sx q[2];
rz(-1.6479658) q[2];
sx q[2];
rz(-2.992291) q[2];
rz(-1.5721372) q[3];
sx q[3];
rz(-0.41437498) q[3];
sx q[3];
rz(-2.6140183) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
