OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711174) q[0];
sx q[0];
rz(-1.820571) q[0];
sx q[0];
rz(-0.063637861) q[0];
x q[1];
rz(-1.9185669) q[2];
sx q[2];
rz(-0.43453056) q[2];
sx q[2];
rz(2.6968616) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.3738457) q[1];
sx q[1];
rz(-2.0031553) q[1];
sx q[1];
rz(0.14111622) q[1];
rz(-pi) q[2];
rz(0.31580117) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(-0.11581126) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(0.23981747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46562425) q[0];
sx q[0];
rz(-1.3249319) q[0];
sx q[0];
rz(-1.4093536) q[0];
rz(2.9302836) q[2];
sx q[2];
rz(-1.276187) q[2];
sx q[2];
rz(-0.78795563) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51057306) q[1];
sx q[1];
rz(-2.3595855) q[1];
sx q[1];
rz(-1.4960947) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2858743) q[3];
sx q[3];
rz(-0.021589605) q[3];
sx q[3];
rz(-1.3947226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.6372765) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2770237) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0720172) q[0];
sx q[0];
rz(-1.5298651) q[0];
sx q[0];
rz(-1.3224365) q[0];
x q[1];
rz(-0.5553603) q[2];
sx q[2];
rz(-1.4484324) q[2];
sx q[2];
rz(2.8107779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2324893) q[1];
sx q[1];
rz(-2.4965582) q[1];
sx q[1];
rz(2.049936) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90956456) q[3];
sx q[3];
rz(-1.665691) q[3];
sx q[3];
rz(-2.6385917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(3.1090453) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(-1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.6548086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0006205) q[0];
sx q[0];
rz(-1.9019433) q[0];
sx q[0];
rz(-1.8676057) q[0];
x q[1];
rz(-0.60947946) q[2];
sx q[2];
rz(-1.8595427) q[2];
sx q[2];
rz(2.2878873) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1292124) q[1];
sx q[1];
rz(-2.817569) q[1];
sx q[1];
rz(-2.654241) q[1];
rz(-pi) q[2];
rz(3.1348226) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(-3.0778459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(-2.0452943) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0356045) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.460357) q[0];
rz(-1.5530855) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-0.016074093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8668629) q[0];
sx q[0];
rz(-2.1511973) q[0];
sx q[0];
rz(-2.4972563) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1400914) q[2];
sx q[2];
rz(-1.5632196) q[2];
sx q[2];
rz(-1.8097144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9863723) q[1];
sx q[1];
rz(-1.1481789) q[1];
sx q[1];
rz(-0.95814725) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1297242) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(0.12479898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(-1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(-0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-2.2348256) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(1.9168568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8794494) q[0];
sx q[0];
rz(-2.1855542) q[0];
sx q[0];
rz(-0.26457796) q[0];
rz(-pi) q[1];
rz(1.9903509) q[2];
sx q[2];
rz(-2.662979) q[2];
sx q[2];
rz(-1.7350369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3983706) q[1];
sx q[1];
rz(-0.9658769) q[1];
sx q[1];
rz(0.066992316) q[1];
x q[2];
rz(1.46035) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(-0.16050592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619103) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(0.58037037) q[0];
rz(-1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-2.4408128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.128292) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(1.1178455) q[0];
rz(-0.3406616) q[2];
sx q[2];
rz(-0.33764687) q[2];
sx q[2];
rz(1.1866736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.339387) q[1];
sx q[1];
rz(-1.4805668) q[1];
sx q[1];
rz(2.1369364) q[1];
rz(0.84007646) q[3];
sx q[3];
rz(-1.6335765) q[3];
sx q[3];
rz(-1.9907794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(2.3968905) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78684029) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(-1.4165075) q[0];
rz(-1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69670024) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-2.9219887) q[0];
rz(0.58015577) q[2];
sx q[2];
rz(-1.1859425) q[2];
sx q[2];
rz(0.58154026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2878694) q[1];
sx q[1];
rz(-0.81765491) q[1];
sx q[1];
rz(0.57662782) q[1];
rz(-pi) q[2];
rz(1.1710839) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(-1.8613929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.7766215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7576335) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(-0.061116771) q[0];
rz(-pi) q[1];
rz(1.1599837) q[2];
sx q[2];
rz(-1.64752) q[2];
sx q[2];
rz(3.0551747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2024467) q[1];
sx q[1];
rz(-0.89679634) q[1];
sx q[1];
rz(2.0161611) q[1];
rz(-pi) q[2];
rz(1.8678719) q[3];
sx q[3];
rz(-1.3798514) q[3];
sx q[3];
rz(1.8553268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.6636794) q[0];
rz(-1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(0.96819425) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68010846) q[0];
sx q[0];
rz(-2.619954) q[0];
sx q[0];
rz(2.3011544) q[0];
rz(-pi) q[1];
rz(0.88629006) q[2];
sx q[2];
rz(-0.61294014) q[2];
sx q[2];
rz(2.431589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8755175) q[1];
sx q[1];
rz(-0.35421687) q[1];
sx q[1];
rz(-0.43550272) q[1];
x q[2];
rz(-1.2486542) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-2.145888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.39682) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(-0.36623476) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(-2.3920849) q[2];
sx q[2];
rz(-1.8029677) q[2];
sx q[2];
rz(0.36063902) q[2];
rz(-2.1063741) q[3];
sx q[3];
rz(-2.513701) q[3];
sx q[3];
rz(2.7660478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
