OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3482398) q[0];
sx q[0];
rz(-1.5982331) q[0];
sx q[0];
rz(-1.5016851) q[0];
rz(1.7465254) q[1];
sx q[1];
rz(-2.7434064) q[1];
sx q[1];
rz(-0.15287486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3216696) q[0];
sx q[0];
rz(-1.4355441) q[0];
sx q[0];
rz(-0.090094264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42793226) q[2];
sx q[2];
rz(-0.38201354) q[2];
sx q[2];
rz(-2.0150507) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6907883) q[1];
sx q[1];
rz(-1.5105472) q[1];
sx q[1];
rz(-1.3817527) q[1];
x q[2];
rz(2.6809815) q[3];
sx q[3];
rz(-0.95129993) q[3];
sx q[3];
rz(1.7509489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0632826) q[2];
sx q[2];
rz(-1.7141043) q[2];
sx q[2];
rz(-0.59986344) q[2];
rz(-2.6317224) q[3];
sx q[3];
rz(-2.8190835) q[3];
sx q[3];
rz(2.9738284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55709368) q[0];
sx q[0];
rz(-2.2826513) q[0];
sx q[0];
rz(2.9993045) q[0];
rz(1.8498869) q[1];
sx q[1];
rz(-2.6085491) q[1];
sx q[1];
rz(-0.60089111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37764022) q[0];
sx q[0];
rz(-0.48177347) q[0];
sx q[0];
rz(-2.0773621) q[0];
x q[1];
rz(-1.1942062) q[2];
sx q[2];
rz(-0.96842679) q[2];
sx q[2];
rz(-0.9488238) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4838769) q[1];
sx q[1];
rz(-2.2536088) q[1];
sx q[1];
rz(-1.3394936) q[1];
x q[2];
rz(-1.7606335) q[3];
sx q[3];
rz(-2.077436) q[3];
sx q[3];
rz(2.02158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1945524) q[2];
sx q[2];
rz(-1.9779454) q[2];
sx q[2];
rz(-2.2770503) q[2];
rz(2.7583127) q[3];
sx q[3];
rz(-2.7195103) q[3];
sx q[3];
rz(-0.88467902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8466723) q[0];
sx q[0];
rz(-0.28691322) q[0];
sx q[0];
rz(-2.3160146) q[0];
rz(0.83746743) q[1];
sx q[1];
rz(-2.1223964) q[1];
sx q[1];
rz(0.79140633) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90368962) q[0];
sx q[0];
rz(-2.3775953) q[0];
sx q[0];
rz(1.4542411) q[0];
rz(-1.4405652) q[2];
sx q[2];
rz(-1.669906) q[2];
sx q[2];
rz(-2.2203022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0703146) q[1];
sx q[1];
rz(-2.5025939) q[1];
sx q[1];
rz(-0.60296873) q[1];
rz(-1.0867001) q[3];
sx q[3];
rz(-1.181385) q[3];
sx q[3];
rz(2.2661346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5241663) q[2];
sx q[2];
rz(-2.7238621) q[2];
sx q[2];
rz(-1.8903271) q[2];
rz(-1.1795801) q[3];
sx q[3];
rz(-2.0056632) q[3];
sx q[3];
rz(-2.484926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11611045) q[0];
sx q[0];
rz(-1.9132834) q[0];
sx q[0];
rz(0.81892282) q[0];
rz(-3.0602835) q[1];
sx q[1];
rz(-0.58650494) q[1];
sx q[1];
rz(0.78048817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9355192) q[0];
sx q[0];
rz(-1.9447826) q[0];
sx q[0];
rz(1.5078578) q[0];
x q[1];
rz(1.6967746) q[2];
sx q[2];
rz(-0.65444817) q[2];
sx q[2];
rz(2.12754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38950237) q[1];
sx q[1];
rz(-2.6679975) q[1];
sx q[1];
rz(-2.447763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1673766) q[3];
sx q[3];
rz(-1.2781004) q[3];
sx q[3];
rz(-0.044805275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9235733) q[2];
sx q[2];
rz(-1.6417445) q[2];
sx q[2];
rz(-1.8531331) q[2];
rz(-1.4675379) q[3];
sx q[3];
rz(-2.5947184) q[3];
sx q[3];
rz(-0.43681496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502289) q[0];
sx q[0];
rz(-0.65058351) q[0];
sx q[0];
rz(-0.22317602) q[0];
rz(-0.079809345) q[1];
sx q[1];
rz(-2.1640919) q[1];
sx q[1];
rz(-2.3056183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2612117) q[0];
sx q[0];
rz(-1.1376842) q[0];
sx q[0];
rz(-2.8773737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4673632) q[2];
sx q[2];
rz(-0.59665702) q[2];
sx q[2];
rz(2.3643189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0093352) q[1];
sx q[1];
rz(-2.461488) q[1];
sx q[1];
rz(-2.2538005) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9014408) q[3];
sx q[3];
rz(-2.504494) q[3];
sx q[3];
rz(1.6589056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9072546) q[2];
sx q[2];
rz(-2.2865488) q[2];
sx q[2];
rz(-1.7503395) q[2];
rz(1.2587345) q[3];
sx q[3];
rz(-1.7573059) q[3];
sx q[3];
rz(0.39906991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2788972) q[0];
sx q[0];
rz(-0.47877043) q[0];
sx q[0];
rz(0.71267772) q[0];
rz(-2.0172987) q[1];
sx q[1];
rz(-1.5713888) q[1];
sx q[1];
rz(-2.1052776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6725715) q[0];
sx q[0];
rz(-2.2037134) q[0];
sx q[0];
rz(-1.4255037) q[0];
rz(-1.2769624) q[2];
sx q[2];
rz(-1.7987408) q[2];
sx q[2];
rz(-2.0541525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6129174) q[1];
sx q[1];
rz(-2.5096748) q[1];
sx q[1];
rz(-1.9252434) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2676663) q[3];
sx q[3];
rz(-1.4269265) q[3];
sx q[3];
rz(1.914639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0816281) q[2];
sx q[2];
rz(-0.57264239) q[2];
sx q[2];
rz(2.6247978) q[2];
rz(1.9948657) q[3];
sx q[3];
rz(-0.99249339) q[3];
sx q[3];
rz(0.049588047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0674755) q[0];
sx q[0];
rz(-2.4570486) q[0];
sx q[0];
rz(1.6790947) q[0];
rz(-2.0430203) q[1];
sx q[1];
rz(-1.699828) q[1];
sx q[1];
rz(0.77902737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.532004) q[0];
sx q[0];
rz(-1.6250984) q[0];
sx q[0];
rz(1.5691343) q[0];
x q[1];
rz(-0.46709664) q[2];
sx q[2];
rz(-1.0425241) q[2];
sx q[2];
rz(1.686765) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0454899) q[1];
sx q[1];
rz(-1.742846) q[1];
sx q[1];
rz(2.9330611) q[1];
rz(-pi) q[2];
rz(-2.7610711) q[3];
sx q[3];
rz(-2.2448934) q[3];
sx q[3];
rz(0.014691513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74759185) q[2];
sx q[2];
rz(-0.70351768) q[2];
sx q[2];
rz(2.7819989) q[2];
rz(0.30073419) q[3];
sx q[3];
rz(-0.81529236) q[3];
sx q[3];
rz(2.1451779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312298) q[0];
sx q[0];
rz(-1.2082986) q[0];
sx q[0];
rz(0.30074686) q[0];
rz(2.6077479) q[1];
sx q[1];
rz(-2.0678935) q[1];
sx q[1];
rz(-0.11375443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14989242) q[0];
sx q[0];
rz(-1.0731369) q[0];
sx q[0];
rz(2.0927621) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0871307) q[2];
sx q[2];
rz(-2.5321333) q[2];
sx q[2];
rz(3.0888588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2054322) q[1];
sx q[1];
rz(-2.0179085) q[1];
sx q[1];
rz(-0.7128678) q[1];
rz(1.9756687) q[3];
sx q[3];
rz(-0.40388122) q[3];
sx q[3];
rz(0.98458588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45812312) q[2];
sx q[2];
rz(-0.14294954) q[2];
sx q[2];
rz(-0.45041034) q[2];
rz(-2.2266375) q[3];
sx q[3];
rz(-2.2136642) q[3];
sx q[3];
rz(-1.3373059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8704855) q[0];
sx q[0];
rz(-0.73942375) q[0];
sx q[0];
rz(-0.41573218) q[0];
rz(0.43139002) q[1];
sx q[1];
rz(-0.58512551) q[1];
sx q[1];
rz(-2.364667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5901075) q[0];
sx q[0];
rz(-0.05478207) q[0];
sx q[0];
rz(-1.3339186) q[0];
rz(-pi) q[1];
rz(0.84085744) q[2];
sx q[2];
rz(-2.1321725) q[2];
sx q[2];
rz(-0.64602588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7207234) q[1];
sx q[1];
rz(-0.72378473) q[1];
sx q[1];
rz(2.4813244) q[1];
rz(-pi) q[2];
rz(0.52962065) q[3];
sx q[3];
rz(-0.52859113) q[3];
sx q[3];
rz(1.9749255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0100157) q[2];
sx q[2];
rz(-0.78993979) q[2];
sx q[2];
rz(1.5124403) q[2];
rz(-2.7029964) q[3];
sx q[3];
rz(-0.83493835) q[3];
sx q[3];
rz(-2.6252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55385357) q[0];
sx q[0];
rz(-2.6022311) q[0];
sx q[0];
rz(1.3158276) q[0];
rz(-3.059803) q[1];
sx q[1];
rz(-1.9606699) q[1];
sx q[1];
rz(-1.2710424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2270541) q[0];
sx q[0];
rz(-2.7167121) q[0];
sx q[0];
rz(1.4999313) q[0];
rz(-pi) q[1];
rz(2.3009335) q[2];
sx q[2];
rz(-1.5224384) q[2];
sx q[2];
rz(-2.4052467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0942022) q[1];
sx q[1];
rz(-1.2823815) q[1];
sx q[1];
rz(-0.80077313) q[1];
x q[2];
rz(0.081689452) q[3];
sx q[3];
rz(-1.5458071) q[3];
sx q[3];
rz(0.38988926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79909331) q[2];
sx q[2];
rz(-1.7010331) q[2];
sx q[2];
rz(2.0889757) q[2];
rz(-1.3275702) q[3];
sx q[3];
rz(-1.9460287) q[3];
sx q[3];
rz(2.6322406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159745) q[0];
sx q[0];
rz(-1.3798036) q[0];
sx q[0];
rz(1.8984541) q[0];
rz(1.3432518) q[1];
sx q[1];
rz(-2.6814798) q[1];
sx q[1];
rz(0.36569256) q[1];
rz(1.0940362) q[2];
sx q[2];
rz(-1.7443716) q[2];
sx q[2];
rz(-2.3732408) q[2];
rz(3.1111191) q[3];
sx q[3];
rz(-2.7298374) q[3];
sx q[3];
rz(0.8047837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
