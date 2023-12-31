OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(3.119757) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(4.6012576) q[1];
sx q[1];
rz(9.6396946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41479933) q[0];
sx q[0];
rz(-0.3814632) q[0];
sx q[0];
rz(2.3261855) q[0];
rz(-pi) q[1];
rz(2.6724042) q[2];
sx q[2];
rz(-1.708963) q[2];
sx q[2];
rz(2.5772622) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3494527) q[1];
sx q[1];
rz(-0.98590241) q[1];
sx q[1];
rz(-0.5485512) q[1];
rz(-pi) q[2];
rz(1.0825726) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(-1.9762951) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78643909) q[0];
sx q[0];
rz(-2.4053898) q[0];
sx q[0];
rz(1.8650706) q[0];
rz(-0.16883822) q[2];
sx q[2];
rz(-1.642792) q[2];
sx q[2];
rz(0.26693401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0577382) q[1];
sx q[1];
rz(-1.4469622) q[1];
sx q[1];
rz(-0.79353516) q[1];
rz(-pi) q[2];
rz(-2.3223022) q[3];
sx q[3];
rz(-1.8835526) q[3];
sx q[3];
rz(0.4526588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4774287) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(2.7880923) q[0];
rz(-pi) q[1];
rz(-0.25031309) q[2];
sx q[2];
rz(-1.7134943) q[2];
sx q[2];
rz(2.755969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0355465) q[1];
sx q[1];
rz(-1.7340845) q[1];
sx q[1];
rz(-2.0590904) q[1];
rz(2.7828091) q[3];
sx q[3];
rz(-2.3453418) q[3];
sx q[3];
rz(2.1735454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-2.8386774) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-2.8767169) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6661975) q[0];
sx q[0];
rz(-0.65984939) q[0];
sx q[0];
rz(-0.049113627) q[0];
rz(2.3189544) q[2];
sx q[2];
rz(-0.3974786) q[2];
sx q[2];
rz(-0.17695225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-2.2717443) q[1];
sx q[1];
rz(-2.1125395) q[1];
rz(-0.78062765) q[3];
sx q[3];
rz(-2.160191) q[3];
sx q[3];
rz(1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(0.95265257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081776) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(1.4334701) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2383934) q[2];
sx q[2];
rz(-1.9661511) q[2];
sx q[2];
rz(-2.7796641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2561803) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(-2.1229565) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5289375) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(-0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(2.9690202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694002) q[0];
sx q[0];
rz(-1.5844371) q[0];
sx q[0];
rz(1.9872679) q[0];
rz(1.1144981) q[2];
sx q[2];
rz(-1.3367532) q[2];
sx q[2];
rz(2.8161088) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3249358) q[1];
sx q[1];
rz(-1.6759733) q[1];
sx q[1];
rz(1.1549969) q[1];
rz(-2.5112126) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(0.58296766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(-0.75659928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3698759) q[0];
sx q[0];
rz(-1.4780095) q[0];
sx q[0];
rz(1.4183527) q[0];
rz(-pi) q[1];
rz(-1.7724178) q[2];
sx q[2];
rz(-1.9747509) q[2];
sx q[2];
rz(-2.902365) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3488256) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(-1.8030333) q[1];
rz(-pi) q[2];
rz(0.19161253) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0044272) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(0.19443092) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(-2.1527122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37799997) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(2.2647122) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16915377) q[2];
sx q[2];
rz(-3.0468468) q[2];
sx q[2];
rz(2.6731238) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58395308) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(0.072903452) q[1];
rz(0.51560651) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(0.069375667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(-2.2176946) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(2.2682155) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0926889) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(-1.2094686) q[0];
rz(-pi) q[1];
rz(2.3960605) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(1.0996639) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59233353) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(1.9418282) q[1];
rz(-2.7018413) q[3];
sx q[3];
rz(-2.1034749) q[3];
sx q[3];
rz(2.7524878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(-2.7673289) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.3964765) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6122702) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(-0.32353185) q[0];
rz(-pi) q[1];
x q[1];
rz(1.868532) q[2];
sx q[2];
rz(-2.2072788) q[2];
sx q[2];
rz(2.0993078) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64618387) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(-2.3133548) q[1];
x q[2];
rz(-0.26225984) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-0.15979016) q[2];
rz(0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044534279) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(-0.13327577) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(0.940154) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(1.521048) q[3];
sx q[3];
rz(-2.1038901) q[3];
sx q[3];
rz(-0.62705561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
