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
rz(-2.1751997) q[0];
sx q[0];
rz(-1.3858495) q[0];
sx q[0];
rz(0.59544271) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(2.5505677) q[1];
sx q[1];
rz(12.443065) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57014556) q[0];
sx q[0];
rz(-2.4295904) q[0];
sx q[0];
rz(-0.90306905) q[0];
x q[1];
rz(-2.3507471) q[2];
sx q[2];
rz(-0.95480761) q[2];
sx q[2];
rz(-3.0012584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5414303) q[1];
sx q[1];
rz(-1.3568391) q[1];
sx q[1];
rz(-3.0453277) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7233347) q[3];
sx q[3];
rz(-1.9328874) q[3];
sx q[3];
rz(-2.3063776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82723242) q[2];
sx q[2];
rz(-1.102697) q[2];
sx q[2];
rz(3.0296791) q[2];
rz(1.3917475) q[3];
sx q[3];
rz(-1.1575969) q[3];
sx q[3];
rz(-0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7258485) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(-0.14990212) q[0];
rz(0.28597486) q[1];
sx q[1];
rz(-1.3949225) q[1];
sx q[1];
rz(-2.012595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647346) q[0];
sx q[0];
rz(-1.6010512) q[0];
sx q[0];
rz(-1.2824351) q[0];
x q[1];
rz(-1.6304134) q[2];
sx q[2];
rz(-1.5562415) q[2];
sx q[2];
rz(0.83257914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4399229) q[1];
sx q[1];
rz(-1.7365125) q[1];
sx q[1];
rz(-2.7782337) q[1];
rz(0.942248) q[3];
sx q[3];
rz(-1.0744922) q[3];
sx q[3];
rz(-1.4013578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67843208) q[2];
sx q[2];
rz(-0.9223991) q[2];
sx q[2];
rz(1.9909667) q[2];
rz(1.9567418) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(-0.38159698) q[0];
rz(1.8474139) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(-0.19827422) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564045) q[0];
sx q[0];
rz(-1.3964147) q[0];
sx q[0];
rz(1.3402433) q[0];
x q[1];
rz(1.6724104) q[2];
sx q[2];
rz(-2.4696015) q[2];
sx q[2];
rz(-0.72173126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4130062) q[1];
sx q[1];
rz(-0.31794128) q[1];
sx q[1];
rz(-2.4012662) q[1];
x q[2];
rz(-1.6842849) q[3];
sx q[3];
rz(-2.9656565) q[3];
sx q[3];
rz(0.709155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76564378) q[2];
sx q[2];
rz(-0.15268923) q[2];
sx q[2];
rz(-2.622733) q[2];
rz(-1.2505924) q[3];
sx q[3];
rz(-2.0542681) q[3];
sx q[3];
rz(-0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64312235) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(0.44922391) q[0];
rz(-1.4462224) q[1];
sx q[1];
rz(-2.4999373) q[1];
sx q[1];
rz(-2.5208688) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45263824) q[0];
sx q[0];
rz(-1.7151378) q[0];
sx q[0];
rz(-2.4792433) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5081257) q[2];
sx q[2];
rz(-1.4711498) q[2];
sx q[2];
rz(-1.5657305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52445946) q[1];
sx q[1];
rz(-1.7726328) q[1];
sx q[1];
rz(2.733888) q[1];
rz(-pi) q[2];
rz(-1.4820547) q[3];
sx q[3];
rz(-2.6013298) q[3];
sx q[3];
rz(-1.0080573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0351403) q[2];
sx q[2];
rz(-1.8702714) q[2];
sx q[2];
rz(2.980496) q[2];
rz(-0.39786878) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(2.9919992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905269) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(-2.8845442) q[0];
rz(-0.88047475) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(-1.2535198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7505676) q[0];
sx q[0];
rz(-1.5526958) q[0];
sx q[0];
rz(1.6068293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3257371) q[2];
sx q[2];
rz(-0.7759717) q[2];
sx q[2];
rz(1.4910335) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1977928) q[1];
sx q[1];
rz(-1.5256923) q[1];
sx q[1];
rz(0.79774858) q[1];
rz(0.34714602) q[3];
sx q[3];
rz(-1.5691363) q[3];
sx q[3];
rz(-1.1356789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72838655) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(1.2665292) q[2];
rz(-2.5201216) q[3];
sx q[3];
rz(-2.5340243) q[3];
sx q[3];
rz(-0.5411886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40015873) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(-3.1165282) q[0];
rz(1.9301682) q[1];
sx q[1];
rz(-1.2592659) q[1];
sx q[1];
rz(2.8866344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0023277442) q[0];
sx q[0];
rz(-0.049590913) q[0];
sx q[0];
rz(0.10164405) q[0];
rz(-pi) q[1];
rz(-2.4676085) q[2];
sx q[2];
rz(-0.70491284) q[2];
sx q[2];
rz(-0.92943905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7595664) q[1];
sx q[1];
rz(-0.81795482) q[1];
sx q[1];
rz(-0.76063971) q[1];
x q[2];
rz(1.2473769) q[3];
sx q[3];
rz(-2.5445523) q[3];
sx q[3];
rz(-3.0003594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8956464) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(0.45905217) q[2];
rz(1.4970655) q[3];
sx q[3];
rz(-3.0247757) q[3];
sx q[3];
rz(-3.0204401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46886214) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(-2.22279) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(-1.930621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0491517) q[0];
sx q[0];
rz(-1.3685797) q[0];
sx q[0];
rz(-2.8114955) q[0];
x q[1];
rz(1.0031021) q[2];
sx q[2];
rz(-1.9493812) q[2];
sx q[2];
rz(-2.1673983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1001491) q[1];
sx q[1];
rz(-2.7584834) q[1];
sx q[1];
rz(-2.1418397) q[1];
x q[2];
rz(1.2065229) q[3];
sx q[3];
rz(-1.9134054) q[3];
sx q[3];
rz(0.87876696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81898895) q[2];
sx q[2];
rz(-0.49590597) q[2];
sx q[2];
rz(0.71713478) q[2];
rz(-2.1142193) q[3];
sx q[3];
rz(-1.4077978) q[3];
sx q[3];
rz(-0.43924847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1502007) q[0];
sx q[0];
rz(-0.99656492) q[0];
sx q[0];
rz(0.014658654) q[0];
rz(2.390059) q[1];
sx q[1];
rz(-1.2208168) q[1];
sx q[1];
rz(-1.6709447) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5293619) q[0];
sx q[0];
rz(-2.9069773) q[0];
sx q[0];
rz(-2.505872) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25663767) q[2];
sx q[2];
rz(-1.6483288) q[2];
sx q[2];
rz(-2.3687349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60622207) q[1];
sx q[1];
rz(-2.722123) q[1];
sx q[1];
rz(1.8505881) q[1];
rz(-pi) q[2];
rz(-1.2257158) q[3];
sx q[3];
rz(-1.5672641) q[3];
sx q[3];
rz(-2.6659903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8810001) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(-0.27349681) q[2];
rz(1.1547487) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1028033) q[0];
sx q[0];
rz(-2.0162855) q[0];
sx q[0];
rz(-1.0754841) q[0];
rz(-0.12818809) q[1];
sx q[1];
rz(-2.2905541) q[1];
sx q[1];
rz(-2.7959965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49508383) q[0];
sx q[0];
rz(-2.2301607) q[0];
sx q[0];
rz(2.93883) q[0];
x q[1];
rz(2.8198411) q[2];
sx q[2];
rz(-2.1158127) q[2];
sx q[2];
rz(-2.7779752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4491475) q[1];
sx q[1];
rz(-1.0303823) q[1];
sx q[1];
rz(1.497333) q[1];
x q[2];
rz(1.4050499) q[3];
sx q[3];
rz(-2.1624544) q[3];
sx q[3];
rz(1.6935284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(1.5209939) q[2];
rz(-1.3711551) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(2.7267406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3138251) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-0.23571043) q[0];
rz(2.56855) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(-1.3689573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9873857) q[0];
sx q[0];
rz(-1.9999749) q[0];
sx q[0];
rz(2.3230419) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70765453) q[2];
sx q[2];
rz(-0.81205149) q[2];
sx q[2];
rz(2.6132513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9006834) q[1];
sx q[1];
rz(-2.5473809) q[1];
sx q[1];
rz(1.2628984) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24457358) q[3];
sx q[3];
rz(-0.37062708) q[3];
sx q[3];
rz(1.2115492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77091757) q[2];
sx q[2];
rz(-1.8149899) q[2];
sx q[2];
rz(0.053000432) q[2];
rz(-2.3897589) q[3];
sx q[3];
rz(-2.1078608) q[3];
sx q[3];
rz(2.2161765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5175405) q[0];
sx q[0];
rz(-2.1255827) q[0];
sx q[0];
rz(-0.95638635) q[0];
rz(1.3507631) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(-1.440329) q[2];
sx q[2];
rz(-0.87776504) q[2];
sx q[2];
rz(-0.19765111) q[2];
rz(0.57784537) q[3];
sx q[3];
rz(-2.2457849) q[3];
sx q[3];
rz(-0.59413322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
