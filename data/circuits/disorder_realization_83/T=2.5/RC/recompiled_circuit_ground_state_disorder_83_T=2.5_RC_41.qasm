OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72252005) q[0];
sx q[0];
rz(2.4604586) q[0];
sx q[0];
rz(13.606986) q[0];
rz(2.430727) q[1];
sx q[1];
rz(-2.4912973) q[1];
sx q[1];
rz(1.8898036) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1052763) q[0];
sx q[0];
rz(-1.9364089) q[0];
sx q[0];
rz(-2.7514639) q[0];
x q[1];
rz(-2.5898143) q[2];
sx q[2];
rz(-1.1349196) q[2];
sx q[2];
rz(0.18407719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7573479) q[1];
sx q[1];
rz(-0.66435821) q[1];
sx q[1];
rz(0.039907736) q[1];
x q[2];
rz(-2.4667418) q[3];
sx q[3];
rz(-2.1938257) q[3];
sx q[3];
rz(2.2365776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17994443) q[2];
sx q[2];
rz(-1.0893931) q[2];
sx q[2];
rz(-1.7845478) q[2];
rz(2.0876136) q[3];
sx q[3];
rz(-1.6097924) q[3];
sx q[3];
rz(1.0668782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2422975) q[0];
sx q[0];
rz(-2.5488148) q[0];
sx q[0];
rz(-1.5574667) q[0];
rz(1.8862995) q[1];
sx q[1];
rz(-2.2785432) q[1];
sx q[1];
rz(-0.042162808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164207) q[0];
sx q[0];
rz(-0.76043425) q[0];
sx q[0];
rz(-2.6359933) q[0];
x q[1];
rz(1.043528) q[2];
sx q[2];
rz(-2.0828491) q[2];
sx q[2];
rz(1.5917503) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0709585) q[1];
sx q[1];
rz(-1.8629304) q[1];
sx q[1];
rz(-0.16865428) q[1];
rz(-0.20037074) q[3];
sx q[3];
rz(-1.0215825) q[3];
sx q[3];
rz(0.39796838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92370799) q[2];
sx q[2];
rz(-1.2486685) q[2];
sx q[2];
rz(-2.5200747) q[2];
rz(0.38226852) q[3];
sx q[3];
rz(-1.1416124) q[3];
sx q[3];
rz(-0.83466616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37654787) q[0];
sx q[0];
rz(-1.6355729) q[0];
sx q[0];
rz(-1.7682834) q[0];
rz(0.37257591) q[1];
sx q[1];
rz(-0.5245477) q[1];
sx q[1];
rz(0.12038825) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62154462) q[0];
sx q[0];
rz(-0.95764521) q[0];
sx q[0];
rz(0.91213062) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4311426) q[2];
sx q[2];
rz(-0.28208986) q[2];
sx q[2];
rz(0.26821801) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87625831) q[1];
sx q[1];
rz(-1.2560532) q[1];
sx q[1];
rz(-2.3708484) q[1];
x q[2];
rz(-2.8684738) q[3];
sx q[3];
rz(-2.0865174) q[3];
sx q[3];
rz(-1.7777594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4176214) q[2];
sx q[2];
rz(-1.3105023) q[2];
sx q[2];
rz(-1.2869147) q[2];
rz(-2.3567965) q[3];
sx q[3];
rz(-2.0311821) q[3];
sx q[3];
rz(-0.77814656) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50412905) q[0];
sx q[0];
rz(-0.4937208) q[0];
sx q[0];
rz(-0.24566393) q[0];
rz(-3.1371112) q[1];
sx q[1];
rz(-2.5878398) q[1];
sx q[1];
rz(-0.063668879) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0086567) q[0];
sx q[0];
rz(-1.3244761) q[0];
sx q[0];
rz(2.4755347) q[0];
rz(2.8174899) q[2];
sx q[2];
rz(-1.1253005) q[2];
sx q[2];
rz(-0.53037077) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50925991) q[1];
sx q[1];
rz(-2.2508988) q[1];
sx q[1];
rz(-0.24056178) q[1];
rz(-pi) q[2];
rz(-1.7710835) q[3];
sx q[3];
rz(-1.1326579) q[3];
sx q[3];
rz(2.6753257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29023805) q[2];
sx q[2];
rz(-1.3685702) q[2];
sx q[2];
rz(1.2151388) q[2];
rz(-2.3859207) q[3];
sx q[3];
rz(-2.1561421) q[3];
sx q[3];
rz(-2.9198666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82472411) q[0];
sx q[0];
rz(-1.2687954) q[0];
sx q[0];
rz(-0.49279898) q[0];
rz(2.1700676) q[1];
sx q[1];
rz(-1.3243472) q[1];
sx q[1];
rz(0.18878254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.550866) q[0];
sx q[0];
rz(-1.8456689) q[0];
sx q[0];
rz(-1.7455533) q[0];
rz(-0.014585693) q[2];
sx q[2];
rz(-1.5684394) q[2];
sx q[2];
rz(-0.49864008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3535315) q[1];
sx q[1];
rz(-2.6279098) q[1];
sx q[1];
rz(-1.4151001) q[1];
rz(-0.32922642) q[3];
sx q[3];
rz(-2.4705973) q[3];
sx q[3];
rz(1.383757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8652953) q[2];
sx q[2];
rz(-1.6388288) q[2];
sx q[2];
rz(0.067528188) q[2];
rz(-1.6826132) q[3];
sx q[3];
rz(-2.429481) q[3];
sx q[3];
rz(2.0290831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825298) q[0];
sx q[0];
rz(-1.9628061) q[0];
sx q[0];
rz(-1.4558526) q[0];
rz(0.23446941) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(0.04105982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9735255) q[0];
sx q[0];
rz(-0.70190185) q[0];
sx q[0];
rz(-3.0977416) q[0];
rz(-1.3439937) q[2];
sx q[2];
rz(-2.5580346) q[2];
sx q[2];
rz(2.9160519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31672656) q[1];
sx q[1];
rz(-1.1371964) q[1];
sx q[1];
rz(1.0528013) q[1];
rz(0.84223522) q[3];
sx q[3];
rz(-1.6819309) q[3];
sx q[3];
rz(-0.59167093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91512758) q[2];
sx q[2];
rz(-2.3751986) q[2];
sx q[2];
rz(2.8506632) q[2];
rz(-0.21864024) q[3];
sx q[3];
rz(-0.69059697) q[3];
sx q[3];
rz(0.87029988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17756322) q[0];
sx q[0];
rz(-0.95528269) q[0];
sx q[0];
rz(2.5780504) q[0];
rz(0.82849416) q[1];
sx q[1];
rz(-0.68873134) q[1];
sx q[1];
rz(0.73961893) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1391382) q[0];
sx q[0];
rz(-2.1147303) q[0];
sx q[0];
rz(-0.55940658) q[0];
rz(2.0881485) q[2];
sx q[2];
rz(-1.7202677) q[2];
sx q[2];
rz(-2.8959413) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1355734) q[1];
sx q[1];
rz(-1.83653) q[1];
sx q[1];
rz(-2.0905079) q[1];
rz(-pi) q[2];
rz(-2.9097203) q[3];
sx q[3];
rz(-0.64157569) q[3];
sx q[3];
rz(0.48504638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9699817) q[2];
sx q[2];
rz(-0.7526528) q[2];
sx q[2];
rz(0.29772154) q[2];
rz(-3.1155078) q[3];
sx q[3];
rz(-1.1931158) q[3];
sx q[3];
rz(-0.74371663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.96935) q[0];
sx q[0];
rz(-1.5650711) q[0];
sx q[0];
rz(-1.7283537) q[0];
rz(-2.6353432) q[1];
sx q[1];
rz(-2.1126316) q[1];
sx q[1];
rz(-1.9821573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0191553) q[0];
sx q[0];
rz(-1.2673237) q[0];
sx q[0];
rz(1.6772683) q[0];
rz(0.069483453) q[2];
sx q[2];
rz(-2.5589856) q[2];
sx q[2];
rz(-2.9779788) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.373542) q[1];
sx q[1];
rz(-2.6061879) q[1];
sx q[1];
rz(2.6635976) q[1];
rz(-0.13661249) q[3];
sx q[3];
rz(-2.3834565) q[3];
sx q[3];
rz(1.5026758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96782288) q[2];
sx q[2];
rz(-0.92542595) q[2];
sx q[2];
rz(-1.3646763) q[2];
rz(2.3260498) q[3];
sx q[3];
rz(-1.848685) q[3];
sx q[3];
rz(1.3968141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5705465) q[0];
sx q[0];
rz(-2.1111574) q[0];
sx q[0];
rz(1.9739157) q[0];
rz(-2.9171464) q[1];
sx q[1];
rz(-2.3851604) q[1];
sx q[1];
rz(2.6969553) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015405795) q[0];
sx q[0];
rz(-3.0660136) q[0];
sx q[0];
rz(2.3919657) q[0];
x q[1];
rz(2.3063956) q[2];
sx q[2];
rz(-0.85934256) q[2];
sx q[2];
rz(0.85853133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4363775) q[1];
sx q[1];
rz(-1.4333998) q[1];
sx q[1];
rz(-1.2295396) q[1];
x q[2];
rz(-0.44944758) q[3];
sx q[3];
rz(-1.6995113) q[3];
sx q[3];
rz(-1.4851324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.726752) q[2];
sx q[2];
rz(-1.8291992) q[2];
sx q[2];
rz(-1.4078338) q[2];
rz(-1.5358745) q[3];
sx q[3];
rz(-1.1955669) q[3];
sx q[3];
rz(-2.2685952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6294412) q[0];
sx q[0];
rz(-0.81708556) q[0];
sx q[0];
rz(1.0907115) q[0];
rz(2.0953983) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(0.88669056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2883023) q[0];
sx q[0];
rz(-1.9294943) q[0];
sx q[0];
rz(-0.38889287) q[0];
x q[1];
rz(-2.1616814) q[2];
sx q[2];
rz(-1.8136029) q[2];
sx q[2];
rz(1.3473912) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1613933) q[1];
sx q[1];
rz(-2.7100013) q[1];
sx q[1];
rz(0.031567911) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7575023) q[3];
sx q[3];
rz(-2.2010097) q[3];
sx q[3];
rz(-1.5648791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0118759) q[2];
sx q[2];
rz(-2.1265714) q[2];
sx q[2];
rz(1.2474308) q[2];
rz(-0.85793197) q[3];
sx q[3];
rz(-1.4638289) q[3];
sx q[3];
rz(-1.5990545) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7958551) q[0];
sx q[0];
rz(-1.7248187) q[0];
sx q[0];
rz(3.0378573) q[0];
rz(-0.70032447) q[1];
sx q[1];
rz(-1.2979957) q[1];
sx q[1];
rz(-1.8630984) q[1];
rz(-0.55194912) q[2];
sx q[2];
rz(-2.911534) q[2];
sx q[2];
rz(-0.96155675) q[2];
rz(-0.92209254) q[3];
sx q[3];
rz(-2.5917883) q[3];
sx q[3];
rz(-2.5592309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
