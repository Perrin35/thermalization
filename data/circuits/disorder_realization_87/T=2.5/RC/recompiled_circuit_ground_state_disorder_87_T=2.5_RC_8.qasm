OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(-0.41271451) q[0];
sx q[0];
rz(-2.3053919) q[0];
rz(0.77464453) q[1];
sx q[1];
rz(-0.062008468) q[1];
sx q[1];
rz(-1.0083415) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22577408) q[0];
sx q[0];
rz(-1.4918707) q[0];
sx q[0];
rz(1.945334) q[0];
x q[1];
rz(1.8062334) q[2];
sx q[2];
rz(-1.9411692) q[2];
sx q[2];
rz(-2.3595904) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0550067) q[1];
sx q[1];
rz(-1.0957967) q[1];
sx q[1];
rz(1.7473298) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4654245) q[3];
sx q[3];
rz(-2.947315) q[3];
sx q[3];
rz(-2.7716178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8334373) q[2];
sx q[2];
rz(-0.94352949) q[2];
sx q[2];
rz(2.4885139) q[2];
rz(-0.11450442) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(-1.0623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36650518) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9089841) q[0];
sx q[0];
rz(-1.2261008) q[0];
sx q[0];
rz(-0.20362754) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1158022) q[2];
sx q[2];
rz(-2.4206851) q[2];
sx q[2];
rz(1.4269331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.008604) q[1];
sx q[1];
rz(-2.3594366) q[1];
sx q[1];
rz(2.6166366) q[1];
rz(-pi) q[2];
rz(2.1866131) q[3];
sx q[3];
rz(-1.338025) q[3];
sx q[3];
rz(-1.7580838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6985942) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(0.2300187) q[2];
rz(1.5864141) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(-2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25552937) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(-2.0607167) q[0];
rz(0.64322645) q[1];
sx q[1];
rz(-2.3422362) q[1];
sx q[1];
rz(-3.0562775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9832343) q[0];
sx q[0];
rz(-2.2081349) q[0];
sx q[0];
rz(-2.4054224) q[0];
rz(-pi) q[1];
rz(0.32344748) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(1.310854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6416723) q[1];
sx q[1];
rz(-2.9292078) q[1];
sx q[1];
rz(2.4714064) q[1];
x q[2];
rz(-3.078104) q[3];
sx q[3];
rz(-0.62219071) q[3];
sx q[3];
rz(2.6486462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2049415) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(0.19634518) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(-1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.008217) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(-0.38544449) q[0];
rz(-1.7281945) q[1];
sx q[1];
rz(-1.2047267) q[1];
sx q[1];
rz(-2.8916496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37873822) q[0];
sx q[0];
rz(-2.1702499) q[0];
sx q[0];
rz(-2.6361639) q[0];
rz(-pi) q[1];
rz(0.087935178) q[2];
sx q[2];
rz(-1.6307978) q[2];
sx q[2];
rz(1.1052002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39482597) q[1];
sx q[1];
rz(-1.4020846) q[1];
sx q[1];
rz(-0.49612237) q[1];
x q[2];
rz(1.3718666) q[3];
sx q[3];
rz(-0.78708157) q[3];
sx q[3];
rz(0.50686344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6803153) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(0.2363905) q[2];
rz(-1.8810898) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(-2.4642956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054955) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(-0.48511037) q[0];
rz(-1.1803892) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1331951) q[0];
sx q[0];
rz(-2.6566165) q[0];
sx q[0];
rz(-2.4147245) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31827688) q[2];
sx q[2];
rz(-1.9564544) q[2];
sx q[2];
rz(1.9463271) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.066557601) q[1];
sx q[1];
rz(-1.6187985) q[1];
sx q[1];
rz(-2.2407299) q[1];
x q[2];
rz(-1.5077293) q[3];
sx q[3];
rz(-1.9708037) q[3];
sx q[3];
rz(1.1253175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89852077) q[2];
sx q[2];
rz(-2.409755) q[2];
sx q[2];
rz(2.3720429) q[2];
rz(2.2122993) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(-0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1997851) q[0];
sx q[0];
rz(-2.6564044) q[0];
sx q[0];
rz(-0.76882452) q[0];
rz(1.7898412) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(2.2519055) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2425491) q[0];
sx q[0];
rz(-2.7022323) q[0];
sx q[0];
rz(-1.0602632) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5264411) q[2];
sx q[2];
rz(-1.5694478) q[2];
sx q[2];
rz(1.9951472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8823679) q[1];
sx q[1];
rz(-2.9881757) q[1];
sx q[1];
rz(1.8795561) q[1];
x q[2];
rz(-2.9075648) q[3];
sx q[3];
rz(-2.6796973) q[3];
sx q[3];
rz(1.6017101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49030226) q[2];
sx q[2];
rz(-1.4375261) q[2];
sx q[2];
rz(1.5550522) q[2];
rz(1.8709315) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(-3.0095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9689869) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(1.9975115) q[0];
rz(-0.83089685) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(0.78713083) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68987304) q[0];
sx q[0];
rz(-1.7053118) q[0];
sx q[0];
rz(0.13264199) q[0];
x q[1];
rz(2.2550341) q[2];
sx q[2];
rz(-1.1349003) q[2];
sx q[2];
rz(0.8800216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8528008) q[1];
sx q[1];
rz(-2.0388985) q[1];
sx q[1];
rz(-1.0352943) q[1];
rz(-0.58148099) q[3];
sx q[3];
rz(-0.6430917) q[3];
sx q[3];
rz(0.067300178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0234915) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(0.033163158) q[2];
rz(1.5432594) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(-1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.24932662) q[0];
sx q[0];
rz(-1.5858269) q[0];
sx q[0];
rz(-1.7484885) q[0];
rz(-0.45685592) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(1.1005864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3063359) q[0];
sx q[0];
rz(-3.0655711) q[0];
sx q[0];
rz(3.118909) q[0];
rz(0.92488168) q[2];
sx q[2];
rz(-2.1585228) q[2];
sx q[2];
rz(2.4921668) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4405358) q[1];
sx q[1];
rz(-0.66796366) q[1];
sx q[1];
rz(-1.6255767) q[1];
rz(2.8131717) q[3];
sx q[3];
rz(-2.7304318) q[3];
sx q[3];
rz(0.31842768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0988203) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(-2.3780195) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(2.1055351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7324657) q[0];
sx q[0];
rz(-0.75508535) q[0];
sx q[0];
rz(-0.29148802) q[0];
rz(-1.905722) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(0.12289563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.204958) q[0];
sx q[0];
rz(-0.35047784) q[0];
sx q[0];
rz(2.9913783) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8529296) q[2];
sx q[2];
rz(-1.7846493) q[2];
sx q[2];
rz(0.58889533) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5220241) q[1];
sx q[1];
rz(-1.777473) q[1];
sx q[1];
rz(3.0289438) q[1];
rz(-pi) q[2];
rz(-3.1355643) q[3];
sx q[3];
rz(-1.42618) q[3];
sx q[3];
rz(0.58352375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5128532) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(-2.0780308) q[2];
rz(0.17744803) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(2.1232846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7040831) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(-0.60072947) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(2.5240555) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3445471) q[0];
sx q[0];
rz(-2.8918242) q[0];
sx q[0];
rz(2.0956371) q[0];
x q[1];
rz(-2.3138758) q[2];
sx q[2];
rz(-1.9736145) q[2];
sx q[2];
rz(0.043443505) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2899298) q[1];
sx q[1];
rz(-1.549822) q[1];
sx q[1];
rz(2.1037654) q[1];
rz(-pi) q[2];
rz(2.6838949) q[3];
sx q[3];
rz(-0.62333306) q[3];
sx q[3];
rz(1.3184354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90126976) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(1.0824925) q[3];
sx q[3];
rz(-0.87477028) q[3];
sx q[3];
rz(3.0860743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33901535) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(-0.045724178) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(1.0411714) q[2];
sx q[2];
rz(-2.0176135) q[2];
sx q[2];
rz(1.8420646) q[2];
rz(2.4847538) q[3];
sx q[3];
rz(-1.4359063) q[3];
sx q[3];
rz(-2.0834854) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
