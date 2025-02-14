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
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(1.0083415) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5429687) q[0];
sx q[0];
rz(-2.7592139) q[0];
sx q[0];
rz(1.7837058) q[0];
x q[1];
rz(0.37990976) q[2];
sx q[2];
rz(-1.7899982) q[2];
sx q[2];
rz(-2.4394112) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.086585933) q[1];
sx q[1];
rz(-2.0457959) q[1];
sx q[1];
rz(-1.7473298) q[1];
rz(-pi) q[2];
rz(1.6761682) q[3];
sx q[3];
rz(-0.19427768) q[3];
sx q[3];
rz(0.36997488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8334373) q[2];
sx q[2];
rz(-0.94352949) q[2];
sx q[2];
rz(-0.65307871) q[2];
rz(-0.11450442) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7750875) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-2.7313045) q[1];
sx q[1];
rz(-0.62058273) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4575093) q[0];
sx q[0];
rz(-0.39827049) q[0];
sx q[0];
rz(1.0578367) q[0];
x q[1];
rz(2.7730915) q[2];
sx q[2];
rz(-0.93612367) q[2];
sx q[2];
rz(-2.2919185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0479792) q[1];
sx q[1];
rz(-1.9318214) q[1];
sx q[1];
rz(-2.4314636) q[1];
rz(-pi) q[2];
rz(-2.1866131) q[3];
sx q[3];
rz(-1.338025) q[3];
sx q[3];
rz(-1.3835088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4429984) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(-0.2300187) q[2];
rz(1.5551785) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25552937) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(2.0607167) q[0];
rz(2.4983662) q[1];
sx q[1];
rz(-2.3422362) q[1];
sx q[1];
rz(-0.08531514) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9935308) q[0];
sx q[0];
rz(-0.9328273) q[0];
sx q[0];
rz(2.30739) q[0];
rz(-pi) q[1];
rz(2.8181452) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(-1.310854) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7300295) q[1];
sx q[1];
rz(-1.4394898) q[1];
sx q[1];
rz(-2.9741785) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5203568) q[3];
sx q[3];
rz(-1.6077822) q[3];
sx q[3];
rz(-2.1153574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9366511) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(-0.19634518) q[2];
rz(2.8593072) q[3];
sx q[3];
rz(-2.0245656) q[3];
sx q[3];
rz(1.0452247) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.008217) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(0.38544449) q[0];
rz(-1.7281945) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(-0.24994303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871368) q[0];
sx q[0];
rz(-0.76351316) q[0];
sx q[0];
rz(-0.95434965) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6310299) q[2];
sx q[2];
rz(-1.4830198) q[2];
sx q[2];
rz(-2.6812832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0562607) q[1];
sx q[1];
rz(-2.0592494) q[1];
sx q[1];
rz(-1.3794829) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3718666) q[3];
sx q[3];
rz(-2.3545111) q[3];
sx q[3];
rz(0.50686344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6803153) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(-0.2363905) q[2];
rz(1.2605028) q[3];
sx q[3];
rz(-0.25006306) q[3];
sx q[3];
rz(2.4642956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054955) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(-0.48511037) q[0];
rz(1.1803892) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(0.83650437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9125875) q[0];
sx q[0];
rz(-1.8857755) q[0];
sx q[0];
rz(-0.37512414) q[0];
x q[1];
rz(-0.31827688) q[2];
sx q[2];
rz(-1.9564544) q[2];
sx q[2];
rz(1.1952656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4662469) q[1];
sx q[1];
rz(-2.2398178) q[1];
sx q[1];
rz(0.061208486) q[1];
rz(-pi) q[2];
rz(-2.9936166) q[3];
sx q[3];
rz(-2.736909) q[3];
sx q[3];
rz(-1.2860822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89852077) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(2.3720429) q[2];
rz(-2.2122993) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(-0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1997851) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(0.76882452) q[0];
rz(1.3517514) q[1];
sx q[1];
rz(-2.6685346) q[1];
sx q[1];
rz(-0.88968712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6883967) q[0];
sx q[0];
rz(-1.9510117) q[0];
sx q[0];
rz(-2.9158457) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0013498505) q[2];
sx q[2];
rz(-1.5264411) q[2];
sx q[2];
rz(-0.42429098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8823679) q[1];
sx q[1];
rz(-0.15341694) q[1];
sx q[1];
rz(-1.2620366) q[1];
x q[2];
rz(2.9075648) q[3];
sx q[3];
rz(-2.6796973) q[3];
sx q[3];
rz(-1.6017101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6512904) q[2];
sx q[2];
rz(-1.4375261) q[2];
sx q[2];
rz(1.5550522) q[2];
rz(1.2706612) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(-0.13203013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9689869) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(1.1440811) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(2.3544618) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68987304) q[0];
sx q[0];
rz(-1.4362808) q[0];
sx q[0];
rz(3.0089507) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93573715) q[2];
sx q[2];
rz(-2.3496186) q[2];
sx q[2];
rz(-2.9284649) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28879189) q[1];
sx q[1];
rz(-1.1026942) q[1];
sx q[1];
rz(1.0352943) q[1];
x q[2];
rz(-1.1803341) q[3];
sx q[3];
rz(-1.045908) q[3];
sx q[3];
rz(-2.386664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1181011) q[2];
sx q[2];
rz(-0.96994895) q[2];
sx q[2];
rz(3.1084295) q[2];
rz(-1.5983332) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892266) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(-1.3931042) q[0];
rz(-2.6847367) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(-1.1005864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2835866) q[0];
sx q[0];
rz(-1.6467983) q[0];
sx q[0];
rz(-1.5690687) q[0];
x q[1];
rz(2.4068836) q[2];
sx q[2];
rz(-2.2976954) q[2];
sx q[2];
rz(-1.5555842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7010569) q[1];
sx q[1];
rz(-0.66796366) q[1];
sx q[1];
rz(-1.5160159) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3914137) q[3];
sx q[3];
rz(-1.7000704) q[3];
sx q[3];
rz(-1.5551274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.042772375) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(0.76357311) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(-1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(0.29148802) q[0];
rz(1.2358707) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(0.12289563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50704573) q[0];
sx q[0];
rz(-1.6222008) q[0];
sx q[0];
rz(2.794751) q[0];
rz(-2.2332623) q[2];
sx q[2];
rz(-0.35229063) q[2];
sx q[2];
rz(0.34991821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0258515) q[1];
sx q[1];
rz(-2.9065955) q[1];
sx q[1];
rz(2.0629289) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0060283383) q[3];
sx q[3];
rz(-1.7154126) q[3];
sx q[3];
rz(-0.58352375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62873944) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(2.0780308) q[2];
rz(-0.17744803) q[3];
sx q[3];
rz(-0.95824233) q[3];
sx q[3];
rz(2.1232846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43750957) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.7105239) q[0];
rz(2.5408632) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(0.61753714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8059313) q[0];
sx q[0];
rz(-1.3552203) q[0];
sx q[0];
rz(-3.0144601) q[0];
x q[1];
rz(2.1328636) q[2];
sx q[2];
rz(-2.3150839) q[2];
sx q[2];
rz(1.1240608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2899298) q[1];
sx q[1];
rz(-1.5917707) q[1];
sx q[1];
rz(1.0378273) q[1];
rz(-pi) q[2];
rz(-1.8784058) q[3];
sx q[3];
rz(-2.121939) q[3];
sx q[3];
rz(1.8637509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2403229) q[2];
sx q[2];
rz(-2.1412886) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(1.0824925) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-3.0860743) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8025773) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(-0.045724178) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(-0.81188079) q[2];
sx q[2];
rz(-0.67888454) q[2];
sx q[2];
rz(0.90711403) q[2];
rz(1.4010728) q[3];
sx q[3];
rz(-2.2206497) q[3];
sx q[3];
rz(2.7322265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
