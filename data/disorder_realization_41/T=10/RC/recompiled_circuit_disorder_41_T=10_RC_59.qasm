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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37786814) q[0];
sx q[0];
rz(-1.8452497) q[0];
sx q[0];
rz(-2.8732357) q[0];
rz(1.4161413) q[2];
sx q[2];
rz(-1.1064331) q[2];
sx q[2];
rz(-2.2048339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3494527) q[1];
sx q[1];
rz(-2.1556902) q[1];
sx q[1];
rz(-2.5930415) q[1];
x q[2];
rz(-2.4258852) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(0.79431278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(1.2692497) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3551536) q[0];
sx q[0];
rz(-2.4053898) q[0];
sx q[0];
rz(-1.2765221) q[0];
rz(-pi) q[1];
rz(2.7361761) q[2];
sx q[2];
rz(-2.9581796) q[2];
sx q[2];
rz(2.237052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0577382) q[1];
sx q[1];
rz(-1.6946304) q[1];
sx q[1];
rz(2.3480575) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41665839) q[3];
sx q[3];
rz(-0.86371242) q[3];
sx q[3];
rz(1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-2.7064586) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86706485) q[0];
sx q[0];
rz(-1.9221677) q[0];
sx q[0];
rz(-1.4562796) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(-1.692786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0355465) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(2.0590904) q[1];
x q[2];
rz(-0.76336236) q[3];
sx q[3];
rz(-1.8244787) q[3];
sx q[3];
rz(0.34624472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(-1.3251925) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6661975) q[0];
sx q[0];
rz(-2.4817433) q[0];
sx q[0];
rz(-3.092479) q[0];
rz(1.2722837) q[2];
sx q[2];
rz(-1.304317) q[2];
sx q[2];
rz(2.4556015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5323822) q[1];
sx q[1];
rz(-2.2717443) q[1];
sx q[1];
rz(2.1125395) q[1];
x q[2];
rz(2.3817252) q[3];
sx q[3];
rz(-0.93899512) q[3];
sx q[3];
rz(2.2894273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081776) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.7081225) q[0];
x q[1];
rz(-2.7258337) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(-1.3410459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(2.8847242) q[1];
rz(-0.080658241) q[3];
sx q[3];
rz(-2.6622052) q[3];
sx q[3];
rz(2.5731034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7197363) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0046376) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(-0.01491551) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8820011) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(1.358658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29218201) q[1];
sx q[1];
rz(-1.9841571) q[1];
sx q[1];
rz(0.11489111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8302912) q[3];
sx q[3];
rz(-0.65400306) q[3];
sx q[3];
rz(1.9037387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-0.66147584) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(0.75659928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7998357) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(1.0210277) q[0];
rz(-pi) q[1];
rz(2.7301894) q[2];
sx q[2];
rz(-1.3855993) q[2];
sx q[2];
rz(-1.411737) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8826897) q[1];
sx q[1];
rz(-1.8001302) q[1];
sx q[1];
rz(-0.16112666) q[1];
rz(-pi) q[2];
rz(0.19161253) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38935223) q[0];
sx q[0];
rz(-1.9703431) q[0];
sx q[0];
rz(2.7826392) q[0];
rz(-pi) q[1];
rz(-3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(1.8708558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5576396) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(0.072903452) q[1];
rz(-pi) q[2];
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
rz(-1.6797592) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(-2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(0.60920238) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-2.2682155) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83335919) q[0];
sx q[0];
rz(-1.6382268) q[0];
sx q[0];
rz(1.7504577) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9940669) q[2];
sx q[2];
rz(-1.7056871) q[2];
sx q[2];
rz(2.877176) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59233353) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(-1.1997644) q[1];
rz(0.94536762) q[3];
sx q[3];
rz(-0.67694596) q[3];
sx q[3];
rz(-0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(0.68230391) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.7451161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6122702) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(-2.8180608) q[0];
rz(-pi) q[1];
rz(-1.2730607) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(1.0422848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1435946) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(0.29923156) q[1];
rz(1.5139919) q[3];
sx q[3];
rz(-1.3089404) q[3];
sx q[3];
rz(1.9522304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-0.15979016) q[2];
rz(-0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(1.8272022) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(-0.53363579) q[3];
sx q[3];
rz(-1.5279557) q[3];
sx q[3];
rz(0.91844311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
