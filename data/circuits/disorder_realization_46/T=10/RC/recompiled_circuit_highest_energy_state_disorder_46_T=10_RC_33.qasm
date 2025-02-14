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
rz(2.1815648) q[0];
sx q[0];
rz(-0.61909827) q[0];
sx q[0];
rz(1.4033138) q[0];
rz(2.4461441) q[1];
sx q[1];
rz(-0.56013501) q[1];
sx q[1];
rz(-2.4897895) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8389002) q[0];
sx q[0];
rz(-2.3639285) q[0];
sx q[0];
rz(1.0833141) q[0];
rz(-1.5648834) q[2];
sx q[2];
rz(-1.8522709) q[2];
sx q[2];
rz(0.23879063) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4308555) q[1];
sx q[1];
rz(-2.9086782) q[1];
sx q[1];
rz(-3.0431299) q[1];
rz(2.321389) q[3];
sx q[3];
rz(-1.3705882) q[3];
sx q[3];
rz(2.8516796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2191676) q[2];
sx q[2];
rz(-1.3895915) q[2];
sx q[2];
rz(-2.8233042) q[2];
rz(2.2394032) q[3];
sx q[3];
rz(-0.92034322) q[3];
sx q[3];
rz(0.86377803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1282463) q[0];
sx q[0];
rz(-2.8885169) q[0];
sx q[0];
rz(2.6213562) q[0];
rz(3.0546313) q[1];
sx q[1];
rz(-2.088701) q[1];
sx q[1];
rz(0.797995) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3083852) q[0];
sx q[0];
rz(-1.7705581) q[0];
sx q[0];
rz(2.8248592) q[0];
rz(-pi) q[1];
rz(-2.3981061) q[2];
sx q[2];
rz(-1.6891589) q[2];
sx q[2];
rz(1.3918882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50988156) q[1];
sx q[1];
rz(-2.7260511) q[1];
sx q[1];
rz(0.35453896) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9738484) q[3];
sx q[3];
rz(-1.155793) q[3];
sx q[3];
rz(0.61674207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8331376) q[2];
sx q[2];
rz(-1.986958) q[2];
sx q[2];
rz(-2.7412565) q[2];
rz(-0.18276246) q[3];
sx q[3];
rz(-0.57923135) q[3];
sx q[3];
rz(-1.2919174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12277814) q[0];
sx q[0];
rz(-3.0570539) q[0];
sx q[0];
rz(-1.9345181) q[0];
rz(1.5030376) q[1];
sx q[1];
rz(-1.259558) q[1];
sx q[1];
rz(0.17301339) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4543145) q[0];
sx q[0];
rz(-1.1862687) q[0];
sx q[0];
rz(-3.049535) q[0];
x q[1];
rz(0.9890828) q[2];
sx q[2];
rz(-1.4787041) q[2];
sx q[2];
rz(-1.3710772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.971774) q[1];
sx q[1];
rz(-2.1425852) q[1];
sx q[1];
rz(-2.1717768) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9152619) q[3];
sx q[3];
rz(-2.3809098) q[3];
sx q[3];
rz(2.5193391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5485237) q[2];
sx q[2];
rz(-2.5544281) q[2];
sx q[2];
rz(2.1068088) q[2];
rz(-0.89094025) q[3];
sx q[3];
rz(-0.8689298) q[3];
sx q[3];
rz(-1.9152036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064875038) q[0];
sx q[0];
rz(-0.95219505) q[0];
sx q[0];
rz(3.0127443) q[0];
rz(2.7128362) q[1];
sx q[1];
rz(-1.4219475) q[1];
sx q[1];
rz(-1.6645924) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3698206) q[0];
sx q[0];
rz(-1.6390529) q[0];
sx q[0];
rz(-0.34567709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3643983) q[2];
sx q[2];
rz(-2.5487196) q[2];
sx q[2];
rz(-0.89250084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0018938) q[1];
sx q[1];
rz(-0.47283462) q[1];
sx q[1];
rz(2.0496778) q[1];
rz(1.5724395) q[3];
sx q[3];
rz(-2.5503359) q[3];
sx q[3];
rz(1.3639243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9513272) q[2];
sx q[2];
rz(-1.724388) q[2];
sx q[2];
rz(1.3235271) q[2];
rz(-1.6773112) q[3];
sx q[3];
rz(-0.5885632) q[3];
sx q[3];
rz(0.67414362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73439634) q[0];
sx q[0];
rz(-0.36235991) q[0];
sx q[0];
rz(1.6998442) q[0];
rz(-2.7184519) q[1];
sx q[1];
rz(-0.95760456) q[1];
sx q[1];
rz(2.8513681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30681243) q[0];
sx q[0];
rz(-2.4270822) q[0];
sx q[0];
rz(0.45551703) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3168081) q[2];
sx q[2];
rz(-0.76003492) q[2];
sx q[2];
rz(-1.7766376) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.90846244) q[1];
sx q[1];
rz(-1.6987579) q[1];
sx q[1];
rz(1.5217693) q[1];
rz(-0.69921826) q[3];
sx q[3];
rz(-2.7479246) q[3];
sx q[3];
rz(0.75692219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71771249) q[2];
sx q[2];
rz(-2.3545357) q[2];
sx q[2];
rz(-0.94669739) q[2];
rz(2.8801019) q[3];
sx q[3];
rz(-1.7777092) q[3];
sx q[3];
rz(2.7546895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.43712) q[0];
sx q[0];
rz(-1.8081212) q[0];
sx q[0];
rz(0.76960808) q[0];
rz(0.17598027) q[1];
sx q[1];
rz(-1.8006005) q[1];
sx q[1];
rz(1.5699068) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0933857) q[0];
sx q[0];
rz(-0.95952672) q[0];
sx q[0];
rz(-0.78637357) q[0];
rz(0.85167517) q[2];
sx q[2];
rz(-1.5181418) q[2];
sx q[2];
rz(-2.5824314) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9573201) q[1];
sx q[1];
rz(-1.9369003) q[1];
sx q[1];
rz(1.7680607) q[1];
rz(0.79738252) q[3];
sx q[3];
rz(-2.4719072) q[3];
sx q[3];
rz(2.0978239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1311329) q[2];
sx q[2];
rz(-0.95644462) q[2];
sx q[2];
rz(-1.825911) q[2];
rz(1.9192421) q[3];
sx q[3];
rz(-1.1244011) q[3];
sx q[3];
rz(-2.840672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2216126) q[0];
sx q[0];
rz(-1.6918809) q[0];
sx q[0];
rz(-1.8020887) q[0];
rz(-1.9446531) q[1];
sx q[1];
rz(-1.0547799) q[1];
sx q[1];
rz(-2.1766677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20112637) q[0];
sx q[0];
rz(-2.622704) q[0];
sx q[0];
rz(-0.19176439) q[0];
rz(-pi) q[1];
x q[1];
rz(2.262124) q[2];
sx q[2];
rz(-2.6973638) q[2];
sx q[2];
rz(-2.6951126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33789148) q[1];
sx q[1];
rz(-1.6543904) q[1];
sx q[1];
rz(2.9210832) q[1];
rz(-pi) q[2];
rz(-1.3960725) q[3];
sx q[3];
rz(-1.0936605) q[3];
sx q[3];
rz(1.9486547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40184608) q[2];
sx q[2];
rz(-1.9332989) q[2];
sx q[2];
rz(1.1811264) q[2];
rz(0.014569672) q[3];
sx q[3];
rz(-0.40112344) q[3];
sx q[3];
rz(-1.1723664) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75083098) q[0];
sx q[0];
rz(-1.2068692) q[0];
sx q[0];
rz(-2.3857351) q[0];
rz(-2.4769056) q[1];
sx q[1];
rz(-1.1988147) q[1];
sx q[1];
rz(-1.2896779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9826848) q[0];
sx q[0];
rz(-2.3450627) q[0];
sx q[0];
rz(3.0523053) q[0];
rz(-pi) q[1];
rz(1.8721183) q[2];
sx q[2];
rz(-2.0908643) q[2];
sx q[2];
rz(0.10094303) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3933223) q[1];
sx q[1];
rz(-0.38667187) q[1];
sx q[1];
rz(2.3063117) q[1];
rz(2.8323402) q[3];
sx q[3];
rz(-1.6669295) q[3];
sx q[3];
rz(0.94667305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1332625) q[2];
sx q[2];
rz(-1.8066758) q[2];
sx q[2];
rz(1.5670212) q[2];
rz(1.745892) q[3];
sx q[3];
rz(-1.402366) q[3];
sx q[3];
rz(2.375026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.72346) q[0];
sx q[0];
rz(-2.8657931) q[0];
sx q[0];
rz(-2.7511399) q[0];
rz(1.1013203) q[1];
sx q[1];
rz(-1.7951671) q[1];
sx q[1];
rz(-2.2094545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2208122) q[0];
sx q[0];
rz(-1.6075862) q[0];
sx q[0];
rz(-3.1139741) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53932346) q[2];
sx q[2];
rz(-1.5289837) q[2];
sx q[2];
rz(-1.5635179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4313163) q[1];
sx q[1];
rz(-1.2936771) q[1];
sx q[1];
rz(1.3701012) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6868709) q[3];
sx q[3];
rz(-2.623436) q[3];
sx q[3];
rz(-1.3738812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0497047) q[2];
sx q[2];
rz(-0.69685093) q[2];
sx q[2];
rz(1.4616802) q[2];
rz(-1.0852496) q[3];
sx q[3];
rz(-2.0094252) q[3];
sx q[3];
rz(0.6404883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28814462) q[0];
sx q[0];
rz(-2.4818821) q[0];
sx q[0];
rz(1.1261384) q[0];
rz(2.0938734) q[1];
sx q[1];
rz(-0.62647096) q[1];
sx q[1];
rz(-1.610021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5384993) q[0];
sx q[0];
rz(-1.5694008) q[0];
sx q[0];
rz(-3.1183232) q[0];
rz(-pi) q[1];
rz(0.45827403) q[2];
sx q[2];
rz(-1.7802139) q[2];
sx q[2];
rz(2.6334514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7992656) q[1];
sx q[1];
rz(-1.5060358) q[1];
sx q[1];
rz(1.4928476) q[1];
rz(-pi) q[2];
rz(1.6592386) q[3];
sx q[3];
rz(-1.0103265) q[3];
sx q[3];
rz(-1.1286398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62022007) q[2];
sx q[2];
rz(-2.3282101) q[2];
sx q[2];
rz(-2.9353471) q[2];
rz(-2.5824052) q[3];
sx q[3];
rz(-1.5234448) q[3];
sx q[3];
rz(2.0172113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6013721) q[0];
sx q[0];
rz(-1.843353) q[0];
sx q[0];
rz(-2.699615) q[0];
rz(1.1454918) q[1];
sx q[1];
rz(-1.3064697) q[1];
sx q[1];
rz(-1.3189955) q[1];
rz(-2.3445208) q[2];
sx q[2];
rz(-2.7596724) q[2];
sx q[2];
rz(1.8308664) q[2];
rz(-1.3973936) q[3];
sx q[3];
rz(-0.85416334) q[3];
sx q[3];
rz(-0.76696542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
