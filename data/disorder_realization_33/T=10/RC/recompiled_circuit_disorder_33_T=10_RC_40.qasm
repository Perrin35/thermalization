OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(-0.27591053) q[0];
sx q[0];
rz(1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966272) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(2.3636723) q[0];
rz(0.76531305) q[2];
sx q[2];
rz(-1.0368477) q[2];
sx q[2];
rz(3.090976) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.03304122) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(-2.7944195) q[1];
rz(-pi) q[2];
rz(0.19574638) q[3];
sx q[3];
rz(-1.9276852) q[3];
sx q[3];
rz(-1.2008047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8685535) q[0];
sx q[0];
rz(-2.2076063) q[0];
sx q[0];
rz(-0.36006948) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1255323) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(2.1330657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.195897) q[1];
sx q[1];
rz(-0.90598124) q[1];
sx q[1];
rz(2.871454) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2651477) q[3];
sx q[3];
rz(-2.090824) q[3];
sx q[3];
rz(-2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.8537834) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(-0.93260971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639097) q[0];
sx q[0];
rz(-1.253486) q[0];
sx q[0];
rz(-0.056563932) q[0];
rz(0.18205299) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(-1.6857266) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0886791) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(-1.4312137) q[1];
x q[2];
rz(2.695735) q[3];
sx q[3];
rz(-1.0599531) q[3];
sx q[3];
rz(-2.0542991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(1.0926584) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(-0.50022593) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(1.6436228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7786176) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(1.0233364) q[0];
rz(-pi) q[1];
x q[1];
rz(1.285032) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(2.6464268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56470358) q[1];
sx q[1];
rz(-0.35208811) q[1];
sx q[1];
rz(2.4114154) q[1];
rz(-pi) q[2];
rz(-0.46521503) q[3];
sx q[3];
rz(-1.1606996) q[3];
sx q[3];
rz(1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-2.4397819) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-1.048208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807555) q[0];
sx q[0];
rz(-1.8341944) q[0];
sx q[0];
rz(-2.8834228) q[0];
rz(1.5585209) q[2];
sx q[2];
rz(-2.2222812) q[2];
sx q[2];
rz(-0.94142454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.198846) q[1];
sx q[1];
rz(-1.8837351) q[1];
sx q[1];
rz(-1.320977) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1225554) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(0.90941959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-3.0117603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31045612) q[0];
sx q[0];
rz(-1.9400915) q[0];
sx q[0];
rz(2.5184758) q[0];
rz(-pi) q[1];
rz(-2.7733299) q[2];
sx q[2];
rz(-2.1149181) q[2];
sx q[2];
rz(-2.1899109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8772014) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(-1.1370204) q[1];
rz(1.767166) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3123902) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(-0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.6947421) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(2.4553305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26353729) q[0];
sx q[0];
rz(-1.1362846) q[0];
sx q[0];
rz(0.52699071) q[0];
rz(-pi) q[1];
rz(-0.94481988) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(-3.0996029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2974907) q[1];
sx q[1];
rz(-0.19980783) q[1];
sx q[1];
rz(1.954345) q[1];
x q[2];
rz(-2.7191914) q[3];
sx q[3];
rz(-1.2761315) q[3];
sx q[3];
rz(-2.9871123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.4461393) q[0];
rz(0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.6400281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9672464) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(1.0780225) q[0];
x q[1];
rz(3.115032) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(2.314687) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20128076) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(-2.8708354) q[1];
rz(-0.99174188) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(0.21608298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.2040899) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4708913) q[0];
sx q[0];
rz(-1.1180709) q[0];
sx q[0];
rz(2.7265413) q[0];
rz(-pi) q[1];
rz(2.4497689) q[2];
sx q[2];
rz(-1.3335388) q[2];
sx q[2];
rz(2.0123864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(-0.07304904) q[1];
x q[2];
rz(-0.44378186) q[3];
sx q[3];
rz(-3.0203331) q[3];
sx q[3];
rz(-1.6463943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(2.0857247) q[3];
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
rz(1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(-0.19432755) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(-1.0338354) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72458306) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(-0.91086046) q[0];
x q[1];
rz(0.50844426) q[2];
sx q[2];
rz(-0.93548453) q[2];
sx q[2];
rz(2.9490162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49627134) q[1];
sx q[1];
rz(-1.6722582) q[1];
sx q[1];
rz(-2.1437777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71318993) q[3];
sx q[3];
rz(-1.7367559) q[3];
sx q[3];
rz(-2.1470269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0620492) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(2.5058084) q[2];
rz(-0.27030269) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-0.96799093) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(-3.071143) q[3];
sx q[3];
rz(-1.9000713) q[3];
sx q[3];
rz(-2.6303359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
