OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0140822) q[0];
sx q[0];
rz(-0.2427225) q[0];
sx q[0];
rz(2.2665562) q[0];
rz(-1.8889282) q[1];
sx q[1];
rz(-1.6369605) q[1];
sx q[1];
rz(2.5445282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.208114) q[0];
sx q[0];
rz(-0.39486082) q[0];
sx q[0];
rz(1.3482679) q[0];
x q[1];
rz(-1.784034) q[2];
sx q[2];
rz(-0.87940826) q[2];
sx q[2];
rz(0.85644507) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41433147) q[1];
sx q[1];
rz(-2.4684069) q[1];
sx q[1];
rz(3.0607231) q[1];
x q[2];
rz(0.050757082) q[3];
sx q[3];
rz(-1.4043706) q[3];
sx q[3];
rz(0.084648477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.732932) q[2];
sx q[2];
rz(-0.084736846) q[2];
sx q[2];
rz(1.6049467) q[2];
rz(-1.5417277) q[3];
sx q[3];
rz(-2.7432975) q[3];
sx q[3];
rz(1.0140243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23564553) q[0];
sx q[0];
rz(-2.5962317) q[0];
sx q[0];
rz(1.1646618) q[0];
rz(0.91854489) q[1];
sx q[1];
rz(-2.6530177) q[1];
sx q[1];
rz(-2.4426544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0526657) q[0];
sx q[0];
rz(-2.1667826) q[0];
sx q[0];
rz(1.862103) q[0];
x q[1];
rz(0.21007819) q[2];
sx q[2];
rz(-0.89001402) q[2];
sx q[2];
rz(-1.4821881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0684425) q[1];
sx q[1];
rz(-3.0038976) q[1];
sx q[1];
rz(0.87611468) q[1];
rz(0.1546448) q[3];
sx q[3];
rz(-0.15081146) q[3];
sx q[3];
rz(0.035816222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32236448) q[2];
sx q[2];
rz(-1.7027723) q[2];
sx q[2];
rz(1.2054075) q[2];
rz(-0.15959218) q[3];
sx q[3];
rz(-1.3887838) q[3];
sx q[3];
rz(-3.0157183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.50652248) q[0];
sx q[0];
rz(-3.0634026) q[0];
sx q[0];
rz(-0.11678326) q[0];
rz(2.0109239) q[1];
sx q[1];
rz(-3.0394381) q[1];
sx q[1];
rz(-0.048361383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65927261) q[0];
sx q[0];
rz(-3.0249615) q[0];
sx q[0];
rz(-0.92088033) q[0];
x q[1];
rz(1.5226256) q[2];
sx q[2];
rz(-1.4496231) q[2];
sx q[2];
rz(-0.19751422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.060052311) q[1];
sx q[1];
rz(-2.2339666) q[1];
sx q[1];
rz(2.7377048) q[1];
rz(-pi) q[2];
rz(-0.039325754) q[3];
sx q[3];
rz(-1.4429386) q[3];
sx q[3];
rz(0.10165299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9174663) q[2];
sx q[2];
rz(-1.5175061) q[2];
sx q[2];
rz(-0.65829128) q[2];
rz(1.0607464) q[3];
sx q[3];
rz(-2.3297533) q[3];
sx q[3];
rz(-2.9909824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.1310405) q[0];
sx q[0];
rz(-0.083715938) q[0];
sx q[0];
rz(-0.87378275) q[0];
rz(2.1281758) q[1];
sx q[1];
rz(-3.1384917) q[1];
sx q[1];
rz(-0.82469624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8270758) q[0];
sx q[0];
rz(-1.5387879) q[0];
sx q[0];
rz(1.5289983) q[0];
rz(-pi) q[1];
rz(-2.1120788) q[2];
sx q[2];
rz(-0.80179399) q[2];
sx q[2];
rz(3.1355592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4383126) q[1];
sx q[1];
rz(-1.6973799) q[1];
sx q[1];
rz(-1.2939318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6225067) q[3];
sx q[3];
rz(-1.9955764) q[3];
sx q[3];
rz(0.82889601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.970845) q[2];
sx q[2];
rz(-2.8837995) q[2];
sx q[2];
rz(2.8802059) q[2];
rz(1.9445253) q[3];
sx q[3];
rz(-1.7944929) q[3];
sx q[3];
rz(-0.66024595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9561387) q[0];
sx q[0];
rz(-3.0281797) q[0];
sx q[0];
rz(-2.043682) q[0];
rz(-2.5838891) q[1];
sx q[1];
rz(-3.1126366) q[1];
sx q[1];
rz(1.8580233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4056617) q[0];
sx q[0];
rz(-0.064583555) q[0];
sx q[0];
rz(-0.92900999) q[0];
rz(0.27622707) q[2];
sx q[2];
rz(-0.27135393) q[2];
sx q[2];
rz(0.63947397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9030212) q[1];
sx q[1];
rz(-1.0728749) q[1];
sx q[1];
rz(1.4378888) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2992283) q[3];
sx q[3];
rz(-2.0443562) q[3];
sx q[3];
rz(2.0061988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6534609) q[2];
sx q[2];
rz(-0.9119432) q[2];
sx q[2];
rz(-0.27257356) q[2];
rz(-1.9778947) q[3];
sx q[3];
rz(-1.8018164) q[3];
sx q[3];
rz(0.94289601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96384752) q[0];
sx q[0];
rz(-3.0680532) q[0];
sx q[0];
rz(-0.22200008) q[0];
rz(1.47413) q[1];
sx q[1];
rz(-0.024802955) q[1];
sx q[1];
rz(2.7892392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8519634) q[0];
sx q[0];
rz(-1.5705938) q[0];
sx q[0];
rz(-1.5727497) q[0];
x q[1];
rz(-0.58587296) q[2];
sx q[2];
rz(-2.2112741) q[2];
sx q[2];
rz(-2.7197995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7166802) q[1];
sx q[1];
rz(-1.0031478) q[1];
sx q[1];
rz(2.6006954) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0025939) q[3];
sx q[3];
rz(-0.96240265) q[3];
sx q[3];
rz(-2.268057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12677975) q[2];
sx q[2];
rz(-1.959602) q[2];
sx q[2];
rz(-0.78753769) q[2];
rz(3.0102503) q[3];
sx q[3];
rz(-0.98024875) q[3];
sx q[3];
rz(0.76287764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545559) q[0];
sx q[0];
rz(-3.1058703) q[0];
sx q[0];
rz(0.63294739) q[0];
rz(0.58730209) q[1];
sx q[1];
rz(-3.1226698) q[1];
sx q[1];
rz(-2.6515554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5929664) q[0];
sx q[0];
rz(-1.6621309) q[0];
sx q[0];
rz(1.4110325) q[0];
x q[1];
rz(2.7964726) q[2];
sx q[2];
rz(-1.39087) q[2];
sx q[2];
rz(-0.82054446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32691218) q[1];
sx q[1];
rz(-0.94409724) q[1];
sx q[1];
rz(-1.9912849) q[1];
rz(0.96411772) q[3];
sx q[3];
rz(-0.58219606) q[3];
sx q[3];
rz(-2.4011998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6106762) q[2];
sx q[2];
rz(-0.55675113) q[2];
sx q[2];
rz(0.69182932) q[2];
rz(-2.1107215) q[3];
sx q[3];
rz(-2.2014047) q[3];
sx q[3];
rz(0.33215365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0523025) q[0];
sx q[0];
rz(-0.060929935) q[0];
sx q[0];
rz(0.45284176) q[0];
rz(-2.7678658) q[1];
sx q[1];
rz(-3.0174297) q[1];
sx q[1];
rz(-2.9492818) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6598337) q[0];
sx q[0];
rz(-3.026986) q[0];
sx q[0];
rz(-3.0069754) q[0];
rz(1.9779124) q[2];
sx q[2];
rz(-1.0958899) q[2];
sx q[2];
rz(-2.6559336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3329703) q[1];
sx q[1];
rz(-2.0870627) q[1];
sx q[1];
rz(-2.9201169) q[1];
x q[2];
rz(1.9662596) q[3];
sx q[3];
rz(-2.7928958) q[3];
sx q[3];
rz(-2.6570081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33334357) q[2];
sx q[2];
rz(-3.0875751) q[2];
sx q[2];
rz(-0.56600189) q[2];
rz(0.69837767) q[3];
sx q[3];
rz(-1.799182) q[3];
sx q[3];
rz(-3.0799871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9653387) q[0];
sx q[0];
rz(-2.7955671) q[0];
sx q[0];
rz(1.6842496) q[0];
rz(-2.1894042) q[1];
sx q[1];
rz(-0.0038650611) q[1];
sx q[1];
rz(1.7239408) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24704449) q[0];
sx q[0];
rz(-1.7590917) q[0];
sx q[0];
rz(-1.2597127) q[0];
rz(-pi) q[1];
x q[1];
rz(0.023073828) q[2];
sx q[2];
rz(-2.0947813) q[2];
sx q[2];
rz(0.11701458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6030857) q[1];
sx q[1];
rz(-2.9630447) q[1];
sx q[1];
rz(-2.958955) q[1];
rz(-pi) q[2];
rz(1.438497) q[3];
sx q[3];
rz(-1.6198564) q[3];
sx q[3];
rz(2.8956312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62966043) q[2];
sx q[2];
rz(-0.28714445) q[2];
sx q[2];
rz(2.6426537) q[2];
rz(2.387909) q[3];
sx q[3];
rz(-0.53871173) q[3];
sx q[3];
rz(-2.5628569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27624398) q[0];
sx q[0];
rz(-2.8911599) q[0];
sx q[0];
rz(-2.8540386) q[0];
rz(-2.4240671) q[1];
sx q[1];
rz(-3.0375807) q[1];
sx q[1];
rz(1.2984553) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8632057) q[0];
sx q[0];
rz(-2.1223377) q[0];
sx q[0];
rz(-2.2831149) q[0];
rz(-pi) q[1];
rz(0.98787843) q[2];
sx q[2];
rz(-2.0955502) q[2];
sx q[2];
rz(-0.74452213) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76261452) q[1];
sx q[1];
rz(-1.3022547) q[1];
sx q[1];
rz(-1.1648127) q[1];
rz(-pi) q[2];
rz(-2.82622) q[3];
sx q[3];
rz(-1.4829691) q[3];
sx q[3];
rz(1.4855291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.31460497) q[2];
sx q[2];
rz(-1.2890451) q[2];
sx q[2];
rz(2.2513466) q[2];
rz(-0.024197772) q[3];
sx q[3];
rz(-0.20166339) q[3];
sx q[3];
rz(2.3307687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496307) q[0];
sx q[0];
rz(-0.55593119) q[0];
sx q[0];
rz(2.2297106) q[0];
rz(-0.98056071) q[1];
sx q[1];
rz(-0.33976561) q[1];
sx q[1];
rz(0.38358546) q[1];
rz(-1.3487074) q[2];
sx q[2];
rz(-0.32094706) q[2];
sx q[2];
rz(1.361346) q[2];
rz(-1.2804563) q[3];
sx q[3];
rz(-1.803259) q[3];
sx q[3];
rz(-2.4358376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
