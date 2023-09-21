OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(4.3447189) q[1];
sx q[1];
rz(0.523518) q[1];
sx q[1];
rz(8.5365774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22904299) q[0];
sx q[0];
rz(-0.17670822) q[0];
sx q[0];
rz(-2.8852709) q[0];
rz(-pi) q[1];
rz(2.1076803) q[2];
sx q[2];
rz(-1.4281338) q[2];
sx q[2];
rz(0.85096525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29477316) q[1];
sx q[1];
rz(-1.5713437) q[1];
sx q[1];
rz(0.16695395) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6039443) q[3];
sx q[3];
rz(-0.87246934) q[3];
sx q[3];
rz(2.2693199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28329864) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(2.0236012) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-3.0956691) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730826) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.1647859) q[1];
sx q[1];
rz(-3.0156946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7315642) q[0];
sx q[0];
rz(-1.0985939) q[0];
sx q[0];
rz(-2.7239951) q[0];
x q[1];
rz(-0.84463859) q[2];
sx q[2];
rz(-2.3181049) q[2];
sx q[2];
rz(2.5708432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44015917) q[1];
sx q[1];
rz(-1.6007746) q[1];
sx q[1];
rz(-1.4440126) q[1];
rz(-pi) q[2];
rz(1.5829854) q[3];
sx q[3];
rz(-2.4180531) q[3];
sx q[3];
rz(2.7140868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(-2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(-1.1598587) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.738049) q[0];
sx q[0];
rz(-1.4883211) q[0];
sx q[0];
rz(1.0595881) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50767501) q[2];
sx q[2];
rz(-2.1360364) q[2];
sx q[2];
rz(0.6914247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0594271) q[1];
sx q[1];
rz(-1.9485928) q[1];
sx q[1];
rz(0.49281812) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4781811) q[3];
sx q[3];
rz(-2.5941656) q[3];
sx q[3];
rz(-2.3401591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7362061) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(0.64374271) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(-2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(1.6595586) q[0];
rz(-2.4404793) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(-0.28465095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4581504) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(3.0577858) q[0];
rz(-pi) q[1];
rz(-2.1842723) q[2];
sx q[2];
rz(-0.34733221) q[2];
sx q[2];
rz(-0.50328244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2440566) q[1];
sx q[1];
rz(-1.4923864) q[1];
sx q[1];
rz(1.4439911) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1553667) q[3];
sx q[3];
rz(-0.75891906) q[3];
sx q[3];
rz(-0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.231679) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.652997) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.8150785) q[0];
rz(1.1524221) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(2.2036536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4834578) q[0];
sx q[0];
rz(-1.5926653) q[0];
sx q[0];
rz(1.6629501) q[0];
x q[1];
rz(0.83154251) q[2];
sx q[2];
rz(-2.4387896) q[2];
sx q[2];
rz(2.0895095) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2069619) q[1];
sx q[1];
rz(-0.043255581) q[1];
sx q[1];
rz(-2.2224777) q[1];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(-0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9662629) q[2];
sx q[2];
rz(-2.9523409) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(2.5842216) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(2.5851137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59069955) q[0];
sx q[0];
rz(-1.7771582) q[0];
sx q[0];
rz(-2.7569689) q[0];
rz(-pi) q[1];
rz(-0.85530497) q[2];
sx q[2];
rz(-1.7024634) q[2];
sx q[2];
rz(-1.635074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52292697) q[1];
sx q[1];
rz(-1.3892801) q[1];
sx q[1];
rz(1.3154485) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3658386) q[3];
sx q[3];
rz(-1.6072825) q[3];
sx q[3];
rz(3.0690103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(-1.9167985) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(-1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0790134) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(3.0986837) q[0];
rz(-2.2242916) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(2.6409805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635054) q[0];
sx q[0];
rz(-0.18112077) q[0];
sx q[0];
rz(1.7630793) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31502864) q[2];
sx q[2];
rz(-1.050596) q[2];
sx q[2];
rz(2.5556285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4088879) q[1];
sx q[1];
rz(-1.171247) q[1];
sx q[1];
rz(-1.4399745) q[1];
x q[2];
rz(-0.75617366) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(2.365436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(2.4659757) q[2];
rz(-0.44089857) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(-2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0764517) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(-1.0674397) q[0];
rz(-2.7087129) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(1.4656461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6976801) q[0];
sx q[0];
rz(-1.4364103) q[0];
sx q[0];
rz(2.0356376) q[0];
x q[1];
rz(-0.18657121) q[2];
sx q[2];
rz(-1.4636883) q[2];
sx q[2];
rz(-0.69498108) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17897478) q[1];
sx q[1];
rz(-2.0638421) q[1];
sx q[1];
rz(0.38260539) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.397714) q[3];
sx q[3];
rz(-0.79259593) q[3];
sx q[3];
rz(2.511123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33891588) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(-0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(-0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-0.73721686) q[0];
rz(-0.018741477) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(-0.92528701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5697445) q[0];
sx q[0];
rz(-1.0751372) q[0];
sx q[0];
rz(0.75832383) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7032353) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(0.027564136) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3146989) q[1];
sx q[1];
rz(-1.3876545) q[1];
sx q[1];
rz(1.7032743) q[1];
x q[2];
rz(-2.0426644) q[3];
sx q[3];
rz(-2.3657551) q[3];
sx q[3];
rz(2.8692506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(-2.1378689) q[2];
rz(-0.090099661) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(1.2596624) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(2.3840747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82523358) q[0];
sx q[0];
rz(-0.23783437) q[0];
sx q[0];
rz(2.1912078) q[0];
rz(0.9805571) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(-2.5622501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6706558) q[1];
sx q[1];
rz(-0.48479143) q[1];
sx q[1];
rz(-2.0623341) q[1];
rz(-pi) q[2];
rz(0.032978756) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(-0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(0.79375664) q[2];
rz(0.63888597) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6486075) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-2.5095148) q[2];
sx q[2];
rz(-0.52543228) q[2];
sx q[2];
rz(0.57731522) q[2];
rz(2.8201841) q[3];
sx q[3];
rz(-2.1264429) q[3];
sx q[3];
rz(0.85254729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
