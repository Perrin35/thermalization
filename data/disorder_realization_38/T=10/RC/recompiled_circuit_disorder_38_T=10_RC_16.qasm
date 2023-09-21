OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(-0.5501774) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43139797) q[0];
sx q[0];
rz(-0.22127998) q[0];
sx q[0];
rz(1.6943323) q[0];
x q[1];
rz(-1.0371738) q[2];
sx q[2];
rz(-1.1871157) q[2];
sx q[2];
rz(-1.2984315) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3783592) q[1];
sx q[1];
rz(-0.86508703) q[1];
sx q[1];
rz(2.6325429) q[1];
rz(-3.0580871) q[3];
sx q[3];
rz(-2.180035) q[3];
sx q[3];
rz(3.1374251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7258519) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(-0.63981167) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(-1.6289904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4166116) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(-2.7532817) q[0];
x q[1];
rz(1.5674044) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(2.0397182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0257033) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(0.89110903) q[1];
x q[2];
rz(-1.8566425) q[3];
sx q[3];
rz(-1.3405521) q[3];
sx q[3];
rz(2.2765991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-0.025618205) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(0.29552466) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94507664) q[0];
sx q[0];
rz(-1.7917504) q[0];
sx q[0];
rz(-0.093284746) q[0];
rz(-pi) q[1];
rz(0.27772851) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(3.1090528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5902139) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(2.6022634) q[1];
rz(-pi) q[2];
rz(-2.6882719) q[3];
sx q[3];
rz(-1.7954149) q[3];
sx q[3];
rz(0.9325222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-2.8682958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1647987) q[0];
sx q[0];
rz(-1.2480038) q[0];
sx q[0];
rz(0.093683634) q[0];
rz(-pi) q[1];
rz(-2.1347413) q[2];
sx q[2];
rz(-1.0014921) q[2];
sx q[2];
rz(3.0727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71972388) q[1];
sx q[1];
rz(-0.27742741) q[1];
sx q[1];
rz(2.1348907) q[1];
x q[2];
rz(1.3917543) q[3];
sx q[3];
rz(-2.2125707) q[3];
sx q[3];
rz(0.87953506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79167241) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(-2.6376574) q[2];
rz(-0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-0.3751522) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50901978) q[0];
sx q[0];
rz(-2.1633254) q[0];
sx q[0];
rz(0.49125262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5802025) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(2.8807092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0070237006) q[1];
sx q[1];
rz(-0.94788523) q[1];
sx q[1];
rz(-1.9593777) q[1];
x q[2];
rz(3.0994814) q[3];
sx q[3];
rz(-1.9633506) q[3];
sx q[3];
rz(-2.7819355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8905028) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(0.791839) q[0];
rz(-2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148221) q[0];
sx q[0];
rz(-1.556067) q[0];
sx q[0];
rz(1.8526088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.860414) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(0.16298018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4931603) q[1];
sx q[1];
rz(-1.5409768) q[1];
sx q[1];
rz(1.0489419) q[1];
rz(-pi) q[2];
x q[2];
rz(1.389723) q[3];
sx q[3];
rz(-1.4317792) q[3];
sx q[3];
rz(0.74336038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(-0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(-0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-3.1380222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9650426) q[0];
sx q[0];
rz(-1.1561484) q[0];
sx q[0];
rz(1.167017) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75242075) q[2];
sx q[2];
rz(-1.8634897) q[2];
sx q[2];
rz(-2.9257286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6047302) q[1];
sx q[1];
rz(-0.68568789) q[1];
sx q[1];
rz(2.0525949) q[1];
rz(-pi) q[2];
rz(-0.076038578) q[3];
sx q[3];
rz(-1.7685316) q[3];
sx q[3];
rz(0.34464371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-2.2214831) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.7238808) q[0];
rz(-0.40813804) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-2.4628941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5355276) q[0];
sx q[0];
rz(-1.3074271) q[0];
sx q[0];
rz(0.33959629) q[0];
rz(-2.3085824) q[2];
sx q[2];
rz(-1.2150803) q[2];
sx q[2];
rz(-2.9290111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39290909) q[1];
sx q[1];
rz(-0.52782413) q[1];
sx q[1];
rz(1.4455568) q[1];
x q[2];
rz(-2.7160866) q[3];
sx q[3];
rz(-1.7153499) q[3];
sx q[3];
rz(-1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1552102) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(0.52694595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53370595) q[0];
sx q[0];
rz(-2.7237646) q[0];
sx q[0];
rz(2.582344) q[0];
x q[1];
rz(1.9138463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(2.3234141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38339112) q[1];
sx q[1];
rz(-1.2609298) q[1];
sx q[1];
rz(-1.2916702) q[1];
rz(-pi) q[2];
rz(-0.44550495) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(1.6642237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(0.79088598) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(-0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0246968) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(-3.1162529) q[0];
x q[1];
rz(1.0656409) q[2];
sx q[2];
rz(-0.2444707) q[2];
sx q[2];
rz(-0.62015647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1211179) q[1];
sx q[1];
rz(-1.7133683) q[1];
sx q[1];
rz(0.74417443) q[1];
x q[2];
rz(0.1694069) q[3];
sx q[3];
rz(-1.7161887) q[3];
sx q[3];
rz(-1.8190847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(1.1311244) q[2];
sx q[2];
rz(-2.9154073) q[2];
sx q[2];
rz(0.67868457) q[2];
rz(0.73137024) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
