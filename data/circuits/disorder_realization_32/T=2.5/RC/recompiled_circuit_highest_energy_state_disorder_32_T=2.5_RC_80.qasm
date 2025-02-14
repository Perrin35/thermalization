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
rz(0.12953144) q[0];
sx q[0];
rz(-0.13106267) q[0];
sx q[0];
rz(-0.60453209) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(-0.25783917) q[1];
sx q[1];
rz(1.9115619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8800921) q[0];
sx q[0];
rz(-1.3208436) q[0];
sx q[0];
rz(0.34323741) q[0];
rz(-pi) q[1];
rz(0.31754874) q[2];
sx q[2];
rz(-0.63830909) q[2];
sx q[2];
rz(-2.2090863) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2369654) q[1];
sx q[1];
rz(-2.2481027) q[1];
sx q[1];
rz(-0.64579247) q[1];
rz(-pi) q[2];
rz(-0.12768605) q[3];
sx q[3];
rz(-2.0632072) q[3];
sx q[3];
rz(-3.1000053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0543542) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(2.6916091) q[3];
sx q[3];
rz(-0.68322244) q[3];
sx q[3];
rz(0.78342485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1012652) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(2.878317) q[0];
rz(2.1030262) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(0.77478772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2939071) q[0];
sx q[0];
rz(-2.4890702) q[0];
sx q[0];
rz(0.34282617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4500174) q[2];
sx q[2];
rz(-1.0333158) q[2];
sx q[2];
rz(-0.36851685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1549708) q[1];
sx q[1];
rz(-0.38119527) q[1];
sx q[1];
rz(1.4710557) q[1];
rz(-pi) q[2];
rz(-0.70330422) q[3];
sx q[3];
rz(-1.9580132) q[3];
sx q[3];
rz(-2.4037698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7210377) q[2];
sx q[2];
rz(-1.5169531) q[2];
sx q[2];
rz(0.59845412) q[2];
rz(1.7789486) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(2.8708598) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8552928) q[0];
sx q[0];
rz(-0.4628276) q[0];
sx q[0];
rz(1.5077952) q[0];
rz(0.15448054) q[1];
sx q[1];
rz(-1.4198317) q[1];
sx q[1];
rz(-2.3522164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76217795) q[0];
sx q[0];
rz(-1.6573849) q[0];
sx q[0];
rz(0.53126727) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28114281) q[2];
sx q[2];
rz(-1.707555) q[2];
sx q[2];
rz(3.1355372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7676413) q[1];
sx q[1];
rz(-2.1726296) q[1];
sx q[1];
rz(-0.37202073) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79929914) q[3];
sx q[3];
rz(-1.2033312) q[3];
sx q[3];
rz(0.41900837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7827451) q[2];
sx q[2];
rz(-0.91155702) q[2];
sx q[2];
rz(-1.4770799) q[2];
rz(-1.2137671) q[3];
sx q[3];
rz(-1.341308) q[3];
sx q[3];
rz(-1.7267797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724991) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(2.5308727) q[0];
rz(-0.67131132) q[1];
sx q[1];
rz(-1.6339615) q[1];
sx q[1];
rz(1.7866561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31025654) q[0];
sx q[0];
rz(-0.68643565) q[0];
sx q[0];
rz(-2.3879334) q[0];
rz(-1.2569095) q[2];
sx q[2];
rz(-2.7531024) q[2];
sx q[2];
rz(1.1780648) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23933168) q[1];
sx q[1];
rz(-0.25333764) q[1];
sx q[1];
rz(-1.317666) q[1];
rz(-pi) q[2];
rz(1.7132961) q[3];
sx q[3];
rz(-0.5915407) q[3];
sx q[3];
rz(1.174389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2063724) q[2];
sx q[2];
rz(-0.78038961) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(2.3937461) q[3];
sx q[3];
rz(-2.3220389) q[3];
sx q[3];
rz(-1.6929172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95626962) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(-2.4638033) q[0];
rz(2.9970844) q[1];
sx q[1];
rz(-0.78737193) q[1];
sx q[1];
rz(1.4871303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6078454) q[0];
sx q[0];
rz(-3.0032039) q[0];
sx q[0];
rz(-1.7615303) q[0];
x q[1];
rz(-1.1688091) q[2];
sx q[2];
rz(-1.6239407) q[2];
sx q[2];
rz(-1.0099908) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2320391) q[1];
sx q[1];
rz(-2.0318446) q[1];
sx q[1];
rz(2.1559245) q[1];
rz(-pi) q[2];
rz(-1.4492338) q[3];
sx q[3];
rz(-0.37827493) q[3];
sx q[3];
rz(-0.88751436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80374485) q[2];
sx q[2];
rz(-2.8779112) q[2];
sx q[2];
rz(-2.7860876) q[2];
rz(1.0454987) q[3];
sx q[3];
rz(-1.9020566) q[3];
sx q[3];
rz(-0.56226468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1657555) q[0];
sx q[0];
rz(-0.58335692) q[0];
sx q[0];
rz(-3.052886) q[0];
rz(-1.6070131) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(3.065899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3189118) q[0];
sx q[0];
rz(-1.3782129) q[0];
sx q[0];
rz(0.67355021) q[0];
x q[1];
rz(-1.6275109) q[2];
sx q[2];
rz(-2.4091588) q[2];
sx q[2];
rz(-0.54679856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9348889) q[1];
sx q[1];
rz(-2.5674501) q[1];
sx q[1];
rz(-1.8341792) q[1];
rz(0.47786062) q[3];
sx q[3];
rz(-0.63242542) q[3];
sx q[3];
rz(2.6101108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82464108) q[2];
sx q[2];
rz(-1.757916) q[2];
sx q[2];
rz(-2.1023882) q[2];
rz(-0.21229395) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067326389) q[0];
sx q[0];
rz(-0.75263158) q[0];
sx q[0];
rz(2.0972032) q[0];
rz(2.3692865) q[1];
sx q[1];
rz(-0.65985313) q[1];
sx q[1];
rz(-0.73371249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70673101) q[0];
sx q[0];
rz(-1.0821663) q[0];
sx q[0];
rz(0.74669331) q[0];
rz(-2.9915221) q[2];
sx q[2];
rz(-1.9378621) q[2];
sx q[2];
rz(2.31524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69960475) q[1];
sx q[1];
rz(-1.2803704) q[1];
sx q[1];
rz(-1.6922349) q[1];
rz(-pi) q[2];
rz(0.10802631) q[3];
sx q[3];
rz(-2.3026534) q[3];
sx q[3];
rz(-2.9895003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5250728) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(-0.60866848) q[2];
rz(1.0673374) q[3];
sx q[3];
rz(-0.87456861) q[3];
sx q[3];
rz(-2.5909891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6390425) q[0];
sx q[0];
rz(-0.35824963) q[0];
sx q[0];
rz(-1.6453561) q[0];
rz(1.7773588) q[1];
sx q[1];
rz(-2.5271006) q[1];
sx q[1];
rz(0.38633698) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7832121) q[0];
sx q[0];
rz(-2.2370501) q[0];
sx q[0];
rz(-1.9080286) q[0];
rz(1.3851829) q[2];
sx q[2];
rz(-1.6395373) q[2];
sx q[2];
rz(-0.5663213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1586502) q[1];
sx q[1];
rz(-1.6931931) q[1];
sx q[1];
rz(-0.33315181) q[1];
x q[2];
rz(1.531485) q[3];
sx q[3];
rz(-0.98061258) q[3];
sx q[3];
rz(3.085768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51314696) q[2];
sx q[2];
rz(-1.7690423) q[2];
sx q[2];
rz(0.15303843) q[2];
rz(-2.2506574) q[3];
sx q[3];
rz(-2.1864083) q[3];
sx q[3];
rz(-2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1560169) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(2.7942221) q[0];
rz(-3.103745) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(-0.52245021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7670531) q[0];
sx q[0];
rz(-1.660083) q[0];
sx q[0];
rz(-1.1566628) q[0];
x q[1];
rz(-0.00084472617) q[2];
sx q[2];
rz(-0.36549308) q[2];
sx q[2];
rz(-1.5216684) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2938483) q[1];
sx q[1];
rz(-1.9894682) q[1];
sx q[1];
rz(-2.5878504) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1303923) q[3];
sx q[3];
rz(-1.2204477) q[3];
sx q[3];
rz(-0.94061414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6309506) q[2];
sx q[2];
rz(-2.6783671) q[2];
sx q[2];
rz(1.9412712) q[2];
rz(0.55780324) q[3];
sx q[3];
rz(-1.6226945) q[3];
sx q[3];
rz(2.7571078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2931622) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(-1.4137319) q[0];
rz(-3.0005455) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(-2.486855) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075874627) q[0];
sx q[0];
rz(-1.7786586) q[0];
sx q[0];
rz(-1.8130568) q[0];
rz(-2.5467196) q[2];
sx q[2];
rz(-2.0030177) q[2];
sx q[2];
rz(3.1405666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.480994) q[1];
sx q[1];
rz(-2.8428322) q[1];
sx q[1];
rz(0.95282747) q[1];
rz(0.61956866) q[3];
sx q[3];
rz(-1.1244785) q[3];
sx q[3];
rz(0.60947641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7901788) q[2];
sx q[2];
rz(-1.1620136) q[2];
sx q[2];
rz(0.72511017) q[2];
rz(2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17726041) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(-2.019885) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(2.4779392) q[2];
sx q[2];
rz(-1.2812231) q[2];
sx q[2];
rz(-0.44322586) q[2];
rz(0.18425758) q[3];
sx q[3];
rz(-0.72132106) q[3];
sx q[3];
rz(-3.0527243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
