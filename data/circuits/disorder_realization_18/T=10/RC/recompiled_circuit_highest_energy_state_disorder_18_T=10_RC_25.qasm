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
rz(1.7914766) q[0];
sx q[0];
rz(-0.67576367) q[0];
sx q[0];
rz(3.1337373) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(1.6176728) q[1];
sx q[1];
rz(9.6674506) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382671) q[0];
sx q[0];
rz(-1.6969862) q[0];
sx q[0];
rz(-0.85889205) q[0];
rz(-1.2062293) q[2];
sx q[2];
rz(-1.1375712) q[2];
sx q[2];
rz(-1.792576) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3625177) q[1];
sx q[1];
rz(-1.9240161) q[1];
sx q[1];
rz(-1.5455957) q[1];
rz(0.66565973) q[3];
sx q[3];
rz(-1.375631) q[3];
sx q[3];
rz(1.631537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(-0.71199065) q[2];
rz(-0.30501929) q[3];
sx q[3];
rz(-2.9192393) q[3];
sx q[3];
rz(-3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63118339) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(1.7741868) q[0];
rz(-3.101688) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(0.59919277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0829701) q[0];
sx q[0];
rz(-0.48958594) q[0];
sx q[0];
rz(-2.8955196) q[0];
rz(0.56053253) q[2];
sx q[2];
rz(-2.2403702) q[2];
sx q[2];
rz(0.59434429) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.42223052) q[1];
sx q[1];
rz(-1.9047202) q[1];
sx q[1];
rz(0.67113282) q[1];
x q[2];
rz(0.33984025) q[3];
sx q[3];
rz(-0.71469864) q[3];
sx q[3];
rz(-0.13091892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0590608) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(1.833029) q[2];
rz(0.023905309) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5628975) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4557274) q[0];
sx q[0];
rz(-0.66615921) q[0];
sx q[0];
rz(0.84079963) q[0];
rz(-2.1127545) q[1];
sx q[1];
rz(-0.40887555) q[1];
sx q[1];
rz(1.4604481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8864635) q[0];
sx q[0];
rz(-2.6570519) q[0];
sx q[0];
rz(0.31130975) q[0];
x q[1];
rz(-0.47346327) q[2];
sx q[2];
rz(-1.9046648) q[2];
sx q[2];
rz(0.75227458) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3884133) q[1];
sx q[1];
rz(-2.3824661) q[1];
sx q[1];
rz(2.6469346) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6275241) q[3];
sx q[3];
rz(-0.9642082) q[3];
sx q[3];
rz(-1.1067672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8847522) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(2.9849226) q[2];
rz(-1.5001851) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0665862) q[0];
sx q[0];
rz(-1.2970507) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(-1.1219885) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(-1.5544308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.156896) q[0];
sx q[0];
rz(-1.5050355) q[0];
sx q[0];
rz(-2.3143011) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6115613) q[2];
sx q[2];
rz(-1.4546548) q[2];
sx q[2];
rz(0.16414205) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5579368) q[1];
sx q[1];
rz(-0.93975818) q[1];
sx q[1];
rz(-2.2028752) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46458475) q[3];
sx q[3];
rz(-0.66028336) q[3];
sx q[3];
rz(0.44251501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2598205) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(2.1176977) q[2];
rz(-0.77670589) q[3];
sx q[3];
rz(-1.1506162) q[3];
sx q[3];
rz(2.1385433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7316932) q[0];
sx q[0];
rz(-0.85407805) q[0];
sx q[0];
rz(-2.5140629) q[0];
rz(-2.4420786) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(-0.90739179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32226792) q[0];
sx q[0];
rz(-1.5307679) q[0];
sx q[0];
rz(2.3884474) q[0];
rz(1.4110231) q[2];
sx q[2];
rz(-1.7587307) q[2];
sx q[2];
rz(-0.86216506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4470565) q[1];
sx q[1];
rz(-1.9354104) q[1];
sx q[1];
rz(1.0178119) q[1];
rz(0.058919546) q[3];
sx q[3];
rz(-1.0353966) q[3];
sx q[3];
rz(2.8297092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0044535) q[2];
sx q[2];
rz(-1.2848102) q[2];
sx q[2];
rz(0.84490204) q[2];
rz(-2.4431303) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(0.79160488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38248211) q[0];
sx q[0];
rz(-1.6805205) q[0];
sx q[0];
rz(-1.2680898) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(2.7472034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1295107) q[0];
sx q[0];
rz(-2.2144268) q[0];
sx q[0];
rz(-0.056931007) q[0];
rz(-pi) q[1];
rz(-2.4066448) q[2];
sx q[2];
rz(-1.6326666) q[2];
sx q[2];
rz(2.8501373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2339237) q[1];
sx q[1];
rz(-0.65429293) q[1];
sx q[1];
rz(-3.0680502) q[1];
x q[2];
rz(1.6485571) q[3];
sx q[3];
rz(-1.9304515) q[3];
sx q[3];
rz(-0.80870562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4322728) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(-2.0105441) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.5442218) q[0];
sx q[0];
rz(-2.9591296) q[0];
sx q[0];
rz(2.9083948) q[0];
rz(-0.11416642) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(-0.17359576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1602992) q[0];
sx q[0];
rz(-0.7795142) q[0];
sx q[0];
rz(2.1008089) q[0];
rz(-2.7596682) q[2];
sx q[2];
rz(-1.1537103) q[2];
sx q[2];
rz(-2.797319) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0376589) q[1];
sx q[1];
rz(-1.3500431) q[1];
sx q[1];
rz(-0.16640618) q[1];
rz(2.9529157) q[3];
sx q[3];
rz(-1.5828307) q[3];
sx q[3];
rz(-2.6493612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1058098) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(-2.2868273) q[2];
rz(3.0586976) q[3];
sx q[3];
rz(-1.1707183) q[3];
sx q[3];
rz(0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4891124) q[0];
sx q[0];
rz(-1.880045) q[0];
sx q[0];
rz(-1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.2636355) q[1];
sx q[1];
rz(2.9875535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48705593) q[0];
sx q[0];
rz(-0.77870071) q[0];
sx q[0];
rz(2.6683183) q[0];
x q[1];
rz(2.3588347) q[2];
sx q[2];
rz(-2.2680757) q[2];
sx q[2];
rz(3.1408666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5402962) q[1];
sx q[1];
rz(-2.5754693) q[1];
sx q[1];
rz(-2.152312) q[1];
rz(-0.2897183) q[3];
sx q[3];
rz(-1.8236092) q[3];
sx q[3];
rz(-0.022965206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4607294) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(2.4737849) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(-0.63647979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77968303) q[0];
sx q[0];
rz(-1.4754262) q[0];
sx q[0];
rz(0.59378004) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(1.2978172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0827328) q[0];
sx q[0];
rz(-1.642859) q[0];
sx q[0];
rz(-0.58499344) q[0];
x q[1];
rz(0.4164575) q[2];
sx q[2];
rz(-2.4707332) q[2];
sx q[2];
rz(2.351298) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62420708) q[1];
sx q[1];
rz(-2.7605857) q[1];
sx q[1];
rz(0.34492774) q[1];
rz(-pi) q[2];
rz(-2.3343349) q[3];
sx q[3];
rz(-1.689925) q[3];
sx q[3];
rz(1.8594683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7456776) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(-1.5094666) q[2];
rz(0.11219003) q[3];
sx q[3];
rz(-2.1144919) q[3];
sx q[3];
rz(-1.3794911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701732) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(-0.26563409) q[0];
rz(0.98948014) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(2.8094453) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0743389) q[0];
sx q[0];
rz(-1.5746207) q[0];
sx q[0];
rz(-1.2373562) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1411601) q[2];
sx q[2];
rz(-2.5186335) q[2];
sx q[2];
rz(2.0596383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89236081) q[1];
sx q[1];
rz(-0.3438102) q[1];
sx q[1];
rz(1.7147786) q[1];
rz(-0.048236851) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0997448) q[2];
sx q[2];
rz(-2.0691278) q[2];
sx q[2];
rz(-2.9912046) q[2];
rz(2.2930875) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0583508) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(-1.3311483) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(-1.4990357) q[2];
sx q[2];
rz(-2.4007779) q[2];
sx q[2];
rz(3.130198) q[2];
rz(-0.030134044) q[3];
sx q[3];
rz(-2.4345955) q[3];
sx q[3];
rz(1.101936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
