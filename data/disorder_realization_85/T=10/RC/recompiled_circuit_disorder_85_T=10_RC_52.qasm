OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24297548) q[0];
sx q[0];
rz(-2.3529422) q[0];
sx q[0];
rz(-2.2126161) q[0];
rz(0.55144989) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(1.431682) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7264759) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(2.9813558) q[1];
rz(-pi) q[2];
rz(-1.6269496) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(-2.5889531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(-2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-2.4172799) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9162826) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(-1.1067953) q[0];
rz(-0.19947796) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(-2.2897838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1867379) q[1];
sx q[1];
rz(-1.2192982) q[1];
sx q[1];
rz(-0.21912205) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23726666) q[3];
sx q[3];
rz(-1.3788584) q[3];
sx q[3];
rz(2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.3954337) q[0];
rz(-0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0116545) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(-2.8610693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29580446) q[2];
sx q[2];
rz(-0.89106262) q[2];
sx q[2];
rz(0.55065292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0137274) q[1];
sx q[1];
rz(-1.3294819) q[1];
sx q[1];
rz(1.024854) q[1];
rz(-pi) q[2];
rz(-3.1261256) q[3];
sx q[3];
rz(-1.3028212) q[3];
sx q[3];
rz(2.1249352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(-1.4455618) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-0.82675654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402156) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(2.4504689) q[0];
rz(-0.13917285) q[2];
sx q[2];
rz(-0.77251245) q[2];
sx q[2];
rz(-3.006209) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0485059) q[1];
sx q[1];
rz(-0.48492453) q[1];
sx q[1];
rz(-2.4946458) q[1];
rz(-pi) q[2];
rz(-2.4171962) q[3];
sx q[3];
rz(-1.3339692) q[3];
sx q[3];
rz(-2.608992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(-3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(0.28516969) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171705) q[0];
sx q[0];
rz(-1.6132014) q[0];
sx q[0];
rz(2.0457343) q[0];
rz(-pi) q[1];
rz(-0.4816149) q[2];
sx q[2];
rz(-0.72313213) q[2];
sx q[2];
rz(-2.8358104) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4080216) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(-1.8838521) q[1];
x q[2];
rz(-0.22484803) q[3];
sx q[3];
rz(-0.41089155) q[3];
sx q[3];
rz(1.2332682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(2.9845797) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(1.7745811) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088928662) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(-2.0220387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81142601) q[2];
sx q[2];
rz(-1.1277414) q[2];
sx q[2];
rz(-1.7048938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(1.8606436) q[1];
x q[2];
rz(-1.6744162) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(-0.87408376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8453025) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(2.2677299) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.9708995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6252977) q[0];
sx q[0];
rz(-1.5646311) q[0];
sx q[0];
rz(3.1185634) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5198176) q[2];
sx q[2];
rz(-1.4547252) q[2];
sx q[2];
rz(2.799084) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2611744) q[1];
sx q[1];
rz(-0.92217991) q[1];
sx q[1];
rz(-0.76678126) q[1];
x q[2];
rz(1.2675769) q[3];
sx q[3];
rz(-1.056864) q[3];
sx q[3];
rz(-0.60790387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532115) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(-1.3636916) q[0];
rz(-pi) q[1];
rz(-0.32142873) q[2];
sx q[2];
rz(-0.49389631) q[2];
sx q[2];
rz(2.0733881) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99107332) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(1.9626161) q[1];
rz(-pi) q[2];
rz(-0.93692245) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(2.6867552) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(-0.79750693) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-0.10841766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2415376) q[0];
sx q[0];
rz(-2.2336322) q[0];
sx q[0];
rz(3.0682949) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7494781) q[2];
sx q[2];
rz(-1.5249426) q[2];
sx q[2];
rz(2.1669441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11549599) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(0.91067578) q[1];
x q[2];
rz(-0.18687825) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0124399) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(-0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-0.055158786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178135) q[0];
sx q[0];
rz(-1.9473416) q[0];
sx q[0];
rz(1.5012653) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55969413) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(0.44911227) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(2.39141) q[1];
rz(-1.6454562) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(-0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(1.6677042) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(2.9604838) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
