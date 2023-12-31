OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(2.7603005) q[1];
sx q[1];
rz(-2.5420904) q[1];
sx q[1];
rz(-1.376027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3289514) q[0];
sx q[0];
rz(-1.1321804) q[0];
sx q[0];
rz(-2.2493275) q[0];
rz(-pi) q[1];
rz(-0.55144989) q[2];
sx q[2];
rz(-2.3460238) q[2];
sx q[2];
rz(-1.7099107) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4151167) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(-0.16023689) q[1];
rz(-pi) q[2];
rz(-1.514643) q[3];
sx q[3];
rz(-0.75222844) q[3];
sx q[3];
rz(-2.5889531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(-2.3108216) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-2.4172799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3121719) q[0];
sx q[0];
rz(-0.90733007) q[0];
sx q[0];
rz(-2.7396766) q[0];
x q[1];
rz(-0.19947796) q[2];
sx q[2];
rz(-1.4885159) q[2];
sx q[2];
rz(2.2897838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.449125) q[1];
sx q[1];
rz(-1.7763224) q[1];
sx q[1];
rz(1.2114026) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69085391) q[3];
sx q[3];
rz(-0.30403954) q[3];
sx q[3];
rz(-0.090304852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90536845) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419462) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(-0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0116545) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(-0.28052335) q[0];
x q[1];
rz(0.8692603) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(-1.9321835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(-2.0134258) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.838802) q[3];
sx q[3];
rz(-1.5857113) q[3];
sx q[3];
rz(-0.55004317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22359426) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0543538) q[0];
sx q[0];
rz(-2.0728489) q[0];
sx q[0];
rz(-2.0421844) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76765676) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(1.3354288) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0485059) q[1];
sx q[1];
rz(-0.48492453) q[1];
sx q[1];
rz(-0.64694689) q[1];
rz(-pi) q[2];
rz(-0.34927807) q[3];
sx q[3];
rz(-0.75540245) q[3];
sx q[3];
rz(-2.3625771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(1.7626804) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681686) q[0];
sx q[0];
rz(-2.0452721) q[0];
sx q[0];
rz(0.04767496) q[0];
rz(-1.1826913) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(-2.2270122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8385182) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(1.5721553) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4739591) q[3];
sx q[3];
rz(-1.9707465) q[3];
sx q[3];
rz(1.6638343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(0.30682492) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(2.9845797) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.7745811) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54129823) q[0];
sx q[0];
rz(-1.5026662) q[0];
sx q[0];
rz(0.033072254) q[0];
x q[1];
rz(-2.3301666) q[2];
sx q[2];
rz(-2.0138513) q[2];
sx q[2];
rz(-1.7048938) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94783084) q[1];
sx q[1];
rz(-1.281922) q[1];
sx q[1];
rz(3.0572901) q[1];
x q[2];
rz(1.4671765) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8453025) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(-0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(-1.9708995) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6252977) q[0];
sx q[0];
rz(-1.5646311) q[0];
sx q[0];
rz(3.1185634) q[0];
rz(3.0253719) q[2];
sx q[2];
rz(-1.6214317) q[2];
sx q[2];
rz(1.9192139) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88041828) q[1];
sx q[1];
rz(-2.2194127) q[1];
sx q[1];
rz(-0.76678126) q[1];
rz(-pi) q[2];
rz(1.2675769) q[3];
sx q[3];
rz(-2.0847287) q[3];
sx q[3];
rz(0.60790387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(-2.4954605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1387678) q[0];
sx q[0];
rz(-1.3647623) q[0];
sx q[0];
rz(-3.0384008) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4023151) q[2];
sx q[2];
rz(-1.104276) q[2];
sx q[2];
rz(-0.70665765) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38478002) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(2.6095819) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93692245) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-0.10841766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71589564) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(-2.2349368) q[0];
rz(-pi) q[1];
rz(-1.7494781) q[2];
sx q[2];
rz(-1.5249426) q[2];
sx q[2];
rz(0.97464857) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81356955) q[1];
sx q[1];
rz(-1.7377825) q[1];
sx q[1];
rz(1.7896673) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7897723) q[3];
sx q[3];
rz(-0.7162381) q[3];
sx q[3];
rz(0.66974528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(2.4160014) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(0.055158786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.178135) q[0];
sx q[0];
rz(-1.194251) q[0];
sx q[0];
rz(1.5012653) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.360838) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(-1.9090261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-2.2844237) q[1];
sx q[1];
rz(0.75018261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6454562) q[3];
sx q[3];
rz(-0.45711043) q[3];
sx q[3];
rz(2.8198869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(2.9329119) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(1.4738884) q[2];
sx q[2];
rz(-1.2599535) q[2];
sx q[2];
rz(-2.8355666) q[2];
rz(1.6155852) q[3];
sx q[3];
rz(-1.8106034) q[3];
sx q[3];
rz(-2.8869224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
