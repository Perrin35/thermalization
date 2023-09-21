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
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(2.9266761) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.874274) q[0];
sx q[0];
rz(-1.8288757) q[0];
sx q[0];
rz(-1.2866856) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29834892) q[2];
sx q[2];
rz(-2.6539408) q[2];
sx q[2];
rz(1.2717441) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6282181) q[1];
sx q[1];
rz(-0.77925032) q[1];
sx q[1];
rz(-0.90374225) q[1];
rz(-2.05902) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.2692497) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.6033283) q[1];
sx q[1];
rz(-0.89675084) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7433886) q[0];
sx q[0];
rz(-2.2687015) q[0];
sx q[0];
rz(-2.8845805) q[0];
rz(-pi) q[1];
rz(2.9727544) q[2];
sx q[2];
rz(-1.642792) q[2];
sx q[2];
rz(0.26693401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36205081) q[1];
sx q[1];
rz(-2.3565787) q[1];
sx q[1];
rz(-1.7464459) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1284626) q[3];
sx q[3];
rz(-0.80207523) q[3];
sx q[3];
rz(-0.80004063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(1.0401475) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(2.7064586) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2745278) q[0];
sx q[0];
rz(-1.2194249) q[0];
sx q[0];
rz(-1.685313) q[0];
rz(-pi) q[1];
rz(1.4235731) q[2];
sx q[2];
rz(-1.3230811) q[2];
sx q[2];
rz(-1.9927646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16776925) q[1];
sx q[1];
rz(-0.5127738) q[1];
sx q[1];
rz(-1.9085401) q[1];
x q[2];
rz(2.7828091) q[3];
sx q[3];
rz(-0.79625087) q[3];
sx q[3];
rz(0.96804726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-2.8386774) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(2.8767169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47539513) q[0];
sx q[0];
rz(-0.65984939) q[0];
sx q[0];
rz(-3.092479) q[0];
rz(-pi) q[1];
rz(2.8633966) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(0.96565914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(-1.0290531) q[1];
rz(0.78062765) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.2300223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.778487) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(2.662861) q[0];
rz(-2.1084673) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(0.95265257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081776) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(-1.7081225) q[0];
x q[1];
rz(2.7258337) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(1.8005467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9597185) q[1];
sx q[1];
rz(-1.3337143) q[1];
sx q[1];
rz(1.9738166) q[1];
rz(-pi) q[2];
rz(-1.5289375) q[3];
sx q[3];
rz(-1.0930982) q[3];
sx q[3];
rz(0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1369551) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(3.1266771) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8820011) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(-1.7829347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8494106) q[1];
sx q[1];
rz(-1.9841571) q[1];
sx q[1];
rz(0.11489111) q[1];
rz(-pi) q[2];
rz(-0.63038007) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(-2.3849934) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7998357) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(2.1205649) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7724178) q[2];
sx q[2];
rz(-1.9747509) q[2];
sx q[2];
rz(2.902365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36180624) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(0.96868412) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19161253) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(-2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-2.1527122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38935223) q[0];
sx q[0];
rz(-1.9703431) q[0];
sx q[0];
rz(2.7826392) q[0];
rz(-pi) q[1];
rz(3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(1.2707368) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58395308) q[1];
sx q[1];
rz(-1.156731) q[1];
sx q[1];
rz(-3.0686892) q[1];
x q[2];
rz(0.51560651) q[3];
sx q[3];
rz(-2.0022087) q[3];
sx q[3];
rz(3.072217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(-0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-0.87337714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3082335) q[0];
sx q[0];
rz(-1.5033659) q[0];
sx q[0];
rz(1.391135) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4344425) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(-1.855195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5492591) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(1.1997644) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7018413) q[3];
sx q[3];
rz(-1.0381178) q[3];
sx q[3];
rz(-2.7524878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1197027) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-2.9746829) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.3964765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5293225) q[0];
sx q[0];
rz(-1.6565874) q[0];
sx q[0];
rz(-2.8180608) q[0];
rz(-pi) q[1];
rz(-0.65812494) q[2];
sx q[2];
rz(-1.3326367) q[2];
sx q[2];
rz(2.4326774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6203306) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(-1.2465338) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20874899) q[3];
sx q[3];
rz(-0.26780805) q[3];
sx q[3];
rz(2.1684614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(2.9818025) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(-2.2014387) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(1.521048) q[3];
sx q[3];
rz(-2.1038901) q[3];
sx q[3];
rz(-0.62705561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
