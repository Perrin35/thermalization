OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.81808972) q[0];
sx q[0];
rz(-0.40677318) q[0];
sx q[0];
rz(9.9343205) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(6.4767467) q[1];
sx q[1];
rz(11.560796) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5307546) q[0];
sx q[0];
rz(-1.4959529) q[0];
sx q[0];
rz(0.2006528) q[0];
x q[1];
rz(-2.7948954) q[2];
sx q[2];
rz(-0.95333734) q[2];
sx q[2];
rz(0.098010691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7889346) q[1];
sx q[1];
rz(-2.1772258) q[1];
sx q[1];
rz(-0.85373099) q[1];
rz(-0.60723234) q[3];
sx q[3];
rz(-2.2650026) q[3];
sx q[3];
rz(1.9429562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2931508) q[2];
sx q[2];
rz(-1.8895431) q[2];
sx q[2];
rz(-1.1085917) q[2];
rz(2.0075924) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(0.081324287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9950681) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(-2.8296237) q[0];
rz(2.9106855) q[1];
sx q[1];
rz(-1.1038019) q[1];
sx q[1];
rz(1.1280967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19418476) q[0];
sx q[0];
rz(-1.1028521) q[0];
sx q[0];
rz(2.9293961) q[0];
rz(-2.6083469) q[2];
sx q[2];
rz(-0.74138481) q[2];
sx q[2];
rz(2.7864151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.384239) q[1];
sx q[1];
rz(-1.4600672) q[1];
sx q[1];
rz(0.44902451) q[1];
x q[2];
rz(2.3771068) q[3];
sx q[3];
rz(-1.2448881) q[3];
sx q[3];
rz(-0.21316646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4565304) q[2];
sx q[2];
rz(-1.5038467) q[2];
sx q[2];
rz(-3.077363) q[2];
rz(-1.3226994) q[3];
sx q[3];
rz(-2.5048246) q[3];
sx q[3];
rz(2.2548811) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927758) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(-1.2491666) q[0];
rz(-2.8104172) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(1.1963199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6702658) q[0];
sx q[0];
rz(-1.5705137) q[0];
sx q[0];
rz(3.1413881) q[0];
rz(-1.565479) q[2];
sx q[2];
rz(-1.8182767) q[2];
sx q[2];
rz(2.4385045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1041682) q[1];
sx q[1];
rz(-1.0598039) q[1];
sx q[1];
rz(2.4591706) q[1];
x q[2];
rz(2.3881641) q[3];
sx q[3];
rz(-1.2810858) q[3];
sx q[3];
rz(-3.0353607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5195878) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(-0.91274846) q[2];
rz(1.7031472) q[3];
sx q[3];
rz(-1.5104975) q[3];
sx q[3];
rz(1.335817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790633) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(0.67034876) q[0];
rz(2.4743075) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(0.60428062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42603618) q[0];
sx q[0];
rz(-0.73764387) q[0];
sx q[0];
rz(-0.064427388) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4154345) q[2];
sx q[2];
rz(-1.1929907) q[2];
sx q[2];
rz(2.8239259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88612823) q[1];
sx q[1];
rz(-1.3166952) q[1];
sx q[1];
rz(0.21826616) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84008645) q[3];
sx q[3];
rz(-1.4021676) q[3];
sx q[3];
rz(-0.11883277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0303354) q[2];
sx q[2];
rz(-1.5680485) q[2];
sx q[2];
rz(2.7643909) q[2];
rz(-3.0180569) q[3];
sx q[3];
rz(-1.7081407) q[3];
sx q[3];
rz(1.7326573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1971624) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(0.29931983) q[0];
rz(-2.0206644) q[1];
sx q[1];
rz(-1.2052373) q[1];
sx q[1];
rz(0.96251208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045557307) q[0];
sx q[0];
rz(-1.0518414) q[0];
sx q[0];
rz(-0.23764289) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2268874) q[2];
sx q[2];
rz(-1.9491299) q[2];
sx q[2];
rz(-1.8112184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.923063) q[1];
sx q[1];
rz(-2.0835365) q[1];
sx q[1];
rz(-0.96455892) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4351746) q[3];
sx q[3];
rz(-0.73143857) q[3];
sx q[3];
rz(2.4789435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0472497) q[2];
sx q[2];
rz(-2.3361358) q[2];
sx q[2];
rz(0.94364014) q[2];
rz(1.0654248) q[3];
sx q[3];
rz(-1.3506972) q[3];
sx q[3];
rz(2.0431199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074325398) q[0];
sx q[0];
rz(-1.704498) q[0];
sx q[0];
rz(3.1332916) q[0];
rz(2.6240194) q[1];
sx q[1];
rz(-2.4868496) q[1];
sx q[1];
rz(-1.5932721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28866523) q[0];
sx q[0];
rz(-1.2199243) q[0];
sx q[0];
rz(-0.025625833) q[0];
rz(-pi) q[1];
rz(-2.6363434) q[2];
sx q[2];
rz(-1.3224241) q[2];
sx q[2];
rz(-0.58282436) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8928194) q[1];
sx q[1];
rz(-1.237251) q[1];
sx q[1];
rz(1.1897586) q[1];
x q[2];
rz(-2.1287789) q[3];
sx q[3];
rz(-2.4468385) q[3];
sx q[3];
rz(1.6266803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8325309) q[2];
sx q[2];
rz(-1.4893724) q[2];
sx q[2];
rz(1.3365411) q[2];
rz(-0.055770326) q[3];
sx q[3];
rz(-0.99168188) q[3];
sx q[3];
rz(0.43558863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5830773) q[0];
sx q[0];
rz(-0.4774839) q[0];
sx q[0];
rz(-0.68786311) q[0];
rz(-1.4631924) q[1];
sx q[1];
rz(-2.4122767) q[1];
sx q[1];
rz(-2.8942143) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98589135) q[0];
sx q[0];
rz(-1.3257605) q[0];
sx q[0];
rz(-1.3831397) q[0];
rz(-pi) q[1];
rz(0.94893564) q[2];
sx q[2];
rz(-1.0100216) q[2];
sx q[2];
rz(-1.4863297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0717031) q[1];
sx q[1];
rz(-0.99537288) q[1];
sx q[1];
rz(-2.6399122) q[1];
rz(-pi) q[2];
rz(-0.12567802) q[3];
sx q[3];
rz(-1.3397371) q[3];
sx q[3];
rz(2.8383858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6972202) q[2];
sx q[2];
rz(-1.8070544) q[2];
sx q[2];
rz(2.7174301) q[2];
rz(-1.0423202) q[3];
sx q[3];
rz(-0.76516953) q[3];
sx q[3];
rz(0.55646363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384077) q[0];
sx q[0];
rz(-2.1557032) q[0];
sx q[0];
rz(-2.5307122) q[0];
rz(2.8914087) q[1];
sx q[1];
rz(-1.4004204) q[1];
sx q[1];
rz(-2.6649323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395916) q[0];
sx q[0];
rz(-1.6098813) q[0];
sx q[0];
rz(2.0530967) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3872434) q[2];
sx q[2];
rz(-1.6590365) q[2];
sx q[2];
rz(2.1004408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5077153) q[1];
sx q[1];
rz(-0.65836473) q[1];
sx q[1];
rz(0.69067278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1564756) q[3];
sx q[3];
rz(-1.9352311) q[3];
sx q[3];
rz(-1.6062615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1499947) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(-0.67982137) q[2];
rz(2.6796135) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(1.6405039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944549) q[0];
sx q[0];
rz(-1.7663904) q[0];
sx q[0];
rz(2.9875901) q[0];
rz(2.9601861) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(3.0214686) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4171158) q[0];
sx q[0];
rz(-1.6301148) q[0];
sx q[0];
rz(1.1671221) q[0];
rz(-2.207802) q[2];
sx q[2];
rz(-0.32378886) q[2];
sx q[2];
rz(1.1188913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57615818) q[1];
sx q[1];
rz(-1.4729958) q[1];
sx q[1];
rz(-1.3937772) q[1];
rz(2.613853) q[3];
sx q[3];
rz(-1.6956009) q[3];
sx q[3];
rz(0.11858701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4274365) q[2];
sx q[2];
rz(-1.3849881) q[2];
sx q[2];
rz(-2.6825421) q[2];
rz(0.51042026) q[3];
sx q[3];
rz(-0.84158689) q[3];
sx q[3];
rz(0.69971219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089461483) q[0];
sx q[0];
rz(-2.1928146) q[0];
sx q[0];
rz(-2.7556038) q[0];
rz(1.5486859) q[1];
sx q[1];
rz(-2.6451151) q[1];
sx q[1];
rz(-1.5923502) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44819866) q[0];
sx q[0];
rz(-1.5992461) q[0];
sx q[0];
rz(1.4017808) q[0];
rz(-pi) q[1];
rz(2.5086918) q[2];
sx q[2];
rz(-1.0833246) q[2];
sx q[2];
rz(-0.097214708) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1116699) q[1];
sx q[1];
rz(-3.0360465) q[1];
sx q[1];
rz(-1.3131857) q[1];
x q[2];
rz(1.4172961) q[3];
sx q[3];
rz(-1.5046538) q[3];
sx q[3];
rz(-1.9754174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3647032) q[2];
sx q[2];
rz(-2.2025755) q[2];
sx q[2];
rz(-1.4466064) q[2];
rz(-1.0363091) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(3.115263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25630367) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(-0.22008315) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(-0.34192495) q[2];
sx q[2];
rz(-0.76818633) q[2];
sx q[2];
rz(1.1053602) q[2];
rz(-0.55752936) q[3];
sx q[3];
rz(-2.5685608) q[3];
sx q[3];
rz(1.3309042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
