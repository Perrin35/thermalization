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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(4.876457) q[1];
sx q[1];
rz(9.193037) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2278175) q[0];
sx q[0];
rz(-1.639606) q[0];
sx q[0];
rz(-0.34389307) q[0];
rz(-3.0475869) q[2];
sx q[2];
rz(-2.5931907) q[2];
sx q[2];
rz(-1.7552055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8084131) q[1];
sx q[1];
rz(-1.9486548) q[1];
sx q[1];
rz(-2.4757132) q[1];
x q[2];
rz(2.8032325) q[3];
sx q[3];
rz(-1.6099811) q[3];
sx q[3];
rz(0.071166347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.090652466) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(-1.862662) q[2];
rz(2.2517962) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(-2.8029627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5559674) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(0.92581785) q[0];
rz(1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(0.085478641) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8166072) q[0];
sx q[0];
rz(-2.7517635) q[0];
sx q[0];
rz(-3.1329324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7932554) q[2];
sx q[2];
rz(-2.1536963) q[2];
sx q[2];
rz(-0.23886853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2457471) q[1];
sx q[1];
rz(-1.8008044) q[1];
sx q[1];
rz(-1.0279015) q[1];
x q[2];
rz(2.5898629) q[3];
sx q[3];
rz(-0.76653102) q[3];
sx q[3];
rz(0.54099247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(-1.9246842) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81095186) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(2.6445342) q[0];
rz(-2.0992384) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(-0.29464468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9884856) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(-1.4979304) q[0];
rz(-pi) q[1];
rz(-2.9109086) q[2];
sx q[2];
rz(-0.65197021) q[2];
sx q[2];
rz(2.1459652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3556488) q[1];
sx q[1];
rz(-2.3501767) q[1];
sx q[1];
rz(-0.9407465) q[1];
x q[2];
rz(2.9975843) q[3];
sx q[3];
rz(-2.1668808) q[3];
sx q[3];
rz(2.2726187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7872539) q[2];
sx q[2];
rz(-0.86091176) q[2];
sx q[2];
rz(-1.5617237) q[2];
rz(2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(-0.92897433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314986) q[0];
sx q[0];
rz(-1.6153233) q[0];
sx q[0];
rz(0.22931799) q[0];
rz(-1.5062821) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(-3.07952) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631043) q[0];
sx q[0];
rz(-1.0764499) q[0];
sx q[0];
rz(0.93017471) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8624195) q[2];
sx q[2];
rz(-1.9491157) q[2];
sx q[2];
rz(-2.407848) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5531299) q[1];
sx q[1];
rz(-2.3847918) q[1];
sx q[1];
rz(-2.7314145) q[1];
x q[2];
rz(-0.43020474) q[3];
sx q[3];
rz(-0.22048002) q[3];
sx q[3];
rz(1.9752432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.092992358) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(1.830706) q[2];
rz(-0.038289573) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(-0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467317) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(-0.89299655) q[0];
rz(-2.112174) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(3.0457048) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9150118) q[0];
sx q[0];
rz(-1.1553197) q[0];
sx q[0];
rz(1.0092495) q[0];
x q[1];
rz(-2.2402359) q[2];
sx q[2];
rz(-3.0283805) q[2];
sx q[2];
rz(0.33153807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3600625) q[1];
sx q[1];
rz(-1.2075829) q[1];
sx q[1];
rz(0.59477721) q[1];
x q[2];
rz(-2.5378102) q[3];
sx q[3];
rz(-0.60294916) q[3];
sx q[3];
rz(-0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9868077) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(0.5031684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7285889) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(3.0244306) q[0];
rz(1.4940184) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(1.0711627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18497047) q[0];
sx q[0];
rz(-1.8529467) q[0];
sx q[0];
rz(2.6376702) q[0];
rz(-pi) q[1];
rz(0.44354673) q[2];
sx q[2];
rz(-0.58813349) q[2];
sx q[2];
rz(-3.0556222) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2342959) q[1];
sx q[1];
rz(-1.3456555) q[1];
sx q[1];
rz(2.2683558) q[1];
x q[2];
rz(-2.6758606) q[3];
sx q[3];
rz(-0.35251401) q[3];
sx q[3];
rz(-0.62437144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57227197) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(-0.039483698) q[2];
rz(0.076400541) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(-2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(1.9901265) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(-1.3075525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9614544) q[0];
sx q[0];
rz(-0.39648065) q[0];
sx q[0];
rz(-0.97361083) q[0];
rz(-0.37613873) q[2];
sx q[2];
rz(-1.6385632) q[2];
sx q[2];
rz(-2.4532831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7286018) q[1];
sx q[1];
rz(-2.6003693) q[1];
sx q[1];
rz(1.9899998) q[1];
rz(0.96486196) q[3];
sx q[3];
rz(-1.4601225) q[3];
sx q[3];
rz(-2.0697274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(-0.55366984) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.6863916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.11947908) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(2.8346862) q[0];
rz(1.8652929) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(1.2409522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6825323) q[0];
sx q[0];
rz(-1.1272361) q[0];
sx q[0];
rz(0.17291594) q[0];
rz(1.9996793) q[2];
sx q[2];
rz(-2.1253573) q[2];
sx q[2];
rz(2.6311324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13682374) q[1];
sx q[1];
rz(-1.3826177) q[1];
sx q[1];
rz(0.31463639) q[1];
rz(-pi) q[2];
rz(-2.5316174) q[3];
sx q[3];
rz(-1.2362567) q[3];
sx q[3];
rz(-1.1776678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3500195) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(1.2903068) q[2];
rz(0.6811412) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(-1.5017989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27555585) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(-1.1736897) q[0];
rz(2.5545919) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(3.0618844) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3143334) q[0];
sx q[0];
rz(-1.200935) q[0];
sx q[0];
rz(-0.24263675) q[0];
rz(-2.9201512) q[2];
sx q[2];
rz(-1.4613073) q[2];
sx q[2];
rz(1.8851282) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.798493) q[1];
sx q[1];
rz(-1.0086372) q[1];
sx q[1];
rz(2.4141099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3354638) q[3];
sx q[3];
rz(-1.3622096) q[3];
sx q[3];
rz(-0.94193469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6012663) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(-1.8076753) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48225668) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.851086) q[0];
rz(-0.036529649) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8528906) q[0];
sx q[0];
rz(-1.1452617) q[0];
sx q[0];
rz(0.3279) q[0];
x q[1];
rz(2.8743083) q[2];
sx q[2];
rz(-1.2517831) q[2];
sx q[2];
rz(2.6397982) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.579291) q[1];
sx q[1];
rz(-2.3962483) q[1];
sx q[1];
rz(-1.7991123) q[1];
x q[2];
rz(-2.815991) q[3];
sx q[3];
rz(-2.7956617) q[3];
sx q[3];
rz(-0.17445645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(-1.3247789) q[2];
rz(-2.3850208) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751547) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(2.6773166) q[2];
sx q[2];
rz(-0.79002476) q[2];
sx q[2];
rz(2.644047) q[2];
rz(0.57831709) q[3];
sx q[3];
rz(-0.74902799) q[3];
sx q[3];
rz(2.0109162) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
