OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.30224213) q[0];
sx q[0];
rz(-1.2737162) q[0];
sx q[0];
rz(-2.1923375) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(2.1646808) q[1];
sx q[1];
rz(10.77471) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5744517) q[0];
sx q[0];
rz(-2.0724839) q[0];
sx q[0];
rz(1.4760963) q[0];
rz(2.737713) q[2];
sx q[2];
rz(-0.35940659) q[2];
sx q[2];
rz(-3.0702116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9266859) q[1];
sx q[1];
rz(-1.6388571) q[1];
sx q[1];
rz(-2.7665274) q[1];
rz(-pi) q[2];
rz(2.0032521) q[3];
sx q[3];
rz(-1.8205595) q[3];
sx q[3];
rz(0.43507155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.533941) q[2];
sx q[2];
rz(-2.0484296) q[2];
sx q[2];
rz(-0.1926113) q[2];
rz(-1.2827778) q[3];
sx q[3];
rz(-2.0359813) q[3];
sx q[3];
rz(-2.5530596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31805661) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(-2.4349924) q[0];
rz(2.508714) q[1];
sx q[1];
rz(-2.0426079) q[1];
sx q[1];
rz(2.852827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59161883) q[0];
sx q[0];
rz(-1.5697877) q[0];
sx q[0];
rz(0.002480025) q[0];
rz(2.9339004) q[2];
sx q[2];
rz(-1.4121659) q[2];
sx q[2];
rz(-2.6191448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0763467) q[1];
sx q[1];
rz(-1.7829121) q[1];
sx q[1];
rz(0.29233039) q[1];
rz(2.3815126) q[3];
sx q[3];
rz(-1.2562498) q[3];
sx q[3];
rz(-3.1414978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67549813) q[2];
sx q[2];
rz(-1.0885295) q[2];
sx q[2];
rz(-1.6509854) q[2];
rz(-2.364482) q[3];
sx q[3];
rz(-1.6831574) q[3];
sx q[3];
rz(1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260088) q[0];
sx q[0];
rz(-1.6827787) q[0];
sx q[0];
rz(-2.7050731) q[0];
rz(0.79477683) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(-2.2025542) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8191231) q[0];
sx q[0];
rz(-1.7717585) q[0];
sx q[0];
rz(-0.51216258) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.375023) q[2];
sx q[2];
rz(-1.1423282) q[2];
sx q[2];
rz(0.049190532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66575501) q[1];
sx q[1];
rz(-1.2547973) q[1];
sx q[1];
rz(2.5412987) q[1];
rz(-pi) q[2];
rz(-2.9211798) q[3];
sx q[3];
rz(-2.1142981) q[3];
sx q[3];
rz(1.8956748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7006435) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(-0.6353333) q[2];
rz(2.3144531) q[3];
sx q[3];
rz(-0.45378903) q[3];
sx q[3];
rz(1.872983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24469911) q[0];
sx q[0];
rz(-3.0785955) q[0];
sx q[0];
rz(2.6196106) q[0];
rz(0.43102795) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(2.4897051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9126221) q[0];
sx q[0];
rz(-2.5755799) q[0];
sx q[0];
rz(0.45155756) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84649936) q[2];
sx q[2];
rz(-1.9545467) q[2];
sx q[2];
rz(2.4946314) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74002162) q[1];
sx q[1];
rz(-1.5844939) q[1];
sx q[1];
rz(-1.1915951) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4383955) q[3];
sx q[3];
rz(-1.0872989) q[3];
sx q[3];
rz(-1.0994224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42644694) q[2];
sx q[2];
rz(-1.2812252) q[2];
sx q[2];
rz(-2.35671) q[2];
rz(2.6134885) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(-1.2877134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89609471) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(-0.78952638) q[0];
rz(0.49742571) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(1.3015889) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808944) q[0];
sx q[0];
rz(-2.5142365) q[0];
sx q[0];
rz(-0.30904667) q[0];
rz(0.33904262) q[2];
sx q[2];
rz(-1.6280988) q[2];
sx q[2];
rz(2.6471241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6005206) q[1];
sx q[1];
rz(-2.3883935) q[1];
sx q[1];
rz(-0.25347565) q[1];
x q[2];
rz(1.6297518) q[3];
sx q[3];
rz(-2.7736933) q[3];
sx q[3];
rz(2.2045362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5975981) q[2];
sx q[2];
rz(-1.602729) q[2];
sx q[2];
rz(2.8653115) q[2];
rz(-2.1918519) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(0.55324078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5211869) q[0];
sx q[0];
rz(-3.1053472) q[0];
sx q[0];
rz(2.1976443) q[0];
rz(-1.2983407) q[1];
sx q[1];
rz(-2.0746168) q[1];
sx q[1];
rz(-2.3640769) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.992618) q[0];
sx q[0];
rz(-1.0273522) q[0];
sx q[0];
rz(1.1538366) q[0];
x q[1];
rz(1.3087511) q[2];
sx q[2];
rz(-2.2876559) q[2];
sx q[2];
rz(-2.1070631) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3537763) q[1];
sx q[1];
rz(-2.750038) q[1];
sx q[1];
rz(-0.33772525) q[1];
x q[2];
rz(1.603807) q[3];
sx q[3];
rz(-1.2668161) q[3];
sx q[3];
rz(-0.98807166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62766084) q[2];
sx q[2];
rz(-1.7646503) q[2];
sx q[2];
rz(2.7505006) q[2];
rz(0.62134653) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(2.3656316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710627) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(-1.712557) q[0];
rz(-2.1757226) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(2.3557854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98322372) q[0];
sx q[0];
rz(-2.2136723) q[0];
sx q[0];
rz(-1.6352562) q[0];
rz(-pi) q[1];
rz(-1.1675535) q[2];
sx q[2];
rz(-0.87783646) q[2];
sx q[2];
rz(-2.5041818) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.475226) q[1];
sx q[1];
rz(-2.6738648) q[1];
sx q[1];
rz(2.9592614) q[1];
rz(-pi) q[2];
rz(-1.2811411) q[3];
sx q[3];
rz(-1.9928586) q[3];
sx q[3];
rz(-2.6817123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69661951) q[2];
sx q[2];
rz(-2.5225621) q[2];
sx q[2];
rz(2.6444198) q[2];
rz(2.2677126) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9484321) q[0];
sx q[0];
rz(-0.30696294) q[0];
sx q[0];
rz(0.69751414) q[0];
rz(-0.3802309) q[1];
sx q[1];
rz(-1.1153328) q[1];
sx q[1];
rz(0.14022216) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33454681) q[0];
sx q[0];
rz(-1.3893034) q[0];
sx q[0];
rz(1.6491778) q[0];
rz(-pi) q[1];
rz(0.44084844) q[2];
sx q[2];
rz(-1.5356283) q[2];
sx q[2];
rz(-0.22415417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6483874) q[1];
sx q[1];
rz(-0.81118656) q[1];
sx q[1];
rz(-0.56583515) q[1];
rz(-1.7626552) q[3];
sx q[3];
rz(-1.1664806) q[3];
sx q[3];
rz(-1.2412925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2908638) q[2];
sx q[2];
rz(-1.9242492) q[2];
sx q[2];
rz(-0.49120894) q[2];
rz(2.9344007) q[3];
sx q[3];
rz(-2.5184293) q[3];
sx q[3];
rz(1.4108968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17289138) q[0];
sx q[0];
rz(-1.4145114) q[0];
sx q[0];
rz(0.15039314) q[0];
rz(-0.70101678) q[1];
sx q[1];
rz(-1.0944518) q[1];
sx q[1];
rz(1.4775803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9532861) q[0];
sx q[0];
rz(-0.64670783) q[0];
sx q[0];
rz(-0.20215277) q[0];
x q[1];
rz(0.16123743) q[2];
sx q[2];
rz(-1.45597) q[2];
sx q[2];
rz(1.6063362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25178441) q[1];
sx q[1];
rz(-1.4191333) q[1];
sx q[1];
rz(-1.4633578) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2712237) q[3];
sx q[3];
rz(-0.75787395) q[3];
sx q[3];
rz(1.2927428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56665862) q[2];
sx q[2];
rz(-2.7329972) q[2];
sx q[2];
rz(1.5236141) q[2];
rz(-0.59897113) q[3];
sx q[3];
rz(-2.6409769) q[3];
sx q[3];
rz(-0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4481675) q[0];
sx q[0];
rz(-2.7126815) q[0];
sx q[0];
rz(-1.6397788) q[0];
rz(0.93694726) q[1];
sx q[1];
rz(-1.0277156) q[1];
sx q[1];
rz(2.1243336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1382683) q[0];
sx q[0];
rz(-1.956245) q[0];
sx q[0];
rz(-1.9786563) q[0];
rz(3.0363778) q[2];
sx q[2];
rz(-1.6561837) q[2];
sx q[2];
rz(-0.24959031) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81472936) q[1];
sx q[1];
rz(-1.6039324) q[1];
sx q[1];
rz(-1.8327946) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3689092) q[3];
sx q[3];
rz(-1.1491778) q[3];
sx q[3];
rz(3.0516171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1169869) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(-1.619722) q[2];
rz(1.4874602) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37987729) q[0];
sx q[0];
rz(-1.3991671) q[0];
sx q[0];
rz(0.73200926) q[0];
rz(-1.4749745) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(2.1521062) q[2];
sx q[2];
rz(-1.9099997) q[2];
sx q[2];
rz(2.9861502) q[2];
rz(0.75763221) q[3];
sx q[3];
rz(-2.7136867) q[3];
sx q[3];
rz(2.4242662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
