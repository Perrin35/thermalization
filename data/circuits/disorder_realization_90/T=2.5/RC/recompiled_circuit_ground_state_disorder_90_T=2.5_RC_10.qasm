OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6381792) q[0];
sx q[0];
rz(-1.8621651) q[0];
sx q[0];
rz(2.3655565) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(0.59901839) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020928362) q[0];
sx q[0];
rz(-2.0428223) q[0];
sx q[0];
rz(-0.13735227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2602735) q[2];
sx q[2];
rz(-1.3506839) q[2];
sx q[2];
rz(-0.51267363) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4091275) q[1];
sx q[1];
rz(-0.62668398) q[1];
sx q[1];
rz(-0.62925689) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1407858) q[3];
sx q[3];
rz(-0.83807349) q[3];
sx q[3];
rz(-3.0504459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1385931) q[2];
sx q[2];
rz(-1.1370167) q[2];
sx q[2];
rz(-0.19621672) q[2];
rz(-1.044322) q[3];
sx q[3];
rz(-0.82572562) q[3];
sx q[3];
rz(-2.830982) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0394548) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(-0.85773221) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(-2.3589755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1038541) q[0];
sx q[0];
rz(-1.4840992) q[0];
sx q[0];
rz(2.187665) q[0];
rz(1.247632) q[2];
sx q[2];
rz(-2.4944381) q[2];
sx q[2];
rz(-2.0314856) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23797135) q[1];
sx q[1];
rz(-0.97717932) q[1];
sx q[1];
rz(1.3376544) q[1];
x q[2];
rz(2.2563062) q[3];
sx q[3];
rz(-1.1210223) q[3];
sx q[3];
rz(1.3441844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9767849) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(-2.5505193) q[2];
rz(-1.0278541) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(-0.13636057) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117821) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(-2.0853364) q[0];
rz(-0.99110574) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(-2.3371005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5493889) q[0];
sx q[0];
rz(-2.443586) q[0];
sx q[0];
rz(2.7207359) q[0];
rz(0.94203888) q[2];
sx q[2];
rz(-1.647718) q[2];
sx q[2];
rz(2.6475737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.04136297) q[1];
sx q[1];
rz(-1.3303262) q[1];
sx q[1];
rz(0.23794707) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8896585) q[3];
sx q[3];
rz(-2.478699) q[3];
sx q[3];
rz(-0.54352647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6911917) q[2];
sx q[2];
rz(-2.2914026) q[2];
sx q[2];
rz(-2.5737393) q[2];
rz(1.9495226) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49517) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(1.1871185) q[0];
rz(0.25539708) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(0.38937169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0856649) q[0];
sx q[0];
rz(-1.4367668) q[0];
sx q[0];
rz(-2.9197951) q[0];
x q[1];
rz(0.03002982) q[2];
sx q[2];
rz(-2.2976934) q[2];
sx q[2];
rz(-0.50315693) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43453211) q[1];
sx q[1];
rz(-1.7599306) q[1];
sx q[1];
rz(2.2347336) q[1];
rz(-3.0087972) q[3];
sx q[3];
rz(-1.7657874) q[3];
sx q[3];
rz(-2.8094069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(1.7741989) q[2];
rz(0.5395475) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(-1.8168943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668199) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(-2.5812126) q[0];
rz(2.9226774) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(0.17094368) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4881011) q[0];
sx q[0];
rz(-2.2224036) q[0];
sx q[0];
rz(2.1733858) q[0];
rz(-pi) q[1];
rz(2.937491) q[2];
sx q[2];
rz(-1.4320254) q[2];
sx q[2];
rz(-0.52817527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39300181) q[1];
sx q[1];
rz(-0.95470631) q[1];
sx q[1];
rz(-2.2115117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1910652) q[3];
sx q[3];
rz(-1.5135362) q[3];
sx q[3];
rz(-1.7050153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35974744) q[2];
sx q[2];
rz(-2.013194) q[2];
sx q[2];
rz(0.95853364) q[2];
rz(-0.16252276) q[3];
sx q[3];
rz(-2.2992117) q[3];
sx q[3];
rz(1.5925647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-3.077226) q[0];
sx q[0];
rz(2.3934613) q[0];
rz(-2.023078) q[1];
sx q[1];
rz(-0.41901127) q[1];
sx q[1];
rz(-2.8245139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1823163) q[0];
sx q[0];
rz(-1.3249363) q[0];
sx q[0];
rz(1.5005174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9002764) q[2];
sx q[2];
rz(-0.20212999) q[2];
sx q[2];
rz(-1.2723107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3534155) q[1];
sx q[1];
rz(-1.6420351) q[1];
sx q[1];
rz(0.80435462) q[1];
rz(-1.6104782) q[3];
sx q[3];
rz(-2.3473266) q[3];
sx q[3];
rz(-0.021707699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76724425) q[2];
sx q[2];
rz(-2.2129462) q[2];
sx q[2];
rz(1.7283776) q[2];
rz(3.0610541) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(-2.9197689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11248511) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(1.8674194) q[0];
rz(1.3451276) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(1.8003731) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2224436) q[0];
sx q[0];
rz(-2.0608927) q[0];
sx q[0];
rz(0.63700139) q[0];
rz(-1.28181) q[2];
sx q[2];
rz(-2.0353051) q[2];
sx q[2];
rz(-1.7211085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6051424) q[1];
sx q[1];
rz(-1.4172232) q[1];
sx q[1];
rz(1.3773514) q[1];
rz(0.11734924) q[3];
sx q[3];
rz(-2.2175991) q[3];
sx q[3];
rz(-2.5600159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93366569) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(-2.6742317) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(-1.8925586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6811328) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(-1.6946174) q[0];
rz(-2.815822) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(1.5209341) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50624079) q[0];
sx q[0];
rz(-1.304848) q[0];
sx q[0];
rz(2.6145934) q[0];
x q[1];
rz(1.5711741) q[2];
sx q[2];
rz(-1.4791569) q[2];
sx q[2];
rz(-2.3204892) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9190678) q[1];
sx q[1];
rz(-1.7456858) q[1];
sx q[1];
rz(-1.9293849) q[1];
rz(-2.2748442) q[3];
sx q[3];
rz(-1.4102077) q[3];
sx q[3];
rz(-0.15204568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96523607) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(-2.526324) q[2];
rz(1.2728914) q[3];
sx q[3];
rz(-0.26510173) q[3];
sx q[3];
rz(-0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1212921) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(-0.85025382) q[0];
rz(-2.0203159) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(-0.77176315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61311713) q[0];
sx q[0];
rz(-2.4474047) q[0];
sx q[0];
rz(1.608169) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1807801) q[2];
sx q[2];
rz(-1.5018936) q[2];
sx q[2];
rz(-3.1071752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1771639) q[1];
sx q[1];
rz(-1.5156974) q[1];
sx q[1];
rz(1.3269618) q[1];
x q[2];
rz(1.8949299) q[3];
sx q[3];
rz(-1.8283024) q[3];
sx q[3];
rz(2.3070564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11543342) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(0.8405295) q[2];
rz(0.080951512) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7131272) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(1.9848829) q[0];
rz(1.0276065) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(-2.3311232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9543957) q[0];
sx q[0];
rz(-2.3141197) q[0];
sx q[0];
rz(-1.8660956) q[0];
rz(-pi) q[1];
rz(-0.54525156) q[2];
sx q[2];
rz(-1.8601189) q[2];
sx q[2];
rz(-0.74301737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15511423) q[1];
sx q[1];
rz(-0.41448516) q[1];
sx q[1];
rz(-2.3398967) q[1];
rz(1.719172) q[3];
sx q[3];
rz(-2.1273489) q[3];
sx q[3];
rz(2.839193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79601866) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(1.6176809) q[2];
rz(-1.1003305) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(-2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0236459) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(-3.0991411) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(1.7003851) q[2];
sx q[2];
rz(-2.7322506) q[2];
sx q[2];
rz(2.6553497) q[2];
rz(1.6067947) q[3];
sx q[3];
rz(-1.1136342) q[3];
sx q[3];
rz(-1.3200091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
