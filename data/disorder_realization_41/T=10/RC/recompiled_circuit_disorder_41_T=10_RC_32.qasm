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
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(4.6012576) q[1];
sx q[1];
rz(9.6396946) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.874274) q[0];
sx q[0];
rz(-1.312717) q[0];
sx q[0];
rz(1.2866856) q[0];
rz(-pi) q[1];
rz(-1.4161413) q[2];
sx q[2];
rz(-1.1064331) q[2];
sx q[2];
rz(2.2048339) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51337459) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(2.2378504) q[1];
rz(-pi) q[2];
rz(2.05902) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(-1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(1.365186) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(0.89675084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1363163) q[0];
sx q[0];
rz(-1.7668084) q[0];
sx q[0];
rz(2.2851903) q[0];
rz(1.6438269) q[2];
sx q[2];
rz(-1.4023997) q[2];
sx q[2];
rz(1.3161236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36205081) q[1];
sx q[1];
rz(-2.3565787) q[1];
sx q[1];
rz(-1.7464459) q[1];
rz(-0.81929042) q[3];
sx q[3];
rz(-1.2580401) q[3];
sx q[3];
rz(-2.6889338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(2.7064586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66416392) q[0];
sx q[0];
rz(-1.4633044) q[0];
sx q[0];
rz(-2.7880923) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4235731) q[2];
sx q[2];
rz(-1.8185116) q[2];
sx q[2];
rz(-1.148828) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5907026) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(-2.9571556) q[1];
rz(-0.35878351) q[3];
sx q[3];
rz(-2.3453418) q[3];
sx q[3];
rz(-0.96804726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(0.91745013) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47539513) q[0];
sx q[0];
rz(-2.4817433) q[0];
sx q[0];
rz(-0.049113627) q[0];
x q[1];
rz(-1.2722837) q[2];
sx q[2];
rz(-1.304317) q[2];
sx q[2];
rz(-2.4556015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5323822) q[1];
sx q[1];
rz(-2.2717443) q[1];
sx q[1];
rz(2.1125395) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81569205) q[3];
sx q[3];
rz(-0.94592735) q[3];
sx q[3];
rz(-2.9790785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.778487) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(-1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(2.1889401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1088516) q[0];
sx q[0];
rz(-1.4534338) q[0];
sx q[0];
rz(2.5928241) q[0];
rz(2.7258337) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(-1.8005467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.488727) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(0.25686849) q[1];
x q[2];
rz(-3.0609344) q[3];
sx q[3];
rz(-2.6622052) q[3];
sx q[3];
rz(-2.5731034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7197363) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(-0.83166844) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96777746) q[0];
sx q[0];
rz(-0.41668188) q[0];
sx q[0];
rz(-1.5370876) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8820011) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(1.358658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3249358) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(1.9865958) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5112126) q[3];
sx q[3];
rz(-1.3833589) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(2.0149569) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28850266) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(0.75659928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3698759) q[0];
sx q[0];
rz(-1.6635832) q[0];
sx q[0];
rz(1.4183527) q[0];
rz(-2.7034764) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(2.9012836) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25890299) q[1];
sx q[1];
rz(-1.8001302) q[1];
sx q[1];
rz(2.980466) q[1];
x q[2];
rz(1.5114622) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(-2.4997366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-2.9471617) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(3.074926) q[0];
rz(0.32456675) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522404) q[0];
sx q[0];
rz(-1.1712495) q[0];
sx q[0];
rz(0.35895343) q[0];
rz(-pi) q[1];
rz(-0.093401508) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(-1.8708558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5576396) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(3.0686892) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.057468) q[3];
sx q[3];
rz(-1.1063965) q[3];
sx q[3];
rz(1.407479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-0.40763339) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-2.5323903) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0489037) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(-1.2094686) q[0];
rz(-pi) q[1];
rz(0.74553211) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(-1.0996639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59233353) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(-1.1997644) q[1];
x q[2];
rz(-0.94536762) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(-0.35811801) q[3];
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
rz(0.68230391) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(-0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.3964765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8499334) q[0];
sx q[0];
rz(-2.8072661) q[0];
sx q[0];
rz(-2.8773984) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2730607) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(-2.0993078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5212621) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(-1.2465338) q[1];
x q[2];
rz(-1.6276007) q[3];
sx q[3];
rz(-1.8326522) q[3];
sx q[3];
rz(1.1893623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0970584) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(-1.8272022) q[2];
sx q[2];
rz(-0.64654965) q[2];
sx q[2];
rz(1.8241573) q[2];
rz(2.6079569) q[3];
sx q[3];
rz(-1.5279557) q[3];
sx q[3];
rz(0.91844311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
