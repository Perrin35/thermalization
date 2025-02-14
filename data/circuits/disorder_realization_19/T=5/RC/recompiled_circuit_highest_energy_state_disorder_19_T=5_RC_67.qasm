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
rz(0.49993604) q[0];
sx q[0];
rz(-1.191782) q[0];
sx q[0];
rz(1.7697822) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(-0.86725441) q[1];
sx q[1];
rz(0.39630085) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4258408) q[0];
sx q[0];
rz(-1.4914728) q[0];
sx q[0];
rz(-0.21004814) q[0];
rz(-pi) q[1];
rz(0.17163817) q[2];
sx q[2];
rz(-1.2002103) q[2];
sx q[2];
rz(2.5570392) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5057136) q[1];
sx q[1];
rz(-1.4766271) q[1];
sx q[1];
rz(2.0999184) q[1];
rz(-pi) q[2];
rz(-0.62403535) q[3];
sx q[3];
rz(-2.6665832) q[3];
sx q[3];
rz(1.6866682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1402011) q[2];
sx q[2];
rz(-2.3346257) q[2];
sx q[2];
rz(-1.3414471) q[2];
rz(2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(-1.3816381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694201) q[0];
sx q[0];
rz(-0.91745806) q[0];
sx q[0];
rz(-2.4826352) q[0];
rz(3.068889) q[1];
sx q[1];
rz(-1.0560938) q[1];
sx q[1];
rz(1.2603849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1457659) q[0];
sx q[0];
rz(-2.6814309) q[0];
sx q[0];
rz(-1.7547248) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2748806) q[2];
sx q[2];
rz(-2.5051691) q[2];
sx q[2];
rz(-1.1918024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5301828) q[1];
sx q[1];
rz(-1.6356704) q[1];
sx q[1];
rz(-0.26574175) q[1];
x q[2];
rz(1.4006413) q[3];
sx q[3];
rz(-2.5239787) q[3];
sx q[3];
rz(-0.028384203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2468804) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(-1.3512208) q[2];
rz(0.61659914) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63224822) q[0];
sx q[0];
rz(-2.6743439) q[0];
sx q[0];
rz(-0.21009357) q[0];
rz(-0.083077438) q[1];
sx q[1];
rz(-1.2385085) q[1];
sx q[1];
rz(0.067213623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7930896) q[0];
sx q[0];
rz(-1.5843035) q[0];
sx q[0];
rz(1.5828787) q[0];
x q[1];
rz(-1.018173) q[2];
sx q[2];
rz(-2.6913096) q[2];
sx q[2];
rz(2.1076815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.018118) q[1];
sx q[1];
rz(-2.2384081) q[1];
sx q[1];
rz(2.7099744) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6609069) q[3];
sx q[3];
rz(-2.2743974) q[3];
sx q[3];
rz(-1.4058324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82789603) q[2];
sx q[2];
rz(-1.4915497) q[2];
sx q[2];
rz(1.776604) q[2];
rz(-0.40417534) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(1.9425758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2760524) q[0];
sx q[0];
rz(-1.7409538) q[0];
sx q[0];
rz(2.6633967) q[0];
rz(-2.1887691) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(2.731954) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.206868) q[0];
sx q[0];
rz(-1.544401) q[0];
sx q[0];
rz(2.4858263) q[0];
rz(-pi) q[1];
rz(-0.19721376) q[2];
sx q[2];
rz(-0.46425113) q[2];
sx q[2];
rz(1.0997538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.86120104) q[1];
sx q[1];
rz(-0.93547677) q[1];
sx q[1];
rz(2.3767571) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8919935) q[3];
sx q[3];
rz(-1.4347335) q[3];
sx q[3];
rz(0.67067671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59565583) q[2];
sx q[2];
rz(-1.4384392) q[2];
sx q[2];
rz(-2.0036009) q[2];
rz(-2.1465837) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.397641) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(-2.8635136) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.2828628) q[1];
sx q[1];
rz(1.5778731) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5681158) q[0];
sx q[0];
rz(-1.2088803) q[0];
sx q[0];
rz(3.0058577) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7817338) q[2];
sx q[2];
rz(-2.3198876) q[2];
sx q[2];
rz(-2.6287959) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1335781) q[1];
sx q[1];
rz(-2.1431794) q[1];
sx q[1];
rz(-1.7126669) q[1];
rz(-1.4117354) q[3];
sx q[3];
rz(-1.4791227) q[3];
sx q[3];
rz(-1.5810458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9331253) q[2];
sx q[2];
rz(-1.6105904) q[2];
sx q[2];
rz(2.7739286) q[2];
rz(1.095088) q[3];
sx q[3];
rz(-2.118066) q[3];
sx q[3];
rz(-2.9663405) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(1.3391986) q[0];
rz(3.0195492) q[1];
sx q[1];
rz(-1.782878) q[1];
sx q[1];
rz(-0.36062127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80194762) q[0];
sx q[0];
rz(-2.3390798) q[0];
sx q[0];
rz(2.5967477) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4819988) q[2];
sx q[2];
rz(-2.7200903) q[2];
sx q[2];
rz(2.5777566) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5721467) q[1];
sx q[1];
rz(-1.1317557) q[1];
sx q[1];
rz(-0.31198685) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.020645647) q[3];
sx q[3];
rz(-0.25293487) q[3];
sx q[3];
rz(-1.3731352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7755166) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(-1.1386846) q[2];
rz(-0.25389296) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(2.3937288) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89895407) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(-0.97578543) q[0];
rz(1.1294533) q[1];
sx q[1];
rz(-2.6934846) q[1];
sx q[1];
rz(-1.7879558) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0669374) q[0];
sx q[0];
rz(-2.2569509) q[0];
sx q[0];
rz(1.9524491) q[0];
rz(-pi) q[1];
rz(1.9949525) q[2];
sx q[2];
rz(-1.4978906) q[2];
sx q[2];
rz(-0.50114252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8541214) q[1];
sx q[1];
rz(-1.439289) q[1];
sx q[1];
rz(2.8033957) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72809345) q[3];
sx q[3];
rz(-0.59641664) q[3];
sx q[3];
rz(2.8567258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5740616) q[2];
sx q[2];
rz(-2.3109544) q[2];
sx q[2];
rz(-2.2243824) q[2];
rz(-1.5489102) q[3];
sx q[3];
rz(-0.33532381) q[3];
sx q[3];
rz(-2.3622021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2061283) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(0.024854831) q[0];
rz(-1.286233) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(1.53055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12537374) q[0];
sx q[0];
rz(-1.9589931) q[0];
sx q[0];
rz(-2.6632705) q[0];
x q[1];
rz(1.3127682) q[2];
sx q[2];
rz(-0.90094968) q[2];
sx q[2];
rz(-2.9473398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81653412) q[1];
sx q[1];
rz(-1.6055425) q[1];
sx q[1];
rz(0.28157061) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3601264) q[3];
sx q[3];
rz(-1.7066741) q[3];
sx q[3];
rz(1.6257515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1207235) q[2];
sx q[2];
rz(-1.0039696) q[2];
sx q[2];
rz(1.8570159) q[2];
rz(-0.32127109) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(-2.9200714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1247509) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(-2.2220213) q[0];
rz(-0.59182566) q[1];
sx q[1];
rz(-2.070919) q[1];
sx q[1];
rz(2.8048973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3172835) q[0];
sx q[0];
rz(-1.4364087) q[0];
sx q[0];
rz(1.3008253) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7783103) q[2];
sx q[2];
rz(-1.4575301) q[2];
sx q[2];
rz(2.670778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7408161) q[1];
sx q[1];
rz(-1.0662406) q[1];
sx q[1];
rz(-1.4228805) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1235663) q[3];
sx q[3];
rz(-2.8447731) q[3];
sx q[3];
rz(0.33949131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1153974) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(-0.78286147) q[2];
rz(0.10910263) q[3];
sx q[3];
rz(-2.1249873) q[3];
sx q[3];
rz(1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03610177) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(-0.92392695) q[0];
rz(0.52664122) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(0.99871666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17820464) q[0];
sx q[0];
rz(-2.9792157) q[0];
sx q[0];
rz(0.33372648) q[0];
rz(-1.8615103) q[2];
sx q[2];
rz(-1.4913342) q[2];
sx q[2];
rz(-3.0018011) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0044788) q[1];
sx q[1];
rz(-0.58885899) q[1];
sx q[1];
rz(-1.5084672) q[1];
rz(-1.2015593) q[3];
sx q[3];
rz(-2.2155361) q[3];
sx q[3];
rz(-0.45756868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0052789) q[2];
sx q[2];
rz(-2.2515991) q[2];
sx q[2];
rz(-0.33373731) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.4849097) q[3];
sx q[3];
rz(-0.0066283289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6776047) q[0];
sx q[0];
rz(-1.9515568) q[0];
sx q[0];
rz(-2.4757181) q[0];
rz(-0.55466501) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(1.8205504) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(-2.4246115) q[3];
sx q[3];
rz(-2.3752799) q[3];
sx q[3];
rz(0.3175288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
