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
rz(-0.63646746) q[0];
sx q[0];
rz(-0.31390733) q[0];
sx q[0];
rz(-1.7933581) q[0];
rz(2.1836166) q[1];
sx q[1];
rz(-1.7254683) q[1];
sx q[1];
rz(0.15813601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30017553) q[0];
sx q[0];
rz(-2.149507) q[0];
sx q[0];
rz(-2.7912223) q[0];
rz(-2.4329107) q[2];
sx q[2];
rz(-2.6593609) q[2];
sx q[2];
rz(0.9672375) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6414538) q[1];
sx q[1];
rz(-1.5654501) q[1];
sx q[1];
rz(3.1405906) q[1];
x q[2];
rz(-1.1920378) q[3];
sx q[3];
rz(-2.023732) q[3];
sx q[3];
rz(-0.68309802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30369514) q[2];
sx q[2];
rz(-0.94472307) q[2];
sx q[2];
rz(-2.5486805) q[2];
rz(0.29218519) q[3];
sx q[3];
rz(-0.018915011) q[3];
sx q[3];
rz(-1.7570447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050046571) q[0];
sx q[0];
rz(-0.36588359) q[0];
sx q[0];
rz(-3.0475317) q[0];
rz(1.7496109) q[1];
sx q[1];
rz(-1.6110907) q[1];
sx q[1];
rz(-1.403341) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7624403) q[0];
sx q[0];
rz(-1.8383674) q[0];
sx q[0];
rz(-1.3570157) q[0];
rz(-pi) q[1];
rz(-1.044079) q[2];
sx q[2];
rz(-1.7657914) q[2];
sx q[2];
rz(-3.1205683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5793094) q[1];
sx q[1];
rz(-2.1726949) q[1];
sx q[1];
rz(-3.137869) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7330328) q[3];
sx q[3];
rz(-1.4848391) q[3];
sx q[3];
rz(-2.3547821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89468181) q[2];
sx q[2];
rz(-0.44591388) q[2];
sx q[2];
rz(1.1723588) q[2];
rz(0.42919484) q[3];
sx q[3];
rz(-2.6515638) q[3];
sx q[3];
rz(1.711285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9046852) q[0];
sx q[0];
rz(-2.560736) q[0];
sx q[0];
rz(2.7823606) q[0];
rz(-1.5776186) q[1];
sx q[1];
rz(-2.3465395) q[1];
sx q[1];
rz(-0.94594812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1079179) q[0];
sx q[0];
rz(-1.6527411) q[0];
sx q[0];
rz(1.7235403) q[0];
rz(1.6851186) q[2];
sx q[2];
rz(-1.6094484) q[2];
sx q[2];
rz(-1.1383411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11752973) q[1];
sx q[1];
rz(-0.27654031) q[1];
sx q[1];
rz(-1.0599972) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69880693) q[3];
sx q[3];
rz(-0.10891373) q[3];
sx q[3];
rz(2.1366936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.642091) q[2];
sx q[2];
rz(-1.6104128) q[2];
sx q[2];
rz(1.1031319) q[2];
rz(-1.057386) q[3];
sx q[3];
rz(-1.5921581) q[3];
sx q[3];
rz(-2.9252083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0353521) q[0];
sx q[0];
rz(-0.092656605) q[0];
sx q[0];
rz(2.4308391) q[0];
rz(-2.4757929) q[1];
sx q[1];
rz(-0.0087105287) q[1];
sx q[1];
rz(-0.3054558) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75765002) q[0];
sx q[0];
rz(-0.096127495) q[0];
sx q[0];
rz(2.686475) q[0];
x q[1];
rz(-2.6389559) q[2];
sx q[2];
rz(-1.8515311) q[2];
sx q[2];
rz(1.0681149) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98887029) q[1];
sx q[1];
rz(-0.45256361) q[1];
sx q[1];
rz(1.7532224) q[1];
rz(-pi) q[2];
rz(1.0287343) q[3];
sx q[3];
rz(-2.4031478) q[3];
sx q[3];
rz(-1.5800103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1315883) q[2];
sx q[2];
rz(-1.5350716) q[2];
sx q[2];
rz(-0.32001495) q[2];
rz(-0.53712505) q[3];
sx q[3];
rz(-2.7607626) q[3];
sx q[3];
rz(0.91184688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.491275) q[0];
sx q[0];
rz(-0.22458751) q[0];
sx q[0];
rz(-3.0564296) q[0];
rz(2.5501309) q[1];
sx q[1];
rz(-3.1384835) q[1];
sx q[1];
rz(1.1726146) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64567287) q[0];
sx q[0];
rz(-1.572019) q[0];
sx q[0];
rz(-1.5652577) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3306562) q[2];
sx q[2];
rz(-1.2954172) q[2];
sx q[2];
rz(-0.20698838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53211305) q[1];
sx q[1];
rz(-2.0572741) q[1];
sx q[1];
rz(2.3064613) q[1];
rz(-pi) q[2];
rz(-2.8211604) q[3];
sx q[3];
rz(-2.6310807) q[3];
sx q[3];
rz(1.3714457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0255967) q[2];
sx q[2];
rz(-1.3469348) q[2];
sx q[2];
rz(-1.4541413) q[2];
rz(-1.2429552) q[3];
sx q[3];
rz(-0.60296139) q[3];
sx q[3];
rz(-2.3292144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1989307) q[0];
sx q[0];
rz(-2.9781065) q[0];
sx q[0];
rz(-1.7401975) q[0];
rz(-0.61680782) q[1];
sx q[1];
rz(-0.016409358) q[1];
sx q[1];
rz(1.1161463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45081954) q[0];
sx q[0];
rz(-1.6661394) q[0];
sx q[0];
rz(1.3006163) q[0];
x q[1];
rz(2.7737439) q[2];
sx q[2];
rz(-0.90613885) q[2];
sx q[2];
rz(-0.14205113) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20780288) q[1];
sx q[1];
rz(-2.3456165) q[1];
sx q[1];
rz(-0.7178623) q[1];
rz(-pi) q[2];
rz(2.6924388) q[3];
sx q[3];
rz(-1.7633668) q[3];
sx q[3];
rz(1.4447007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8650032) q[2];
sx q[2];
rz(-1.208655) q[2];
sx q[2];
rz(1.8763982) q[2];
rz(-2.4943634) q[3];
sx q[3];
rz(-0.45014683) q[3];
sx q[3];
rz(1.885672) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7728421) q[0];
sx q[0];
rz(-2.5970646) q[0];
sx q[0];
rz(-0.37753373) q[0];
rz(-0.15097161) q[1];
sx q[1];
rz(-0.010846373) q[1];
sx q[1];
rz(-0.46802256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6418544) q[0];
sx q[0];
rz(-1.5473135) q[0];
sx q[0];
rz(-2.6905943) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0156844) q[2];
sx q[2];
rz(-0.85239886) q[2];
sx q[2];
rz(0.56294051) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49761729) q[1];
sx q[1];
rz(-2.4598498) q[1];
sx q[1];
rz(1.6830744) q[1];
x q[2];
rz(-0.19429732) q[3];
sx q[3];
rz(-1.0602942) q[3];
sx q[3];
rz(2.0795151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2630792) q[2];
sx q[2];
rz(-0.24363467) q[2];
sx q[2];
rz(1.1070975) q[2];
rz(2.257972) q[3];
sx q[3];
rz(-2.3487909) q[3];
sx q[3];
rz(2.4293374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368619) q[0];
sx q[0];
rz(-0.87918133) q[0];
sx q[0];
rz(-2.9469446) q[0];
rz(1.4600935) q[1];
sx q[1];
rz(-3.1250521) q[1];
sx q[1];
rz(0.3009235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096357249) q[0];
sx q[0];
rz(-2.3060252) q[0];
sx q[0];
rz(0.51765963) q[0];
x q[1];
rz(-0.76809068) q[2];
sx q[2];
rz(-0.67736191) q[2];
sx q[2];
rz(-1.0528262) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4408988) q[1];
sx q[1];
rz(-1.3446302) q[1];
sx q[1];
rz(-2.6146099) q[1];
x q[2];
rz(-1.5125723) q[3];
sx q[3];
rz(-1.749862) q[3];
sx q[3];
rz(1.5032081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.459317) q[2];
sx q[2];
rz(-1.6271017) q[2];
sx q[2];
rz(0.21919361) q[2];
rz(1.7949665) q[3];
sx q[3];
rz(-3.0248088) q[3];
sx q[3];
rz(2.35675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.858736) q[0];
sx q[0];
rz(-2.5114926) q[0];
sx q[0];
rz(0.0081188763) q[0];
rz(-0.66932622) q[1];
sx q[1];
rz(-0.012834276) q[1];
sx q[1];
rz(-2.377811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8421461) q[0];
sx q[0];
rz(-1.1538635) q[0];
sx q[0];
rz(-1.3119158) q[0];
rz(3.0349651) q[2];
sx q[2];
rz(-2.5322057) q[2];
sx q[2];
rz(-2.0963017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8593401) q[1];
sx q[1];
rz(-2.1164923) q[1];
sx q[1];
rz(0.061435862) q[1];
rz(-2.2586601) q[3];
sx q[3];
rz(-2.2730423) q[3];
sx q[3];
rz(-1.8851999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1179489) q[2];
sx q[2];
rz(-0.014160784) q[2];
sx q[2];
rz(0.99948779) q[2];
rz(0.1855447) q[3];
sx q[3];
rz(-1.0368291) q[3];
sx q[3];
rz(3.1256092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213181) q[0];
sx q[0];
rz(-0.11225926) q[0];
sx q[0];
rz(2.4768594) q[0];
rz(0.96066535) q[1];
sx q[1];
rz(-3.101109) q[1];
sx q[1];
rz(-1.6368846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5793631) q[0];
sx q[0];
rz(-2.7873171) q[0];
sx q[0];
rz(2.1112403) q[0];
rz(-pi) q[1];
rz(-3.0111752) q[2];
sx q[2];
rz(-2.7777024) q[2];
sx q[2];
rz(-0.39375776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0814652) q[1];
sx q[1];
rz(-1.0614396) q[1];
sx q[1];
rz(2.1101923) q[1];
x q[2];
rz(0.89791132) q[3];
sx q[3];
rz(-2.5012883) q[3];
sx q[3];
rz(1.4434356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30161834) q[2];
sx q[2];
rz(-3.1349389) q[2];
sx q[2];
rz(1.9807695) q[2];
rz(-0.079856722) q[3];
sx q[3];
rz(-3.1294332) q[3];
sx q[3];
rz(-1.9399835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53497159) q[0];
sx q[0];
rz(-1.4749682) q[0];
sx q[0];
rz(-1.5780021) q[0];
rz(3.1040991) q[1];
sx q[1];
rz(-2.9480724) q[1];
sx q[1];
rz(-2.9185157) q[1];
rz(-1.1983331) q[2];
sx q[2];
rz(-1.6732775) q[2];
sx q[2];
rz(0.35107935) q[2];
rz(-0.014122176) q[3];
sx q[3];
rz(-0.89992186) q[3];
sx q[3];
rz(1.810771) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
