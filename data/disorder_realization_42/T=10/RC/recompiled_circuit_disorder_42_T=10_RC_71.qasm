OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6289829) q[0];
sx q[0];
rz(4.1689685) q[0];
sx q[0];
rz(12.203759) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5255614) q[0];
sx q[0];
rz(-1.4481059) q[0];
sx q[0];
rz(-2.7659155) q[0];
rz(-pi) q[1];
rz(-2.5506637) q[2];
sx q[2];
rz(-1.8701631) q[2];
sx q[2];
rz(-2.9654944) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8475854) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(-1.0311544) q[1];
rz(-pi) q[2];
rz(-2.9176641) q[3];
sx q[3];
rz(-1.9609946) q[3];
sx q[3];
rz(-1.8412561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(-2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0618806) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.99615) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(-0.72584814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84859914) q[0];
sx q[0];
rz(-1.2147875) q[0];
sx q[0];
rz(2.5550935) q[0];
rz(-2.5088596) q[2];
sx q[2];
rz(-0.69721141) q[2];
sx q[2];
rz(-0.34174191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1314288) q[1];
sx q[1];
rz(-1.5457075) q[1];
sx q[1];
rz(-1.2618834) q[1];
rz(-2.8511091) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(-3.052352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-0.99728161) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230351) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(-2.0130656) q[0];
rz(1.3035125) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(0.9054786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826686) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(2.8959136) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.682914) q[2];
sx q[2];
rz(-2.3346666) q[2];
sx q[2];
rz(-2.5708831) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9943774) q[1];
sx q[1];
rz(-0.71502393) q[1];
sx q[1];
rz(-1.4579525) q[1];
x q[2];
rz(2.2842992) q[3];
sx q[3];
rz(-1.6033353) q[3];
sx q[3];
rz(-1.1710222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(3.0958214) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0320597) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(0.12606829) q[0];
rz(0.18445045) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(-1.9741845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8465189) q[0];
sx q[0];
rz(-2.3210038) q[0];
sx q[0];
rz(1.8473162) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8235769) q[2];
sx q[2];
rz(-0.52581767) q[2];
sx q[2];
rz(-1.9981245) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0427525) q[1];
sx q[1];
rz(-2.4272857) q[1];
sx q[1];
rz(0.40978281) q[1];
rz(-pi) q[2];
rz(-0.91058369) q[3];
sx q[3];
rz(-1.0331717) q[3];
sx q[3];
rz(2.0349353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2771153) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(-2.946092) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1642078) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-2.6745093) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5669117) q[0];
sx q[0];
rz(-1.3710638) q[0];
sx q[0];
rz(1.1916222) q[0];
x q[1];
rz(-0.37093016) q[2];
sx q[2];
rz(-2.5226454) q[2];
sx q[2];
rz(1.3528454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0568697) q[1];
sx q[1];
rz(-1.1493756) q[1];
sx q[1];
rz(-2.9028068) q[1];
x q[2];
rz(2.2561982) q[3];
sx q[3];
rz(-1.6939031) q[3];
sx q[3];
rz(2.3791594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-0.49218991) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(1.6451947) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841227) q[0];
sx q[0];
rz(-0.5623445) q[0];
sx q[0];
rz(-1.9447295) q[0];
x q[1];
rz(2.1093844) q[2];
sx q[2];
rz(-1.7787873) q[2];
sx q[2];
rz(1.6230735) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3277668) q[1];
sx q[1];
rz(-1.5470978) q[1];
sx q[1];
rz(-1.2162672) q[1];
rz(-pi) q[2];
rz(2.4286527) q[3];
sx q[3];
rz(-1.5804277) q[3];
sx q[3];
rz(-0.18129098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4795586) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(3.1299652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893338) q[0];
sx q[0];
rz(-2.6703435) q[0];
sx q[0];
rz(-0.61347368) q[0];
rz(-pi) q[1];
rz(-2.3064012) q[2];
sx q[2];
rz(-0.037926849) q[2];
sx q[2];
rz(1.4284301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7736518) q[1];
sx q[1];
rz(-2.5950187) q[1];
sx q[1];
rz(1.6163338) q[1];
x q[2];
rz(-0.67354789) q[3];
sx q[3];
rz(-1.827497) q[3];
sx q[3];
rz(1.7208769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(-1.9912432) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(-1.2121375) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(1.4935965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9858915) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(2.5168688) q[0];
rz(-1.1475032) q[2];
sx q[2];
rz(-2.9967542) q[2];
sx q[2];
rz(-0.3651948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32840604) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(2.3996668) q[1];
x q[2];
rz(0.11843189) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-1.1008319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(-1.2517694) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-0.1782724) q[0];
rz(0.072862236) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7992135) q[0];
sx q[0];
rz(-1.2782492) q[0];
sx q[0];
rz(0.53036687) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91782848) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(-0.81673056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0852016) q[1];
sx q[1];
rz(-1.8279652) q[1];
sx q[1];
rz(0.53253865) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0719028) q[3];
sx q[3];
rz(-1.2341098) q[3];
sx q[3];
rz(-1.6285553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(2.1098095) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(-2.0526989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63554791) q[0];
sx q[0];
rz(-1.5061) q[0];
sx q[0];
rz(-0.062775469) q[0];
rz(-pi) q[1];
rz(-0.18386545) q[2];
sx q[2];
rz(-0.92391787) q[2];
sx q[2];
rz(-1.136214) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1354462) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(1.0165434) q[1];
rz(1.9128996) q[3];
sx q[3];
rz(-1.7497239) q[3];
sx q[3];
rz(0.10781328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(0.87674117) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4927647) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(1.6246673) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(2.7355315) q[2];
sx q[2];
rz(-0.95352298) q[2];
sx q[2];
rz(2.7609115) q[2];
rz(1.6886961) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
