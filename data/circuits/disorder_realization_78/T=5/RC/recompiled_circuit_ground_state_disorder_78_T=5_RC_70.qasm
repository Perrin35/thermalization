OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(-2.3999441) q[1];
sx q[1];
rz(2.9902966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97340527) q[0];
sx q[0];
rz(-0.12790132) q[0];
sx q[0];
rz(1.3692877) q[0];
rz(-3.0252671) q[2];
sx q[2];
rz(-2.1172197) q[2];
sx q[2];
rz(-0.31508581) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1946434) q[1];
sx q[1];
rz(-0.8425172) q[1];
sx q[1];
rz(-0.76232736) q[1];
rz(-2.862418) q[3];
sx q[3];
rz(-0.82133365) q[3];
sx q[3];
rz(2.4031635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1260881) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(-0.31910953) q[2];
rz(-2.0305521) q[3];
sx q[3];
rz(-0.55912656) q[3];
sx q[3];
rz(1.8266953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023271712) q[0];
sx q[0];
rz(-0.51902223) q[0];
sx q[0];
rz(2.5530489) q[0];
rz(2.5449246) q[1];
sx q[1];
rz(-1.3304973) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7241192) q[0];
sx q[0];
rz(-2.3406174) q[0];
sx q[0];
rz(0.95306122) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.071429169) q[2];
sx q[2];
rz(-1.0712581) q[2];
sx q[2];
rz(0.75354924) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6728404) q[1];
sx q[1];
rz(-0.44768718) q[1];
sx q[1];
rz(-1.3393558) q[1];
x q[2];
rz(-0.94542687) q[3];
sx q[3];
rz(-1.4875879) q[3];
sx q[3];
rz(-2.0997467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0236464) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(3.0493951) q[2];
rz(2.2495031) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(1.0649072) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104765) q[0];
sx q[0];
rz(-1.6350063) q[0];
sx q[0];
rz(-0.098467501) q[0];
rz(2.5473728) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(2.1509511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38483175) q[0];
sx q[0];
rz(-0.38479003) q[0];
sx q[0];
rz(-1.8033474) q[0];
x q[1];
rz(0.76108257) q[2];
sx q[2];
rz(-1.154261) q[2];
sx q[2];
rz(1.0786982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4842867) q[1];
sx q[1];
rz(-2.4542744) q[1];
sx q[1];
rz(1.1794075) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1029937) q[3];
sx q[3];
rz(-1.3672223) q[3];
sx q[3];
rz(2.9969203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46382612) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(-0.80219913) q[2];
rz(1.5944611) q[3];
sx q[3];
rz(-1.0651257) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39740729) q[0];
sx q[0];
rz(-2.7825401) q[0];
sx q[0];
rz(1.3209976) q[0];
rz(-1.5253223) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(-0.1604518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9994151) q[0];
sx q[0];
rz(-0.83926187) q[0];
sx q[0];
rz(2.8150303) q[0];
rz(1.4549667) q[2];
sx q[2];
rz(-1.7421075) q[2];
sx q[2];
rz(2.0320323) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8453809) q[1];
sx q[1];
rz(-1.7239162) q[1];
sx q[1];
rz(-2.2232242) q[1];
rz(0.64063756) q[3];
sx q[3];
rz(-2.5036219) q[3];
sx q[3];
rz(0.20221329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3204331) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(1.2218366) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4959167) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(-2.3107279) q[0];
rz(1.9310541) q[1];
sx q[1];
rz(-1.6485063) q[1];
sx q[1];
rz(-0.99162203) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9308335) q[0];
sx q[0];
rz(-0.88085876) q[0];
sx q[0];
rz(-2.7372975) q[0];
rz(-pi) q[1];
rz(2.2575602) q[2];
sx q[2];
rz(-1.3515389) q[2];
sx q[2];
rz(-0.21603661) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8960921) q[1];
sx q[1];
rz(-0.84907167) q[1];
sx q[1];
rz(2.9338899) q[1];
x q[2];
rz(0.87144676) q[3];
sx q[3];
rz(-1.1394403) q[3];
sx q[3];
rz(-0.68701216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3811938) q[2];
sx q[2];
rz(-0.92279592) q[2];
sx q[2];
rz(2.7823616) q[2];
rz(2.4257816) q[3];
sx q[3];
rz(-2.3575213) q[3];
sx q[3];
rz(1.1904967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0642218) q[0];
sx q[0];
rz(-1.8697898) q[0];
sx q[0];
rz(-0.52128681) q[0];
rz(-2.298666) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(0.11046031) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2409637) q[0];
sx q[0];
rz(-0.74375737) q[0];
sx q[0];
rz(2.7636011) q[0];
x q[1];
rz(3.1119124) q[2];
sx q[2];
rz(-0.92845193) q[2];
sx q[2];
rz(-2.6386767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1350491) q[1];
sx q[1];
rz(-1.1056756) q[1];
sx q[1];
rz(-2.9198482) q[1];
rz(-pi) q[2];
x q[2];
rz(1.13917) q[3];
sx q[3];
rz(-2.5768498) q[3];
sx q[3];
rz(-0.016591681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9342186) q[2];
sx q[2];
rz(-1.5997581) q[2];
sx q[2];
rz(0.99384394) q[2];
rz(-2.5908616) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(-1.5229092) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2391424) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(0.19110876) q[0];
rz(0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-2.5083127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0201418) q[0];
sx q[0];
rz(-0.25853863) q[0];
sx q[0];
rz(-2.7516737) q[0];
rz(-2.6197144) q[2];
sx q[2];
rz(-2.3503135) q[2];
sx q[2];
rz(-0.14497862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7821473) q[1];
sx q[1];
rz(-1.7759062) q[1];
sx q[1];
rz(-0.19711776) q[1];
rz(-pi) q[2];
rz(2.9464989) q[3];
sx q[3];
rz(-2.4332402) q[3];
sx q[3];
rz(0.1699902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9621027) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(-0.38086677) q[2];
rz(1.3880091) q[3];
sx q[3];
rz(-2.1151147) q[3];
sx q[3];
rz(1.4306205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4717344) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(-1.5420472) q[0];
rz(-1.3767287) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(-2.210604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0536097) q[0];
sx q[0];
rz(-1.6585104) q[0];
sx q[0];
rz(-2.5045519) q[0];
x q[1];
rz(-2.687876) q[2];
sx q[2];
rz(-1.798011) q[2];
sx q[2];
rz(1.5308282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9500095) q[1];
sx q[1];
rz(-1.5693136) q[1];
sx q[1];
rz(-2.3133469) q[1];
rz(0.52102153) q[3];
sx q[3];
rz(-2.1590194) q[3];
sx q[3];
rz(0.38834342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57637438) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(3.0726688) q[2];
rz(-1.8779514) q[3];
sx q[3];
rz(-2.7694323) q[3];
sx q[3];
rz(2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7428335) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(-0.051890705) q[0];
rz(-2.9715111) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-2.7896519) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1387685) q[0];
sx q[0];
rz(-1.616298) q[0];
sx q[0];
rz(2.7753434) q[0];
x q[1];
rz(-2.5150256) q[2];
sx q[2];
rz(-2.0800903) q[2];
sx q[2];
rz(1.2207104) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70507732) q[1];
sx q[1];
rz(-1.962858) q[1];
sx q[1];
rz(-0.99832557) q[1];
rz(-1.7991583) q[3];
sx q[3];
rz(-2.8274837) q[3];
sx q[3];
rz(2.4921162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(-0.43668288) q[2];
rz(-0.39143482) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(-2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142414) q[0];
sx q[0];
rz(-0.7069718) q[0];
sx q[0];
rz(-0.11908764) q[0];
rz(-1.8424312) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(1.3577168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085310809) q[0];
sx q[0];
rz(-1.9765903) q[0];
sx q[0];
rz(0.88078518) q[0];
x q[1];
rz(-2.7804271) q[2];
sx q[2];
rz(-2.2884011) q[2];
sx q[2];
rz(-1.4360365) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57798856) q[1];
sx q[1];
rz(-1.23659) q[1];
sx q[1];
rz(-0.6329221) q[1];
rz(-2.1351027) q[3];
sx q[3];
rz(-0.83632937) q[3];
sx q[3];
rz(-1.8759954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9327717) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(-1.5116073) q[3];
sx q[3];
rz(-0.67174086) q[3];
sx q[3];
rz(3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83196249) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(3.0713785) q[2];
sx q[2];
rz(-1.4660942) q[2];
sx q[2];
rz(1.6906665) q[2];
rz(2.3336505) q[3];
sx q[3];
rz(-1.0026889) q[3];
sx q[3];
rz(-1.9849594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
