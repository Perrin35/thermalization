OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(-1.0426961) q[0];
sx q[0];
rz(0.27867499) q[0];
rz(1.1468118) q[1];
sx q[1];
rz(-0.52462259) q[1];
sx q[1];
rz(1.8705179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69158254) q[0];
sx q[0];
rz(-1.5890536) q[0];
sx q[0];
rz(-2.8810347) q[0];
rz(-pi) q[1];
rz(0.15058168) q[2];
sx q[2];
rz(-1.7007593) q[2];
sx q[2];
rz(2.5559363) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9624176) q[1];
sx q[1];
rz(-1.6500874) q[1];
sx q[1];
rz(-0.85722629) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1500312) q[3];
sx q[3];
rz(-2.09822) q[3];
sx q[3];
rz(0.19471951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.110454) q[2];
sx q[2];
rz(-1.217507) q[2];
sx q[2];
rz(0.26220599) q[2];
rz(-2.0860705) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(-0.087337703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156886) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(-2.9600034) q[0];
rz(-1.7407821) q[1];
sx q[1];
rz(-1.1927651) q[1];
sx q[1];
rz(2.3006732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5530295) q[0];
sx q[0];
rz(-0.32051495) q[0];
sx q[0];
rz(-3.0856087) q[0];
rz(-pi) q[1];
rz(2.6662331) q[2];
sx q[2];
rz(-0.98612758) q[2];
sx q[2];
rz(-0.66051403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5121756) q[1];
sx q[1];
rz(-1.7390774) q[1];
sx q[1];
rz(-0.03668935) q[1];
rz(-pi) q[2];
rz(2.0729468) q[3];
sx q[3];
rz(-0.6752033) q[3];
sx q[3];
rz(1.5218376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.759364) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(-2.0767029) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.8062785) q[3];
sx q[3];
rz(1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29838022) q[0];
sx q[0];
rz(-1.1792264) q[0];
sx q[0];
rz(-2.9425353) q[0];
rz(-0.82636034) q[1];
sx q[1];
rz(-0.91824707) q[1];
sx q[1];
rz(-0.77702776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8733002) q[0];
sx q[0];
rz(-1.4140633) q[0];
sx q[0];
rz(0.066989338) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7119646) q[2];
sx q[2];
rz(-2.0734534) q[2];
sx q[2];
rz(0.70870295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18622929) q[1];
sx q[1];
rz(-1.0252153) q[1];
sx q[1];
rz(-0.30382352) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8641508) q[3];
sx q[3];
rz(-2.146477) q[3];
sx q[3];
rz(2.4572008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0394773) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(-0.51304212) q[2];
rz(-2.3593864) q[3];
sx q[3];
rz(-1.9358044) q[3];
sx q[3];
rz(-1.7647083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.2143329) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(1.4904892) q[0];
rz(0.53228846) q[1];
sx q[1];
rz(-0.63096255) q[1];
sx q[1];
rz(-1.6040364) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89808116) q[0];
sx q[0];
rz(-1.6771206) q[0];
sx q[0];
rz(-0.16379078) q[0];
x q[1];
rz(-1.7236716) q[2];
sx q[2];
rz(-0.55143967) q[2];
sx q[2];
rz(1.3496292) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58504471) q[1];
sx q[1];
rz(-2.2050207) q[1];
sx q[1];
rz(-3.0536431) q[1];
rz(-pi) q[2];
rz(2.7521637) q[3];
sx q[3];
rz(-0.5618605) q[3];
sx q[3];
rz(-1.495468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.535546) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(0.01288506) q[2];
rz(1.9384725) q[3];
sx q[3];
rz(-1.67484) q[3];
sx q[3];
rz(-1.0217246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1383698) q[0];
sx q[0];
rz(-2.5165181) q[0];
sx q[0];
rz(1.3997929) q[0];
rz(-0.35203448) q[1];
sx q[1];
rz(-1.0540009) q[1];
sx q[1];
rz(-2.3037516) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96124803) q[0];
sx q[0];
rz(-2.367428) q[0];
sx q[0];
rz(1.7889687) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3532392) q[2];
sx q[2];
rz(-0.78561312) q[2];
sx q[2];
rz(2.2503302) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6988236) q[1];
sx q[1];
rz(-1.6976446) q[1];
sx q[1];
rz(0.5684828) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75148186) q[3];
sx q[3];
rz(-0.22322907) q[3];
sx q[3];
rz(-3.1187559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1167404) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.1836729) q[2];
rz(0.27070326) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(-1.5326327) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026641332) q[0];
sx q[0];
rz(-0.6237492) q[0];
sx q[0];
rz(2.5624516) q[0];
rz(1.9193513) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(2.0319669) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0298808) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(-0.16686186) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8366361) q[2];
sx q[2];
rz(-1.5169608) q[2];
sx q[2];
rz(0.75537813) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5803924) q[1];
sx q[1];
rz(-1.8186186) q[1];
sx q[1];
rz(2.7996254) q[1];
rz(-pi) q[2];
rz(0.81619451) q[3];
sx q[3];
rz(-0.69218721) q[3];
sx q[3];
rz(0.029267197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9455202) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(-0.35745364) q[2];
rz(-1.1899905) q[3];
sx q[3];
rz(-1.2001195) q[3];
sx q[3];
rz(-1.801871) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59915197) q[0];
sx q[0];
rz(-2.1245133) q[0];
sx q[0];
rz(-0.79208148) q[0];
rz(-0.7041086) q[1];
sx q[1];
rz(-2.6857565) q[1];
sx q[1];
rz(-1.5830931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38403758) q[0];
sx q[0];
rz(-0.4641986) q[0];
sx q[0];
rz(-3.1350101) q[0];
rz(-pi) q[1];
rz(0.81659813) q[2];
sx q[2];
rz(-2.3629521) q[2];
sx q[2];
rz(-1.1731847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4729135) q[1];
sx q[1];
rz(-1.9967419) q[1];
sx q[1];
rz(-1.5329453) q[1];
rz(2.9005403) q[3];
sx q[3];
rz(-2.3755382) q[3];
sx q[3];
rz(-2.0360006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.030297) q[2];
sx q[2];
rz(-2.3534677) q[2];
sx q[2];
rz(3.1374068) q[2];
rz(0.26675102) q[3];
sx q[3];
rz(-1.6175852) q[3];
sx q[3];
rz(-2.4193616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5608212) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(-0.15604493) q[0];
rz(-1.5849628) q[1];
sx q[1];
rz(-0.86235756) q[1];
sx q[1];
rz(0.59476605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9879919) q[0];
sx q[0];
rz(-0.56454851) q[0];
sx q[0];
rz(-0.49086824) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7960816) q[2];
sx q[2];
rz(-2.4191769) q[2];
sx q[2];
rz(1.142921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3799499) q[1];
sx q[1];
rz(-2.159286) q[1];
sx q[1];
rz(0.87169991) q[1];
rz(-2.4202926) q[3];
sx q[3];
rz(-2.2709322) q[3];
sx q[3];
rz(2.457778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5742089) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(-2.1136005) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-0.65516156) q[3];
sx q[3];
rz(-1.8067693) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26355711) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(-2.3433319) q[0];
rz(-2.3410666) q[1];
sx q[1];
rz(-0.48897484) q[1];
sx q[1];
rz(-0.54623234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9437708) q[0];
sx q[0];
rz(-1.6795625) q[0];
sx q[0];
rz(2.1740211) q[0];
rz(-pi) q[1];
rz(1.9990218) q[2];
sx q[2];
rz(-2.1167123) q[2];
sx q[2];
rz(-2.2466618) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6896587) q[1];
sx q[1];
rz(-2.0469189) q[1];
sx q[1];
rz(2.1228288) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30719884) q[3];
sx q[3];
rz(-0.23715487) q[3];
sx q[3];
rz(-3.0320771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.74786782) q[2];
sx q[2];
rz(-2.0515029) q[2];
sx q[2];
rz(-3.0785576) q[2];
rz(2.4437599) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49552712) q[0];
sx q[0];
rz(-2.263948) q[0];
sx q[0];
rz(-0.416042) q[0];
rz(1.7795732) q[1];
sx q[1];
rz(-1.1241309) q[1];
sx q[1];
rz(1.8488041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36609367) q[0];
sx q[0];
rz(-1.8387357) q[0];
sx q[0];
rz(-1.8084722) q[0];
rz(-pi) q[1];
rz(-1.2913537) q[2];
sx q[2];
rz(-2.0115174) q[2];
sx q[2];
rz(-2.8677577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67329183) q[1];
sx q[1];
rz(-2.3929962) q[1];
sx q[1];
rz(2.2842201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.148801) q[3];
sx q[3];
rz(-1.8099603) q[3];
sx q[3];
rz(-2.0668427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4585939) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(-0.90751737) q[2];
rz(1.2609743) q[3];
sx q[3];
rz(-2.5670299) q[3];
sx q[3];
rz(1.449409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64611971) q[0];
sx q[0];
rz(-1.6874122) q[0];
sx q[0];
rz(0.73943403) q[0];
rz(2.6932035) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(0.35198718) q[2];
sx q[2];
rz(-1.2766196) q[2];
sx q[2];
rz(2.5310493) q[2];
rz(-0.84239324) q[3];
sx q[3];
rz(-1.471023) q[3];
sx q[3];
rz(-1.4820549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
