OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57269078) q[0];
sx q[0];
rz(4.1299835) q[0];
sx q[0];
rz(9.873793) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(0.13154498) q[1];
sx q[1];
rz(8.2934525) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97206632) q[0];
sx q[0];
rz(-0.6837877) q[0];
sx q[0];
rz(-1.0975361) q[0];
rz(-pi) q[1];
rz(0.63458459) q[2];
sx q[2];
rz(-1.1401083) q[2];
sx q[2];
rz(-0.3542977) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4428569) q[1];
sx q[1];
rz(-1.7444102) q[1];
sx q[1];
rz(-0.97877494) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9692124) q[3];
sx q[3];
rz(-1.0339435) q[3];
sx q[3];
rz(-0.44042021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1863056) q[2];
sx q[2];
rz(-2.7489642) q[2];
sx q[2];
rz(0.10360959) q[2];
rz(-2.7747532) q[3];
sx q[3];
rz(-1.47374) q[3];
sx q[3];
rz(1.8240671) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1845301) q[0];
sx q[0];
rz(-0.4758895) q[0];
sx q[0];
rz(0.5262419) q[0];
rz(-1.2435675) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(-0.19613656) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58914069) q[0];
sx q[0];
rz(-1.6375443) q[0];
sx q[0];
rz(1.1132973) q[0];
rz(-pi) q[1];
rz(1.4043442) q[2];
sx q[2];
rz(-2.2734518) q[2];
sx q[2];
rz(-2.8748517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6013697) q[1];
sx q[1];
rz(-1.1246846) q[1];
sx q[1];
rz(-2.9838954) q[1];
rz(-0.011845592) q[3];
sx q[3];
rz(-1.5603754) q[3];
sx q[3];
rz(1.1529779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5474995) q[2];
sx q[2];
rz(-0.27442351) q[2];
sx q[2];
rz(-0.37626949) q[2];
rz(-2.2606692) q[3];
sx q[3];
rz(-1.3353525) q[3];
sx q[3];
rz(-0.81203619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45102099) q[0];
sx q[0];
rz(-0.32830992) q[0];
sx q[0];
rz(-2.1300533) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.1625544) q[1];
sx q[1];
rz(0.54723251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5315192) q[0];
sx q[0];
rz(-1.7721869) q[0];
sx q[0];
rz(-1.5904443) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7906495) q[2];
sx q[2];
rz(-1.3185274) q[2];
sx q[2];
rz(-0.44544912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7113674) q[1];
sx q[1];
rz(-2.2177296) q[1];
sx q[1];
rz(-0.30374668) q[1];
rz(-pi) q[2];
rz(0.8143592) q[3];
sx q[3];
rz(-1.7118771) q[3];
sx q[3];
rz(0.94179487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4554567) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(1.5020465) q[2];
rz(2.5655668) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(2.7801133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308924) q[0];
sx q[0];
rz(-2.3402813) q[0];
sx q[0];
rz(-0.33962387) q[0];
rz(0.54061186) q[1];
sx q[1];
rz(-2.4410591) q[1];
sx q[1];
rz(0.9300173) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7226573) q[0];
sx q[0];
rz(-1.6412853) q[0];
sx q[0];
rz(1.427729) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9564863) q[2];
sx q[2];
rz(-2.1197332) q[2];
sx q[2];
rz(1.5654636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2423289) q[1];
sx q[1];
rz(-1.5771958) q[1];
sx q[1];
rz(1.0914299) q[1];
rz(-pi) q[2];
rz(0.93149473) q[3];
sx q[3];
rz(-0.75041295) q[3];
sx q[3];
rz(1.9281333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12863079) q[2];
sx q[2];
rz(-1.7065115) q[2];
sx q[2];
rz(1.5677412) q[2];
rz(-0.74357998) q[3];
sx q[3];
rz(-2.3502374) q[3];
sx q[3];
rz(-0.96405205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4670694) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(-0.056644406) q[0];
rz(-1.665834) q[1];
sx q[1];
rz(-1.1328127) q[1];
sx q[1];
rz(-1.7830361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9525653) q[0];
sx q[0];
rz(-2.0533105) q[0];
sx q[0];
rz(2.5607462) q[0];
x q[1];
rz(-0.65746324) q[2];
sx q[2];
rz(-1.7320447) q[2];
sx q[2];
rz(-2.8828414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1093692) q[1];
sx q[1];
rz(-1.0894766) q[1];
sx q[1];
rz(2.4774423) q[1];
x q[2];
rz(-1.6886466) q[3];
sx q[3];
rz(-1.0123583) q[3];
sx q[3];
rz(-1.8359566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.27611) q[2];
sx q[2];
rz(-1.411974) q[2];
sx q[2];
rz(-1.126368) q[2];
rz(1.8051091) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(2.0520463) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4514076) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(0.13352808) q[0];
rz(-0.97767699) q[1];
sx q[1];
rz(-2.1656499) q[1];
sx q[1];
rz(-0.28688637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27659076) q[0];
sx q[0];
rz(-1.6185068) q[0];
sx q[0];
rz(-2.2722831) q[0];
rz(1.2426161) q[2];
sx q[2];
rz(-1.5300473) q[2];
sx q[2];
rz(2.2475257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9911954) q[1];
sx q[1];
rz(-0.6145454) q[1];
sx q[1];
rz(1.1769562) q[1];
rz(1.3135776) q[3];
sx q[3];
rz(-0.99204274) q[3];
sx q[3];
rz(2.677315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4213244) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.486091) q[2];
rz(-2.7739575) q[3];
sx q[3];
rz(-1.0272107) q[3];
sx q[3];
rz(-1.0953085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7169645) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(0.47384438) q[0];
rz(-2.4773856) q[1];
sx q[1];
rz(-1.9240446) q[1];
sx q[1];
rz(-1.7582105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6598005) q[0];
sx q[0];
rz(-2.0542025) q[0];
sx q[0];
rz(-2.3792335) q[0];
rz(1.5084615) q[2];
sx q[2];
rz(-0.41878715) q[2];
sx q[2];
rz(-0.58277786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7935087) q[1];
sx q[1];
rz(-1.7941107) q[1];
sx q[1];
rz(0.23621724) q[1];
rz(1.3455799) q[3];
sx q[3];
rz(-1.7449505) q[3];
sx q[3];
rz(2.1586777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15709269) q[2];
sx q[2];
rz(-1.4089156) q[2];
sx q[2];
rz(-0.3248997) q[2];
rz(-0.38069185) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(2.0438173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070504524) q[0];
sx q[0];
rz(-2.4779713) q[0];
sx q[0];
rz(2.2027503) q[0];
rz(-2.4414869) q[1];
sx q[1];
rz(-2.5036) q[1];
sx q[1];
rz(2.8525888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8696006) q[0];
sx q[0];
rz(-1.4933407) q[0];
sx q[0];
rz(-1.800022) q[0];
rz(-pi) q[1];
rz(2.1949057) q[2];
sx q[2];
rz(-2.3923229) q[2];
sx q[2];
rz(2.5890337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49456143) q[1];
sx q[1];
rz(-0.53003487) q[1];
sx q[1];
rz(-1.7575592) q[1];
rz(-1.6427755) q[3];
sx q[3];
rz(-1.1907309) q[3];
sx q[3];
rz(-0.38004181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8812022) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(-2.728906) q[2];
rz(2.2212501) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(1.9119561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.5145787) q[0];
sx q[0];
rz(-2.161442) q[0];
sx q[0];
rz(2.8421616) q[0];
rz(-2.0191655) q[1];
sx q[1];
rz(-1.9405245) q[1];
sx q[1];
rz(1.0312414) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2107687) q[0];
sx q[0];
rz(-1.5600191) q[0];
sx q[0];
rz(1.204654) q[0];
x q[1];
rz(-0.38487969) q[2];
sx q[2];
rz(-0.49321929) q[2];
sx q[2];
rz(-0.2156336) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6935307) q[1];
sx q[1];
rz(-1.281257) q[1];
sx q[1];
rz(-2.2953643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58166418) q[3];
sx q[3];
rz(-2.1278893) q[3];
sx q[3];
rz(2.795199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1624182) q[2];
sx q[2];
rz(-1.7981497) q[2];
sx q[2];
rz(2.8988885) q[2];
rz(-1.8414712) q[3];
sx q[3];
rz(-1.0588812) q[3];
sx q[3];
rz(2.9259031) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31371394) q[0];
sx q[0];
rz(-0.21610459) q[0];
sx q[0];
rz(-1.7207654) q[0];
rz(-2.8315262) q[1];
sx q[1];
rz(-0.87691751) q[1];
sx q[1];
rz(1.5975331) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102073) q[0];
sx q[0];
rz(-2.4500896) q[0];
sx q[0];
rz(-2.4639936) q[0];
x q[1];
rz(-1.3256959) q[2];
sx q[2];
rz(-0.35159207) q[2];
sx q[2];
rz(-1.0838255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1088562) q[1];
sx q[1];
rz(-1.3110771) q[1];
sx q[1];
rz(-1.3459567) q[1];
rz(-1.4306817) q[3];
sx q[3];
rz(-1.4129253) q[3];
sx q[3];
rz(2.637459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0073504) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(1.3664112) q[2];
rz(2.2640696) q[3];
sx q[3];
rz(-1.4656504) q[3];
sx q[3];
rz(-0.88959488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919745) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(-0.86434518) q[1];
sx q[1];
rz(-2.1280011) q[1];
sx q[1];
rz(-0.43307532) q[1];
rz(2.6411459) q[2];
sx q[2];
rz(-1.3896349) q[2];
sx q[2];
rz(-1.096772) q[2];
rz(3.0861985) q[3];
sx q[3];
rz(-2.0492036) q[3];
sx q[3];
rz(2.6072469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
