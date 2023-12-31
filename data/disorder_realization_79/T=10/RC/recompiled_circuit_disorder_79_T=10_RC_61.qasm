OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0177512) q[0];
sx q[0];
rz(-3.0525644) q[0];
sx q[0];
rz(0.71319367) q[0];
rz(-pi) q[1];
rz(-0.97857742) q[2];
sx q[2];
rz(-2.7263612) q[2];
sx q[2];
rz(0.65939553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.623917) q[1];
sx q[1];
rz(-2.3530934) q[1];
sx q[1];
rz(-2.4967525) q[1];
rz(-pi) q[2];
rz(-2.7636823) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.6764486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10509051) q[0];
sx q[0];
rz(-0.72421342) q[0];
sx q[0];
rz(3.1378531) q[0];
rz(-pi) q[1];
rz(0.25603489) q[2];
sx q[2];
rz(-2.0453774) q[2];
sx q[2];
rz(1.7325967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0959024) q[1];
sx q[1];
rz(-1.1046788) q[1];
sx q[1];
rz(-2.5065266) q[1];
rz(2.39605) q[3];
sx q[3];
rz(-1.8288444) q[3];
sx q[3];
rz(0.52449709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-2.3410472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24414177) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(-1.98154) q[0];
x q[1];
rz(2.4750701) q[2];
sx q[2];
rz(-0.98072532) q[2];
sx q[2];
rz(2.9556264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87677466) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(-3.1289711) q[1];
x q[2];
rz(1.21739) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67733726) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-2.0565313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5497919) q[0];
sx q[0];
rz(-1.8013957) q[0];
sx q[0];
rz(2.6016298) q[0];
rz(1.9374574) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(1.0166849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76772056) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(-2.9994681) q[1];
x q[2];
rz(-1.163108) q[3];
sx q[3];
rz(-2.1424322) q[3];
sx q[3];
rz(-2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(-0.51182169) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(0.98168215) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3571346) q[0];
sx q[0];
rz(-2.031209) q[0];
sx q[0];
rz(-0.61607342) q[0];
rz(0.15820299) q[2];
sx q[2];
rz(-1.8100097) q[2];
sx q[2];
rz(1.6895134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9202068) q[1];
sx q[1];
rz(-1.1077987) q[1];
sx q[1];
rz(1.4057926) q[1];
x q[2];
rz(0.074023789) q[3];
sx q[3];
rz(-0.97462666) q[3];
sx q[3];
rz(-1.0130458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(2.9763124) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.107347) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(-2.8427567) q[0];
x q[1];
rz(1.500962) q[2];
sx q[2];
rz(-0.88431057) q[2];
sx q[2];
rz(-1.3697461) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11152553) q[1];
sx q[1];
rz(-1.0377874) q[1];
sx q[1];
rz(2.0433321) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23541707) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1897141) q[0];
sx q[0];
rz(-1.5218381) q[0];
sx q[0];
rz(-0.017107054) q[0];
rz(-pi) q[1];
rz(-1.2320802) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(0.92600694) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4941102) q[1];
sx q[1];
rz(-2.3586015) q[1];
sx q[1];
rz(2.8857735) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50076671) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-2.4367743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0050469) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4066276) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3085262) q[2];
sx q[2];
rz(-1.3563915) q[2];
sx q[2];
rz(-3.0467141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8340048) q[1];
sx q[1];
rz(-1.3655791) q[1];
sx q[1];
rz(2.8819041) q[1];
x q[2];
rz(1.7796302) q[3];
sx q[3];
rz(-0.71912557) q[3];
sx q[3];
rz(-1.5987087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-2.4618861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4363791) q[0];
sx q[0];
rz(-2.0853015) q[0];
sx q[0];
rz(-1.3374469) q[0];
rz(-pi) q[1];
rz(-1.7858511) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(-3.0422473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9437127) q[1];
sx q[1];
rz(-2.3730952) q[1];
sx q[1];
rz(-2.6430623) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9548095) q[3];
sx q[3];
rz(-0.98815742) q[3];
sx q[3];
rz(0.59347502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(-2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(-1.1766599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51417527) q[0];
sx q[0];
rz(-1.6192993) q[0];
sx q[0];
rz(-1.4072627) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54729692) q[2];
sx q[2];
rz(-1.0101057) q[2];
sx q[2];
rz(-1.0123569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1365876) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(-3.0535166) q[1];
rz(-pi) q[2];
rz(-1.482974) q[3];
sx q[3];
rz(-1.0613958) q[3];
sx q[3];
rz(-1.0305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(2.0991142) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-0.30504967) q[2];
sx q[2];
rz(-0.76022824) q[2];
sx q[2];
rz(-0.40679731) q[2];
rz(-3.1191961) q[3];
sx q[3];
rz(-0.35084421) q[3];
sx q[3];
rz(2.4435333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
