OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(-0.21259354) q[0];
sx q[0];
rz(-2.4074182) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(2.7956378) q[1];
sx q[1];
rz(12.591127) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2233362) q[0];
sx q[0];
rz(-1.5084711) q[0];
sx q[0];
rz(1.5274836) q[0];
rz(-1.3369629) q[2];
sx q[2];
rz(-1.4604092) q[2];
sx q[2];
rz(-1.8821007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5879257) q[1];
sx q[1];
rz(-0.30372639) q[1];
sx q[1];
rz(2.9656253) q[1];
rz(-pi) q[2];
rz(1.3583899) q[3];
sx q[3];
rz(-1.5161361) q[3];
sx q[3];
rz(-2.2509991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88432246) q[2];
sx q[2];
rz(-1.6486282) q[2];
sx q[2];
rz(-2.8327668) q[2];
rz(-2.3158) q[3];
sx q[3];
rz(-1.0079931) q[3];
sx q[3];
rz(-2.3776313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.9708213) q[0];
sx q[0];
rz(-1.9124557) q[0];
sx q[0];
rz(2.1139076) q[0];
rz(0.81614256) q[1];
sx q[1];
rz(-1.1183389) q[1];
sx q[1];
rz(-0.81248409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4554169) q[0];
sx q[0];
rz(-2.7551965) q[0];
sx q[0];
rz(-2.0040345) q[0];
rz(2.7682253) q[2];
sx q[2];
rz(-1.5840216) q[2];
sx q[2];
rz(0.91886273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3025488) q[1];
sx q[1];
rz(-1.7086281) q[1];
sx q[1];
rz(1.0562357) q[1];
x q[2];
rz(-0.3489209) q[3];
sx q[3];
rz(-1.9066208) q[3];
sx q[3];
rz(0.15062697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5673148) q[2];
sx q[2];
rz(-1.4894166) q[2];
sx q[2];
rz(-2.2226649) q[2];
rz(2.610176) q[3];
sx q[3];
rz(-2.0960977) q[3];
sx q[3];
rz(-3.1004068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.70319217) q[0];
sx q[0];
rz(-2.0719318) q[0];
sx q[0];
rz(-2.0042787) q[0];
rz(1.7510341) q[1];
sx q[1];
rz(-0.5568234) q[1];
sx q[1];
rz(-2.4086319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61444401) q[0];
sx q[0];
rz(-1.9839459) q[0];
sx q[0];
rz(-1.2507417) q[0];
x q[1];
rz(-2.175019) q[2];
sx q[2];
rz(-1.086965) q[2];
sx q[2];
rz(-1.0072198) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.809124) q[1];
sx q[1];
rz(-2.4753503) q[1];
sx q[1];
rz(-0.82102832) q[1];
x q[2];
rz(1.5844581) q[3];
sx q[3];
rz(-0.44876305) q[3];
sx q[3];
rz(0.45873935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1375894) q[2];
sx q[2];
rz(-1.4121476) q[2];
sx q[2];
rz(1.0388733) q[2];
rz(-1.002257) q[3];
sx q[3];
rz(-1.8398617) q[3];
sx q[3];
rz(2.8514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20659474) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(-2.2344053) q[0];
rz(-2.9838003) q[1];
sx q[1];
rz(-1.0148427) q[1];
sx q[1];
rz(1.8603604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8195551) q[0];
sx q[0];
rz(-1.0519618) q[0];
sx q[0];
rz(-2.4253393) q[0];
rz(2.2864986) q[2];
sx q[2];
rz(-1.2629328) q[2];
sx q[2];
rz(1.5474873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80824796) q[1];
sx q[1];
rz(-0.95511694) q[1];
sx q[1];
rz(2.0991825) q[1];
rz(-2.5077742) q[3];
sx q[3];
rz(-1.0019284) q[3];
sx q[3];
rz(-0.19396362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8944051) q[2];
sx q[2];
rz(-0.93623585) q[2];
sx q[2];
rz(-1.8518764) q[2];
rz(-1.3442518) q[3];
sx q[3];
rz(-1.4569837) q[3];
sx q[3];
rz(-2.4414506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6808788) q[0];
sx q[0];
rz(-0.49837708) q[0];
sx q[0];
rz(-2.1233249) q[0];
rz(1.1030039) q[1];
sx q[1];
rz(-0.7119199) q[1];
sx q[1];
rz(1.3404554) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0809652) q[0];
sx q[0];
rz(-0.61765352) q[0];
sx q[0];
rz(0.65184848) q[0];
rz(-1.4730155) q[2];
sx q[2];
rz(-1.9488397) q[2];
sx q[2];
rz(2.5387272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75234883) q[1];
sx q[1];
rz(-1.443814) q[1];
sx q[1];
rz(-1.4402706) q[1];
rz(1.4748853) q[3];
sx q[3];
rz(-0.54236327) q[3];
sx q[3];
rz(-0.23254843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5995784) q[2];
sx q[2];
rz(-2.2530563) q[2];
sx q[2];
rz(-2.1395903) q[2];
rz(-2.4501948) q[3];
sx q[3];
rz(-1.1804429) q[3];
sx q[3];
rz(1.6208167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5387251) q[0];
sx q[0];
rz(-2.0724917) q[0];
sx q[0];
rz(-2.5163203) q[0];
rz(-1.9484776) q[1];
sx q[1];
rz(-0.68980491) q[1];
sx q[1];
rz(-0.56732059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.389844) q[0];
sx q[0];
rz(-1.746467) q[0];
sx q[0];
rz(-2.3542464) q[0];
rz(-pi) q[1];
rz(-1.5626603) q[2];
sx q[2];
rz(-2.3971791) q[2];
sx q[2];
rz(-1.2071963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.76657) q[1];
sx q[1];
rz(-2.5225897) q[1];
sx q[1];
rz(1.6155711) q[1];
rz(-pi) q[2];
rz(-2.3472957) q[3];
sx q[3];
rz(-2.8927589) q[3];
sx q[3];
rz(0.90422309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.304473) q[2];
sx q[2];
rz(-1.8581055) q[2];
sx q[2];
rz(-1.5597957) q[2];
rz(2.5290153) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(2.9858203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0462069) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(2.0147391) q[0];
rz(1.2696179) q[1];
sx q[1];
rz(-1.9536628) q[1];
sx q[1];
rz(-2.6205305) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714037) q[0];
sx q[0];
rz(-2.2260283) q[0];
sx q[0];
rz(-2.6960424) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1087569) q[2];
sx q[2];
rz(-1.2419789) q[2];
sx q[2];
rz(1.1177899) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46494486) q[1];
sx q[1];
rz(-1.0527234) q[1];
sx q[1];
rz(0.89333109) q[1];
x q[2];
rz(-1.295213) q[3];
sx q[3];
rz(-1.0465989) q[3];
sx q[3];
rz(-0.39488246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0236987) q[2];
sx q[2];
rz(-1.7774899) q[2];
sx q[2];
rz(-1.2269616) q[2];
rz(-1.6795233) q[3];
sx q[3];
rz(-1.6242124) q[3];
sx q[3];
rz(-2.7000361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284477) q[0];
sx q[0];
rz(-1.0297091) q[0];
sx q[0];
rz(2.5944769) q[0];
rz(3.0534577) q[1];
sx q[1];
rz(-1.3476177) q[1];
sx q[1];
rz(2.7190582) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8582701) q[0];
sx q[0];
rz(-1.2898603) q[0];
sx q[0];
rz(0.022780943) q[0];
x q[1];
rz(-0.59178517) q[2];
sx q[2];
rz(-2.5349504) q[2];
sx q[2];
rz(3.0559412) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48848885) q[1];
sx q[1];
rz(-3.0526027) q[1];
sx q[1];
rz(2.2527534) q[1];
rz(-pi) q[2];
x q[2];
rz(3.001207) q[3];
sx q[3];
rz(-1.3891798) q[3];
sx q[3];
rz(2.3395305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9162468) q[2];
sx q[2];
rz(-1.5044745) q[2];
sx q[2];
rz(2.8543191) q[2];
rz(-0.54287994) q[3];
sx q[3];
rz(-0.77016872) q[3];
sx q[3];
rz(3.0834901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420237) q[0];
sx q[0];
rz(-2.4574807) q[0];
sx q[0];
rz(-2.3642819) q[0];
rz(2.2151392) q[1];
sx q[1];
rz(-1.1421721) q[1];
sx q[1];
rz(1.3949589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762601) q[0];
sx q[0];
rz(-1.9585573) q[0];
sx q[0];
rz(-0.34404018) q[0];
rz(-pi) q[1];
rz(1.821076) q[2];
sx q[2];
rz(-2.1553073) q[2];
sx q[2];
rz(-2.3601041) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3923126) q[1];
sx q[1];
rz(-2.5651076) q[1];
sx q[1];
rz(-2.2063401) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5245669) q[3];
sx q[3];
rz(-1.0279547) q[3];
sx q[3];
rz(1.9844685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42561439) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(2.051029) q[2];
rz(0.96450949) q[3];
sx q[3];
rz(-1.0114074) q[3];
sx q[3];
rz(-1.9264268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5665117) q[0];
sx q[0];
rz(-2.8031741) q[0];
sx q[0];
rz(1.6315208) q[0];
rz(0.034320023) q[1];
sx q[1];
rz(-1.8020554) q[1];
sx q[1];
rz(0.43201772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1572947) q[0];
sx q[0];
rz(-2.108886) q[0];
sx q[0];
rz(0.068416083) q[0];
x q[1];
rz(-1.9218742) q[2];
sx q[2];
rz(-1.6475999) q[2];
sx q[2];
rz(1.5157645) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0610785) q[1];
sx q[1];
rz(-0.17560683) q[1];
sx q[1];
rz(-2.5422515) q[1];
rz(-pi) q[2];
rz(-1.179078) q[3];
sx q[3];
rz(-0.61810247) q[3];
sx q[3];
rz(-2.6276692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67046514) q[2];
sx q[2];
rz(-1.9894783) q[2];
sx q[2];
rz(1.0738037) q[2];
rz(-2.9527169) q[3];
sx q[3];
rz(-0.60868588) q[3];
sx q[3];
rz(2.7591738) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043924532) q[0];
sx q[0];
rz(-0.32079874) q[0];
sx q[0];
rz(-0.70377845) q[0];
rz(-1.7882998) q[1];
sx q[1];
rz(-0.67519338) q[1];
sx q[1];
rz(2.304362) q[1];
rz(-0.7424217) q[2];
sx q[2];
rz(-1.0992285) q[2];
sx q[2];
rz(0.64115094) q[2];
rz(0.13808098) q[3];
sx q[3];
rz(-0.46317536) q[3];
sx q[3];
rz(-2.4899766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
