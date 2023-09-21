OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407335) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(-2.9564234) q[0];
rz(-0.5483746) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(-0.89474364) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6686033) q[1];
sx q[1];
rz(-2.0358634) q[1];
sx q[1];
rz(2.4238062) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1080997) q[3];
sx q[3];
rz(-2.8994312) q[3];
sx q[3];
rz(-0.36377454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(-1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(-2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-2.0181296) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7682122) q[0];
sx q[0];
rz(-2.0313615) q[0];
sx q[0];
rz(3.0773613) q[0];
rz(-1.5078817) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(-0.5069678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9532721) q[1];
sx q[1];
rz(-0.76474944) q[1];
sx q[1];
rz(0.79337593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0081604) q[3];
sx q[3];
rz(-0.99346549) q[3];
sx q[3];
rz(-1.8562223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(2.4675026) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.6154217) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-2.0085874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4988574) q[0];
sx q[0];
rz(-2.9992636) q[0];
sx q[0];
rz(-1.5940773) q[0];
rz(-1.0622382) q[2];
sx q[2];
rz(-2.2937751) q[2];
sx q[2];
rz(1.522097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23302151) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(0.33645333) q[1];
rz(-pi) q[2];
rz(2.1902309) q[3];
sx q[3];
rz(-1.0361443) q[3];
sx q[3];
rz(1.149328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8901849) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(-1.2934925) q[2];
rz(-0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(-2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(0.13555759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7134705) q[0];
sx q[0];
rz(-2.2566416) q[0];
sx q[0];
rz(1.1629521) q[0];
x q[1];
rz(-0.074684871) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(-0.89786868) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5038824) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(-1.3824944) q[1];
rz(-0.1578219) q[3];
sx q[3];
rz(-1.8421679) q[3];
sx q[3];
rz(1.8872758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(-0.044163477) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(1.0838881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7679374) q[0];
sx q[0];
rz(-0.3070139) q[0];
sx q[0];
rz(0.92371773) q[0];
rz(-0.25123698) q[2];
sx q[2];
rz(-1.305797) q[2];
sx q[2];
rz(-2.3064409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0368082) q[1];
sx q[1];
rz(-2.0953062) q[1];
sx q[1];
rz(-0.12496897) q[1];
rz(-pi) q[2];
rz(-1.5203939) q[3];
sx q[3];
rz(-1.0574697) q[3];
sx q[3];
rz(-0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-3.0474512) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0953513) q[0];
sx q[0];
rz(-1.5329251) q[0];
sx q[0];
rz(0.34237679) q[0];
rz(-1.3307829) q[2];
sx q[2];
rz(-2.1905106) q[2];
sx q[2];
rz(2.0680708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9858866) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(-1.8272912) q[1];
rz(-pi) q[2];
x q[2];
rz(1.496109) q[3];
sx q[3];
rz(-1.5731305) q[3];
sx q[3];
rz(-0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(0.41263321) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.2566459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518258) q[0];
sx q[0];
rz(-1.3194114) q[0];
sx q[0];
rz(3.0999523) q[0];
x q[1];
rz(2.0357382) q[2];
sx q[2];
rz(-1.3824116) q[2];
sx q[2];
rz(-0.18907324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17649594) q[1];
sx q[1];
rz(-0.28897044) q[1];
sx q[1];
rz(0.22775905) q[1];
x q[2];
rz(-0.7092181) q[3];
sx q[3];
rz(-1.8172482) q[3];
sx q[3];
rz(2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(-1.9576661) q[0];
rz(-pi) q[1];
rz(2.6110821) q[2];
sx q[2];
rz(-0.94420099) q[2];
sx q[2];
rz(-2.8811787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(3.0216316) q[1];
rz(-0.39051315) q[3];
sx q[3];
rz(-2.6115341) q[3];
sx q[3];
rz(0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-2.7015838) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14116645) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(-2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(3.0922906) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5987199) q[0];
sx q[0];
rz(-0.96519404) q[0];
sx q[0];
rz(0.70830317) q[0];
rz(-pi) q[1];
rz(2.2531613) q[2];
sx q[2];
rz(-1.1738136) q[2];
sx q[2];
rz(-1.9156485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2652492) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(2.6190119) q[1];
rz(-pi) q[2];
rz(-0.74069549) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(1.5965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(-2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.7396897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7125177) q[0];
sx q[0];
rz(-1.6970465) q[0];
sx q[0];
rz(1.35668) q[0];
rz(-pi) q[1];
rz(-2.9981668) q[2];
sx q[2];
rz(-0.18866814) q[2];
sx q[2];
rz(2.8667237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8733752) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(2.7482277) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.510342) q[3];
sx q[3];
rz(-1.7389986) q[3];
sx q[3];
rz(-1.2390618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(-1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(2.4441392) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.993492) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(0.79402906) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(-1.6021452) q[3];
sx q[3];
rz(-0.51732224) q[3];
sx q[3];
rz(2.8696685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
