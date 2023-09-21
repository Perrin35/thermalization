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
rz(-3.1414519) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407335) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(-0.1851693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46618669) q[2];
sx q[2];
rz(-2.5417915) q[2];
sx q[2];
rz(-2.8592062) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27543435) q[1];
sx q[1];
rz(-0.94238102) q[1];
sx q[1];
rz(-2.1584312) q[1];
x q[2];
rz(1.3532577) q[3];
sx q[3];
rz(-1.6780403) q[3];
sx q[3];
rz(1.4835964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-1.123463) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2259953) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(2.0321839) q[0];
rz(-3.0779755) q[2];
sx q[2];
rz(-0.78084313) q[2];
sx q[2];
rz(2.7240679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9065735) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(0.97096918) q[1];
rz(2.1334322) q[3];
sx q[3];
rz(-2.1481272) q[3];
sx q[3];
rz(-1.8562223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27750257) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.8664237) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(2.0085874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4988574) q[0];
sx q[0];
rz(-0.14232902) q[0];
sx q[0];
rz(1.5940773) q[0];
rz(-pi) q[1];
rz(0.50425597) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(-2.2222663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6269835) q[1];
sx q[1];
rz(-1.2830462) q[1];
sx q[1];
rz(2.1327553) q[1];
rz(-pi) q[2];
rz(-0.95136178) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(-1.149328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2599729) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(0.89598957) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(0.13555759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7134705) q[0];
sx q[0];
rz(-0.8849511) q[0];
sx q[0];
rz(1.1629521) q[0];
rz(1.9549471) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(-2.4406976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.63771026) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(1.3824944) q[1];
rz(-pi) q[2];
rz(2.0849864) q[3];
sx q[3];
rz(-0.31294542) q[3];
sx q[3];
rz(-1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9049412) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(0.044163477) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-1.0132382) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(1.0838881) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.097682) q[0];
sx q[0];
rz(-1.8143166) q[0];
sx q[0];
rz(-0.18885141) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8439699) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(-2.4730686) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.52884) q[1];
sx q[1];
rz(-1.6788947) q[1];
sx q[1];
rz(-1.0428863) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6211987) q[3];
sx q[3];
rz(-1.0574697) q[3];
sx q[3];
rz(-2.7798775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23285398) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58127922) q[0];
sx q[0];
rz(-2.7972097) q[0];
sx q[0];
rz(3.0292105) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8199811) q[2];
sx q[2];
rz(-0.65882896) q[2];
sx q[2];
rz(-0.67473251) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6678644) q[1];
sx q[1];
rz(-1.6147991) q[1];
sx q[1];
rz(-1.4021137) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1392519) q[3];
sx q[3];
rz(-1.4961092) q[3];
sx q[3];
rz(-0.60126388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.8849467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.885868) q[0];
sx q[0];
rz(-0.2547383) q[0];
sx q[0];
rz(-1.4101009) q[0];
rz(-pi) q[1];
rz(-1.9728327) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(1.4025584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9650967) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(-0.22775905) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2506966) q[3];
sx q[3];
rz(-2.2543636) q[3];
sx q[3];
rz(-0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(1.3640277) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(-1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-0.38696188) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68650866) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(1.9576661) q[0];
rz(2.6110821) q[2];
sx q[2];
rz(-2.1973917) q[2];
sx q[2];
rz(-0.26041398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-0.67589251) q[1];
sx q[1];
rz(-3.0216316) q[1];
rz(-2.7510795) q[3];
sx q[3];
rz(-0.53005855) q[3];
sx q[3];
rz(-2.183225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14116645) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6599777) q[0];
sx q[0];
rz(-1.0064631) q[0];
sx q[0];
rz(-2.3100287) q[0];
rz(-2.1575035) q[2];
sx q[2];
rz(-2.3684635) q[2];
sx q[2];
rz(-3.0423321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8763435) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(2.6190119) q[1];
rz(-pi) q[2];
rz(0.52600577) q[3];
sx q[3];
rz(-2.3283539) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.07842841) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(2.06185) q[0];
rz(-1.059277) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.7396897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3831543) q[0];
sx q[0];
rz(-2.8935195) q[0];
sx q[0];
rz(1.0323348) q[0];
rz(-2.9981668) q[2];
sx q[2];
rz(-0.18866814) q[2];
sx q[2];
rz(2.8667237) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55587308) q[1];
sx q[1];
rz(-1.8750637) q[1];
sx q[1];
rz(0.85822206) q[1];
rz(0.63125061) q[3];
sx q[3];
rz(-1.7389986) q[3];
sx q[3];
rz(1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.520291) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14810066) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(-0.91611721) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-2.312071) q[2];
sx q[2];
rz(-2.1567232) q[2];
sx q[2];
rz(0.17270252) q[2];
rz(-1.0536853) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];