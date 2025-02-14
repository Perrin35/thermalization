OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3990134) q[0];
sx q[0];
rz(1.9785545) q[0];
sx q[0];
rz(7.4217441) q[0];
rz(-0.63602716) q[1];
sx q[1];
rz(-0.50995246) q[1];
sx q[1];
rz(0.62503254) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33390663) q[0];
sx q[0];
rz(-1.7117097) q[0];
sx q[0];
rz(-1.5265092) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1209645) q[2];
sx q[2];
rz(-1.0307023) q[2];
sx q[2];
rz(2.2079225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4506065) q[1];
sx q[1];
rz(-2.4391973) q[1];
sx q[1];
rz(-2.2390963) q[1];
rz(2.3287541) q[3];
sx q[3];
rz(-1.715797) q[3];
sx q[3];
rz(-2.6036711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2193489) q[2];
sx q[2];
rz(-0.85180989) q[2];
sx q[2];
rz(2.7395978) q[2];
rz(0.829202) q[3];
sx q[3];
rz(-1.3196557) q[3];
sx q[3];
rz(1.7792262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336693) q[0];
sx q[0];
rz(-2.5897554) q[0];
sx q[0];
rz(-2.0358987) q[0];
rz(-2.3528174) q[1];
sx q[1];
rz(-0.62216798) q[1];
sx q[1];
rz(0.50599352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087990847) q[0];
sx q[0];
rz(-1.8058386) q[0];
sx q[0];
rz(1.6490235) q[0];
rz(-2.3680192) q[2];
sx q[2];
rz(-0.94412747) q[2];
sx q[2];
rz(0.0093912436) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6044652) q[1];
sx q[1];
rz(-1.2187276) q[1];
sx q[1];
rz(0.2381937) q[1];
rz(-2.6384505) q[3];
sx q[3];
rz(-2.6330559) q[3];
sx q[3];
rz(-0.10756216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.904423) q[2];
sx q[2];
rz(-1.0639031) q[2];
sx q[2];
rz(0.87192956) q[2];
rz(-2.0875841) q[3];
sx q[3];
rz(-1.0796615) q[3];
sx q[3];
rz(-1.4409298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089040861) q[0];
sx q[0];
rz(-2.4854923) q[0];
sx q[0];
rz(1.22714) q[0];
rz(0.50651208) q[1];
sx q[1];
rz(-0.7917234) q[1];
sx q[1];
rz(-0.92481771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4445405) q[0];
sx q[0];
rz(-2.9612975) q[0];
sx q[0];
rz(0.96526115) q[0];
rz(-pi) q[1];
rz(0.13972762) q[2];
sx q[2];
rz(-2.681571) q[2];
sx q[2];
rz(-0.42528986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9021685) q[1];
sx q[1];
rz(-1.8730867) q[1];
sx q[1];
rz(2.2151674) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9496578) q[3];
sx q[3];
rz(-0.94516813) q[3];
sx q[3];
rz(-1.4257087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0623124) q[2];
sx q[2];
rz(-1.541905) q[2];
sx q[2];
rz(2.6463032) q[2];
rz(-1.371572) q[3];
sx q[3];
rz(-2.3973231) q[3];
sx q[3];
rz(-1.8194958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0168734) q[0];
sx q[0];
rz(-1.9833516) q[0];
sx q[0];
rz(-0.88791263) q[0];
rz(-1.3814231) q[1];
sx q[1];
rz(-2.1235178) q[1];
sx q[1];
rz(2.6240614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84758112) q[0];
sx q[0];
rz(-1.309762) q[0];
sx q[0];
rz(-3.1270967) q[0];
rz(1.04884) q[2];
sx q[2];
rz(-2.219355) q[2];
sx q[2];
rz(3.0364325) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0486517) q[1];
sx q[1];
rz(-3.042964) q[1];
sx q[1];
rz(1.6706549) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46948423) q[3];
sx q[3];
rz(-1.9133798) q[3];
sx q[3];
rz(2.1367578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.069328221) q[2];
sx q[2];
rz(-1.1257409) q[2];
sx q[2];
rz(-2.1605055) q[2];
rz(-0.4782933) q[3];
sx q[3];
rz(-2.7893453) q[3];
sx q[3];
rz(-0.52811629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9454055) q[0];
sx q[0];
rz(-1.7178752) q[0];
sx q[0];
rz(-0.032935306) q[0];
rz(1.6163274) q[1];
sx q[1];
rz(-2.567629) q[1];
sx q[1];
rz(2.2059435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55452915) q[0];
sx q[0];
rz(-1.1008223) q[0];
sx q[0];
rz(-1.5558467) q[0];
rz(-pi) q[1];
rz(-2.3439336) q[2];
sx q[2];
rz(-1.9211384) q[2];
sx q[2];
rz(-1.2454741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.074240597) q[1];
sx q[1];
rz(-0.79533333) q[1];
sx q[1];
rz(1.4834845) q[1];
rz(-pi) q[2];
rz(-0.58459063) q[3];
sx q[3];
rz(-1.3673615) q[3];
sx q[3];
rz(3.0082321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4719438) q[2];
sx q[2];
rz(-0.38413298) q[2];
sx q[2];
rz(1.498339) q[2];
rz(0.60254997) q[3];
sx q[3];
rz(-0.82585255) q[3];
sx q[3];
rz(1.92441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5947386) q[0];
sx q[0];
rz(-2.8475519) q[0];
sx q[0];
rz(-3.102741) q[0];
rz(-2.7815869) q[1];
sx q[1];
rz(-1.4119166) q[1];
sx q[1];
rz(-2.9927599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15880824) q[0];
sx q[0];
rz(-2.9624334) q[0];
sx q[0];
rz(0.78252234) q[0];
x q[1];
rz(-2.932933) q[2];
sx q[2];
rz(-0.57950117) q[2];
sx q[2];
rz(0.3267056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31735709) q[1];
sx q[1];
rz(-0.83401331) q[1];
sx q[1];
rz(-0.68781091) q[1];
x q[2];
rz(-2.5555796) q[3];
sx q[3];
rz(-1.8253606) q[3];
sx q[3];
rz(-2.3814122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0570602) q[2];
sx q[2];
rz(-0.82814211) q[2];
sx q[2];
rz(-2.7866936) q[2];
rz(1.0585632) q[3];
sx q[3];
rz(-0.94616977) q[3];
sx q[3];
rz(-0.53275776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8232089) q[0];
sx q[0];
rz(-1.4223149) q[0];
sx q[0];
rz(2.3329155) q[0];
rz(0.093756229) q[1];
sx q[1];
rz(-0.84051991) q[1];
sx q[1];
rz(2.7309928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79390271) q[0];
sx q[0];
rz(-1.9436033) q[0];
sx q[0];
rz(1.1469141) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8316531) q[2];
sx q[2];
rz(-0.93376389) q[2];
sx q[2];
rz(-2.2230679) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73170529) q[1];
sx q[1];
rz(-2.6123939) q[1];
sx q[1];
rz(0.7284109) q[1];
rz(-pi) q[2];
rz(-1.6924573) q[3];
sx q[3];
rz(-2.231519) q[3];
sx q[3];
rz(0.66189754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9984596) q[2];
sx q[2];
rz(-2.748558) q[2];
sx q[2];
rz(-0.03579363) q[2];
rz(-1.3068457) q[3];
sx q[3];
rz(-1.9948317) q[3];
sx q[3];
rz(-0.80865639) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82658139) q[0];
sx q[0];
rz(-2.1666574) q[0];
sx q[0];
rz(-2.5052729) q[0];
rz(0.66337216) q[1];
sx q[1];
rz(-1.1095108) q[1];
sx q[1];
rz(0.62796193) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3805035) q[0];
sx q[0];
rz(-3.0596943) q[0];
sx q[0];
rz(2.7837672) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6586893) q[2];
sx q[2];
rz(-1.1750487) q[2];
sx q[2];
rz(-1.7858693) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20013576) q[1];
sx q[1];
rz(-1.0147328) q[1];
sx q[1];
rz(-1.0764313) q[1];
x q[2];
rz(-0.18171715) q[3];
sx q[3];
rz(-2.9510806) q[3];
sx q[3];
rz(-2.8439034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5593354) q[2];
sx q[2];
rz(-0.11056837) q[2];
sx q[2];
rz(-1.7250693) q[2];
rz(2.0420117) q[3];
sx q[3];
rz(-1.0397592) q[3];
sx q[3];
rz(0.85272461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002554) q[0];
sx q[0];
rz(-0.89109963) q[0];
sx q[0];
rz(0.6148327) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-0.95242396) q[1];
sx q[1];
rz(0.27613786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9671738) q[0];
sx q[0];
rz(-1.9959269) q[0];
sx q[0];
rz(-2.9993254) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3036518) q[2];
sx q[2];
rz(-0.54590271) q[2];
sx q[2];
rz(-2.5747908) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8135173) q[1];
sx q[1];
rz(-2.349252) q[1];
sx q[1];
rz(-2.4787419) q[1];
rz(-pi) q[2];
rz(-2.7998447) q[3];
sx q[3];
rz(-2.2015988) q[3];
sx q[3];
rz(-0.83021008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.189956) q[2];
sx q[2];
rz(-2.1026244) q[2];
sx q[2];
rz(3.0432126) q[2];
rz(1.9010057) q[3];
sx q[3];
rz(-1.0616579) q[3];
sx q[3];
rz(3.0978751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8561309) q[0];
sx q[0];
rz(-2.0933445) q[0];
sx q[0];
rz(2.7833126) q[0];
rz(1.0048535) q[1];
sx q[1];
rz(-0.88468164) q[1];
sx q[1];
rz(0.34214941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44624871) q[0];
sx q[0];
rz(-0.98606743) q[0];
sx q[0];
rz(1.8323932) q[0];
x q[1];
rz(-2.1120295) q[2];
sx q[2];
rz(-2.4839253) q[2];
sx q[2];
rz(-2.5785411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.085386203) q[1];
sx q[1];
rz(-1.4634988) q[1];
sx q[1];
rz(0.57319586) q[1];
rz(-pi) q[2];
rz(1.2827286) q[3];
sx q[3];
rz(-0.34980845) q[3];
sx q[3];
rz(1.6988848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7973914) q[2];
sx q[2];
rz(-2.6579393) q[2];
sx q[2];
rz(-2.7698216) q[2];
rz(-0.19112912) q[3];
sx q[3];
rz(-1.976795) q[3];
sx q[3];
rz(2.718486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38901781) q[0];
sx q[0];
rz(-2.6448463) q[0];
sx q[0];
rz(-0.65506558) q[0];
rz(0.3599421) q[1];
sx q[1];
rz(-1.5725726) q[1];
sx q[1];
rz(-1.6307065) q[1];
rz(-3.0790569) q[2];
sx q[2];
rz(-1.8725431) q[2];
sx q[2];
rz(0.098081577) q[2];
rz(2.3287931) q[3];
sx q[3];
rz(-1.5817196) q[3];
sx q[3];
rz(-2.055837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
