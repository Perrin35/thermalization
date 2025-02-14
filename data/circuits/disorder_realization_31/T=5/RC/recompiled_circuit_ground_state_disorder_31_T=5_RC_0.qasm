OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6346729) q[0];
sx q[0];
rz(-1.7393751) q[0];
sx q[0];
rz(2.3911067) q[0];
rz(0.063272417) q[1];
sx q[1];
rz(7.1019389) q[1];
sx q[1];
rz(11.547312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1602331) q[0];
sx q[0];
rz(-0.59688127) q[0];
sx q[0];
rz(2.802538) q[0];
rz(2.8739863) q[2];
sx q[2];
rz(-2.2626167) q[2];
sx q[2];
rz(-2.0694422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.99996725) q[1];
sx q[1];
rz(-0.5988752) q[1];
sx q[1];
rz(-0.24354045) q[1];
rz(-pi) q[2];
rz(0.54368162) q[3];
sx q[3];
rz(-1.8988109) q[3];
sx q[3];
rz(-1.8193434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84746209) q[2];
sx q[2];
rz(-2.1475466) q[2];
sx q[2];
rz(1.4685941) q[2];
rz(0.20644203) q[3];
sx q[3];
rz(-1.3279746) q[3];
sx q[3];
rz(-2.4895721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.9701397) q[0];
sx q[0];
rz(-1.7296706) q[0];
sx q[0];
rz(3.0385802) q[0];
rz(0.39762321) q[1];
sx q[1];
rz(-2.7151974) q[1];
sx q[1];
rz(-2.2893589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6838339) q[0];
sx q[0];
rz(-1.6310098) q[0];
sx q[0];
rz(1.0369861) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2322517) q[2];
sx q[2];
rz(-1.8391092) q[2];
sx q[2];
rz(-1.2635096) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6247762) q[1];
sx q[1];
rz(-2.0943536) q[1];
sx q[1];
rz(-0.015990301) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2198803) q[3];
sx q[3];
rz(-0.19769719) q[3];
sx q[3];
rz(1.3330595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50920606) q[2];
sx q[2];
rz(-2.1337324) q[2];
sx q[2];
rz(-0.81727916) q[2];
rz(-2.9362074) q[3];
sx q[3];
rz(-2.9257264) q[3];
sx q[3];
rz(0.5425905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83194724) q[0];
sx q[0];
rz(-1.0878071) q[0];
sx q[0];
rz(-0.52016869) q[0];
rz(-2.4179516) q[1];
sx q[1];
rz(-1.2836722) q[1];
sx q[1];
rz(-2.9626194) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1083281) q[0];
sx q[0];
rz(-1.1305362) q[0];
sx q[0];
rz(-0.56301337) q[0];
x q[1];
rz(1.7759982) q[2];
sx q[2];
rz(-1.2276863) q[2];
sx q[2];
rz(0.24977906) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3368028) q[1];
sx q[1];
rz(-1.8095142) q[1];
sx q[1];
rz(-1.6126339) q[1];
rz(-0.34647219) q[3];
sx q[3];
rz(-2.4541509) q[3];
sx q[3];
rz(-2.8137852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2554864) q[2];
sx q[2];
rz(-0.28033689) q[2];
sx q[2];
rz(-2.172016) q[2];
rz(-2.7995836) q[3];
sx q[3];
rz(-1.0229144) q[3];
sx q[3];
rz(1.8146993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.904357) q[0];
sx q[0];
rz(-2.3396753) q[0];
sx q[0];
rz(3.0661769) q[0];
rz(2.8341809) q[1];
sx q[1];
rz(-1.1630029) q[1];
sx q[1];
rz(-2.7154162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8318752) q[0];
sx q[0];
rz(-2.2906036) q[0];
sx q[0];
rz(-0.72711522) q[0];
rz(-0.51898848) q[2];
sx q[2];
rz(-0.52757571) q[2];
sx q[2];
rz(1.0949539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4170917) q[1];
sx q[1];
rz(-2.1414653) q[1];
sx q[1];
rz(0.23391926) q[1];
rz(-pi) q[2];
rz(3.0846023) q[3];
sx q[3];
rz(-2.7354623) q[3];
sx q[3];
rz(2.1907091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23224607) q[2];
sx q[2];
rz(-1.987395) q[2];
sx q[2];
rz(-2.2171891) q[2];
rz(2.038548) q[3];
sx q[3];
rz(-1.1140991) q[3];
sx q[3];
rz(-1.8614785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1396609) q[0];
sx q[0];
rz(-0.49134555) q[0];
sx q[0];
rz(-0.77144462) q[0];
rz(1.3123243) q[1];
sx q[1];
rz(-1.7701365) q[1];
sx q[1];
rz(-1.2244474) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77030047) q[0];
sx q[0];
rz(-1.5988358) q[0];
sx q[0];
rz(1.4337149) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5688875) q[2];
sx q[2];
rz(-0.95237007) q[2];
sx q[2];
rz(-1.2313953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9414373) q[1];
sx q[1];
rz(-1.1863882) q[1];
sx q[1];
rz(-3.0430493) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78449747) q[3];
sx q[3];
rz(-1.906708) q[3];
sx q[3];
rz(1.0264736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9839342) q[2];
sx q[2];
rz(-2.5658786) q[2];
sx q[2];
rz(2.0470587) q[2];
rz(0.076722773) q[3];
sx q[3];
rz(-0.50944296) q[3];
sx q[3];
rz(-2.4491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43129608) q[0];
sx q[0];
rz(-2.4831979) q[0];
sx q[0];
rz(-0.21507138) q[0];
rz(0.98945016) q[1];
sx q[1];
rz(-0.96885252) q[1];
sx q[1];
rz(-1.4873803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7152255) q[0];
sx q[0];
rz(-2.0837113) q[0];
sx q[0];
rz(-1.7576393) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39535661) q[2];
sx q[2];
rz(-1.2985419) q[2];
sx q[2];
rz(-2.572042) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3645059) q[1];
sx q[1];
rz(-1.1544656) q[1];
sx q[1];
rz(-2.7440666) q[1];
x q[2];
rz(1.0808252) q[3];
sx q[3];
rz(-2.1741011) q[3];
sx q[3];
rz(1.5321466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4789077) q[2];
sx q[2];
rz(-1.4143133) q[2];
sx q[2];
rz(-1.5642081) q[2];
rz(-0.94791895) q[3];
sx q[3];
rz(-1.0736059) q[3];
sx q[3];
rz(1.4028367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941512) q[0];
sx q[0];
rz(-1.489137) q[0];
sx q[0];
rz(0.58962756) q[0];
rz(1.6472752) q[1];
sx q[1];
rz(-1.3811771) q[1];
sx q[1];
rz(-3.0628915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8183647) q[0];
sx q[0];
rz(-1.282489) q[0];
sx q[0];
rz(-2.7298729) q[0];
x q[1];
rz(-0.45273089) q[2];
sx q[2];
rz(-2.6331278) q[2];
sx q[2];
rz(-1.1217211) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0682875) q[1];
sx q[1];
rz(-1.3228497) q[1];
sx q[1];
rz(-1.7854879) q[1];
x q[2];
rz(-2.614903) q[3];
sx q[3];
rz(-1.1136276) q[3];
sx q[3];
rz(1.7383733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2176723) q[2];
sx q[2];
rz(-2.2675026) q[2];
sx q[2];
rz(-0.11202845) q[2];
rz(-0.40515408) q[3];
sx q[3];
rz(-1.3978037) q[3];
sx q[3];
rz(0.20598327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21796313) q[0];
sx q[0];
rz(-2.1681652) q[0];
sx q[0];
rz(1.7226891) q[0];
rz(-1.0796684) q[1];
sx q[1];
rz(-2.7579312) q[1];
sx q[1];
rz(1.7344281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16033835) q[0];
sx q[0];
rz(-1.8869074) q[0];
sx q[0];
rz(-2.714732) q[0];
x q[1];
rz(2.9687644) q[2];
sx q[2];
rz(-2.7044074) q[2];
sx q[2];
rz(-2.3035216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8720334) q[1];
sx q[1];
rz(-0.57231325) q[1];
sx q[1];
rz(-0.78929269) q[1];
rz(-pi) q[2];
rz(2.892835) q[3];
sx q[3];
rz(-1.837656) q[3];
sx q[3];
rz(-2.998179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63274038) q[2];
sx q[2];
rz(-0.036157046) q[2];
sx q[2];
rz(2.4031694) q[2];
rz(-0.19668713) q[3];
sx q[3];
rz(-0.6179215) q[3];
sx q[3];
rz(0.64016199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7554373) q[0];
sx q[0];
rz(-2.8897987) q[0];
sx q[0];
rz(-0.039903076) q[0];
rz(-2.9207322) q[1];
sx q[1];
rz(-1.5233636) q[1];
sx q[1];
rz(-1.1562645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62558371) q[0];
sx q[0];
rz(-1.7157641) q[0];
sx q[0];
rz(-1.8377393) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73157633) q[2];
sx q[2];
rz(-2.8551364) q[2];
sx q[2];
rz(-0.67026223) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77462536) q[1];
sx q[1];
rz(-1.0242265) q[1];
sx q[1];
rz(0.39440407) q[1];
rz(3.0262104) q[3];
sx q[3];
rz(-1.4658864) q[3];
sx q[3];
rz(-0.87593381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7463189) q[2];
sx q[2];
rz(-2.7776182) q[2];
sx q[2];
rz(3.086997) q[2];
rz(0.78628457) q[3];
sx q[3];
rz(-1.1966642) q[3];
sx q[3];
rz(-1.981885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7603124) q[0];
sx q[0];
rz(-2.3535643) q[0];
sx q[0];
rz(-2.330761) q[0];
rz(2.6122818) q[1];
sx q[1];
rz(-1.7050754) q[1];
sx q[1];
rz(-2.0307821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7082695) q[0];
sx q[0];
rz(-0.71107159) q[0];
sx q[0];
rz(0.61490402) q[0];
rz(-pi) q[1];
rz(2.6535437) q[2];
sx q[2];
rz(-2.077861) q[2];
sx q[2];
rz(1.4442854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0089493) q[1];
sx q[1];
rz(-2.3248782) q[1];
sx q[1];
rz(-1.3470178) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96145351) q[3];
sx q[3];
rz(-1.8224484) q[3];
sx q[3];
rz(0.52667743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3436053) q[2];
sx q[2];
rz(-0.56861773) q[2];
sx q[2];
rz(-2.9000751) q[2];
rz(0.24164116) q[3];
sx q[3];
rz(-1.6499949) q[3];
sx q[3];
rz(-2.9901166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80736137) q[0];
sx q[0];
rz(-1.3483028) q[0];
sx q[0];
rz(-2.469113) q[0];
rz(-0.45202759) q[1];
sx q[1];
rz(-0.94257911) q[1];
sx q[1];
rz(-0.51530757) q[1];
rz(-0.059941344) q[2];
sx q[2];
rz(-1.1287545) q[2];
sx q[2];
rz(-1.9311369) q[2];
rz(2.3005121) q[3];
sx q[3];
rz(-0.65541808) q[3];
sx q[3];
rz(-2.2578166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
