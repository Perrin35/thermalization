OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(1.1507432) q[0];
rz(-0.28252217) q[1];
sx q[1];
rz(-2.0158267) q[1];
sx q[1];
rz(1.1123302) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939358) q[0];
sx q[0];
rz(-0.2875244) q[0];
sx q[0];
rz(-2.7178746) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98872008) q[2];
sx q[2];
rz(-1.4766105) q[2];
sx q[2];
rz(-2.1605952) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.46953) q[1];
sx q[1];
rz(-2.1066252) q[1];
sx q[1];
rz(2.5754926) q[1];
rz(-0.86701534) q[3];
sx q[3];
rz(-0.91546446) q[3];
sx q[3];
rz(-2.1239779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7925966) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(1.3953155) q[2];
rz(-3.0830749) q[3];
sx q[3];
rz(-1.7250215) q[3];
sx q[3];
rz(1.831656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.42404744) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(-0.44977093) q[0];
rz(-0.46478081) q[1];
sx q[1];
rz(-1.3100781) q[1];
sx q[1];
rz(-0.25258499) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.775894) q[0];
sx q[0];
rz(-1.5868385) q[0];
sx q[0];
rz(1.9220065) q[0];
rz(2.4812883) q[2];
sx q[2];
rz(-1.2129158) q[2];
sx q[2];
rz(1.7972657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8014268) q[1];
sx q[1];
rz(-2.8599842) q[1];
sx q[1];
rz(-1.5810709) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59201805) q[3];
sx q[3];
rz(-2.238174) q[3];
sx q[3];
rz(-0.84718466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8978591) q[2];
sx q[2];
rz(-2.2458138) q[2];
sx q[2];
rz(0.016156999) q[2];
rz(-2.2195418) q[3];
sx q[3];
rz(-0.99370304) q[3];
sx q[3];
rz(-3.0097358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.09318) q[0];
sx q[0];
rz(-1.484363) q[0];
sx q[0];
rz(-0.68159252) q[0];
rz(-2.5438578) q[1];
sx q[1];
rz(-1.6860516) q[1];
sx q[1];
rz(-0.32424232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77406835) q[0];
sx q[0];
rz(-1.2307967) q[0];
sx q[0];
rz(-2.6794479) q[0];
rz(-pi) q[1];
rz(1.6382255) q[2];
sx q[2];
rz(-2.0766269) q[2];
sx q[2];
rz(2.1492599) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2940087) q[1];
sx q[1];
rz(-1.4022744) q[1];
sx q[1];
rz(0.042673341) q[1];
x q[2];
rz(0.16894582) q[3];
sx q[3];
rz(-0.96833951) q[3];
sx q[3];
rz(1.4698616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9165667) q[2];
sx q[2];
rz(-0.39174199) q[2];
sx q[2];
rz(-2.0908835) q[2];
rz(2.4010036) q[3];
sx q[3];
rz(-1.2935484) q[3];
sx q[3];
rz(-3.0802687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447613) q[0];
sx q[0];
rz(-2.3065541) q[0];
sx q[0];
rz(0.9541676) q[0];
rz(-1.1369811) q[1];
sx q[1];
rz(-1.6258207) q[1];
sx q[1];
rz(0.7581996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0643142) q[0];
sx q[0];
rz(-1.2998733) q[0];
sx q[0];
rz(3.0999557) q[0];
rz(-2.3978176) q[2];
sx q[2];
rz(-1.4267046) q[2];
sx q[2];
rz(-0.38915044) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2705225) q[1];
sx q[1];
rz(-0.70065537) q[1];
sx q[1];
rz(-1.7228222) q[1];
rz(-0.011298374) q[3];
sx q[3];
rz(-0.73888981) q[3];
sx q[3];
rz(2.9113462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.053146426) q[2];
sx q[2];
rz(-2.5909178) q[2];
sx q[2];
rz(-1.3147563) q[2];
rz(-2.8407319) q[3];
sx q[3];
rz(-1.340539) q[3];
sx q[3];
rz(-0.68527591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.143173) q[0];
sx q[0];
rz(-1.5564065) q[0];
sx q[0];
rz(1.5013303) q[0];
rz(-2.7550664) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(-1.5949257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0536755) q[0];
sx q[0];
rz(-2.2639416) q[0];
sx q[0];
rz(-1.7960834) q[0];
rz(-0.96225454) q[2];
sx q[2];
rz(-1.3013162) q[2];
sx q[2];
rz(0.4294006) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.23473492) q[1];
sx q[1];
rz(-1.4868007) q[1];
sx q[1];
rz(0.43387167) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5677979) q[3];
sx q[3];
rz(-1.6125896) q[3];
sx q[3];
rz(2.6190663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3950562) q[2];
sx q[2];
rz(-2.655513) q[2];
sx q[2];
rz(-1.0699832) q[2];
rz(0.88548958) q[3];
sx q[3];
rz(-1.0642137) q[3];
sx q[3];
rz(1.1522256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162061) q[0];
sx q[0];
rz(-2.8031271) q[0];
sx q[0];
rz(-1.0534519) q[0];
rz(0.18928754) q[1];
sx q[1];
rz(-2.3129251) q[1];
sx q[1];
rz(-1.6493571) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5439592) q[0];
sx q[0];
rz(-1.9372371) q[0];
sx q[0];
rz(-2.9935915) q[0];
rz(1.8771421) q[2];
sx q[2];
rz(-1.3679081) q[2];
sx q[2];
rz(-2.5032147) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1188394) q[1];
sx q[1];
rz(-0.51278881) q[1];
sx q[1];
rz(-0.20365479) q[1];
rz(-pi) q[2];
rz(0.80811852) q[3];
sx q[3];
rz(-1.5991773) q[3];
sx q[3];
rz(-0.24695268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80387992) q[2];
sx q[2];
rz(-1.0151851) q[2];
sx q[2];
rz(1.052915) q[2];
rz(3.0247011) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(1.4611999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148024) q[0];
sx q[0];
rz(-2.5616665) q[0];
sx q[0];
rz(-0.54452288) q[0];
rz(0.25360423) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(2.8869218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47557467) q[0];
sx q[0];
rz(-3.0694002) q[0];
sx q[0];
rz(-1.484314) q[0];
rz(-2.7484863) q[2];
sx q[2];
rz(-1.6507033) q[2];
sx q[2];
rz(-2.7167632) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0510301) q[1];
sx q[1];
rz(-1.5080308) q[1];
sx q[1];
rz(-0.81845567) q[1];
rz(-pi) q[2];
rz(2.3319844) q[3];
sx q[3];
rz(-2.0008464) q[3];
sx q[3];
rz(2.3187217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2250681) q[2];
sx q[2];
rz(-2.2038286) q[2];
sx q[2];
rz(2.7226105) q[2];
rz(-0.032912832) q[3];
sx q[3];
rz(-1.2337647) q[3];
sx q[3];
rz(2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58881775) q[0];
sx q[0];
rz(-2.3989615) q[0];
sx q[0];
rz(1.8137929) q[0];
rz(-2.1090419) q[1];
sx q[1];
rz(-1.7951218) q[1];
sx q[1];
rz(2.5520777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1731271) q[0];
sx q[0];
rz(-2.0636981) q[0];
sx q[0];
rz(2.1987134) q[0];
rz(0.2915129) q[2];
sx q[2];
rz(-1.8290724) q[2];
sx q[2];
rz(-1.8669389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3302916) q[1];
sx q[1];
rz(-1.9988235) q[1];
sx q[1];
rz(0.93702646) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0699851) q[3];
sx q[3];
rz(-1.2195865) q[3];
sx q[3];
rz(0.90027383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0699658) q[2];
sx q[2];
rz(-2.4667141) q[2];
sx q[2];
rz(2.2281632) q[2];
rz(2.6895798) q[3];
sx q[3];
rz(-1.2830696) q[3];
sx q[3];
rz(-1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2695059) q[0];
sx q[0];
rz(-2.1156613) q[0];
sx q[0];
rz(-1.5103229) q[0];
rz(-1.3144846) q[1];
sx q[1];
rz(-2.1035106) q[1];
sx q[1];
rz(0.41183919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3519234) q[0];
sx q[0];
rz(-2.4487307) q[0];
sx q[0];
rz(-2.0122276) q[0];
rz(-pi) q[1];
rz(-0.14526315) q[2];
sx q[2];
rz(-0.32215873) q[2];
sx q[2];
rz(-1.8069428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5631521) q[1];
sx q[1];
rz(-2.4870589) q[1];
sx q[1];
rz(0.5963948) q[1];
x q[2];
rz(1.8985073) q[3];
sx q[3];
rz(-1.2043038) q[3];
sx q[3];
rz(1.2372563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1371548) q[2];
sx q[2];
rz(-1.4919446) q[2];
sx q[2];
rz(0.23492661) q[2];
rz(0.91819417) q[3];
sx q[3];
rz(-2.4166959) q[3];
sx q[3];
rz(-2.0535859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10201564) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(-1.8812195) q[0];
rz(-1.6873987) q[1];
sx q[1];
rz(-1.6834384) q[1];
sx q[1];
rz(1.6967324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6337109) q[0];
sx q[0];
rz(-2.0145881) q[0];
sx q[0];
rz(-1.1133975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6131975) q[2];
sx q[2];
rz(-0.95167347) q[2];
sx q[2];
rz(2.2592659) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6548938) q[1];
sx q[1];
rz(-1.6954665) q[1];
sx q[1];
rz(2.423684) q[1];
x q[2];
rz(-2.7637611) q[3];
sx q[3];
rz(-1.8820541) q[3];
sx q[3];
rz(1.6938083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5626278) q[2];
sx q[2];
rz(-0.29890385) q[2];
sx q[2];
rz(-0.12358269) q[2];
rz(-0.28892162) q[3];
sx q[3];
rz(-1.2762504) q[3];
sx q[3];
rz(0.45946521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73744437) q[0];
sx q[0];
rz(-1.7973719) q[0];
sx q[0];
rz(1.6112882) q[0];
rz(2.1743446) q[1];
sx q[1];
rz(-0.9728685) q[1];
sx q[1];
rz(-2.2546993) q[1];
rz(-2.5423302) q[2];
sx q[2];
rz(-1.6159372) q[2];
sx q[2];
rz(1.898832) q[2];
rz(-2.5206) q[3];
sx q[3];
rz(-0.98747333) q[3];
sx q[3];
rz(2.3458521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
