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
rz(-3.0566293) q[0];
sx q[0];
rz(-0.30240107) q[0];
sx q[0];
rz(-3.1095355) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(-2.6377331) q[1];
sx q[1];
rz(-2.388968) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0085474) q[0];
sx q[0];
rz(-0.38131443) q[0];
sx q[0];
rz(-2.5757936) q[0];
rz(-pi) q[1];
rz(-2.339509) q[2];
sx q[2];
rz(-0.49979106) q[2];
sx q[2];
rz(-2.0715203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17129414) q[1];
sx q[1];
rz(-2.7134905) q[1];
sx q[1];
rz(-2.5004205) q[1];
rz(-pi) q[2];
rz(2.0977375) q[3];
sx q[3];
rz(-0.40832061) q[3];
sx q[3];
rz(0.645831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(2.7005633) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012861982) q[0];
sx q[0];
rz(-0.91330376) q[0];
sx q[0];
rz(0.41020694) q[0];
rz(-0.38093105) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.358323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.298734) q[0];
sx q[0];
rz(-2.4141387) q[0];
sx q[0];
rz(2.1841315) q[0];
rz(-pi) q[1];
rz(1.0274506) q[2];
sx q[2];
rz(-1.1305692) q[2];
sx q[2];
rz(1.5981975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0002901) q[1];
sx q[1];
rz(-1.4030255) q[1];
sx q[1];
rz(-3.0649158) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7776599) q[3];
sx q[3];
rz(-0.65012041) q[3];
sx q[3];
rz(-2.2083685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7831948) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(-2.6386293) q[2];
rz(1.2991692) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(2.9097596) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(-0.9758392) q[0];
rz(-0.93636912) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(3.0036614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9693404) q[0];
sx q[0];
rz(-1.2069877) q[0];
sx q[0];
rz(2.0813638) q[0];
rz(-0.39769002) q[2];
sx q[2];
rz(-1.4639336) q[2];
sx q[2];
rz(-1.4390838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.960818) q[1];
sx q[1];
rz(-1.4394234) q[1];
sx q[1];
rz(0.8701433) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6343752) q[3];
sx q[3];
rz(-0.52195575) q[3];
sx q[3];
rz(-2.0824661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58933538) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(-1.662558) q[2];
rz(2.3432483) q[3];
sx q[3];
rz(-0.84404293) q[3];
sx q[3];
rz(-3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524356) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(-1.5929219) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(0.41935316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28946149) q[0];
sx q[0];
rz(-1.7999819) q[0];
sx q[0];
rz(2.9782692) q[0];
rz(-pi) q[1];
rz(-2.9871171) q[2];
sx q[2];
rz(-0.53336582) q[2];
sx q[2];
rz(-1.6090924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49653445) q[1];
sx q[1];
rz(-0.36318159) q[1];
sx q[1];
rz(1.7523604) q[1];
rz(-pi) q[2];
rz(0.97030117) q[3];
sx q[3];
rz(-0.56823778) q[3];
sx q[3];
rz(-2.0214391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7833917) q[2];
sx q[2];
rz(-2.6322067) q[2];
sx q[2];
rz(1.5436714) q[2];
rz(-1.915043) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(1.290087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.876038) q[0];
sx q[0];
rz(-0.72172481) q[0];
sx q[0];
rz(-2.3853886) q[0];
rz(1.2528231) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(-1.4160215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9276816) q[0];
sx q[0];
rz(-0.98386231) q[0];
sx q[0];
rz(-1.8326493) q[0];
x q[1];
rz(0.33966392) q[2];
sx q[2];
rz(-2.77861) q[2];
sx q[2];
rz(-1.5659005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49633138) q[1];
sx q[1];
rz(-2.1095536) q[1];
sx q[1];
rz(1.5625728) q[1];
rz(-pi) q[2];
rz(-2.7783423) q[3];
sx q[3];
rz(-1.3923858) q[3];
sx q[3];
rz(1.6476064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34117928) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(2.9608012) q[2];
rz(1.4888658) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(-1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7086696) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(2.3413626) q[0];
rz(0.284614) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-0.12399331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02888814) q[0];
sx q[0];
rz(-1.7005973) q[0];
sx q[0];
rz(3.1129254) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8154549) q[2];
sx q[2];
rz(-0.24194939) q[2];
sx q[2];
rz(-1.9242147) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2468894) q[1];
sx q[1];
rz(-1.8754301) q[1];
sx q[1];
rz(1.4851634) q[1];
rz(-pi) q[2];
x q[2];
rz(1.049253) q[3];
sx q[3];
rz(-0.66415706) q[3];
sx q[3];
rz(0.60598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2064712) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(3.0653811) q[2];
rz(1.291409) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(-2.5230303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16977075) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(-0.75702697) q[0];
rz(-0.79611671) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-0.98181358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9011089) q[0];
sx q[0];
rz(-2.3496858) q[0];
sx q[0];
rz(-1.9740941) q[0];
x q[1];
rz(-1.4635529) q[2];
sx q[2];
rz(-1.5011884) q[2];
sx q[2];
rz(2.7977242) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.031739) q[1];
sx q[1];
rz(-0.64424911) q[1];
sx q[1];
rz(0.8532614) q[1];
rz(-pi) q[2];
rz(1.8514093) q[3];
sx q[3];
rz(-1.4320489) q[3];
sx q[3];
rz(2.8900103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41165274) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(2.5977503) q[2];
rz(-3.1206711) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(-0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2127317) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(0.62335706) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-0.61505452) q[1];
sx q[1];
rz(-2.8772433) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2580483) q[0];
sx q[0];
rz(-2.2730581) q[0];
sx q[0];
rz(-0.27900286) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4521763) q[2];
sx q[2];
rz(-0.33782321) q[2];
sx q[2];
rz(1.092697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0496638) q[1];
sx q[1];
rz(-0.85159066) q[1];
sx q[1];
rz(-2.767241) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1846011) q[3];
sx q[3];
rz(-0.34117801) q[3];
sx q[3];
rz(-0.5285078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0048206) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(0.47478673) q[2];
rz(0.36561203) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30597618) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(0.48582745) q[0];
rz(1.0659418) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(-0.33445439) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4518648) q[0];
sx q[0];
rz(-0.34669995) q[0];
sx q[0];
rz(0.50677104) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4304305) q[2];
sx q[2];
rz(-2.0054501) q[2];
sx q[2];
rz(0.19186867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85035861) q[1];
sx q[1];
rz(-2.6667074) q[1];
sx q[1];
rz(0.33109457) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1074839) q[3];
sx q[3];
rz(-2.3114738) q[3];
sx q[3];
rz(-2.9024359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36755422) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(1.4878081) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(-3.0706792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.3286572) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(-2.8712414) q[0];
rz(-1.6290172) q[1];
sx q[1];
rz(-2.3051579) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.644076) q[0];
sx q[0];
rz(-0.77267161) q[0];
sx q[0];
rz(1.1091883) q[0];
rz(0.43105189) q[2];
sx q[2];
rz(-2.1630175) q[2];
sx q[2];
rz(-2.2811802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8423564) q[1];
sx q[1];
rz(-1.6116214) q[1];
sx q[1];
rz(-3.051332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7553431) q[3];
sx q[3];
rz(-1.874141) q[3];
sx q[3];
rz(-1.4331499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7865929) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(-1.5691441) q[2];
rz(-2.3908424) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(2.2641838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56562051) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(1.3427973) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-3.1411896) q[2];
sx q[2];
rz(-0.12231356) q[2];
sx q[2];
rz(-0.077153645) q[2];
rz(-0.20082898) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
