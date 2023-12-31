OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1146381) q[0];
sx q[0];
rz(-1.4517598) q[0];
sx q[0];
rz(2.4960158) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6490373) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(-1.210968) q[0];
rz(-1.6662007) q[2];
sx q[2];
rz(-2.5089052) q[2];
sx q[2];
rz(-0.22104095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9420535) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(2.8295423) q[1];
rz(-pi) q[2];
rz(0.50819355) q[3];
sx q[3];
rz(-2.4262706) q[3];
sx q[3];
rz(-3.1014266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(-2.3085964) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.3566646) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67265636) q[0];
sx q[0];
rz(-2.7514973) q[0];
sx q[0];
rz(0.54498589) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(2.3166336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4175407) q[0];
sx q[0];
rz(-1.9736971) q[0];
sx q[0];
rz(2.9379662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5812133) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(1.8889697) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7488885) q[1];
sx q[1];
rz(-2.5671133) q[1];
sx q[1];
rz(2.9267465) q[1];
rz(-pi) q[2];
rz(-2.9460658) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(-0.60002458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(0.17641243) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-0.038539561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3241203) q[0];
sx q[0];
rz(-1.5701553) q[0];
sx q[0];
rz(1.5630432) q[0];
rz(-1.4171717) q[2];
sx q[2];
rz(-1.3426174) q[2];
sx q[2];
rz(1.4128026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4290532) q[1];
sx q[1];
rz(-2.6479811) q[1];
sx q[1];
rz(1.7008971) q[1];
rz(0.58979804) q[3];
sx q[3];
rz(-2.0876948) q[3];
sx q[3];
rz(-0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5588351) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-0.19515881) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.450481) q[0];
sx q[0];
rz(0.093219482) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0267369) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(1.3865711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8198204) q[1];
sx q[1];
rz(-2.0211377) q[1];
sx q[1];
rz(-0.26178534) q[1];
rz(2.5024662) q[3];
sx q[3];
rz(-1.5838606) q[3];
sx q[3];
rz(0.55707219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-0.8403362) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(-0.37809125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3867214) q[0];
sx q[0];
rz(-0.62749642) q[0];
sx q[0];
rz(-2.0621215) q[0];
rz(-pi) q[1];
rz(0.80412229) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(-0.59876195) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0542478) q[1];
sx q[1];
rz(-1.1451045) q[1];
sx q[1];
rz(-0.40360968) q[1];
rz(-pi) q[2];
rz(1.7686754) q[3];
sx q[3];
rz(-1.1117001) q[3];
sx q[3];
rz(1.0517694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75509214) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(-2.2264218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.299303) q[0];
sx q[0];
rz(-1.8931086) q[0];
sx q[0];
rz(-1.6567985) q[0];
rz(1.0038063) q[2];
sx q[2];
rz(-2.0096471) q[2];
sx q[2];
rz(-1.7457419) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5988016) q[1];
sx q[1];
rz(-2.8208477) q[1];
sx q[1];
rz(-0.21943211) q[1];
rz(-pi) q[2];
rz(-1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(-2.3510695) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73615605) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(-1.7162011) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910949) q[0];
sx q[0];
rz(-1.1384581) q[0];
sx q[0];
rz(1.3652486) q[0];
rz(-0.44600365) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(2.8446978) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.62896171) q[1];
sx q[1];
rz(-2.1753864) q[1];
sx q[1];
rz(-2.0357893) q[1];
rz(1.0670788) q[3];
sx q[3];
rz(-1.5590612) q[3];
sx q[3];
rz(2.3624682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(3.0650744) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9936179) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(-2.5698635) q[0];
x q[1];
rz(-0.92213995) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(-0.61308544) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68554316) q[1];
sx q[1];
rz(-1.4577663) q[1];
sx q[1];
rz(-3.0526524) q[1];
rz(2.5520677) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(2.7713283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.8267652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1413404) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(1.89639) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36395145) q[2];
sx q[2];
rz(-1.4462573) q[2];
sx q[2];
rz(3.0968015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6445551) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(2.6130996) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5626838) q[3];
sx q[3];
rz(-2.3981961) q[3];
sx q[3];
rz(-0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(2.8961199) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(0.28668177) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(2.749696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4927917) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(3.1317943) q[0];
rz(-pi) q[1];
rz(-0.91498615) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(-2.8142625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3956986) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.3032773) q[1];
rz(-pi) q[2];
rz(-0.55962555) q[3];
sx q[3];
rz(-1.8436173) q[3];
sx q[3];
rz(0.19459693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3863603) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(2.9677283) q[2];
sx q[2];
rz(-2.3542913) q[2];
sx q[2];
rz(-2.3402294) q[2];
rz(2.1203534) q[3];
sx q[3];
rz(-1.869535) q[3];
sx q[3];
rz(0.46365999) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
