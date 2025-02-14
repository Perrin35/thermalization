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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(1.6562847) q[0];
rz(-2.0909042) q[1];
sx q[1];
rz(-1.3991855) q[1];
sx q[1];
rz(-2.4032226) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2353819) q[0];
sx q[0];
rz(-0.28774187) q[0];
sx q[0];
rz(0.55858992) q[0];
x q[1];
rz(1.9264439) q[2];
sx q[2];
rz(-2.7251149) q[2];
sx q[2];
rz(2.6082325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.76285) q[1];
sx q[1];
rz(-1.1302916) q[1];
sx q[1];
rz(-1.7251013) q[1];
rz(-2.5220736) q[3];
sx q[3];
rz(-1.8219007) q[3];
sx q[3];
rz(0.34003231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-0.28306511) q[2];
sx q[2];
rz(-0.4134678) q[2];
rz(-0.56420285) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(1.9570501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7269932) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(0.10398908) q[0];
rz(-0.14101401) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-2.7242421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0275729) q[0];
sx q[0];
rz(-2.0947523) q[0];
sx q[0];
rz(-2.7922676) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.012243791) q[2];
sx q[2];
rz(-1.9345539) q[2];
sx q[2];
rz(-3.0573483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1097393) q[1];
sx q[1];
rz(-1.8735587) q[1];
sx q[1];
rz(3.0979373) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6127698) q[3];
sx q[3];
rz(-1.0695056) q[3];
sx q[3];
rz(-2.5426585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.523681) q[2];
sx q[2];
rz(-1.0289501) q[2];
sx q[2];
rz(2.9376612) q[2];
rz(1.5150874) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6025036) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(-2.0643945) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(-2.5881252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578149) q[0];
sx q[0];
rz(-2.4071065) q[0];
sx q[0];
rz(-0.015588394) q[0];
rz(-pi) q[1];
rz(0.71142324) q[2];
sx q[2];
rz(-1.7345979) q[2];
sx q[2];
rz(-1.244551) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1225153) q[1];
sx q[1];
rz(-1.919849) q[1];
sx q[1];
rz(2.529789) q[1];
rz(2.9514246) q[3];
sx q[3];
rz(-1.3983418) q[3];
sx q[3];
rz(1.8784539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.074097721) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(-1.925776) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(1.2522987) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85553402) q[0];
sx q[0];
rz(-1.0581886) q[0];
sx q[0];
rz(-1.4372987) q[0];
rz(-1.8057436) q[1];
sx q[1];
rz(-0.62594405) q[1];
sx q[1];
rz(-0.45509532) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6949161) q[0];
sx q[0];
rz(-1.6507663) q[0];
sx q[0];
rz(-0.6178426) q[0];
x q[1];
rz(2.0680769) q[2];
sx q[2];
rz(-2.7233363) q[2];
sx q[2];
rz(-2.4474395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0245106) q[1];
sx q[1];
rz(-0.41943892) q[1];
sx q[1];
rz(-2.504435) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6654799) q[3];
sx q[3];
rz(-1.3631184) q[3];
sx q[3];
rz(-2.7439347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(-2.3095798) q[2];
rz(0.653382) q[3];
sx q[3];
rz(-0.91975206) q[3];
sx q[3];
rz(2.466989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5117383) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(0.6889371) q[0];
rz(0.2116994) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(-0.45417085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31541889) q[0];
sx q[0];
rz(-2.4103904) q[0];
sx q[0];
rz(0.8765799) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7709097) q[2];
sx q[2];
rz(-2.1739568) q[2];
sx q[2];
rz(-0.93450817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3836174) q[1];
sx q[1];
rz(-1.0838008) q[1];
sx q[1];
rz(-2.3112554) q[1];
x q[2];
rz(2.1319785) q[3];
sx q[3];
rz(-0.30721634) q[3];
sx q[3];
rz(0.12061435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5357369) q[2];
sx q[2];
rz(-1.7822632) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(1.2153252) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(-0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257492) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(0.51344839) q[0];
rz(0.91649857) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(1.7880012) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5149263) q[0];
sx q[0];
rz(-2.0139222) q[0];
sx q[0];
rz(1.6624381) q[0];
x q[1];
rz(-2.2960179) q[2];
sx q[2];
rz(-0.52547272) q[2];
sx q[2];
rz(2.2487244) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.533355) q[1];
sx q[1];
rz(-1.5974733) q[1];
sx q[1];
rz(-0.4877301) q[1];
rz(-1.6597802) q[3];
sx q[3];
rz(-1.5060194) q[3];
sx q[3];
rz(-1.806206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6650271) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(-2.636886) q[3];
sx q[3];
rz(-1.4387771) q[3];
sx q[3];
rz(0.75468841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018709239) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(1.9710185) q[0];
rz(2.9508044) q[1];
sx q[1];
rz(-2.6759594) q[1];
sx q[1];
rz(0.64291397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.64775) q[0];
sx q[0];
rz(-2.8558585) q[0];
sx q[0];
rz(-1.2005376) q[0];
x q[1];
rz(-1.9011074) q[2];
sx q[2];
rz(-1.8165556) q[2];
sx q[2];
rz(1.2643697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9289405) q[1];
sx q[1];
rz(-1.3319322) q[1];
sx q[1];
rz(0.3122621) q[1];
rz(-pi) q[2];
rz(1.4516109) q[3];
sx q[3];
rz(-1.168712) q[3];
sx q[3];
rz(1.6329671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99870318) q[2];
sx q[2];
rz(-3.0316752) q[2];
sx q[2];
rz(-1.9996803) q[2];
rz(-0.029684639) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(2.3619385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340045) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(1.3080066) q[0];
rz(1.1169149) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(2.9294779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7362999) q[0];
sx q[0];
rz(-2.2515902) q[0];
sx q[0];
rz(0.41053562) q[0];
x q[1];
rz(-0.30435698) q[2];
sx q[2];
rz(-1.7946912) q[2];
sx q[2];
rz(-2.9606113) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5221716) q[1];
sx q[1];
rz(-0.7070578) q[1];
sx q[1];
rz(0.85925428) q[1];
x q[2];
rz(2.2905825) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(-1.855576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1954631) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(1.2809666) q[2];
rz(1.403275) q[3];
sx q[3];
rz(-2.2303228) q[3];
sx q[3];
rz(-2.4173071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69287777) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(2.2591059) q[0];
rz(-0.29300434) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(2.015347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6544271) q[0];
sx q[0];
rz(-1.2058655) q[0];
sx q[0];
rz(-2.4449744) q[0];
x q[1];
rz(1.5669021) q[2];
sx q[2];
rz(-0.97087395) q[2];
sx q[2];
rz(-1.1614161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80690224) q[1];
sx q[1];
rz(-0.77065361) q[1];
sx q[1];
rz(1.6477963) q[1];
x q[2];
rz(-1.8873439) q[3];
sx q[3];
rz(-1.1831565) q[3];
sx q[3];
rz(0.46771184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5249411) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(0.16723995) q[2];
rz(-3.1104769) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6328218) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(2.3311145) q[0];
rz(1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(-0.097361758) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3980078) q[0];
sx q[0];
rz(-0.79774374) q[0];
sx q[0];
rz(-1.8037075) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4989158) q[2];
sx q[2];
rz(-0.38205636) q[2];
sx q[2];
rz(0.17554131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88579455) q[1];
sx q[1];
rz(-2.5558639) q[1];
sx q[1];
rz(0.08759193) q[1];
rz(3.024289) q[3];
sx q[3];
rz(-1.8364834) q[3];
sx q[3];
rz(-0.22646389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8436766) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(-1.8664912) q[2];
rz(2.7013333) q[3];
sx q[3];
rz(-2.2831254) q[3];
sx q[3];
rz(2.8662203) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-2.5195925) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(0.95876454) q[2];
sx q[2];
rz(-3.0059881) q[2];
sx q[2];
rz(-0.81632951) q[2];
rz(2.4980656) q[3];
sx q[3];
rz(-1.7089927) q[3];
sx q[3];
rz(-0.24571936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
