OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7140952) q[0];
sx q[0];
rz(-2.5795955) q[0];
sx q[0];
rz(-0.23101097) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(-0.45431554) q[1];
sx q[1];
rz(-1.8543724) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79345771) q[0];
sx q[0];
rz(-1.0890111) q[0];
sx q[0];
rz(-2.8346377) q[0];
rz(-pi) q[1];
rz(1.8051992) q[2];
sx q[2];
rz(-1.3318828) q[2];
sx q[2];
rz(-0.51942458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8652628) q[1];
sx q[1];
rz(-1.6693322) q[1];
sx q[1];
rz(0.16698412) q[1];
rz(3.0662905) q[3];
sx q[3];
rz(-1.7455532) q[3];
sx q[3];
rz(-2.99461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7370558) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(2.7810968) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-2.093061) q[3];
sx q[3];
rz(1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26038134) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(2.4847109) q[0];
rz(0.33292133) q[1];
sx q[1];
rz(-1.0924783) q[1];
sx q[1];
rz(-0.93516707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.864185) q[0];
sx q[0];
rz(-1.5404697) q[0];
sx q[0];
rz(3.0369989) q[0];
rz(-pi) q[1];
rz(-0.33719535) q[2];
sx q[2];
rz(-2.8697578) q[2];
sx q[2];
rz(1.4120917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3352901) q[1];
sx q[1];
rz(-0.15220255) q[1];
sx q[1];
rz(2.8782513) q[1];
rz(-pi) q[2];
rz(1.7788497) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(1.6121685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3071345) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(1.9251941) q[2];
rz(-0.52792102) q[3];
sx q[3];
rz(-2.1025751) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5288178) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(0.58498996) q[0];
rz(-3.0168369) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(-1.4488719) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0199658) q[0];
sx q[0];
rz(-0.0040071132) q[0];
sx q[0];
rz(-0.83929707) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9563006) q[2];
sx q[2];
rz(-0.073436471) q[2];
sx q[2];
rz(2.9543608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.587882) q[1];
sx q[1];
rz(-2.0016016) q[1];
sx q[1];
rz(-0.84970857) q[1];
x q[2];
rz(2.8040941) q[3];
sx q[3];
rz(-1.322116) q[3];
sx q[3];
rz(1.6832222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8831732) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(1.6955356) q[2];
rz(-2.3840733) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(-0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3676753) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(2.0539334) q[0];
rz(-0.37519535) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(-2.1813724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4488569) q[0];
sx q[0];
rz(-1.5310643) q[0];
sx q[0];
rz(-1.8360774) q[0];
rz(1.4713418) q[2];
sx q[2];
rz(-1.5410454) q[2];
sx q[2];
rz(-3.0029675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.580048) q[1];
sx q[1];
rz(-0.64148884) q[1];
sx q[1];
rz(-0.80342355) q[1];
x q[2];
rz(-2.8071515) q[3];
sx q[3];
rz(-2.8570647) q[3];
sx q[3];
rz(0.65438017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9356161) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(-1.4975632) q[2];
rz(-0.86132541) q[3];
sx q[3];
rz(-0.70278168) q[3];
sx q[3];
rz(-0.45708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1645666) q[0];
sx q[0];
rz(-2.470546) q[0];
sx q[0];
rz(-2.5373051) q[0];
rz(-1.9970278) q[1];
sx q[1];
rz(-1.4860169) q[1];
sx q[1];
rz(-0.42246517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6910962) q[0];
sx q[0];
rz(-1.7393149) q[0];
sx q[0];
rz(-3.0590579) q[0];
rz(2.1429135) q[2];
sx q[2];
rz(-2.0586176) q[2];
sx q[2];
rz(-1.3582548) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4854086) q[1];
sx q[1];
rz(-2.2586169) q[1];
sx q[1];
rz(-1.6641892) q[1];
rz(2.7606008) q[3];
sx q[3];
rz(-1.2514827) q[3];
sx q[3];
rz(1.8007985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3182688) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(-2.4853415) q[2];
rz(1.5308135) q[3];
sx q[3];
rz(-1.8717513) q[3];
sx q[3];
rz(2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4055279) q[0];
sx q[0];
rz(-2.1298213) q[0];
sx q[0];
rz(0.62498012) q[0];
rz(0.75195733) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(1.5350852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67095516) q[0];
sx q[0];
rz(-2.8482901) q[0];
sx q[0];
rz(2.9706012) q[0];
rz(-0.67120303) q[2];
sx q[2];
rz(-2.3103788) q[2];
sx q[2];
rz(2.9321456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62404406) q[1];
sx q[1];
rz(-0.66429783) q[1];
sx q[1];
rz(0.37599998) q[1];
rz(-pi) q[2];
rz(2.296359) q[3];
sx q[3];
rz(-0.74016011) q[3];
sx q[3];
rz(2.6773767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.1890761) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(0.9291741) q[2];
rz(0.82516986) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(-1.1024124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5850942) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(1.7671385) q[0];
rz(-2.6311686) q[1];
sx q[1];
rz(-2.3511032) q[1];
sx q[1];
rz(0.62320954) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0334367) q[0];
sx q[0];
rz(-1.9029473) q[0];
sx q[0];
rz(-0.16031127) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3750149) q[2];
sx q[2];
rz(-2.0718241) q[2];
sx q[2];
rz(-2.632338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12909266) q[1];
sx q[1];
rz(-2.9575037) q[1];
sx q[1];
rz(0.91839183) q[1];
rz(1.3680063) q[3];
sx q[3];
rz(-2.30808) q[3];
sx q[3];
rz(-1.1356419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0906543) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(-1.9742924) q[2];
rz(-0.022631571) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(2.0126655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8959494) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(2.1242712) q[0];
rz(1.7474878) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(2.5140433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16495569) q[0];
sx q[0];
rz(-0.73574388) q[0];
sx q[0];
rz(2.1829905) q[0];
rz(-pi) q[1];
rz(3.0369736) q[2];
sx q[2];
rz(-2.4596408) q[2];
sx q[2];
rz(1.6624818) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7146475) q[1];
sx q[1];
rz(-1.8083739) q[1];
sx q[1];
rz(-1.3690884) q[1];
x q[2];
rz(-1.9110319) q[3];
sx q[3];
rz(-1.1200532) q[3];
sx q[3];
rz(0.72117248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5158186) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-0.41075692) q[2];
rz(-3.0913894) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(0.075411782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.5589767) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(2.8726752) q[0];
rz(-0.31002632) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(3.0115829) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5101679) q[0];
sx q[0];
rz(-1.893259) q[0];
sx q[0];
rz(2.9084951) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79370462) q[2];
sx q[2];
rz(-0.48173258) q[2];
sx q[2];
rz(-2.6793753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7183154) q[1];
sx q[1];
rz(-2.1287103) q[1];
sx q[1];
rz(1.7608791) q[1];
rz(-pi) q[2];
rz(0.49874108) q[3];
sx q[3];
rz(-2.1864656) q[3];
sx q[3];
rz(1.9174674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93398634) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(-3.0450191) q[2];
rz(-1.6377595) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0132975) q[0];
sx q[0];
rz(-2.9533563) q[0];
sx q[0];
rz(-0.53303322) q[0];
rz(-0.036272613) q[1];
sx q[1];
rz(-0.78250042) q[1];
sx q[1];
rz(1.689555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9617758) q[0];
sx q[0];
rz(-2.4357033) q[0];
sx q[0];
rz(0.96700661) q[0];
rz(-pi) q[1];
rz(-2.8223561) q[2];
sx q[2];
rz(-1.6743273) q[2];
sx q[2];
rz(1.5543303) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2382792) q[1];
sx q[1];
rz(-0.36064816) q[1];
sx q[1];
rz(-3.0566932) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6598654) q[3];
sx q[3];
rz(-1.1067821) q[3];
sx q[3];
rz(0.17645141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4286246) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(-0.0326322) q[2];
rz(2.2964358) q[3];
sx q[3];
rz(-1.2266351) q[3];
sx q[3];
rz(2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7964771) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(-0.79700094) q[1];
sx q[1];
rz(-2.2140257) q[1];
sx q[1];
rz(0.066233403) q[1];
rz(1.8664411) q[2];
sx q[2];
rz(-2.0870123) q[2];
sx q[2];
rz(-0.43475702) q[2];
rz(-0.81099323) q[3];
sx q[3];
rz(-1.4832433) q[3];
sx q[3];
rz(1.087838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
