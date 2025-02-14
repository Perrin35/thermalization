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
rz(-1.929317) q[0];
sx q[0];
rz(-1.1317929) q[0];
sx q[0];
rz(1.7687067) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(-1.7665266) q[1];
sx q[1];
rz(2.3553203) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776913) q[0];
sx q[0];
rz(-2.7161975) q[0];
sx q[0];
rz(1.0546272) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.588618) q[2];
sx q[2];
rz(-0.99054256) q[2];
sx q[2];
rz(-1.0667104) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7934301) q[1];
sx q[1];
rz(-1.600666) q[1];
sx q[1];
rz(1.8601599) q[1];
x q[2];
rz(2.1350098) q[3];
sx q[3];
rz(-0.67970961) q[3];
sx q[3];
rz(2.2077479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0218411) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(1.8915668) q[2];
rz(0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982872) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(3.0702718) q[0];
rz(-1.9885063) q[1];
sx q[1];
rz(-2.3801897) q[1];
sx q[1];
rz(-2.6128795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3246177) q[0];
sx q[0];
rz(-1.8864838) q[0];
sx q[0];
rz(-1.0782918) q[0];
rz(-pi) q[1];
rz(-2.9527093) q[2];
sx q[2];
rz(-1.802465) q[2];
sx q[2];
rz(0.29555368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0866249) q[1];
sx q[1];
rz(-1.0180915) q[1];
sx q[1];
rz(2.5235559) q[1];
rz(-0.22457476) q[3];
sx q[3];
rz(-1.9764426) q[3];
sx q[3];
rz(-2.7499475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4802287) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(-0.028623494) q[2];
rz(-2.7875767) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(0.55907512) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(2.2502374) q[0];
rz(-2.640653) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(3.0827789) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0100493) q[0];
sx q[0];
rz(-1.5230522) q[0];
sx q[0];
rz(0.086011925) q[0];
rz(-pi) q[1];
rz(2.0979314) q[2];
sx q[2];
rz(-1.322515) q[2];
sx q[2];
rz(2.6206526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8938287) q[1];
sx q[1];
rz(-2.2732452) q[1];
sx q[1];
rz(0.6599627) q[1];
rz(-pi) q[2];
rz(0.18043442) q[3];
sx q[3];
rz(-1.5431607) q[3];
sx q[3];
rz(-2.7829426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(0.2365665) q[2];
rz(2.2208354) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(-1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22223602) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(-2.2663569) q[0];
rz(-1.5454166) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(-3.0962871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0939498) q[0];
sx q[0];
rz(-1.6730621) q[0];
sx q[0];
rz(-3.0759252) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2384134) q[2];
sx q[2];
rz(-1.1588237) q[2];
sx q[2];
rz(-2.2087086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.873535) q[1];
sx q[1];
rz(-0.90374871) q[1];
sx q[1];
rz(-0.11063834) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.281146) q[3];
sx q[3];
rz(-2.5389606) q[3];
sx q[3];
rz(-1.8528191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8594325) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(-1.1939987) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.9186391) q[3];
sx q[3];
rz(2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.10876656) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(-0.68242514) q[0];
rz(-1.8531063) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(1.8103745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82560655) q[0];
sx q[0];
rz(-2.2974112) q[0];
sx q[0];
rz(3.1196703) q[0];
x q[1];
rz(-0.26010988) q[2];
sx q[2];
rz(-1.9528198) q[2];
sx q[2];
rz(-1.1191238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4466616) q[1];
sx q[1];
rz(-0.38685054) q[1];
sx q[1];
rz(2.7617161) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87792895) q[3];
sx q[3];
rz(-2.4932541) q[3];
sx q[3];
rz(3.0832269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0266626) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(2.4746573) q[2];
rz(2.5866348) q[3];
sx q[3];
rz(-2.4542377) q[3];
sx q[3];
rz(-2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32894593) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(-0.70372787) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(2.9827548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9861529) q[0];
sx q[0];
rz(-1.4503398) q[0];
sx q[0];
rz(-1.7140165) q[0];
x q[1];
rz(2.4689697) q[2];
sx q[2];
rz(-2.9458698) q[2];
sx q[2];
rz(1.7253523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82032404) q[1];
sx q[1];
rz(-1.0918198) q[1];
sx q[1];
rz(1.7198573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63780906) q[3];
sx q[3];
rz(-1.7553698) q[3];
sx q[3];
rz(-2.0710433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3682897) q[2];
sx q[2];
rz(-2.1998019) q[2];
sx q[2];
rz(-0.85757315) q[2];
rz(-1.9722021) q[3];
sx q[3];
rz(-1.2811067) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9661949) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(0.71863693) q[0];
rz(0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(-2.4729572) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.74092) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(2.7223865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.093542) q[2];
sx q[2];
rz(-0.68533939) q[2];
sx q[2];
rz(0.46721855) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9978511) q[1];
sx q[1];
rz(-0.5023379) q[1];
sx q[1];
rz(2.3010002) q[1];
rz(-pi) q[2];
rz(-3.1304006) q[3];
sx q[3];
rz(-2.5103223) q[3];
sx q[3];
rz(-2.5281363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64688993) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(2.1051712) q[2];
rz(-2.5070665) q[3];
sx q[3];
rz(-0.98792881) q[3];
sx q[3];
rz(-0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(-2.4155937) q[0];
rz(-1.1622102) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(1.9451709) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1370173) q[0];
sx q[0];
rz(-1.5992303) q[0];
sx q[0];
rz(-1.7445376) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4216719) q[2];
sx q[2];
rz(-2.8310378) q[2];
sx q[2];
rz(-0.20937411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70564044) q[1];
sx q[1];
rz(-1.1060103) q[1];
sx q[1];
rz(-2.57354) q[1];
rz(2.7128588) q[3];
sx q[3];
rz(-0.7168684) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2339345) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(0.93351239) q[2];
rz(-0.61698169) q[3];
sx q[3];
rz(-1.0022997) q[3];
sx q[3];
rz(0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(1.9762978) q[0];
sx q[0];
rz(-0.10460654) q[0];
sx q[0];
rz(-3.0883375) q[0];
rz(-0.28757295) q[1];
sx q[1];
rz(-2.2122999) q[1];
sx q[1];
rz(-2.8823749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142562) q[0];
sx q[0];
rz(-1.5467318) q[0];
sx q[0];
rz(1.6238449) q[0];
x q[1];
rz(3.0087148) q[2];
sx q[2];
rz(-1.3406375) q[2];
sx q[2];
rz(-1.2990189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99154918) q[1];
sx q[1];
rz(-0.75356149) q[1];
sx q[1];
rz(1.6607453) q[1];
rz(-pi) q[2];
rz(0.15786981) q[3];
sx q[3];
rz(-2.3693858) q[3];
sx q[3];
rz(1.3881573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58664924) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-2.9602236) q[2];
rz(0.3950611) q[3];
sx q[3];
rz(-0.22160465) q[3];
sx q[3];
rz(-2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7419389) q[0];
sx q[0];
rz(-1.5927277) q[0];
sx q[0];
rz(-0.3821061) q[0];
rz(2.8451879) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(-1.872725) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60703245) q[0];
sx q[0];
rz(-1.9968613) q[0];
sx q[0];
rz(2.0924278) q[0];
rz(-0.13258719) q[2];
sx q[2];
rz(-2.2397723) q[2];
sx q[2];
rz(-2.7821409) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.813432) q[1];
sx q[1];
rz(-2.8012271) q[1];
sx q[1];
rz(-1.634738) q[1];
x q[2];
rz(2.8954072) q[3];
sx q[3];
rz(-0.46588141) q[3];
sx q[3];
rz(1.4108489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2962013) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(0.38983795) q[2];
rz(1.2975533) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(-0.84387422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1525477) q[0];
sx q[0];
rz(-1.7392409) q[0];
sx q[0];
rz(1.637511) q[0];
rz(1.3564431) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(-2.155921) q[2];
sx q[2];
rz(-2.4092842) q[2];
sx q[2];
rz(0.53436919) q[2];
rz(0.33792321) q[3];
sx q[3];
rz(-2.3281964) q[3];
sx q[3];
rz(-2.7184814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
