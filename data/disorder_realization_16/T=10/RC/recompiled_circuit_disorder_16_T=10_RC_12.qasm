OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(4.1339388) q[1];
sx q[1];
rz(9.0864656) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6217123) q[0];
sx q[0];
rz(-1.8147239) q[0];
sx q[0];
rz(1.8983311) q[0];
rz(-pi) q[1];
rz(2.7726735) q[2];
sx q[2];
rz(-2.2152165) q[2];
sx q[2];
rz(-1.8298139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0677883) q[1];
sx q[1];
rz(-1.5177625) q[1];
sx q[1];
rz(0.28633134) q[1];
rz(-1.0584352) q[3];
sx q[3];
rz(-1.5570939) q[3];
sx q[3];
rz(-2.0126359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8006111) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(-0.064963438) q[0];
rz(0.57463542) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(-1.2423135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500538) q[0];
sx q[0];
rz(-2.8613052) q[0];
sx q[0];
rz(2.6602402) q[0];
rz(-pi) q[1];
rz(-0.1588891) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(-1.6275258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2989267) q[1];
sx q[1];
rz(-2.0999523) q[1];
sx q[1];
rz(-1.3859205) q[1];
x q[2];
rz(-1.0817238) q[3];
sx q[3];
rz(-0.92182577) q[3];
sx q[3];
rz(1.058941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(-2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(3.0103502) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-2.4131391) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(2.1247991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166312) q[0];
sx q[0];
rz(-1.6065856) q[0];
sx q[0];
rz(-0.081598452) q[0];
rz(-2.0605893) q[2];
sx q[2];
rz(-2.042633) q[2];
sx q[2];
rz(-0.29957289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46360717) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(0.82171085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8524283) q[3];
sx q[3];
rz(-1.5236119) q[3];
sx q[3];
rz(0.55164528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.1857741) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8124354) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(0.24681117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033337489) q[0];
sx q[0];
rz(-0.70972432) q[0];
sx q[0];
rz(1.3413341) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7814126) q[2];
sx q[2];
rz(-1.1014551) q[2];
sx q[2];
rz(-0.59553713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4775131) q[1];
sx q[1];
rz(-1.941136) q[1];
sx q[1];
rz(-1.5143637) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28624268) q[3];
sx q[3];
rz(-0.97385588) q[3];
sx q[3];
rz(0.7569353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(-2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33070579) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.3274308) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7800956) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(1.6105152) q[0];
rz(0.70154538) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(1.0813576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(2.8413248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8304843) q[3];
sx q[3];
rz(-0.6725544) q[3];
sx q[3];
rz(2.6274519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(-0.38875368) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(3.016901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9799177) q[0];
sx q[0];
rz(-1.7470164) q[0];
sx q[0];
rz(1.6582703) q[0];
x q[1];
rz(1.2776676) q[2];
sx q[2];
rz(-2.1595402) q[2];
sx q[2];
rz(2.6063906) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9166959) q[1];
sx q[1];
rz(-1.6273013) q[1];
sx q[1];
rz(0.7373666) q[1];
rz(-1.3521306) q[3];
sx q[3];
rz(-1.339774) q[3];
sx q[3];
rz(-0.71393379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-2.2119904) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(1.0429617) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-2.0708864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844855) q[0];
sx q[0];
rz(-2.6575408) q[0];
sx q[0];
rz(0.0012782106) q[0];
rz(1.7037017) q[2];
sx q[2];
rz(-1.9153567) q[2];
sx q[2];
rz(2.8466356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0740944) q[1];
sx q[1];
rz(-0.80103445) q[1];
sx q[1];
rz(1.4317516) q[1];
rz(-pi) q[2];
rz(1.5504863) q[3];
sx q[3];
rz(-1.9029402) q[3];
sx q[3];
rz(2.9796245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(2.8179742) q[2];
rz(0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(-1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(-1.2063684) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(2.1941197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652126) q[0];
sx q[0];
rz(-1.0351666) q[0];
sx q[0];
rz(0.11760786) q[0];
rz(-pi) q[1];
rz(-0.1768441) q[2];
sx q[2];
rz(-1.4143922) q[2];
sx q[2];
rz(0.96207372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35093388) q[1];
sx q[1];
rz(-0.62218636) q[1];
sx q[1];
rz(-1.15637) q[1];
rz(1.2097589) q[3];
sx q[3];
rz(-0.79663888) q[3];
sx q[3];
rz(2.3494997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8895421) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(1.8126194) q[0];
rz(0.7912311) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(-2.9387617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79265362) q[0];
sx q[0];
rz(-1.5447504) q[0];
sx q[0];
rz(-3.1025725) q[0];
rz(-0.22325309) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(0.60165652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.030414) q[1];
sx q[1];
rz(-1.6200388) q[1];
sx q[1];
rz(1.4700252) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8392302) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(-2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41708502) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(0.28082401) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1998394) q[0];
sx q[0];
rz(-2.0240677) q[0];
sx q[0];
rz(-2.0487294) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7634723) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(3.0378621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22051375) q[1];
sx q[1];
rz(-1.6069357) q[1];
sx q[1];
rz(-1.5131348) q[1];
rz(-pi) q[2];
rz(1.9032352) q[3];
sx q[3];
rz(-1.3849764) q[3];
sx q[3];
rz(0.18248724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8250371) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(-2.1394219) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-2.5429824) q[2];
sx q[2];
rz(-2.241588) q[2];
sx q[2];
rz(-0.66551756) q[2];
rz(-0.35691805) q[3];
sx q[3];
rz(-1.1834984) q[3];
sx q[3];
rz(-0.055565861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
