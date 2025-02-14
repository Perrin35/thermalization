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
rz(2.4438357) q[0];
sx q[0];
rz(4.3351686) q[0];
sx q[0];
rz(9.6293443) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(-1.3890356) q[1];
sx q[1];
rz(-1.7420446) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8176048) q[0];
sx q[0];
rz(-1.7919375) q[0];
sx q[0];
rz(3.0482685) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7367107) q[2];
sx q[2];
rz(-1.153115) q[2];
sx q[2];
rz(1.4428802) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0006932) q[1];
sx q[1];
rz(-1.6581891) q[1];
sx q[1];
rz(0.42229514) q[1];
rz(-pi) q[2];
rz(-2.4321626) q[3];
sx q[3];
rz(-1.3660079) q[3];
sx q[3];
rz(-0.11573175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7064887) q[2];
sx q[2];
rz(-3.1092643) q[2];
sx q[2];
rz(2.6522563) q[2];
rz(1.0232183) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(1.1255012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0959051) q[0];
sx q[0];
rz(-0.65653312) q[0];
sx q[0];
rz(0.80279654) q[0];
rz(0.071391694) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(-3.0838222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18081576) q[0];
sx q[0];
rz(-1.0313927) q[0];
sx q[0];
rz(-2.8404854) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34681706) q[2];
sx q[2];
rz(-2.8379734) q[2];
sx q[2];
rz(-0.70250073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7747468) q[1];
sx q[1];
rz(-0.54301942) q[1];
sx q[1];
rz(1.1096891) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9936419) q[3];
sx q[3];
rz(-0.86455621) q[3];
sx q[3];
rz(1.0195635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-2.0464996) q[2];
sx q[2];
rz(1.2954953) q[2];
rz(0.96674353) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(0.75743341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029723786) q[0];
sx q[0];
rz(-1.3628549) q[0];
sx q[0];
rz(-1.7080074) q[0];
rz(-3.0729821) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(-2.5624018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96087181) q[0];
sx q[0];
rz(-1.6287043) q[0];
sx q[0];
rz(-1.6625893) q[0];
rz(1.5726818) q[2];
sx q[2];
rz(-2.0465188) q[2];
sx q[2];
rz(1.204551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3855558) q[1];
sx q[1];
rz(-1.3446523) q[1];
sx q[1];
rz(1.6075587) q[1];
x q[2];
rz(1.0924358) q[3];
sx q[3];
rz(-1.1598827) q[3];
sx q[3];
rz(2.876407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27089831) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(0.18343055) q[2];
rz(0.80426788) q[3];
sx q[3];
rz(-2.1427514) q[3];
sx q[3];
rz(0.53406322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19876984) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(2.542069) q[0];
rz(0.41246688) q[1];
sx q[1];
rz(-0.020655276) q[1];
sx q[1];
rz(2.1309158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0507999) q[0];
sx q[0];
rz(-2.7677892) q[0];
sx q[0];
rz(-2.7888377) q[0];
rz(-pi) q[1];
rz(-0.93866556) q[2];
sx q[2];
rz(-1.5297339) q[2];
sx q[2];
rz(-1.2031885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0914826) q[1];
sx q[1];
rz(-2.1569139) q[1];
sx q[1];
rz(-2.8149603) q[1];
rz(-pi) q[2];
rz(-2.385731) q[3];
sx q[3];
rz(-1.6718277) q[3];
sx q[3];
rz(1.7441235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3707054) q[2];
sx q[2];
rz(-2.7949896) q[2];
sx q[2];
rz(0.31751219) q[2];
rz(-2.5192449) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7790826) q[0];
sx q[0];
rz(-0.91479397) q[0];
sx q[0];
rz(-2.1146178) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(1.9245573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.352223) q[0];
sx q[0];
rz(-0.93607157) q[0];
sx q[0];
rz(2.288649) q[0];
x q[1];
rz(0.53834277) q[2];
sx q[2];
rz(-1.6377047) q[2];
sx q[2];
rz(0.72170335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0266307) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(0.86905713) q[1];
x q[2];
rz(-1.1547791) q[3];
sx q[3];
rz(-2.3868594) q[3];
sx q[3];
rz(-2.5861135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43651849) q[2];
sx q[2];
rz(-1.7784092) q[2];
sx q[2];
rz(-2.3684033) q[2];
rz(0.1117205) q[3];
sx q[3];
rz(-1.3216647) q[3];
sx q[3];
rz(0.93938655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47950995) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(-0.42698419) q[0];
rz(2.2110979) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(0.46447909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1982959) q[0];
sx q[0];
rz(-1.3751327) q[0];
sx q[0];
rz(-2.8477459) q[0];
rz(-pi) q[1];
rz(2.1709178) q[2];
sx q[2];
rz(-1.8541186) q[2];
sx q[2];
rz(1.4221869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94333108) q[1];
sx q[1];
rz(-1.2020018) q[1];
sx q[1];
rz(1.53207) q[1];
rz(-pi) q[2];
rz(1.6397912) q[3];
sx q[3];
rz(-2.7164) q[3];
sx q[3];
rz(1.0591648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.53133196) q[2];
sx q[2];
rz(-1.320763) q[2];
sx q[2];
rz(0.28826928) q[2];
rz(1.0432976) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(-0.73879761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6601324) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(0.75755358) q[0];
rz(-3.0689012) q[1];
sx q[1];
rz(-0.025765954) q[1];
sx q[1];
rz(-0.048197897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91201997) q[0];
sx q[0];
rz(-0.98668146) q[0];
sx q[0];
rz(0.035365625) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23400314) q[2];
sx q[2];
rz(-2.1562139) q[2];
sx q[2];
rz(-0.3096748) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1485702) q[1];
sx q[1];
rz(-0.76495586) q[1];
sx q[1];
rz(2.2123442) q[1];
rz(-pi) q[2];
rz(-2.0406575) q[3];
sx q[3];
rz(-1.5239232) q[3];
sx q[3];
rz(-0.27134233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4797719) q[2];
sx q[2];
rz(-1.6419819) q[2];
sx q[2];
rz(-3.0721967) q[2];
rz(1.5540468) q[3];
sx q[3];
rz(-0.78726751) q[3];
sx q[3];
rz(2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7540392) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(-1.7283424) q[0];
rz(2.3130401) q[1];
sx q[1];
rz(-0.041752432) q[1];
sx q[1];
rz(-0.54263306) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12801816) q[0];
sx q[0];
rz(-1.4262222) q[0];
sx q[0];
rz(-1.6909063) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7650039) q[2];
sx q[2];
rz(-1.9053725) q[2];
sx q[2];
rz(1.1194816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9787776) q[1];
sx q[1];
rz(-0.12773027) q[1];
sx q[1];
rz(-0.43600299) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7067634) q[3];
sx q[3];
rz(-1.9169589) q[3];
sx q[3];
rz(1.0258254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98216206) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(2.8816667) q[2];
rz(2.121117) q[3];
sx q[3];
rz(-0.25769886) q[3];
sx q[3];
rz(-2.3475588) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771773) q[0];
sx q[0];
rz(-0.18615119) q[0];
sx q[0];
rz(-1.6625241) q[0];
rz(1.6089449) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(0.74302465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6830269) q[0];
sx q[0];
rz(-1.4560149) q[0];
sx q[0];
rz(1.0558788) q[0];
x q[1];
rz(-0.9833588) q[2];
sx q[2];
rz(-2.5134216) q[2];
sx q[2];
rz(0.52631015) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31752095) q[1];
sx q[1];
rz(-1.5050384) q[1];
sx q[1];
rz(-0.055067896) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0366692) q[3];
sx q[3];
rz(-1.5399974) q[3];
sx q[3];
rz(-1.2237751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6741901) q[2];
sx q[2];
rz(-2.3154066) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(1.449466) q[3];
sx q[3];
rz(-1.2286681) q[3];
sx q[3];
rz(0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529124) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(-2.3895277) q[0];
rz(1.9883142) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(-2.8582252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13334668) q[0];
sx q[0];
rz(-1.0336242) q[0];
sx q[0];
rz(-1.6744012) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81895639) q[2];
sx q[2];
rz(-1.6791846) q[2];
sx q[2];
rz(0.10486952) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9787406) q[1];
sx q[1];
rz(-1.8998441) q[1];
sx q[1];
rz(2.5025337) q[1];
rz(-1.9612736) q[3];
sx q[3];
rz(-0.9538528) q[3];
sx q[3];
rz(-0.65333593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2470384) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(-1.7029597) q[2];
rz(-2.8476207) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(-2.1132052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1080078) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(-2.602018) q[1];
sx q[1];
rz(-0.78782606) q[1];
sx q[1];
rz(0.16001564) q[1];
rz(0.14847427) q[2];
sx q[2];
rz(-1.9273026) q[2];
sx q[2];
rz(-2.9782563) q[2];
rz(0.60565518) q[3];
sx q[3];
rz(-0.45481053) q[3];
sx q[3];
rz(-0.74301471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
