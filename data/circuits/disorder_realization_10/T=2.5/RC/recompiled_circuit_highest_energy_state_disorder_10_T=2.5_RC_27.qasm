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
rz(-0.69775692) q[0];
sx q[0];
rz(-1.193576) q[0];
sx q[0];
rz(2.9370263) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(1.7525571) q[1];
sx q[1];
rz(11.166823) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7273063) q[0];
sx q[0];
rz(-2.901863) q[0];
sx q[0];
rz(1.9637462) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40488196) q[2];
sx q[2];
rz(-1.153115) q[2];
sx q[2];
rz(-1.4428802) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.39069893) q[1];
sx q[1];
rz(-1.9913773) q[1];
sx q[1];
rz(1.4750359) q[1];
x q[2];
rz(1.8379962) q[3];
sx q[3];
rz(-2.2624216) q[3];
sx q[3];
rz(-1.8593781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(0.48933634) q[2];
rz(-2.1183744) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(1.1255012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045687549) q[0];
sx q[0];
rz(-2.4850595) q[0];
sx q[0];
rz(2.3387961) q[0];
rz(-3.070201) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(-3.0838222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7785955) q[0];
sx q[0];
rz(-0.61041202) q[0];
sx q[0];
rz(-1.1108062) q[0];
x q[1];
rz(-0.28654307) q[2];
sx q[2];
rz(-1.6725958) q[2];
sx q[2];
rz(-1.2004146) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80177414) q[1];
sx q[1];
rz(-1.8027824) q[1];
sx q[1];
rz(1.07527) q[1];
rz(-pi) q[2];
rz(1.7419192) q[3];
sx q[3];
rz(-2.4226396) q[3];
sx q[3];
rz(-2.3477682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-1.095093) q[2];
sx q[2];
rz(-1.2954953) q[2];
rz(-2.1748491) q[3];
sx q[3];
rz(-2.3714378) q[3];
sx q[3];
rz(2.3841592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.029723786) q[0];
sx q[0];
rz(-1.3628549) q[0];
sx q[0];
rz(1.7080074) q[0];
rz(-0.068610527) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(-0.57919085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5369956) q[0];
sx q[0];
rz(-1.4791577) q[0];
sx q[0];
rz(3.0834404) q[0];
rz(-2.6658695) q[2];
sx q[2];
rz(-1.5724725) q[2];
sx q[2];
rz(-2.7762108) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2229742) q[1];
sx q[1];
rz(-0.22906216) q[1];
sx q[1];
rz(-0.15840662) q[1];
rz(-pi) q[2];
rz(2.6853232) q[3];
sx q[3];
rz(-1.1351403) q[3];
sx q[3];
rz(1.6317451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.27089831) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(2.9581621) q[2];
rz(-0.80426788) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(-2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.19876984) q[0];
sx q[0];
rz(-3.0267921) q[0];
sx q[0];
rz(-2.542069) q[0];
rz(0.41246688) q[1];
sx q[1];
rz(-0.020655276) q[1];
sx q[1];
rz(-1.0106769) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0507999) q[0];
sx q[0];
rz(-0.37380344) q[0];
sx q[0];
rz(-2.7888377) q[0];
rz(1.5013736) q[2];
sx q[2];
rz(-0.63328082) q[2];
sx q[2];
rz(-0.4236003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4356723) q[1];
sx q[1];
rz(-1.3002035) q[1];
sx q[1];
rz(-2.1822004) q[1];
rz(-pi) q[2];
rz(0.75586163) q[3];
sx q[3];
rz(-1.469765) q[3];
sx q[3];
rz(-1.7441235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3707054) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(-0.31751219) q[2];
rz(-0.62234771) q[3];
sx q[3];
rz(-0.93572664) q[3];
sx q[3];
rz(0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36251003) q[0];
sx q[0];
rz(-0.91479397) q[0];
sx q[0];
rz(-1.0269748) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-0.06409476) q[1];
sx q[1];
rz(1.2170353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.352223) q[0];
sx q[0];
rz(-0.93607157) q[0];
sx q[0];
rz(-0.85294368) q[0];
x q[1];
rz(-3.0116357) q[2];
sx q[2];
rz(-0.54207793) q[2];
sx q[2];
rz(0.96058577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1149619) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(0.86905713) q[1];
rz(1.1547791) q[3];
sx q[3];
rz(-2.3868594) q[3];
sx q[3];
rz(-0.55547914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7050742) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(0.77318937) q[2];
rz(-0.1117205) q[3];
sx q[3];
rz(-1.3216647) q[3];
sx q[3];
rz(-0.93938655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47950995) q[0];
sx q[0];
rz(-2.918512) q[0];
sx q[0];
rz(2.7146085) q[0];
rz(-0.93049479) q[1];
sx q[1];
rz(-3.1246779) q[1];
sx q[1];
rz(2.6771136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94356774) q[0];
sx q[0];
rz(-2.7901398) q[0];
sx q[0];
rz(0.60011421) q[0];
rz(-2.802425) q[2];
sx q[2];
rz(-2.1438823) q[2];
sx q[2];
rz(3.1011875) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0911922) q[1];
sx q[1];
rz(-0.37072966) q[1];
sx q[1];
rz(-3.0417473) q[1];
rz(-0.031207009) q[3];
sx q[3];
rz(-1.1466807) q[3];
sx q[3];
rz(-1.1348789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6102607) q[2];
sx q[2];
rz(-1.8208296) q[2];
sx q[2];
rz(-0.28826928) q[2];
rz(-1.0432976) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(0.73879761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4814602) q[0];
sx q[0];
rz(-1.8539424) q[0];
sx q[0];
rz(-2.3840391) q[0];
rz(3.0689012) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(3.0933948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5023247) q[0];
sx q[0];
rz(-1.5412962) q[0];
sx q[0];
rz(0.98639368) q[0];
x q[1];
rz(-0.9725857) q[2];
sx q[2];
rz(-1.3763104) q[2];
sx q[2];
rz(-1.1301646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.072025) q[1];
sx q[1];
rz(-1.9981019) q[1];
sx q[1];
rz(-0.91520379) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0890325) q[3];
sx q[3];
rz(-2.0400999) q[3];
sx q[3];
rz(1.3232425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4797719) q[2];
sx q[2];
rz(-1.6419819) q[2];
sx q[2];
rz(-3.0721967) q[2];
rz(1.5540468) q[3];
sx q[3];
rz(-2.3543251) q[3];
sx q[3];
rz(-2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3875535) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(1.4132502) q[0];
rz(0.82855254) q[1];
sx q[1];
rz(-0.041752432) q[1];
sx q[1];
rz(0.54263306) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5722902) q[0];
sx q[0];
rz(-0.18768947) q[0];
sx q[0];
rz(2.45298) q[0];
x q[1];
rz(1.9285525) q[2];
sx q[2];
rz(-1.925549) q[2];
sx q[2];
rz(0.32217978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1664913) q[1];
sx q[1];
rz(-1.5169739) q[1];
sx q[1];
rz(-0.11589284) q[1];
x q[2];
rz(2.4335008) q[3];
sx q[3];
rz(-2.5927784) q[3];
sx q[3];
rz(0.08591692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1594306) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(0.25992599) q[2];
rz(-1.0204756) q[3];
sx q[3];
rz(-0.25769886) q[3];
sx q[3];
rz(-2.3475588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771773) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(1.4790685) q[0];
rz(1.6089449) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(0.74302465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047516454) q[0];
sx q[0];
rz(-1.0595982) q[0];
sx q[0];
rz(0.1316977) q[0];
rz(-2.114562) q[2];
sx q[2];
rz(-1.2390422) q[2];
sx q[2];
rz(-0.55032718) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0155269) q[1];
sx q[1];
rz(-0.085745009) q[1];
sx q[1];
rz(0.87460204) q[1];
x q[2];
rz(-0.10492341) q[3];
sx q[3];
rz(-1.6015953) q[3];
sx q[3];
rz(1.9178176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4674025) q[2];
sx q[2];
rz(-2.3154066) q[2];
sx q[2];
rz(-2.3361333) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.2286681) q[3];
sx q[3];
rz(-2.3384371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886803) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(-2.3895277) q[0];
rz(1.9883142) q[1];
sx q[1];
rz(-0.88429943) q[1];
sx q[1];
rz(-0.28336743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8077791) q[0];
sx q[0];
rz(-2.5954843) q[0];
sx q[0];
rz(0.17192745) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7287909) q[2];
sx q[2];
rz(-0.758095) q[2];
sx q[2];
rz(-1.5604863) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1438469) q[1];
sx q[1];
rz(-2.4334868) q[1];
sx q[1];
rz(2.6216402) q[1];
rz(1.180319) q[3];
sx q[3];
rz(-2.1877398) q[3];
sx q[3];
rz(0.65333593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89455426) q[2];
sx q[2];
rz(-3.0587695) q[2];
sx q[2];
rz(-1.438633) q[2];
rz(2.8476207) q[3];
sx q[3];
rz(-0.014466244) q[3];
sx q[3];
rz(1.0283874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.033584874) q[0];
sx q[0];
rz(-1.7243732) q[0];
sx q[0];
rz(1.617817) q[0];
rz(-0.53957466) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(1.930936) q[2];
sx q[2];
rz(-1.4317206) q[2];
sx q[2];
rz(1.7862873) q[2];
rz(2.7593437) q[3];
sx q[3];
rz(-1.8235689) q[3];
sx q[3];
rz(0.27128661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
