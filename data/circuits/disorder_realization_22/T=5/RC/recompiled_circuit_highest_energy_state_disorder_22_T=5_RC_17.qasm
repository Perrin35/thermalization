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
rz(-0.67996112) q[0];
sx q[0];
rz(-0.90587076) q[0];
sx q[0];
rz(-1.0933956) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(-1.1727762) q[1];
sx q[1];
rz(2.7170031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49219201) q[0];
sx q[0];
rz(-1.369242) q[0];
sx q[0];
rz(-3.1026476) q[0];
rz(-3.0908747) q[2];
sx q[2];
rz(-1.4169622) q[2];
sx q[2];
rz(-0.62062973) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0478617) q[1];
sx q[1];
rz(-1.8163396) q[1];
sx q[1];
rz(3.0603916) q[1];
rz(-3.0077259) q[3];
sx q[3];
rz(-2.2765719) q[3];
sx q[3];
rz(5.2701252e-05) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37497416) q[2];
sx q[2];
rz(-1.5593636) q[2];
sx q[2];
rz(1.9826822) q[2];
rz(-2.9186987) q[3];
sx q[3];
rz(-1.4054207) q[3];
sx q[3];
rz(-2.8515653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952154) q[0];
sx q[0];
rz(-1.5257436) q[0];
sx q[0];
rz(-2.0624397) q[0];
rz(1.7201299) q[1];
sx q[1];
rz(-2.454897) q[1];
sx q[1];
rz(-0.064528331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9031398) q[0];
sx q[0];
rz(-1.6232326) q[0];
sx q[0];
rz(-1.658123) q[0];
rz(2.6653397) q[2];
sx q[2];
rz(-1.5459527) q[2];
sx q[2];
rz(-2.3551539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3694683) q[1];
sx q[1];
rz(-1.3288979) q[1];
sx q[1];
rz(2.3568704) q[1];
rz(-pi) q[2];
rz(1.7713624) q[3];
sx q[3];
rz(-0.75955694) q[3];
sx q[3];
rz(-0.49155363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9885538) q[2];
sx q[2];
rz(-1.7861853) q[2];
sx q[2];
rz(-2.0241731) q[2];
rz(-2.2828263) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(0.43706885) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7607464) q[0];
sx q[0];
rz(-2.8570211) q[0];
sx q[0];
rz(2.9685156) q[0];
rz(2.1848047) q[1];
sx q[1];
rz(-2.1134977) q[1];
sx q[1];
rz(-2.079336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040471023) q[0];
sx q[0];
rz(-0.72325828) q[0];
sx q[0];
rz(-1.5803807) q[0];
x q[1];
rz(-2.4281279) q[2];
sx q[2];
rz(-0.77859813) q[2];
sx q[2];
rz(0.32003357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7159285) q[1];
sx q[1];
rz(-0.93837184) q[1];
sx q[1];
rz(-2.4584437) q[1];
x q[2];
rz(2.297256) q[3];
sx q[3];
rz(-0.82714836) q[3];
sx q[3];
rz(-1.8238123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.411285) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(-2.3968706) q[2];
rz(-2.5005285) q[3];
sx q[3];
rz(-1.791626) q[3];
sx q[3];
rz(0.49066576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500126) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(-0.84521729) q[0];
rz(-0.40329626) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(2.8135615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9739537) q[0];
sx q[0];
rz(-0.22952794) q[0];
sx q[0];
rz(2.2295802) q[0];
x q[1];
rz(0.26706605) q[2];
sx q[2];
rz(-0.24644463) q[2];
sx q[2];
rz(1.7676958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9765342) q[1];
sx q[1];
rz(-1.4529902) q[1];
sx q[1];
rz(-0.56748135) q[1];
x q[2];
rz(-1.1370522) q[3];
sx q[3];
rz(-2.5528926) q[3];
sx q[3];
rz(-0.91982809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8959117) q[2];
sx q[2];
rz(-2.3978105) q[2];
sx q[2];
rz(0.15920676) q[2];
rz(3.1253452) q[3];
sx q[3];
rz(-2.1135606) q[3];
sx q[3];
rz(0.73133674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943587) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(1.0900981) q[0];
rz(0.12995003) q[1];
sx q[1];
rz(-0.81472412) q[1];
sx q[1];
rz(1.5257588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1795883) q[0];
sx q[0];
rz(-1.1561014) q[0];
sx q[0];
rz(2.5432822) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4065625) q[2];
sx q[2];
rz(-1.4964074) q[2];
sx q[2];
rz(2.1321572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8860717) q[1];
sx q[1];
rz(-2.204245) q[1];
sx q[1];
rz(1.6048163) q[1];
rz(-pi) q[2];
rz(1.3899743) q[3];
sx q[3];
rz(-1.5360334) q[3];
sx q[3];
rz(1.3270448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7502363) q[2];
sx q[2];
rz(-0.82166925) q[2];
sx q[2];
rz(-0.4772805) q[2];
rz(2.9376049) q[3];
sx q[3];
rz(-2.9450649) q[3];
sx q[3];
rz(-2.5175214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9452962) q[0];
sx q[0];
rz(-2.3260703) q[0];
sx q[0];
rz(-0.77769172) q[0];
rz(0.6174736) q[1];
sx q[1];
rz(-1.6505046) q[1];
sx q[1];
rz(-0.9309353) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6361178) q[0];
sx q[0];
rz(-1.1606998) q[0];
sx q[0];
rz(2.0652886) q[0];
rz(-pi) q[1];
rz(-2.2799643) q[2];
sx q[2];
rz(-0.59979931) q[2];
sx q[2];
rz(1.3280363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41250788) q[1];
sx q[1];
rz(-1.2326483) q[1];
sx q[1];
rz(-0.04674712) q[1];
rz(2.191675) q[3];
sx q[3];
rz(-1.2395879) q[3];
sx q[3];
rz(-0.036605926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65333873) q[2];
sx q[2];
rz(-1.808017) q[2];
sx q[2];
rz(2.874157) q[2];
rz(-0.93287647) q[3];
sx q[3];
rz(-1.8578015) q[3];
sx q[3];
rz(-0.4944087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47144181) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(-0.66584051) q[0];
rz(-2.937607) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(1.1873672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0627619) q[0];
sx q[0];
rz(-1.4785066) q[0];
sx q[0];
rz(-0.4706443) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86249528) q[2];
sx q[2];
rz(-1.1581026) q[2];
sx q[2];
rz(0.19315456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5449808) q[1];
sx q[1];
rz(-2.2447733) q[1];
sx q[1];
rz(1.9121714) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.051249226) q[3];
sx q[3];
rz(-1.3173242) q[3];
sx q[3];
rz(-0.26904256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5953956) q[2];
sx q[2];
rz(-2.6232145) q[2];
sx q[2];
rz(0.67016822) q[2];
rz(-2.8403122) q[3];
sx q[3];
rz(-1.8471085) q[3];
sx q[3];
rz(-0.41845751) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8381074) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(-1.0656892) q[0];
rz(-0.65713716) q[1];
sx q[1];
rz(-0.70825759) q[1];
sx q[1];
rz(1.7117975) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8891212) q[0];
sx q[0];
rz(-1.234113) q[0];
sx q[0];
rz(-2.1648429) q[0];
rz(2.6515342) q[2];
sx q[2];
rz(-2.769751) q[2];
sx q[2];
rz(-1.3272367) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.011466) q[1];
sx q[1];
rz(-1.3523755) q[1];
sx q[1];
rz(-0.41832025) q[1];
x q[2];
rz(-0.64770697) q[3];
sx q[3];
rz(-2.9053024) q[3];
sx q[3];
rz(-1.7648197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8810205) q[2];
sx q[2];
rz(-1.2994956) q[2];
sx q[2];
rz(2.1232429) q[2];
rz(-1.5923422) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(2.3840267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757979) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(0.96762586) q[0];
rz(-0.34010092) q[1];
sx q[1];
rz(-0.83931559) q[1];
sx q[1];
rz(-2.1077154) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32717413) q[0];
sx q[0];
rz(-2.5670739) q[0];
sx q[0];
rz(-1.8273201) q[0];
x q[1];
rz(-0.56615717) q[2];
sx q[2];
rz(-0.70933178) q[2];
sx q[2];
rz(-2.4475803) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.25701538) q[1];
sx q[1];
rz(-2.3126763) q[1];
sx q[1];
rz(-2.1121426) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7179862) q[3];
sx q[3];
rz(-1.3491674) q[3];
sx q[3];
rz(-2.2972884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0391482) q[2];
sx q[2];
rz(-1.6951963) q[2];
sx q[2];
rz(0.038912494) q[2];
rz(0.14032042) q[3];
sx q[3];
rz(-3.0571627) q[3];
sx q[3];
rz(-1.8261568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4843531) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(1.8835618) q[0];
rz(-0.21615061) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(-0.36453882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70750694) q[0];
sx q[0];
rz(-2.0179406) q[0];
sx q[0];
rz(2.8672877) q[0];
rz(2.3380525) q[2];
sx q[2];
rz(-1.8992956) q[2];
sx q[2];
rz(-0.9847509) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8748624) q[1];
sx q[1];
rz(-1.5577496) q[1];
sx q[1];
rz(0.055582837) q[1];
x q[2];
rz(-0.057897827) q[3];
sx q[3];
rz(-1.5385043) q[3];
sx q[3];
rz(0.81626371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8517427) q[2];
sx q[2];
rz(-1.2052636) q[2];
sx q[2];
rz(2.5002948) q[2];
rz(1.6614206) q[3];
sx q[3];
rz(-1.2484173) q[3];
sx q[3];
rz(-2.4680468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4079473) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(0.98450487) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(2.099809) q[2];
sx q[2];
rz(-0.22780252) q[2];
sx q[2];
rz(-3.0115423) q[2];
rz(-2.051864) q[3];
sx q[3];
rz(-1.0952598) q[3];
sx q[3];
rz(2.918603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
