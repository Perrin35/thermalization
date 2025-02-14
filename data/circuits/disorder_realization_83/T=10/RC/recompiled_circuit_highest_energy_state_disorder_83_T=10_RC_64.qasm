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
rz(-0.42521617) q[0];
sx q[0];
rz(4.4758237) q[0];
sx q[0];
rz(9.0444179) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(0.92460728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3596852) q[0];
sx q[0];
rz(-1.9115907) q[0];
sx q[0];
rz(-0.73990344) q[0];
rz(-pi) q[1];
rz(1.2819498) q[2];
sx q[2];
rz(-0.53448662) q[2];
sx q[2];
rz(-1.6030902) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8262337) q[1];
sx q[1];
rz(-0.72972882) q[1];
sx q[1];
rz(-1.6297831) q[1];
rz(-3.0136631) q[3];
sx q[3];
rz(-1.5618213) q[3];
sx q[3];
rz(0.25024807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(2.1535786) q[2];
rz(2.9347349) q[3];
sx q[3];
rz(-1.4068973) q[3];
sx q[3];
rz(-0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67966953) q[0];
sx q[0];
rz(-1.672687) q[0];
sx q[0];
rz(-2.8894506) q[0];
rz(-3.120046) q[1];
sx q[1];
rz(-2.4485059) q[1];
sx q[1];
rz(0.85539877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54411201) q[0];
sx q[0];
rz(-1.5978088) q[0];
sx q[0];
rz(-3.1400694) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1507262) q[2];
sx q[2];
rz(-2.0850811) q[2];
sx q[2];
rz(-0.95555238) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2829078) q[1];
sx q[1];
rz(-2.5577684) q[1];
sx q[1];
rz(-1.1633384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40327252) q[3];
sx q[3];
rz(-2.4174815) q[3];
sx q[3];
rz(-1.4693174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1713193) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-2.8805736) q[2];
rz(0.039693443) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(-0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.71565851) q[0];
sx q[0];
rz(-1.6682699) q[0];
sx q[0];
rz(-2.4470827) q[0];
rz(-2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(-0.99004254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8111525) q[0];
sx q[0];
rz(-1.9188723) q[0];
sx q[0];
rz(-1.6578107) q[0];
rz(-pi) q[1];
rz(-1.2743397) q[2];
sx q[2];
rz(-0.60288376) q[2];
sx q[2];
rz(2.4566513) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1676869) q[1];
sx q[1];
rz(-1.1264633) q[1];
sx q[1];
rz(1.1717887) q[1];
rz(-0.95616266) q[3];
sx q[3];
rz(-2.3376138) q[3];
sx q[3];
rz(1.1756736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29828829) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(-2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0855584) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(0.84247843) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(-0.78559771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89891035) q[0];
sx q[0];
rz(-1.8225095) q[0];
sx q[0];
rz(-2.4221735) q[0];
rz(-0.71707861) q[2];
sx q[2];
rz(-1.684965) q[2];
sx q[2];
rz(-1.9463429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.848104) q[1];
sx q[1];
rz(-1.2624718) q[1];
sx q[1];
rz(-0.56667324) q[1];
rz(-1.3520283) q[3];
sx q[3];
rz(-0.76319088) q[3];
sx q[3];
rz(1.4783875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0164612) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(2.1330736) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763024) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(2.0710131) q[0];
rz(-2.3640682) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(-2.1308965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99791894) q[0];
sx q[0];
rz(-2.5287313) q[0];
sx q[0];
rz(-0.65632239) q[0];
x q[1];
rz(-0.065344409) q[2];
sx q[2];
rz(-1.3962708) q[2];
sx q[2];
rz(-0.068152817) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4324716) q[1];
sx q[1];
rz(-0.84408954) q[1];
sx q[1];
rz(2.5120887) q[1];
x q[2];
rz(2.5794677) q[3];
sx q[3];
rz(-1.486612) q[3];
sx q[3];
rz(0.69417324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(1.5555752) q[2];
rz(-1.0157061) q[3];
sx q[3];
rz(-2.000122) q[3];
sx q[3];
rz(-1.7293845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79065901) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(-0.80107981) q[0];
rz(-0.13310295) q[1];
sx q[1];
rz(-1.3214654) q[1];
sx q[1];
rz(-0.4020234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659365) q[0];
sx q[0];
rz(-0.84235993) q[0];
sx q[0];
rz(0.45439646) q[0];
rz(0.0021462321) q[2];
sx q[2];
rz(-1.2800001) q[2];
sx q[2];
rz(-0.69974297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.133144) q[1];
sx q[1];
rz(-0.56479543) q[1];
sx q[1];
rz(-2.1022335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65816718) q[3];
sx q[3];
rz(-1.4716545) q[3];
sx q[3];
rz(-2.1226573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63783995) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(0.77581882) q[2];
rz(1.4555629) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0966454) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(0.67895472) q[0];
rz(-1.8567765) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(-1.7624034) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75684568) q[0];
sx q[0];
rz(-0.79233903) q[0];
sx q[0];
rz(-1.8154665) q[0];
x q[1];
rz(3.0581362) q[2];
sx q[2];
rz(-1.2922568) q[2];
sx q[2];
rz(-2.4065774) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1759291) q[1];
sx q[1];
rz(-1.129985) q[1];
sx q[1];
rz(-0.55056527) q[1];
rz(-3.0857063) q[3];
sx q[3];
rz(-1.3709785) q[3];
sx q[3];
rz(0.36147396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3369559) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(1.9026559) q[2];
rz(1.8094481) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(3.0031257) q[0];
rz(-2.9578517) q[1];
sx q[1];
rz(-0.67444363) q[1];
sx q[1];
rz(-0.26652452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0040695) q[0];
sx q[0];
rz(-2.7752987) q[0];
sx q[0];
rz(2.0384602) q[0];
rz(-1.0037854) q[2];
sx q[2];
rz(-1.8006174) q[2];
sx q[2];
rz(-0.15574317) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9723818) q[1];
sx q[1];
rz(-1.5145497) q[1];
sx q[1];
rz(1.9906269) q[1];
x q[2];
rz(1.5620392) q[3];
sx q[3];
rz(-2.2976795) q[3];
sx q[3];
rz(1.7576287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0515685) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(-0.23452342) q[2];
rz(1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(2.3852111) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33673564) q[0];
sx q[0];
rz(-0.8661626) q[0];
sx q[0];
rz(2.3418703) q[0];
rz(0.58894482) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(-2.934093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28593913) q[0];
sx q[0];
rz(-0.41623273) q[0];
sx q[0];
rz(-2.8134384) q[0];
x q[1];
rz(-2.8325808) q[2];
sx q[2];
rz(-0.88243077) q[2];
sx q[2];
rz(-2.7613044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55940699) q[1];
sx q[1];
rz(-1.9684125) q[1];
sx q[1];
rz(0.84487265) q[1];
rz(0.028771632) q[3];
sx q[3];
rz(-1.1011916) q[3];
sx q[3];
rz(2.4503277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80692446) q[2];
sx q[2];
rz(-2.25756) q[2];
sx q[2];
rz(-0.36726382) q[2];
rz(-1.7043097) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(-1.9678763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2631898) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(2.3314085) q[0];
rz(-2.0267678) q[1];
sx q[1];
rz(-0.38433847) q[1];
sx q[1];
rz(-2.5883163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081757717) q[0];
sx q[0];
rz(-2.5011269) q[0];
sx q[0];
rz(1.4778796) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94925868) q[2];
sx q[2];
rz(-2.261603) q[2];
sx q[2];
rz(1.1101646) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9844891) q[1];
sx q[1];
rz(-1.2890639) q[1];
sx q[1];
rz(-0.19428044) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35240473) q[3];
sx q[3];
rz(-0.81768546) q[3];
sx q[3];
rz(-1.529983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1222003) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(-0.23966399) q[2];
rz(-1.3278809) q[3];
sx q[3];
rz(-1.0646822) q[3];
sx q[3];
rz(-2.6740668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0620621) q[0];
sx q[0];
rz(-1.3539599) q[0];
sx q[0];
rz(0.39394105) q[0];
rz(2.486034) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(-2.4782933) q[2];
sx q[2];
rz(-0.51021432) q[2];
sx q[2];
rz(-2.9205657) q[2];
rz(-0.58047337) q[3];
sx q[3];
rz(-1.7844641) q[3];
sx q[3];
rz(-2.6337558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
