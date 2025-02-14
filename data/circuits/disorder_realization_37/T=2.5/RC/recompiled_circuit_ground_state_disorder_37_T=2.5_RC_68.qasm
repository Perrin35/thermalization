OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78703824) q[0];
sx q[0];
rz(-0.7631425) q[0];
sx q[0];
rz(0.95844498) q[0];
rz(0.42674843) q[1];
sx q[1];
rz(-0.42438212) q[1];
sx q[1];
rz(2.1688865) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866731) q[0];
sx q[0];
rz(-1.3698319) q[0];
sx q[0];
rz(1.475564) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3742674) q[2];
sx q[2];
rz(-2.9568045) q[2];
sx q[2];
rz(-1.5306461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8805566) q[1];
sx q[1];
rz(-2.2734005) q[1];
sx q[1];
rz(-1.9126585) q[1];
rz(1.7250502) q[3];
sx q[3];
rz(-1.8498382) q[3];
sx q[3];
rz(-2.0786301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.080231754) q[2];
sx q[2];
rz(-1.1876567) q[2];
sx q[2];
rz(-2.3415671) q[2];
rz(-0.13036615) q[3];
sx q[3];
rz(-0.76821199) q[3];
sx q[3];
rz(0.31726328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80113634) q[0];
sx q[0];
rz(-1.9213333) q[0];
sx q[0];
rz(2.1521547) q[0];
rz(2.2438352) q[1];
sx q[1];
rz(-2.0866626) q[1];
sx q[1];
rz(-0.78136939) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2563574) q[0];
sx q[0];
rz(-1.7634749) q[0];
sx q[0];
rz(-1.0568134) q[0];
rz(0.35542458) q[2];
sx q[2];
rz(-2.8925507) q[2];
sx q[2];
rz(0.058573478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20618379) q[1];
sx q[1];
rz(-2.0061135) q[1];
sx q[1];
rz(-1.0401506) q[1];
rz(-pi) q[2];
rz(2.492401) q[3];
sx q[3];
rz(-0.4684557) q[3];
sx q[3];
rz(0.89391232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.117131) q[2];
sx q[2];
rz(-1.7639561) q[2];
sx q[2];
rz(0.77829877) q[2];
rz(0.69027573) q[3];
sx q[3];
rz(-2.2199151) q[3];
sx q[3];
rz(1.5811623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86163259) q[0];
sx q[0];
rz(-0.83502382) q[0];
sx q[0];
rz(-2.7296208) q[0];
rz(-2.1906134) q[1];
sx q[1];
rz(-2.488766) q[1];
sx q[1];
rz(1.9297809) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6854343) q[0];
sx q[0];
rz(-3.0904909) q[0];
sx q[0];
rz(2.2888378) q[0];
rz(-1.3479718) q[2];
sx q[2];
rz(-1.6671902) q[2];
sx q[2];
rz(0.69144648) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6784994) q[1];
sx q[1];
rz(-2.2804567) q[1];
sx q[1];
rz(-1.6072431) q[1];
x q[2];
rz(2.7762967) q[3];
sx q[3];
rz(-1.4202227) q[3];
sx q[3];
rz(0.07380658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95436207) q[2];
sx q[2];
rz(-0.53202859) q[2];
sx q[2];
rz(-1.5032035) q[2];
rz(0.040180834) q[3];
sx q[3];
rz(-0.92622042) q[3];
sx q[3];
rz(-0.23475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16315854) q[0];
sx q[0];
rz(-1.6272767) q[0];
sx q[0];
rz(-3.0117595) q[0];
rz(-1.2854598) q[1];
sx q[1];
rz(-1.1760271) q[1];
sx q[1];
rz(-1.0349549) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87405768) q[0];
sx q[0];
rz(-0.77299905) q[0];
sx q[0];
rz(-0.92138793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1088684) q[2];
sx q[2];
rz(-2.2316859) q[2];
sx q[2];
rz(-2.3488597) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.49832) q[1];
sx q[1];
rz(-0.26277143) q[1];
sx q[1];
rz(-0.33518016) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47607039) q[3];
sx q[3];
rz(-2.2048773) q[3];
sx q[3];
rz(2.7345524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3639565) q[2];
sx q[2];
rz(-1.3866321) q[2];
sx q[2];
rz(-0.094495471) q[2];
rz(-2.7206521) q[3];
sx q[3];
rz(-1.4237483) q[3];
sx q[3];
rz(1.4360992) q[3];
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
rz(2.7218551) q[0];
sx q[0];
rz(-0.57191816) q[0];
sx q[0];
rz(-0.016059248) q[0];
rz(0.84469604) q[1];
sx q[1];
rz(-0.41929308) q[1];
sx q[1];
rz(2.2339581) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.006165531) q[0];
sx q[0];
rz(-0.79815292) q[0];
sx q[0];
rz(0.31240518) q[0];
rz(-pi) q[1];
rz(-2.0727022) q[2];
sx q[2];
rz(-2.6703983) q[2];
sx q[2];
rz(2.1387177) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.065103) q[1];
sx q[1];
rz(-1.2416844) q[1];
sx q[1];
rz(-2.0338716) q[1];
x q[2];
rz(2.7974878) q[3];
sx q[3];
rz(-1.5457166) q[3];
sx q[3];
rz(0.6911975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5661261) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(-3.0779085) q[2];
rz(-2.5746386) q[3];
sx q[3];
rz(-2.3072115) q[3];
sx q[3];
rz(-0.59757346) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31396922) q[0];
sx q[0];
rz(-0.14179985) q[0];
sx q[0];
rz(1.1578479) q[0];
rz(1.4843548) q[1];
sx q[1];
rz(-1.6969705) q[1];
sx q[1];
rz(-0.2690014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3753178) q[0];
sx q[0];
rz(-0.14785375) q[0];
sx q[0];
rz(-2.8327441) q[0];
rz(1.2397175) q[2];
sx q[2];
rz(-0.28713687) q[2];
sx q[2];
rz(1.7652709) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.288876) q[1];
sx q[1];
rz(-1.7501131) q[1];
sx q[1];
rz(-3.0406171) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41407164) q[3];
sx q[3];
rz(-1.3706638) q[3];
sx q[3];
rz(-2.4614863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3645662) q[2];
sx q[2];
rz(-0.51743162) q[2];
sx q[2];
rz(-0.88097921) q[2];
rz(0.12428728) q[3];
sx q[3];
rz(-1.4275987) q[3];
sx q[3];
rz(-2.4833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2800901) q[0];
sx q[0];
rz(-0.26743356) q[0];
sx q[0];
rz(2.4524443) q[0];
rz(2.6742477) q[1];
sx q[1];
rz(-0.57599774) q[1];
sx q[1];
rz(-1.0980094) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9257346) q[0];
sx q[0];
rz(-1.4924865) q[0];
sx q[0];
rz(-2.9758478) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46862342) q[2];
sx q[2];
rz(-0.59983095) q[2];
sx q[2];
rz(-3.1214903) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5913691) q[1];
sx q[1];
rz(-2.6003631) q[1];
sx q[1];
rz(3.079097) q[1];
rz(-pi) q[2];
rz(-0.92236788) q[3];
sx q[3];
rz(-2.4019314) q[3];
sx q[3];
rz(1.9112223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5757307) q[2];
sx q[2];
rz(-0.87811676) q[2];
sx q[2];
rz(2.511054) q[2];
rz(0.60728836) q[3];
sx q[3];
rz(-0.90689617) q[3];
sx q[3];
rz(-1.0926532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31834114) q[0];
sx q[0];
rz(-2.9956151) q[0];
sx q[0];
rz(-1.7952221) q[0];
rz(1.775555) q[1];
sx q[1];
rz(-0.64907688) q[1];
sx q[1];
rz(-0.62873658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36489428) q[0];
sx q[0];
rz(-1.3140251) q[0];
sx q[0];
rz(0.55112324) q[0];
rz(-0.95405719) q[2];
sx q[2];
rz(-2.3696757) q[2];
sx q[2];
rz(1.7454912) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4688411) q[1];
sx q[1];
rz(-2.3769925) q[1];
sx q[1];
rz(-2.4048664) q[1];
rz(-1.2090861) q[3];
sx q[3];
rz(-2.9291398) q[3];
sx q[3];
rz(-0.3750876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8667355) q[2];
sx q[2];
rz(-1.8161512) q[2];
sx q[2];
rz(-0.87230116) q[2];
rz(-2.9288779) q[3];
sx q[3];
rz(-1.3958967) q[3];
sx q[3];
rz(0.58342903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23671959) q[0];
sx q[0];
rz(-1.5308335) q[0];
sx q[0];
rz(-1.2637631) q[0];
rz(-2.1243375) q[1];
sx q[1];
rz(-0.89824289) q[1];
sx q[1];
rz(-0.57473007) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59422661) q[0];
sx q[0];
rz(-1.5061139) q[0];
sx q[0];
rz(-2.6101739) q[0];
rz(-pi) q[1];
rz(-1.310964) q[2];
sx q[2];
rz(-0.75922478) q[2];
sx q[2];
rz(2.7936943) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3176665) q[1];
sx q[1];
rz(-1.9374018) q[1];
sx q[1];
rz(1.7714785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25463661) q[3];
sx q[3];
rz(-2.0395425) q[3];
sx q[3];
rz(2.180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42402521) q[2];
sx q[2];
rz(-1.307345) q[2];
sx q[2];
rz(3.1094816) q[2];
rz(1.9854246) q[3];
sx q[3];
rz(-0.38438946) q[3];
sx q[3];
rz(1.9305852) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1011937) q[0];
sx q[0];
rz(-0.54553425) q[0];
sx q[0];
rz(-2.0174761) q[0];
rz(-2.595937) q[1];
sx q[1];
rz(-0.35174313) q[1];
sx q[1];
rz(2.5901332) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1301981) q[0];
sx q[0];
rz(-1.4470248) q[0];
sx q[0];
rz(-1.3817203) q[0];
rz(-pi) q[1];
rz(-1.2966424) q[2];
sx q[2];
rz(-1.621843) q[2];
sx q[2];
rz(2.5921043) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6956967) q[1];
sx q[1];
rz(-2.4425382) q[1];
sx q[1];
rz(0.76837825) q[1];
rz(-1.7343246) q[3];
sx q[3];
rz(-1.8732605) q[3];
sx q[3];
rz(1.0176942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5896899) q[2];
sx q[2];
rz(-2.6305113) q[2];
sx q[2];
rz(-1.4116633) q[2];
rz(2.3949413) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(0.97486973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31155561) q[0];
sx q[0];
rz(-1.5604326) q[0];
sx q[0];
rz(-1.5695705) q[0];
rz(-0.39881067) q[1];
sx q[1];
rz(-1.8067982) q[1];
sx q[1];
rz(-1.6511818) q[1];
rz(-1.5560935) q[2];
sx q[2];
rz(-0.48853816) q[2];
sx q[2];
rz(1.9449816) q[2];
rz(0.29295425) q[3];
sx q[3];
rz(-2.3200547) q[3];
sx q[3];
rz(-3.0153081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
