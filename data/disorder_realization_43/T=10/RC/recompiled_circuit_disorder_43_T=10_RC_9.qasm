OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2551978) q[0];
sx q[0];
rz(-1.891341) q[0];
sx q[0];
rz(1.3347081) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(0.97595739) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0358409) q[0];
sx q[0];
rz(-0.63397206) q[0];
sx q[0];
rz(-0.41497725) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1337778) q[2];
sx q[2];
rz(-1.8554167) q[2];
sx q[2];
rz(2.9993338) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.576697) q[1];
sx q[1];
rz(-1.7045867) q[1];
sx q[1];
rz(-2.1810075) q[1];
rz(-pi) q[2];
rz(-1.9777771) q[3];
sx q[3];
rz(-2.3325936) q[3];
sx q[3];
rz(-0.33363261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(-1.3403085) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(-0.006342412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9316677) q[0];
sx q[0];
rz(-1.7010265) q[0];
sx q[0];
rz(1.3757214) q[0];
x q[1];
rz(-0.96887529) q[2];
sx q[2];
rz(-2.0191779) q[2];
sx q[2];
rz(0.3647764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1576707) q[1];
sx q[1];
rz(-1.6323286) q[1];
sx q[1];
rz(2.8308949) q[1];
rz(0.38624318) q[3];
sx q[3];
rz(-2.5442903) q[3];
sx q[3];
rz(2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(-3.045851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11287963) q[0];
sx q[0];
rz(-2.9071147) q[0];
sx q[0];
rz(2.2857091) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6355866) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(0.72088748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84903753) q[1];
sx q[1];
rz(-1.109086) q[1];
sx q[1];
rz(-0.95878102) q[1];
rz(-pi) q[2];
rz(1.6060353) q[3];
sx q[3];
rz(-2.1217151) q[3];
sx q[3];
rz(2.2743724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-2.3068008) q[2];
rz(-0.21162027) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-2.1077164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0945275) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(2.4497776) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.526555) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(-1.9499792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2305206) q[1];
sx q[1];
rz(-1.5057179) q[1];
sx q[1];
rz(-1.7622403) q[1];
rz(-pi) q[2];
rz(0.86446188) q[3];
sx q[3];
rz(-2.7221788) q[3];
sx q[3];
rz(-0.97233397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(1.7948077) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(0.41473266) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-0.57410747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60526472) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(1.7954134) q[0];
x q[1];
rz(1.5933883) q[2];
sx q[2];
rz(-0.85561692) q[2];
sx q[2];
rz(1.1720282) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9088604) q[1];
sx q[1];
rz(-0.4545916) q[1];
sx q[1];
rz(-0.67635398) q[1];
rz(1.8168713) q[3];
sx q[3];
rz(-0.68021357) q[3];
sx q[3];
rz(2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85270143) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(1.627702) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-0.69818991) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-2.4093157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66757827) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(1.408512) q[0];
x q[1];
rz(-0.90768355) q[2];
sx q[2];
rz(-0.92539061) q[2];
sx q[2];
rz(1.4484608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1440891) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(1.287582) q[1];
rz(-1.7888277) q[3];
sx q[3];
rz(-1.0112959) q[3];
sx q[3];
rz(0.41527173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5891588) q[0];
sx q[0];
rz(-1.7434412) q[0];
sx q[0];
rz(-2.5347559) q[0];
rz(-2.1595575) q[2];
sx q[2];
rz(-1.2274449) q[2];
sx q[2];
rz(2.8729168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8416482) q[1];
sx q[1];
rz(-0.83924676) q[1];
sx q[1];
rz(-2.4323835) q[1];
x q[2];
rz(-0.50290147) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(2.6754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(-2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-0.36639211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057356) q[0];
sx q[0];
rz(-0.85548399) q[0];
sx q[0];
rz(2.4321796) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.512152) q[2];
sx q[2];
rz(-1.2120314) q[2];
sx q[2];
rz(1.4251054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3185127) q[1];
sx q[1];
rz(-1.4050583) q[1];
sx q[1];
rz(-0.70648273) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0189692) q[3];
sx q[3];
rz(-1.5295715) q[3];
sx q[3];
rz(0.036389694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1239132) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-0.79137897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489483) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(2.2788458) q[0];
rz(-0.4653761) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(-0.79681764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8225122) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(0.49088571) q[1];
x q[2];
rz(0.16718849) q[3];
sx q[3];
rz(-1.0382004) q[3];
sx q[3];
rz(0.54824588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.7609319) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1191694) q[0];
sx q[0];
rz(-1.4403617) q[0];
sx q[0];
rz(-0.93956691) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1157687) q[2];
sx q[2];
rz(-2.0061473) q[2];
sx q[2];
rz(0.045407427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7541411) q[1];
sx q[1];
rz(-1.8879461) q[1];
sx q[1];
rz(-2.0003358) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4139666) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(-2.6550456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(-1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(2.8019194) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(1.2119157) q[3];
sx q[3];
rz(-0.59960312) q[3];
sx q[3];
rz(-3.0740769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
