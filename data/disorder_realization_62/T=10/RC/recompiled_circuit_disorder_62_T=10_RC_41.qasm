OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21270574) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(-1.5972208) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3637949) q[2];
sx q[2];
rz(-1.8282837) q[2];
sx q[2];
rz(2.4553026) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.692688) q[1];
sx q[1];
rz(-1.9069888) q[1];
sx q[1];
rz(2.8640139) q[1];
x q[2];
rz(1.9840368) q[3];
sx q[3];
rz(-0.26502702) q[3];
sx q[3];
rz(0.36039263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5548271) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.9155496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37218371) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(0.43222897) q[0];
rz(-pi) q[1];
rz(-0.77483564) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(1.6006084) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.023149816) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(1.5688194) q[1];
x q[2];
rz(-1.5561043) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79919672) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.249041) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-2.1121315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(2.1973781) q[0];
x q[1];
rz(-0.37000334) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(1.5918658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50476915) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.7117281) q[1];
rz(-pi) q[2];
rz(-0.13266487) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(2.6877407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-2.6339445) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.078331) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(-0.22380933) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1298725) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(-0.59414547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2559291) q[1];
sx q[1];
rz(-2.4293578) q[1];
sx q[1];
rz(2.7182012) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58873119) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(2.246726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2762404) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(-2.4549237) q[0];
rz(-pi) q[1];
rz(2.1331482) q[2];
sx q[2];
rz(-2.214606) q[2];
sx q[2];
rz(1.84562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85256165) q[1];
sx q[1];
rz(-2.3173996) q[1];
sx q[1];
rz(2.186071) q[1];
x q[2];
rz(-2.8551293) q[3];
sx q[3];
rz(-1.5781919) q[3];
sx q[3];
rz(1.8252107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.557945) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.7640132) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-2.6766052) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2842448) q[0];
sx q[0];
rz(-1.9042958) q[0];
sx q[0];
rz(-1.1112945) q[0];
rz(-pi) q[1];
rz(1.497252) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(-1.7433012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28593674) q[1];
sx q[1];
rz(-0.93614139) q[1];
sx q[1];
rz(0.059307701) q[1];
rz(-pi) q[2];
rz(-0.24352169) q[3];
sx q[3];
rz(-0.72391073) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(-0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.4200462) q[0];
rz(-3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50990803) q[0];
sx q[0];
rz(-1.317306) q[0];
sx q[0];
rz(1.7476728) q[0];
x q[1];
rz(2.4620352) q[2];
sx q[2];
rz(-0.62629269) q[2];
sx q[2];
rz(1.2223513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8487726) q[1];
sx q[1];
rz(-2.9584868) q[1];
sx q[1];
rz(1.8715026) q[1];
rz(-2.6050623) q[3];
sx q[3];
rz(-1.1324258) q[3];
sx q[3];
rz(2.8560672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(-0.32564751) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(0.24857323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6876858) q[0];
sx q[0];
rz(-1.0605863) q[0];
sx q[0];
rz(1.8574255) q[0];
x q[1];
rz(-2.3209004) q[2];
sx q[2];
rz(-0.65260115) q[2];
sx q[2];
rz(1.5479659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2182902) q[1];
sx q[1];
rz(-0.81070886) q[1];
sx q[1];
rz(-0.87358012) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.076898889) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(0.83991915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.823002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628172) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(2.0128987) q[0];
x q[1];
rz(-0.96005) q[2];
sx q[2];
rz(-1.0078444) q[2];
sx q[2];
rz(-0.27573953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1773771) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(1.3436951) q[1];
rz(-pi) q[2];
rz(0.41140822) q[3];
sx q[3];
rz(-2.3477926) q[3];
sx q[3];
rz(-2.0421162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.3783962) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.261196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6529918) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(2.7306042) q[0];
x q[1];
rz(-3.10896) q[2];
sx q[2];
rz(-0.59851461) q[2];
sx q[2];
rz(1.061071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2640472) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(-1.4977786) q[1];
rz(1.4114755) q[3];
sx q[3];
rz(-2.2485002) q[3];
sx q[3];
rz(0.89980984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(0.56979257) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-3.1153395) q[2];
sx q[2];
rz(-2.142341) q[2];
sx q[2];
rz(-0.061948902) q[2];
rz(1.431987) q[3];
sx q[3];
rz(-2.1574253) q[3];
sx q[3];
rz(0.16711259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
