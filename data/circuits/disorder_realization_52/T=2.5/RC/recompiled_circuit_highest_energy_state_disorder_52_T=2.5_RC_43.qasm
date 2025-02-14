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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(3.1290913) q[1];
sx q[1];
rz(-1.8379509) q[1];
sx q[1];
rz(-1.5703896) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9725295) q[0];
sx q[0];
rz(-0.95989812) q[0];
sx q[0];
rz(1.1283895) q[0];
rz(-3.0742253) q[2];
sx q[2];
rz(-1.9695896) q[2];
sx q[2];
rz(0.68902868) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.058152288) q[1];
sx q[1];
rz(-1.5136763) q[1];
sx q[1];
rz(-1.2099464) q[1];
x q[2];
rz(3.0812991) q[3];
sx q[3];
rz(-1.8795085) q[3];
sx q[3];
rz(-1.2293881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5300753) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(0.44069904) q[2];
rz(0.079744451) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(2.0174446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1752862) q[0];
sx q[0];
rz(-3.0847302) q[0];
sx q[0];
rz(2.9735907) q[0];
rz(3.1213144) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(-1.6049989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0609866) q[0];
sx q[0];
rz(-1.3362552) q[0];
sx q[0];
rz(-1.7936262) q[0];
rz(1.3290908) q[2];
sx q[2];
rz(-3.1311718) q[2];
sx q[2];
rz(-1.3543606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3602996) q[1];
sx q[1];
rz(-1.5661245) q[1];
sx q[1];
rz(0.0010990573) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5617989) q[3];
sx q[3];
rz(-1.5285042) q[3];
sx q[3];
rz(2.472714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34540471) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(1.7339285) q[2];
rz(2.0965072) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(-2.869587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8318091) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(0.56104863) q[0];
rz(-0.27684119) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(-1.3078088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40807276) q[0];
sx q[0];
rz(-1.2740943) q[0];
sx q[0];
rz(-3.0987334) q[0];
rz(-0.00066999992) q[2];
sx q[2];
rz(-3.1342034) q[2];
sx q[2];
rz(2.0787897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6456994) q[1];
sx q[1];
rz(-0.99826854) q[1];
sx q[1];
rz(3.0719724) q[1];
rz(-2.6827319) q[3];
sx q[3];
rz(-2.1956177) q[3];
sx q[3];
rz(0.030908728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7967367) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(0.5893839) q[2];
rz(-1.2423337) q[3];
sx q[3];
rz(-3.1292249) q[3];
sx q[3];
rz(1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(3.1347443) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(-1.3486598) q[0];
rz(-3.1349365) q[1];
sx q[1];
rz(-1.8204047) q[1];
sx q[1];
rz(3.1080918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703465) q[0];
sx q[0];
rz(-1.7435562) q[0];
sx q[0];
rz(-2.2281649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.037390402) q[2];
sx q[2];
rz(-3.0218736) q[2];
sx q[2];
rz(2.7604288) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9814947) q[1];
sx q[1];
rz(-1.3019053) q[1];
sx q[1];
rz(-3.1310496) q[1];
rz(-pi) q[2];
rz(0.93422814) q[3];
sx q[3];
rz(-0.2033955) q[3];
sx q[3];
rz(-2.0980841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.032430705) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(3.0884009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.771516) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(-2.694743) q[0];
rz(-0.18613786) q[1];
sx q[1];
rz(-3.0797854) q[1];
sx q[1];
rz(-1.4131379) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9139149) q[0];
sx q[0];
rz(-2.9083038) q[0];
sx q[0];
rz(1.4248217) q[0];
x q[1];
rz(-0.21633001) q[2];
sx q[2];
rz(-1.1573175) q[2];
sx q[2];
rz(-1.0461996) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.062033) q[1];
sx q[1];
rz(-1.6158982) q[1];
sx q[1];
rz(1.5212039) q[1];
rz(-pi) q[2];
rz(0.22362169) q[3];
sx q[3];
rz(-1.2960172) q[3];
sx q[3];
rz(1.4449121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83429217) q[2];
sx q[2];
rz(-1.5520381) q[2];
sx q[2];
rz(0.49868047) q[2];
rz(0.57791609) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(-2.5448866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805098) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(0.34641308) q[0];
rz(-0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(-2.3880889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44418535) q[0];
sx q[0];
rz(-1.7427398) q[0];
sx q[0];
rz(-0.19449046) q[0];
x q[1];
rz(3.0664938) q[2];
sx q[2];
rz(-1.4595928) q[2];
sx q[2];
rz(2.5203343) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0776163) q[1];
sx q[1];
rz(-1.156563) q[1];
sx q[1];
rz(-2.2760681) q[1];
rz(-pi) q[2];
rz(3.0493204) q[3];
sx q[3];
rz(-1.408395) q[3];
sx q[3];
rz(-0.55264651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5715013) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(1.6174779) q[2];
rz(3.0213455) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(-2.6038468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8728802) q[0];
sx q[0];
rz(-0.92395067) q[0];
sx q[0];
rz(2.9203316) q[0];
rz(1.4606754) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(3.0642919) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.59931) q[0];
sx q[0];
rz(-3.0998383) q[0];
sx q[0];
rz(1.6058654) q[0];
rz(-pi) q[1];
rz(-1.5796698) q[2];
sx q[2];
rz(-1.5660145) q[2];
sx q[2];
rz(1.7824796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0127924) q[1];
sx q[1];
rz(-0.18928738) q[1];
sx q[1];
rz(-0.38893338) q[1];
rz(-pi) q[2];
rz(-3.0139913) q[3];
sx q[3];
rz(-1.7831047) q[3];
sx q[3];
rz(-1.177976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7920502) q[2];
sx q[2];
rz(-3.1304066) q[2];
sx q[2];
rz(-0.95996094) q[2];
rz(2.8121484) q[3];
sx q[3];
rz(-0.0080778413) q[3];
sx q[3];
rz(-2.2913057) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.195381) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(3.0323113) q[0];
rz(2.7682313) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(-1.9111309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7007992) q[0];
sx q[0];
rz(-0.95958086) q[0];
sx q[0];
rz(0.81328765) q[0];
x q[1];
rz(2.9443342) q[2];
sx q[2];
rz(-1.3840809) q[2];
sx q[2];
rz(-0.017053617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2956063) q[1];
sx q[1];
rz(-1.5593658) q[1];
sx q[1];
rz(1.5039526) q[1];
rz(-pi) q[2];
rz(-2.8318895) q[3];
sx q[3];
rz(-0.94325698) q[3];
sx q[3];
rz(0.58615875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5750778) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(1.3285948) q[2];
rz(1.7447507) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0777271) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(-2.5685837) q[0];
rz(2.833448) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(2.134197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44237374) q[0];
sx q[0];
rz(-1.6773407) q[0];
sx q[0];
rz(-0.005597896) q[0];
rz(1.7083659) q[2];
sx q[2];
rz(-1.614288) q[2];
sx q[2];
rz(1.5303591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3904265) q[1];
sx q[1];
rz(-1.6114283) q[1];
sx q[1];
rz(-1.6540113) q[1];
x q[2];
rz(1.151772) q[3];
sx q[3];
rz(-0.033009987) q[3];
sx q[3];
rz(2.9275683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.31743) q[2];
sx q[2];
rz(-0.62676668) q[2];
sx q[2];
rz(2.7516348) q[2];
rz(-3.0725078) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(2.8100584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.1214509) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(-2.2700229) q[1];
sx q[1];
rz(-1.8337245) q[1];
sx q[1];
rz(-1.647324) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6914929) q[0];
sx q[0];
rz(-2.1039824) q[0];
sx q[0];
rz(1.2144809) q[0];
rz(3.0790555) q[2];
sx q[2];
rz(-0.61574575) q[2];
sx q[2];
rz(3.1341022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20032665) q[1];
sx q[1];
rz(-1.2709874) q[1];
sx q[1];
rz(0.34383641) q[1];
x q[2];
rz(-1.5288197) q[3];
sx q[3];
rz(-1.4577565) q[3];
sx q[3];
rz(0.50769317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-0.031878397) q[2];
rz(0.77830642) q[3];
sx q[3];
rz(-3.1347771) q[3];
sx q[3];
rz(-2.8460898) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185709) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(0.12693916) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(-3.1359966) q[2];
sx q[2];
rz(-1.4311309) q[2];
sx q[2];
rz(-2.8972814) q[2];
rz(-2.245458) q[3];
sx q[3];
rz(-1.6150766) q[3];
sx q[3];
rz(-1.0684039) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
