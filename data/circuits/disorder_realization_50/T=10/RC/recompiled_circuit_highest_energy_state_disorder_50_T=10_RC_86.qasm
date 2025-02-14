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
rz(0.4515689) q[0];
sx q[0];
rz(-1.7234001) q[0];
sx q[0];
rz(-2.7127142) q[0];
rz(-0.14313993) q[1];
sx q[1];
rz(-2.0506471) q[1];
sx q[1];
rz(0.080088869) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1983436) q[0];
sx q[0];
rz(-0.08028537) q[0];
sx q[0];
rz(2.9757418) q[0];
rz(3.0386904) q[2];
sx q[2];
rz(-1.2116625) q[2];
sx q[2];
rz(0.54190901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0804008) q[1];
sx q[1];
rz(-0.75729232) q[1];
sx q[1];
rz(1.1822834) q[1];
rz(-pi) q[2];
rz(-0.24899068) q[3];
sx q[3];
rz(-2.5496229) q[3];
sx q[3];
rz(-2.0651061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0668138) q[2];
sx q[2];
rz(-1.0513693) q[2];
sx q[2];
rz(1.4220413) q[2];
rz(1.413013) q[3];
sx q[3];
rz(-1.6001817) q[3];
sx q[3];
rz(-1.448267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-0.3438147) q[1];
sx q[1];
rz(2.3894892) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1307169) q[0];
sx q[0];
rz(-2.1701522) q[0];
sx q[0];
rz(-0.54979558) q[0];
x q[1];
rz(-3.0235748) q[2];
sx q[2];
rz(-0.22023295) q[2];
sx q[2];
rz(1.756246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.960452) q[1];
sx q[1];
rz(-1.6799095) q[1];
sx q[1];
rz(-1.0427598) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26176287) q[3];
sx q[3];
rz(-1.7883375) q[3];
sx q[3];
rz(-2.1037773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0636474) q[2];
sx q[2];
rz(-2.2293978) q[2];
sx q[2];
rz(-2.4845541) q[2];
rz(-0.71197236) q[3];
sx q[3];
rz(-2.4896121) q[3];
sx q[3];
rz(-0.74487346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7065358) q[0];
sx q[0];
rz(-2.5223795) q[0];
sx q[0];
rz(-1.0414498) q[0];
rz(2.7667747) q[1];
sx q[1];
rz(-1.5033787) q[1];
sx q[1];
rz(-0.55694881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3303012) q[0];
sx q[0];
rz(-2.3200703) q[0];
sx q[0];
rz(-0.4901476) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0224503) q[2];
sx q[2];
rz(-1.5030393) q[2];
sx q[2];
rz(-0.51167831) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3279325) q[1];
sx q[1];
rz(-0.90444649) q[1];
sx q[1];
rz(1.440669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65105533) q[3];
sx q[3];
rz(-1.7604894) q[3];
sx q[3];
rz(0.22471353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53691429) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(-2.7064145) q[2];
rz(-0.52470454) q[3];
sx q[3];
rz(-0.33354959) q[3];
sx q[3];
rz(-0.13478002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0404469) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(-1.6271628) q[0];
rz(-2.0928404) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(1.4357766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7266975) q[0];
sx q[0];
rz(-1.613662) q[0];
sx q[0];
rz(2.9662762) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5975512) q[2];
sx q[2];
rz(-0.61464192) q[2];
sx q[2];
rz(2.7934472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1757112) q[1];
sx q[1];
rz(-2.615965) q[1];
sx q[1];
rz(-1.6986584) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7955608) q[3];
sx q[3];
rz(-1.1308987) q[3];
sx q[3];
rz(-1.4113246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7413896) q[2];
sx q[2];
rz(-1.3578537) q[2];
sx q[2];
rz(-3.1164361) q[2];
rz(-1.4870421) q[3];
sx q[3];
rz(-2.7436723) q[3];
sx q[3];
rz(0.81620836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095489278) q[0];
sx q[0];
rz(-0.64952055) q[0];
sx q[0];
rz(-2.982614) q[0];
rz(2.036463) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(0.22430688) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.255757) q[0];
sx q[0];
rz(-2.1191264) q[0];
sx q[0];
rz(1.095039) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1087628) q[2];
sx q[2];
rz(-1.0107892) q[2];
sx q[2];
rz(-0.0074040961) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.080706656) q[1];
sx q[1];
rz(-1.4554441) q[1];
sx q[1];
rz(1.0447698) q[1];
rz(-pi) q[2];
rz(-0.57711011) q[3];
sx q[3];
rz(-1.6770067) q[3];
sx q[3];
rz(-0.67536992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20532456) q[2];
sx q[2];
rz(-1.779413) q[2];
sx q[2];
rz(2.3929907) q[2];
rz(-0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(-2.7301679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1614138) q[0];
sx q[0];
rz(-0.36853376) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(1.5018564) q[1];
sx q[1];
rz(-1.390099) q[1];
sx q[1];
rz(0.23434848) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4880331) q[0];
sx q[0];
rz(-2.3071978) q[0];
sx q[0];
rz(2.5279765) q[0];
rz(-0.83291556) q[2];
sx q[2];
rz(-2.3489522) q[2];
sx q[2];
rz(-1.8120676) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6987353) q[1];
sx q[1];
rz(-0.1695098) q[1];
sx q[1];
rz(0.28556602) q[1];
rz(-pi) q[2];
rz(-1.1612915) q[3];
sx q[3];
rz(-2.9886768) q[3];
sx q[3];
rz(2.5482938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.948287) q[2];
sx q[2];
rz(-2.327658) q[2];
sx q[2];
rz(-1.8602271) q[2];
rz(-0.60503259) q[3];
sx q[3];
rz(-2.0421959) q[3];
sx q[3];
rz(-1.9090778) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.879377) q[0];
sx q[0];
rz(-2.5285517) q[0];
rz(-2.4609861) q[1];
sx q[1];
rz(-1.1111518) q[1];
sx q[1];
rz(-1.3253164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82241908) q[0];
sx q[0];
rz(-1.028601) q[0];
sx q[0];
rz(1.2577357) q[0];
rz(-0.36779495) q[2];
sx q[2];
rz(-0.22794524) q[2];
sx q[2];
rz(1.6942555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0874859) q[1];
sx q[1];
rz(-1.0102059) q[1];
sx q[1];
rz(-1.4212379) q[1];
x q[2];
rz(1.9716827) q[3];
sx q[3];
rz(-1.0260149) q[3];
sx q[3];
rz(1.137103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1943872) q[2];
sx q[2];
rz(-1.2014061) q[2];
sx q[2];
rz(2.3335333) q[2];
rz(-0.64588532) q[3];
sx q[3];
rz(-0.26791993) q[3];
sx q[3];
rz(-2.8382235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5658257) q[0];
sx q[0];
rz(-0.72661) q[0];
sx q[0];
rz(0.44276825) q[0];
rz(-2.5495461) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(-2.9659042) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239227) q[0];
sx q[0];
rz(-0.45638535) q[0];
sx q[0];
rz(0.66500591) q[0];
rz(0.4515516) q[2];
sx q[2];
rz(-1.6066666) q[2];
sx q[2];
rz(-2.4810239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1877039) q[1];
sx q[1];
rz(-0.64262455) q[1];
sx q[1];
rz(0.072037176) q[1];
x q[2];
rz(-2.2288228) q[3];
sx q[3];
rz(-2.4173749) q[3];
sx q[3];
rz(-0.66373827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5438133) q[2];
sx q[2];
rz(-0.74395776) q[2];
sx q[2];
rz(0.90710863) q[2];
rz(-2.234327) q[3];
sx q[3];
rz(-1.7472569) q[3];
sx q[3];
rz(0.11769122) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1842781) q[0];
sx q[0];
rz(-2.2397581) q[0];
sx q[0];
rz(-0.19004518) q[0];
rz(1.5589145) q[1];
sx q[1];
rz(-1.546944) q[1];
sx q[1];
rz(-1.83439) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652649) q[0];
sx q[0];
rz(-1.4272262) q[0];
sx q[0];
rz(1.2031892) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87086579) q[2];
sx q[2];
rz(-1.3076837) q[2];
sx q[2];
rz(1.7990636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41220442) q[1];
sx q[1];
rz(-1.6784759) q[1];
sx q[1];
rz(3.1086224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7070624) q[3];
sx q[3];
rz(-1.4659681) q[3];
sx q[3];
rz(-1.9275582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0731395) q[2];
sx q[2];
rz(-1.0871202) q[2];
sx q[2];
rz(-0.20115176) q[2];
rz(2.1189832) q[3];
sx q[3];
rz(-0.68251959) q[3];
sx q[3];
rz(1.539591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96711838) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(0.87646595) q[0];
rz(-2.5190952) q[1];
sx q[1];
rz(-0.21177706) q[1];
sx q[1];
rz(-0.34271398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7066318) q[0];
sx q[0];
rz(-1.5049269) q[0];
sx q[0];
rz(-2.1996798) q[0];
rz(-pi) q[1];
rz(-1.9075526) q[2];
sx q[2];
rz(-1.4554664) q[2];
sx q[2];
rz(-0.28240909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5805768) q[1];
sx q[1];
rz(-2.3388712) q[1];
sx q[1];
rz(-0.64508857) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38767918) q[3];
sx q[3];
rz(-1.8053384) q[3];
sx q[3];
rz(-2.925433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1891302) q[2];
sx q[2];
rz(-1.4418944) q[2];
sx q[2];
rz(0.44378898) q[2];
rz(-1.1187547) q[3];
sx q[3];
rz(-1.5464562) q[3];
sx q[3];
rz(2.5877623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(1.2832058) q[0];
sx q[0];
rz(-1.5844185) q[0];
sx q[0];
rz(1.6208741) q[0];
rz(0.063477909) q[1];
sx q[1];
rz(-1.3451481) q[1];
sx q[1];
rz(1.4048911) q[1];
rz(-0.061894682) q[2];
sx q[2];
rz(-1.2409004) q[2];
sx q[2];
rz(0.093766669) q[2];
rz(-2.8982481) q[3];
sx q[3];
rz(-2.1235597) q[3];
sx q[3];
rz(-1.7490993) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
