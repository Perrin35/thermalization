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
rz(-1.834637) q[0];
sx q[0];
rz(-0.87943465) q[0];
sx q[0];
rz(2.6479794) q[0];
rz(1.8618795) q[1];
sx q[1];
rz(-0.6266098) q[1];
sx q[1];
rz(-1.6630747) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8956229) q[0];
sx q[0];
rz(-0.9446656) q[0];
sx q[0];
rz(-2.5139721) q[0];
rz(1.8952719) q[2];
sx q[2];
rz(-2.4518584) q[2];
sx q[2];
rz(1.5943499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6329816) q[1];
sx q[1];
rz(-2.2618083) q[1];
sx q[1];
rz(1.2862842) q[1];
rz(-0.36840393) q[3];
sx q[3];
rz(-2.1737636) q[3];
sx q[3];
rz(-1.4908229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0066068) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(-1.1084278) q[2];
rz(-1.3125575) q[3];
sx q[3];
rz(-0.24459845) q[3];
sx q[3];
rz(2.826622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.7600064) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(-3.1257358) q[0];
rz(2.1511757) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(2.2784065) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28672934) q[0];
sx q[0];
rz(-0.83896381) q[0];
sx q[0];
rz(1.4409723) q[0];
rz(-pi) q[1];
x q[1];
rz(1.241774) q[2];
sx q[2];
rz(-1.8890427) q[2];
sx q[2];
rz(-1.1877738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6241315) q[1];
sx q[1];
rz(-1.7449813) q[1];
sx q[1];
rz(-0.49887533) q[1];
rz(-pi) q[2];
rz(-2.2926278) q[3];
sx q[3];
rz(-1.6833653) q[3];
sx q[3];
rz(-1.4813042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.70011675) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(2.6727943) q[2];
rz(-2.4539454) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-2.8414753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1402635) q[0];
sx q[0];
rz(-2.2183473) q[0];
sx q[0];
rz(0.73626751) q[0];
rz(0.34321347) q[1];
sx q[1];
rz(-2.2540269) q[1];
sx q[1];
rz(-2.9578178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7385173) q[0];
sx q[0];
rz(-0.74664298) q[0];
sx q[0];
rz(-0.097642032) q[0];
x q[1];
rz(0.69664069) q[2];
sx q[2];
rz(-1.9844779) q[2];
sx q[2];
rz(2.6282642) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.80189542) q[1];
sx q[1];
rz(-1.8690799) q[1];
sx q[1];
rz(-2.4043179) q[1];
rz(-pi) q[2];
rz(2.1322429) q[3];
sx q[3];
rz(-0.61547503) q[3];
sx q[3];
rz(-1.4765679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78202128) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(2.6256631) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-0.55086946) q[3];
sx q[3];
rz(1.1512604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9786523) q[0];
sx q[0];
rz(-1.2115703) q[0];
sx q[0];
rz(-2.6287855) q[0];
rz(-0.45979744) q[1];
sx q[1];
rz(-1.5492946) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490246) q[0];
sx q[0];
rz(-1.4563174) q[0];
sx q[0];
rz(-3.0106198) q[0];
rz(-pi) q[1];
rz(-0.41641367) q[2];
sx q[2];
rz(-1.2730619) q[2];
sx q[2];
rz(-0.99978757) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0280007) q[1];
sx q[1];
rz(-1.290844) q[1];
sx q[1];
rz(1.9242313) q[1];
x q[2];
rz(2.7971036) q[3];
sx q[3];
rz(-1.7055344) q[3];
sx q[3];
rz(-1.3400934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1226471) q[2];
sx q[2];
rz(-0.15572369) q[2];
sx q[2];
rz(0.32749185) q[2];
rz(-2.859595) q[3];
sx q[3];
rz(-1.7433386) q[3];
sx q[3];
rz(1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9649488) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(0.97736812) q[0];
rz(2.7727959) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(1.7288953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4810886) q[0];
sx q[0];
rz(-1.3810643) q[0];
sx q[0];
rz(-2.0697748) q[0];
x q[1];
rz(1.5951889) q[2];
sx q[2];
rz(-2.6246723) q[2];
sx q[2];
rz(-2.2395397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2569132) q[1];
sx q[1];
rz(-1.8533145) q[1];
sx q[1];
rz(0.032508548) q[1];
x q[2];
rz(0.47980185) q[3];
sx q[3];
rz(-2.5069624) q[3];
sx q[3];
rz(3.0632085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29350027) q[2];
sx q[2];
rz(-2.3257181) q[2];
sx q[2];
rz(-0.14981848) q[2];
rz(0.39854974) q[3];
sx q[3];
rz(-1.6179251) q[3];
sx q[3];
rz(2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3715816) q[0];
sx q[0];
rz(-3.0122029) q[0];
sx q[0];
rz(0.55321252) q[0];
rz(2.5999056) q[1];
sx q[1];
rz(-2.0952416) q[1];
sx q[1];
rz(-0.4479301) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174367) q[0];
sx q[0];
rz(-1.7154335) q[0];
sx q[0];
rz(-2.7516624) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7121353) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(-2.7384659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74400157) q[1];
sx q[1];
rz(-0.41644704) q[1];
sx q[1];
rz(1.9201502) q[1];
rz(-pi) q[2];
rz(0.33282307) q[3];
sx q[3];
rz(-1.5170005) q[3];
sx q[3];
rz(1.7032361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64122671) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(2.3424303) q[2];
rz(-2.7568119) q[3];
sx q[3];
rz(-2.81541) q[3];
sx q[3];
rz(-1.0001812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36050972) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(2.5497896) q[0];
rz(-0.38161713) q[1];
sx q[1];
rz(-2.7615669) q[1];
sx q[1];
rz(-2.9891678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7397933) q[0];
sx q[0];
rz(-1.2619979) q[0];
sx q[0];
rz(-1.0303506) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.746114) q[2];
sx q[2];
rz(-1.0994871) q[2];
sx q[2];
rz(2.5702554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85411863) q[1];
sx q[1];
rz(-1.9175314) q[1];
sx q[1];
rz(-1.5513913) q[1];
rz(-2.3197738) q[3];
sx q[3];
rz(-2.3082409) q[3];
sx q[3];
rz(0.59274537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3064208) q[2];
sx q[2];
rz(-0.99205899) q[2];
sx q[2];
rz(1.6888118) q[2];
rz(-0.49550223) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(-0.56408322) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38967663) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(2.9801242) q[0];
rz(-0.031919315) q[1];
sx q[1];
rz(-2.4999764) q[1];
sx q[1];
rz(1.8643103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.204435) q[0];
sx q[0];
rz(-1.8291408) q[0];
sx q[0];
rz(0.72465557) q[0];
rz(-0.23060282) q[2];
sx q[2];
rz(-1.1586894) q[2];
sx q[2];
rz(-1.3889564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2252075) q[1];
sx q[1];
rz(-2.8182833) q[1];
sx q[1];
rz(-2.7955452) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6971436) q[3];
sx q[3];
rz(-2.845361) q[3];
sx q[3];
rz(-2.031523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6892467) q[2];
sx q[2];
rz(-0.9114868) q[2];
sx q[2];
rz(-2.7455043) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-2.5466099) q[3];
sx q[3];
rz(-2.2649435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24092995) q[0];
sx q[0];
rz(-3.1234968) q[0];
sx q[0];
rz(0.15116365) q[0];
rz(-0.97686544) q[1];
sx q[1];
rz(-2.7604389) q[1];
sx q[1];
rz(2.3968598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1158651) q[0];
sx q[0];
rz(-2.7926499) q[0];
sx q[0];
rz(2.2141488) q[0];
x q[1];
rz(-1.9836203) q[2];
sx q[2];
rz(-1.4399488) q[2];
sx q[2];
rz(-1.3589588) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7661644) q[1];
sx q[1];
rz(-0.64830983) q[1];
sx q[1];
rz(0.0023771087) q[1];
x q[2];
rz(0.35208296) q[3];
sx q[3];
rz(-1.3439281) q[3];
sx q[3];
rz(2.9405103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(-0.94804478) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.4761997) q[3];
sx q[3];
rz(1.074056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0483765) q[0];
sx q[0];
rz(-0.26234782) q[0];
sx q[0];
rz(-0.23396215) q[0];
rz(1.1853064) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(1.8918461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5025185) q[0];
sx q[0];
rz(-2.5022129) q[0];
sx q[0];
rz(2.6431497) q[0];
rz(-pi) q[1];
rz(-0.62057497) q[2];
sx q[2];
rz(-1.4016376) q[2];
sx q[2];
rz(-2.754654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054259006) q[1];
sx q[1];
rz(-1.0339843) q[1];
sx q[1];
rz(0.76294239) q[1];
x q[2];
rz(-0.21548157) q[3];
sx q[3];
rz(-0.69284791) q[3];
sx q[3];
rz(0.97164916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.027111) q[2];
sx q[2];
rz(-1.1610843) q[2];
sx q[2];
rz(-1.4410045) q[2];
rz(-0.13127413) q[3];
sx q[3];
rz(-2.6797397) q[3];
sx q[3];
rz(2.89768) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048653614) q[0];
sx q[0];
rz(-2.1760512) q[0];
sx q[0];
rz(-2.0078134) q[0];
rz(-2.1009905) q[1];
sx q[1];
rz(-0.84747172) q[1];
sx q[1];
rz(-0.52451959) q[1];
rz(-2.8362989) q[2];
sx q[2];
rz(-2.1758781) q[2];
sx q[2];
rz(0.97347409) q[2];
rz(2.5083422) q[3];
sx q[3];
rz(-1.4489104) q[3];
sx q[3];
rz(-0.62319402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
