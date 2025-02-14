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
rz(0.69819063) q[0];
sx q[0];
rz(-0.31261045) q[0];
sx q[0];
rz(2.1460549) q[0];
rz(-8.3182316) q[1];
sx q[1];
rz(3.9099524) q[1];
sx q[1];
rz(16.04019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8790008) q[0];
sx q[0];
rz(-0.79918062) q[0];
sx q[0];
rz(-0.57882092) q[0];
x q[1];
rz(-1.840749) q[2];
sx q[2];
rz(-0.47856646) q[2];
sx q[2];
rz(1.9437444) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.600647) q[1];
sx q[1];
rz(-1.9917734) q[1];
sx q[1];
rz(1.0599815) q[1];
x q[2];
rz(-2.6237021) q[3];
sx q[3];
rz(-2.0841408) q[3];
sx q[3];
rz(2.1012517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.535061) q[2];
sx q[2];
rz(-2.1273095) q[2];
sx q[2];
rz(-0.052058546) q[2];
rz(-2.0590797) q[3];
sx q[3];
rz(-0.94729298) q[3];
sx q[3];
rz(-1.1516217) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5169736) q[0];
sx q[0];
rz(-3.1095412) q[0];
sx q[0];
rz(0.33729851) q[0];
rz(2.9889122) q[1];
sx q[1];
rz(-0.76714194) q[1];
sx q[1];
rz(-1.5229567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6628159) q[0];
sx q[0];
rz(-2.1838476) q[0];
sx q[0];
rz(1.3663892) q[0];
rz(-1.3682034) q[2];
sx q[2];
rz(-0.94300044) q[2];
sx q[2];
rz(-1.0204878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2238521) q[1];
sx q[1];
rz(-1.187993) q[1];
sx q[1];
rz(0.21767016) q[1];
rz(-pi) q[2];
rz(1.7256152) q[3];
sx q[3];
rz(-0.972675) q[3];
sx q[3];
rz(3.104676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0308257) q[2];
sx q[2];
rz(-0.36588565) q[2];
sx q[2];
rz(-0.3863253) q[2];
rz(1.857916) q[3];
sx q[3];
rz(-1.9648896) q[3];
sx q[3];
rz(1.9234689) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074742643) q[0];
sx q[0];
rz(-2.061494) q[0];
sx q[0];
rz(-1.0190438) q[0];
rz(-2.3661803) q[1];
sx q[1];
rz(-1.8653899) q[1];
sx q[1];
rz(-0.48616854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4751525) q[0];
sx q[0];
rz(-1.5678582) q[0];
sx q[0];
rz(3.1102528) q[0];
rz(-pi) q[1];
rz(2.6576651) q[2];
sx q[2];
rz(-0.84887767) q[2];
sx q[2];
rz(-2.1579822) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8067161) q[1];
sx q[1];
rz(-1.7893605) q[1];
sx q[1];
rz(-2.6456635) q[1];
x q[2];
rz(-0.73681946) q[3];
sx q[3];
rz(-2.7180258) q[3];
sx q[3];
rz(1.119348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7458618) q[2];
sx q[2];
rz(-1.9755325) q[2];
sx q[2];
rz(-2.3365848) q[2];
rz(2.4896367) q[3];
sx q[3];
rz(-2.0396353) q[3];
sx q[3];
rz(-2.2074047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790799) q[0];
sx q[0];
rz(-1.3211687) q[0];
sx q[0];
rz(1.3385734) q[0];
rz(1.592912) q[1];
sx q[1];
rz(-2.4376696) q[1];
sx q[1];
rz(0.42565638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5558472) q[0];
sx q[0];
rz(-0.092215538) q[0];
sx q[0];
rz(-0.088949843) q[0];
rz(-pi) q[1];
rz(-1.2629444) q[2];
sx q[2];
rz(-1.5378219) q[2];
sx q[2];
rz(-2.3127382) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8989694) q[1];
sx q[1];
rz(-2.2950566) q[1];
sx q[1];
rz(2.5893455) q[1];
x q[2];
rz(1.852774) q[3];
sx q[3];
rz(-2.2681142) q[3];
sx q[3];
rz(3.092749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.272133) q[2];
sx q[2];
rz(-1.3674066) q[2];
sx q[2];
rz(-2.7222705) q[2];
rz(-1.1767293) q[3];
sx q[3];
rz(-2.4265225) q[3];
sx q[3];
rz(0.56593219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6856573) q[0];
sx q[0];
rz(-0.74860191) q[0];
sx q[0];
rz(1.5022044) q[0];
rz(-1.4337076) q[1];
sx q[1];
rz(-1.2805484) q[1];
sx q[1];
rz(-0.40275231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6667867) q[0];
sx q[0];
rz(-2.1763069) q[0];
sx q[0];
rz(0.31374991) q[0];
x q[1];
rz(0.5647955) q[2];
sx q[2];
rz(-1.5207054) q[2];
sx q[2];
rz(1.3897071) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13632475) q[1];
sx q[1];
rz(-1.439716) q[1];
sx q[1];
rz(1.6682427) q[1];
rz(1.4536269) q[3];
sx q[3];
rz(-0.90725431) q[3];
sx q[3];
rz(0.63267665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5002284) q[2];
sx q[2];
rz(-1.0470752) q[2];
sx q[2];
rz(0.25263986) q[2];
rz(2.3371475) q[3];
sx q[3];
rz(-0.52144709) q[3];
sx q[3];
rz(-2.6430602) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280403) q[0];
sx q[0];
rz(-0.10843065) q[0];
sx q[0];
rz(-2.257708) q[0];
rz(0.49993316) q[1];
sx q[1];
rz(-1.4686613) q[1];
sx q[1];
rz(-0.51441851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.894569) q[0];
sx q[0];
rz(-1.2366364) q[0];
sx q[0];
rz(-0.96118013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53192682) q[2];
sx q[2];
rz(-1.3585603) q[2];
sx q[2];
rz(-1.3440107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1757092) q[1];
sx q[1];
rz(-2.6789832) q[1];
sx q[1];
rz(1.0641896) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78935577) q[3];
sx q[3];
rz(-1.5522) q[3];
sx q[3];
rz(2.9531053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9875235) q[2];
sx q[2];
rz(-1.6776626) q[2];
sx q[2];
rz(-3.0164111) q[2];
rz(2.3212738) q[3];
sx q[3];
rz(-1.1743952) q[3];
sx q[3];
rz(2.7952747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5888551) q[0];
sx q[0];
rz(-0.6468361) q[0];
sx q[0];
rz(-3.0297025) q[0];
rz(2.0018068) q[1];
sx q[1];
rz(-0.21509376) q[1];
sx q[1];
rz(1.2824167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49408052) q[0];
sx q[0];
rz(-1.843113) q[0];
sx q[0];
rz(2.168505) q[0];
x q[1];
rz(0.23156802) q[2];
sx q[2];
rz(-1.7121669) q[2];
sx q[2];
rz(1.3249152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7805182) q[1];
sx q[1];
rz(-2.2782234) q[1];
sx q[1];
rz(2.646628) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6285283) q[3];
sx q[3];
rz(-2.2542076) q[3];
sx q[3];
rz(2.6635955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45184267) q[2];
sx q[2];
rz(-2.8359783) q[2];
sx q[2];
rz(1.8243194) q[2];
rz(0.02056038) q[3];
sx q[3];
rz(-2.8097184) q[3];
sx q[3];
rz(-2.1338972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2059712) q[0];
sx q[0];
rz(-0.99140778) q[0];
sx q[0];
rz(-0.29888612) q[0];
rz(-1.0609974) q[1];
sx q[1];
rz(-1.7540878) q[1];
sx q[1];
rz(-1.2334197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376801) q[0];
sx q[0];
rz(-2.7938346) q[0];
sx q[0];
rz(1.2832252) q[0];
rz(-pi) q[1];
rz(1.685254) q[2];
sx q[2];
rz(-1.5298163) q[2];
sx q[2];
rz(-1.6270527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1202482) q[1];
sx q[1];
rz(-1.8665736) q[1];
sx q[1];
rz(0.89521726) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0954433) q[3];
sx q[3];
rz(-0.97844175) q[3];
sx q[3];
rz(2.2762959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3491106) q[2];
sx q[2];
rz(-2.0439549) q[2];
sx q[2];
rz(0.20723542) q[2];
rz(-2.4810897) q[3];
sx q[3];
rz(-2.1410172) q[3];
sx q[3];
rz(1.8628666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5527363) q[0];
sx q[0];
rz(-1.739946) q[0];
sx q[0];
rz(1.2218342) q[0];
rz(-0.51512042) q[1];
sx q[1];
rz(-0.11038596) q[1];
sx q[1];
rz(-2.5882904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80712718) q[0];
sx q[0];
rz(-2.5688524) q[0];
sx q[0];
rz(1.80869) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4261999) q[2];
sx q[2];
rz(-2.7167277) q[2];
sx q[2];
rz(-2.1230842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8043704) q[1];
sx q[1];
rz(-2.9644199) q[1];
sx q[1];
rz(-1.6044264) q[1];
rz(-pi) q[2];
rz(3.1114486) q[3];
sx q[3];
rz(-0.6693535) q[3];
sx q[3];
rz(-3.1140259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8052266) q[2];
sx q[2];
rz(-2.0311425) q[2];
sx q[2];
rz(-0.52555788) q[2];
rz(-1.4793652) q[3];
sx q[3];
rz(-2.1879304) q[3];
sx q[3];
rz(-0.76013887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084759921) q[0];
sx q[0];
rz(-2.3469717) q[0];
sx q[0];
rz(-2.7123465) q[0];
rz(-1.6191354) q[1];
sx q[1];
rz(-1.8712964) q[1];
sx q[1];
rz(2.2116908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97668682) q[0];
sx q[0];
rz(-1.9082883) q[0];
sx q[0];
rz(-0.22708186) q[0];
rz(-pi) q[1];
rz(1.3512813) q[2];
sx q[2];
rz(-1.6782614) q[2];
sx q[2];
rz(1.4023413) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.279544) q[1];
sx q[1];
rz(-2.8144208) q[1];
sx q[1];
rz(-1.0485093) q[1];
rz(-pi) q[2];
rz(0.096654722) q[3];
sx q[3];
rz(-2.5864263) q[3];
sx q[3];
rz(-2.6688109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.736019) q[2];
sx q[2];
rz(-0.47601998) q[2];
sx q[2];
rz(-0.19301566) q[2];
rz(0.080282601) q[3];
sx q[3];
rz(-0.86997) q[3];
sx q[3];
rz(1.75753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3748462) q[0];
sx q[0];
rz(-1.2046879) q[0];
sx q[0];
rz(-1.1483703) q[0];
rz(2.171352) q[1];
sx q[1];
rz(-1.2087676) q[1];
sx q[1];
rz(-1.2420775) q[1];
rz(2.8497981) q[2];
sx q[2];
rz(-0.84635432) q[2];
sx q[2];
rz(1.9683471) q[2];
rz(1.5946078) q[3];
sx q[3];
rz(-1.6649369) q[3];
sx q[3];
rz(0.39149951) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
