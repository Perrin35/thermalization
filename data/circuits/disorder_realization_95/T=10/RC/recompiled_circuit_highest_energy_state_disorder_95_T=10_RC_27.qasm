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
rz(0.79701841) q[0];
sx q[0];
rz(-2.2198644) q[0];
sx q[0];
rz(1.2603238) q[0];
rz(0.64646512) q[1];
sx q[1];
rz(4.0151172) q[1];
sx q[1];
rz(9.092796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9063569) q[0];
sx q[0];
rz(-0.42354326) q[0];
sx q[0];
rz(-1.5111501) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69692287) q[2];
sx q[2];
rz(-1.2771318) q[2];
sx q[2];
rz(2.1611093) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6011966) q[1];
sx q[1];
rz(-2.3292377) q[1];
sx q[1];
rz(-1.5271714) q[1];
rz(-0.1370955) q[3];
sx q[3];
rz(-2.6957316) q[3];
sx q[3];
rz(-0.23030989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6648286) q[2];
sx q[2];
rz(-0.9330743) q[2];
sx q[2];
rz(-2.942371) q[2];
rz(-2.9557989) q[3];
sx q[3];
rz(-1.7265065) q[3];
sx q[3];
rz(-1.4411521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.324447) q[0];
sx q[0];
rz(-2.6528093) q[0];
sx q[0];
rz(0.73844886) q[0];
rz(-1.4456519) q[1];
sx q[1];
rz(-1.4001458) q[1];
sx q[1];
rz(-2.6587528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49119887) q[0];
sx q[0];
rz(-0.20696124) q[0];
sx q[0];
rz(-0.82024337) q[0];
x q[1];
rz(-0.17530967) q[2];
sx q[2];
rz(-1.2273437) q[2];
sx q[2];
rz(-0.94409787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.184471) q[1];
sx q[1];
rz(-1.0252684) q[1];
sx q[1];
rz(1.0218744) q[1];
rz(-pi) q[2];
rz(-0.16359292) q[3];
sx q[3];
rz(-1.0181659) q[3];
sx q[3];
rz(-1.1887417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8740251) q[2];
sx q[2];
rz(-1.4635307) q[2];
sx q[2];
rz(-1.2919424) q[2];
rz(1.1235631) q[3];
sx q[3];
rz(-0.70190391) q[3];
sx q[3];
rz(1.7263713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8181151) q[0];
sx q[0];
rz(-0.98387843) q[0];
sx q[0];
rz(1.4418607) q[0];
rz(1.9135176) q[1];
sx q[1];
rz(-1.9107198) q[1];
sx q[1];
rz(1.0736116) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5283877) q[0];
sx q[0];
rz(-2.3484485) q[0];
sx q[0];
rz(2.6499477) q[0];
rz(0.91543897) q[2];
sx q[2];
rz(-1.4129593) q[2];
sx q[2];
rz(-2.137568) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50987989) q[1];
sx q[1];
rz(-2.2728765) q[1];
sx q[1];
rz(3.0982262) q[1];
rz(-pi) q[2];
rz(-1.4074247) q[3];
sx q[3];
rz(-1.1187807) q[3];
sx q[3];
rz(-0.79240914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29232612) q[2];
sx q[2];
rz(-0.30305114) q[2];
sx q[2];
rz(2.1185875) q[2];
rz(3.0813713) q[3];
sx q[3];
rz(-1.9950208) q[3];
sx q[3];
rz(-2.4165418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4802454) q[0];
sx q[0];
rz(-0.1460954) q[0];
sx q[0];
rz(2.6222141) q[0];
rz(-2.5410779) q[1];
sx q[1];
rz(-2.0527716) q[1];
sx q[1];
rz(-0.75469887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31904991) q[0];
sx q[0];
rz(-0.62915914) q[0];
sx q[0];
rz(0.99790093) q[0];
x q[1];
rz(-2.2738931) q[2];
sx q[2];
rz(-1.6804196) q[2];
sx q[2];
rz(-2.2330957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13184758) q[1];
sx q[1];
rz(-2.5452217) q[1];
sx q[1];
rz(0.98747298) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58784501) q[3];
sx q[3];
rz(-2.4949772) q[3];
sx q[3];
rz(0.69348303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34919136) q[2];
sx q[2];
rz(-0.79734048) q[2];
sx q[2];
rz(-2.965773) q[2];
rz(1.1261806) q[3];
sx q[3];
rz(-1.3577941) q[3];
sx q[3];
rz(-0.064182909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46614161) q[0];
sx q[0];
rz(-1.4956681) q[0];
sx q[0];
rz(0.58188907) q[0];
rz(1.1833082) q[1];
sx q[1];
rz(-0.71154037) q[1];
sx q[1];
rz(1.2540832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7447356) q[0];
sx q[0];
rz(-1.893147) q[0];
sx q[0];
rz(-0.8650604) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0462363) q[2];
sx q[2];
rz(-2.0031906) q[2];
sx q[2];
rz(0.43579416) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9068992) q[1];
sx q[1];
rz(-0.36415283) q[1];
sx q[1];
rz(2.1910648) q[1];
rz(-pi) q[2];
rz(-0.50458377) q[3];
sx q[3];
rz(-1.4047897) q[3];
sx q[3];
rz(-2.5291251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1648272) q[2];
sx q[2];
rz(-0.81563121) q[2];
sx q[2];
rz(2.7867479) q[2];
rz(-1.1497078) q[3];
sx q[3];
rz(-1.0747654) q[3];
sx q[3];
rz(2.3000075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717764) q[0];
sx q[0];
rz(-2.8700097) q[0];
sx q[0];
rz(-3.1014882) q[0];
rz(1.4647723) q[1];
sx q[1];
rz(-0.96950871) q[1];
sx q[1];
rz(0.47223314) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6567447) q[0];
sx q[0];
rz(-1.7539331) q[0];
sx q[0];
rz(-1.8369863) q[0];
rz(-pi) q[1];
rz(-0.91798616) q[2];
sx q[2];
rz(-1.5734104) q[2];
sx q[2];
rz(-1.6592342) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80010092) q[1];
sx q[1];
rz(-1.2849548) q[1];
sx q[1];
rz(-1.915201) q[1];
rz(-pi) q[2];
rz(-2.84255) q[3];
sx q[3];
rz(-2.3573993) q[3];
sx q[3];
rz(-1.2687781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.27452305) q[2];
sx q[2];
rz(-2.3752866) q[2];
sx q[2];
rz(1.4806032) q[2];
rz(3.1001422) q[3];
sx q[3];
rz(-2.4292414) q[3];
sx q[3];
rz(0.97431549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90254766) q[0];
sx q[0];
rz(-0.79896611) q[0];
sx q[0];
rz(0.7582742) q[0];
rz(0.43765086) q[1];
sx q[1];
rz(-1.2145019) q[1];
sx q[1];
rz(0.95380107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3054558) q[0];
sx q[0];
rz(-1.1887595) q[0];
sx q[0];
rz(-3.0925095) q[0];
rz(-pi) q[1];
rz(-0.96289159) q[2];
sx q[2];
rz(-0.94568077) q[2];
sx q[2];
rz(-1.6776379) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0702677) q[1];
sx q[1];
rz(-2.4355781) q[1];
sx q[1];
rz(-1.474375) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1498044) q[3];
sx q[3];
rz(-2.6564993) q[3];
sx q[3];
rz(2.4320784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4147676) q[2];
sx q[2];
rz(-2.813952) q[2];
sx q[2];
rz(-0.3978351) q[2];
rz(-1.2658524) q[3];
sx q[3];
rz(-1.7222907) q[3];
sx q[3];
rz(-0.46441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9809113) q[0];
sx q[0];
rz(-2.2137764) q[0];
sx q[0];
rz(-2.8371147) q[0];
rz(1.132384) q[1];
sx q[1];
rz(-2.4696923) q[1];
sx q[1];
rz(-2.0679881) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96655) q[0];
sx q[0];
rz(-1.5538408) q[0];
sx q[0];
rz(-1.5938363) q[0];
x q[1];
rz(0.4173293) q[2];
sx q[2];
rz(-1.9021735) q[2];
sx q[2];
rz(1.1294236) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4661387) q[1];
sx q[1];
rz(-2.4917648) q[1];
sx q[1];
rz(1.9790348) q[1];
rz(-pi) q[2];
rz(1.520548) q[3];
sx q[3];
rz(-2.2138688) q[3];
sx q[3];
rz(-2.1724043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73211804) q[2];
sx q[2];
rz(-1.2827337) q[2];
sx q[2];
rz(-0.049962433) q[2];
rz(0.91160715) q[3];
sx q[3];
rz(-2.215569) q[3];
sx q[3];
rz(0.91134206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53495812) q[0];
sx q[0];
rz(-1.7462523) q[0];
sx q[0];
rz(-2.7789136) q[0];
rz(-2.4210988) q[1];
sx q[1];
rz(-2.3299005) q[1];
sx q[1];
rz(-1.9627176) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3772081) q[0];
sx q[0];
rz(-0.81184719) q[0];
sx q[0];
rz(0.7332515) q[0];
rz(-2.8209053) q[2];
sx q[2];
rz(-1.6455407) q[2];
sx q[2];
rz(2.8976669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8803823) q[1];
sx q[1];
rz(-2.9453813) q[1];
sx q[1];
rz(1.0537149) q[1];
rz(-pi) q[2];
rz(2.5163609) q[3];
sx q[3];
rz(-0.64230157) q[3];
sx q[3];
rz(1.5361749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1371586) q[2];
sx q[2];
rz(-0.68778554) q[2];
sx q[2];
rz(-1.2031817) q[2];
rz(3.1249937) q[3];
sx q[3];
rz(-1.4998452) q[3];
sx q[3];
rz(-1.2703007) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5430629) q[0];
sx q[0];
rz(-0.5235343) q[0];
sx q[0];
rz(-1.6981) q[0];
rz(0.47830018) q[1];
sx q[1];
rz(-2.4595478) q[1];
sx q[1];
rz(-1.2102478) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2624274) q[0];
sx q[0];
rz(-1.0813923) q[0];
sx q[0];
rz(-0.64206815) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9094142) q[2];
sx q[2];
rz(-2.4781151) q[2];
sx q[2];
rz(0.56319204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12809556) q[1];
sx q[1];
rz(-1.6229318) q[1];
sx q[1];
rz(2.8166822) q[1];
rz(-0.5860255) q[3];
sx q[3];
rz(-2.7507493) q[3];
sx q[3];
rz(-1.9236717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9532507) q[2];
sx q[2];
rz(-2.5276999) q[2];
sx q[2];
rz(-0.71315145) q[2];
rz(0.6330511) q[3];
sx q[3];
rz(-2.1745067) q[3];
sx q[3];
rz(0.87575325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0093832) q[0];
sx q[0];
rz(-1.3306946) q[0];
sx q[0];
rz(2.7015986) q[0];
rz(0.72883365) q[1];
sx q[1];
rz(-0.76580096) q[1];
sx q[1];
rz(-2.0892807) q[1];
rz(-1.296879) q[2];
sx q[2];
rz(-2.7379681) q[2];
sx q[2];
rz(-2.2899514) q[2];
rz(1.9464372) q[3];
sx q[3];
rz(-0.79938625) q[3];
sx q[3];
rz(1.1342794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
