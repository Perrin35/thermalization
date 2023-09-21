OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(-1.1480968) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(-2.2009067) q[1];
sx q[1];
rz(1.3936477) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6205665) q[0];
sx q[0];
rz(-1.5687843) q[0];
sx q[0];
rz(0.30962551) q[0];
rz(-pi) q[1];
rz(-0.68140985) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(-1.6207221) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4204971) q[1];
sx q[1];
rz(-1.8407397) q[1];
sx q[1];
rz(-0.47081486) q[1];
rz(-pi) q[2];
rz(2.0066891) q[3];
sx q[3];
rz(-2.5385058) q[3];
sx q[3];
rz(0.0090422612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(-0.17928784) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-2.9860935) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(-0.24138385) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.1855159) q[0];
sx q[0];
rz(-0.38498621) q[0];
rz(-2.3953305) q[2];
sx q[2];
rz(-2.4733739) q[2];
sx q[2];
rz(1.8909188) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1387716) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(-0.9072733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4003795) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(-1.1982329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(-1.2319516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(2.4011491) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(1.0528475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42441472) q[0];
sx q[0];
rz(-0.6783692) q[0];
sx q[0];
rz(0.60027392) q[0];
rz(0.89715965) q[2];
sx q[2];
rz(-1.8067915) q[2];
sx q[2];
rz(1.8434075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7317253) q[1];
sx q[1];
rz(-1.4508529) q[1];
sx q[1];
rz(-2.9883595) q[1];
rz(-pi) q[2];
rz(1.5925519) q[3];
sx q[3];
rz(-2.7113911) q[3];
sx q[3];
rz(0.83963001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8911002) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-2.5946674) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(-0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96994394) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(1.3522211) q[0];
rz(0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(-2.9052177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0683595) q[0];
sx q[0];
rz(-1.7422361) q[0];
sx q[0];
rz(1.839848) q[0];
x q[1];
rz(-2.0868446) q[2];
sx q[2];
rz(-1.2630672) q[2];
sx q[2];
rz(-2.9204521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.826556) q[1];
sx q[1];
rz(-0.16780014) q[1];
sx q[1];
rz(-1.6102953) q[1];
rz(-pi) q[2];
x q[2];
rz(2.866719) q[3];
sx q[3];
rz(-0.50516869) q[3];
sx q[3];
rz(1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(2.5881055) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(-0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.47167641) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(2.5456837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51673698) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(3.0481824) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4679568) q[2];
sx q[2];
rz(-1.4623702) q[2];
sx q[2];
rz(-0.078439586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51922819) q[1];
sx q[1];
rz(-1.9740826) q[1];
sx q[1];
rz(1.7720023) q[1];
x q[2];
rz(2.7704352) q[3];
sx q[3];
rz(-1.3988929) q[3];
sx q[3];
rz(2.2558444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(-0.55111432) q[2];
rz(-0.20714949) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(-2.8463083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2518894) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(0.15051145) q[0];
rz(-pi) q[1];
rz(-0.31111181) q[2];
sx q[2];
rz(-1.7179486) q[2];
sx q[2];
rz(2.1152903) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5438248) q[1];
sx q[1];
rz(-2.4104767) q[1];
sx q[1];
rz(2.3901229) q[1];
rz(3.0764334) q[3];
sx q[3];
rz(-2.0805196) q[3];
sx q[3];
rz(0.62419696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(1.6850083) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0657601) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(2.4767955) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-0.40329969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6751854) q[0];
sx q[0];
rz(-1.9403606) q[0];
sx q[0];
rz(-2.7778366) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67862582) q[2];
sx q[2];
rz(-1.7947065) q[2];
sx q[2];
rz(2.9976171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62337263) q[1];
sx q[1];
rz(-1.5726591) q[1];
sx q[1];
rz(-1.7809479) q[1];
x q[2];
rz(-3.1163252) q[3];
sx q[3];
rz(-1.7606887) q[3];
sx q[3];
rz(2.5437298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(2.1408391) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-2.5411434) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52727115) q[0];
sx q[0];
rz(-1.1724171) q[0];
sx q[0];
rz(2.6398753) q[0];
x q[1];
rz(1.6161643) q[2];
sx q[2];
rz(-2.6160935) q[2];
sx q[2];
rz(-0.85277992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8200092) q[1];
sx q[1];
rz(-1.6899741) q[1];
sx q[1];
rz(2.9585341) q[1];
rz(-pi) q[2];
rz(1.7793711) q[3];
sx q[3];
rz(-1.2731291) q[3];
sx q[3];
rz(1.0514333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(0.93961811) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(2.7240662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5719713) q[0];
sx q[0];
rz(-0.53335359) q[0];
sx q[0];
rz(-1.4825975) q[0];
rz(-1.3067901) q[2];
sx q[2];
rz(-1.4679113) q[2];
sx q[2];
rz(-1.4357944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7716277) q[1];
sx q[1];
rz(-1.7376448) q[1];
sx q[1];
rz(-0.19373993) q[1];
x q[2];
rz(1.1674676) q[3];
sx q[3];
rz(-0.46357337) q[3];
sx q[3];
rz(1.4232672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(0.49368668) q[2];
rz(-2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(2.0102274) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3216074) q[0];
sx q[0];
rz(-1.9157584) q[0];
sx q[0];
rz(2.5489775) q[0];
rz(-pi) q[1];
rz(-2.4360043) q[2];
sx q[2];
rz(-2.0863279) q[2];
sx q[2];
rz(-2.3956092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3494306) q[1];
sx q[1];
rz(-1.3467448) q[1];
sx q[1];
rz(0.38778023) q[1];
x q[2];
rz(2.4990436) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(2.97314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(-2.3317544) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(0.99954188) q[2];
sx q[2];
rz(-1.0456325) q[2];
sx q[2];
rz(0.13620302) q[2];
rz(2.5645732) q[3];
sx q[3];
rz(-2.1885625) q[3];
sx q[3];
rz(2.267005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
