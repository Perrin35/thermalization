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
rz(-0.96002785) q[0];
sx q[0];
rz(-2.5224944) q[0];
sx q[0];
rz(1.7382789) q[0];
rz(2.4461441) q[1];
sx q[1];
rz(-0.56013501) q[1];
sx q[1];
rz(0.6518031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2346828) q[0];
sx q[0];
rz(-1.2359338) q[0];
sx q[0];
rz(0.85483179) q[0];
x q[1];
rz(3.1211463) q[2];
sx q[2];
rz(-2.8600577) q[2];
sx q[2];
rz(0.26007465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0974724) q[1];
sx q[1];
rz(-1.5481045) q[1];
sx q[1];
rz(-2.9097662) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8709098) q[3];
sx q[3];
rz(-0.83867618) q[3];
sx q[3];
rz(-2.0441149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92242509) q[2];
sx q[2];
rz(-1.3895915) q[2];
sx q[2];
rz(-0.31828848) q[2];
rz(-0.90218941) q[3];
sx q[3];
rz(-0.92034322) q[3];
sx q[3];
rz(-2.2778146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.1282463) q[0];
sx q[0];
rz(-2.8885169) q[0];
sx q[0];
rz(2.6213562) q[0];
rz(3.0546313) q[1];
sx q[1];
rz(-2.088701) q[1];
sx q[1];
rz(0.797995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3083852) q[0];
sx q[0];
rz(-1.7705581) q[0];
sx q[0];
rz(-0.31673348) q[0];
rz(-0.74348656) q[2];
sx q[2];
rz(-1.6891589) q[2];
sx q[2];
rz(1.7497045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73436208) q[1];
sx q[1];
rz(-1.4301908) q[1];
sx q[1];
rz(0.39234061) q[1];
rz(1.1677443) q[3];
sx q[3];
rz(-1.9857996) q[3];
sx q[3];
rz(-2.5248506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.308455) q[2];
sx q[2];
rz(-1.1546346) q[2];
sx q[2];
rz(-0.40033611) q[2];
rz(-0.18276246) q[3];
sx q[3];
rz(-0.57923135) q[3];
sx q[3];
rz(1.8496752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12277814) q[0];
sx q[0];
rz(-0.084538758) q[0];
sx q[0];
rz(-1.2070745) q[0];
rz(-1.6385551) q[1];
sx q[1];
rz(-1.8820347) q[1];
sx q[1];
rz(2.9685793) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68727815) q[0];
sx q[0];
rz(-1.1862687) q[0];
sx q[0];
rz(0.092057629) q[0];
rz(-pi) q[1];
rz(1.7373213) q[2];
sx q[2];
rz(-0.58812599) q[2];
sx q[2];
rz(3.0808133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26721482) q[1];
sx q[1];
rz(-2.3372531) q[1];
sx q[1];
rz(-0.72093236) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9152619) q[3];
sx q[3];
rz(-0.76068288) q[3];
sx q[3];
rz(0.62225354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5485237) q[2];
sx q[2];
rz(-2.5544281) q[2];
sx q[2];
rz(-1.0347838) q[2];
rz(0.89094025) q[3];
sx q[3];
rz(-0.8689298) q[3];
sx q[3];
rz(-1.2263891) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0767176) q[0];
sx q[0];
rz(-2.1893976) q[0];
sx q[0];
rz(-3.0127443) q[0];
rz(0.42875641) q[1];
sx q[1];
rz(-1.7196451) q[1];
sx q[1];
rz(1.4770003) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98618988) q[0];
sx q[0];
rz(-2.7895067) q[0];
sx q[0];
rz(-2.9425063) q[0];
x q[1];
rz(0.44775362) q[2];
sx q[2];
rz(-1.9734255) q[2];
sx q[2];
rz(-1.778971) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99812388) q[1];
sx q[1];
rz(-1.7822155) q[1];
sx q[1];
rz(-1.9969673) q[1];
x q[2];
rz(2.1620524) q[3];
sx q[3];
rz(-1.5717122) q[3];
sx q[3];
rz(2.9360848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9513272) q[2];
sx q[2];
rz(-1.4172047) q[2];
sx q[2];
rz(-1.8180656) q[2];
rz(1.4642814) q[3];
sx q[3];
rz(-2.5530294) q[3];
sx q[3];
rz(2.467449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73439634) q[0];
sx q[0];
rz(-0.36235991) q[0];
sx q[0];
rz(1.6998442) q[0];
rz(0.4231407) q[1];
sx q[1];
rz(-0.95760456) q[1];
sx q[1];
rz(-0.29022455) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6184171) q[0];
sx q[0];
rz(-1.8632065) q[0];
sx q[0];
rz(0.66177701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82478451) q[2];
sx q[2];
rz(-2.3815577) q[2];
sx q[2];
rz(1.7766376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90846244) q[1];
sx q[1];
rz(-1.6987579) q[1];
sx q[1];
rz(-1.5217693) q[1];
rz(1.832015) q[3];
sx q[3];
rz(-1.272837) q[3];
sx q[3];
rz(3.1233149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4238802) q[2];
sx q[2];
rz(-2.3545357) q[2];
sx q[2];
rz(-0.94669739) q[2];
rz(-2.8801019) q[3];
sx q[3];
rz(-1.7777092) q[3];
sx q[3];
rz(-2.7546895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.43712) q[0];
sx q[0];
rz(-1.8081212) q[0];
sx q[0];
rz(2.3719846) q[0];
rz(-2.9656124) q[1];
sx q[1];
rz(-1.3409921) q[1];
sx q[1];
rz(-1.5699068) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0024182323) q[0];
sx q[0];
rz(-0.95390924) q[0];
sx q[0];
rz(2.3611538) q[0];
rz(3.0716586) q[2];
sx q[2];
rz(-0.85288793) q[2];
sx q[2];
rz(-2.1760035) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1842725) q[1];
sx q[1];
rz(-1.2046923) q[1];
sx q[1];
rz(-1.7680607) q[1];
rz(-0.50521781) q[3];
sx q[3];
rz(-2.0310302) q[3];
sx q[3];
rz(1.9379747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0104597) q[2];
sx q[2];
rz(-2.185148) q[2];
sx q[2];
rz(-1.3156816) q[2];
rz(1.2223505) q[3];
sx q[3];
rz(-2.0171916) q[3];
sx q[3];
rz(-2.840672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9199801) q[0];
sx q[0];
rz(-1.6918809) q[0];
sx q[0];
rz(1.8020887) q[0];
rz(1.9446531) q[1];
sx q[1];
rz(-1.0547799) q[1];
sx q[1];
rz(-0.96492499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6048991) q[0];
sx q[0];
rz(-1.6654547) q[0];
sx q[0];
rz(-0.51095922) q[0];
rz(-pi) q[1];
rz(-2.262124) q[2];
sx q[2];
rz(-0.44422883) q[2];
sx q[2];
rz(-2.6951126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2141909) q[1];
sx q[1];
rz(-1.7905231) q[1];
sx q[1];
rz(1.6564547) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8172142) q[3];
sx q[3];
rz(-2.6357963) q[3];
sx q[3];
rz(1.5599111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40184608) q[2];
sx q[2];
rz(-1.9332989) q[2];
sx q[2];
rz(1.1811264) q[2];
rz(-0.014569672) q[3];
sx q[3];
rz(-2.7404692) q[3];
sx q[3];
rz(-1.1723664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3907617) q[0];
sx q[0];
rz(-1.9347235) q[0];
sx q[0];
rz(-2.3857351) q[0];
rz(-0.66468704) q[1];
sx q[1];
rz(-1.942778) q[1];
sx q[1];
rz(1.8519148) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15890789) q[0];
sx q[0];
rz(-2.3450627) q[0];
sx q[0];
rz(-0.089287306) q[0];
rz(-pi) q[1];
rz(2.6014156) q[2];
sx q[2];
rz(-1.8312757) q[2];
sx q[2];
rz(1.8249701) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.522053) q[1];
sx q[1];
rz(-1.8541946) q[1];
sx q[1];
rz(0.26669054) q[1];
x q[2];
rz(-1.6716854) q[3];
sx q[3];
rz(-1.8785739) q[3];
sx q[3];
rz(0.65478126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1332625) q[2];
sx q[2];
rz(-1.3349168) q[2];
sx q[2];
rz(1.5745715) q[2];
rz(1.3957006) q[3];
sx q[3];
rz(-1.402366) q[3];
sx q[3];
rz(0.76656669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-1.72346) q[0];
sx q[0];
rz(-2.8657931) q[0];
sx q[0];
rz(-0.3904528) q[0];
rz(2.0402724) q[1];
sx q[1];
rz(-1.3464255) q[1];
sx q[1];
rz(0.93213814) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86499914) q[0];
sx q[0];
rz(-3.0955934) q[0];
sx q[0];
rz(-0.92711751) q[0];
x q[1];
rz(-0.081285921) q[2];
sx q[2];
rz(-2.6008106) q[2];
sx q[2];
rz(0.062449156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7919346) q[1];
sx q[1];
rz(-2.8009543) q[1];
sx q[1];
rz(0.61130166) q[1];
x q[2];
rz(-2.6868709) q[3];
sx q[3];
rz(-2.623436) q[3];
sx q[3];
rz(-1.3738812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.091888) q[2];
sx q[2];
rz(-2.4447417) q[2];
sx q[2];
rz(-1.6799124) q[2];
rz(-1.0852496) q[3];
sx q[3];
rz(-1.1321675) q[3];
sx q[3];
rz(-0.6404883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853448) q[0];
sx q[0];
rz(-2.4818821) q[0];
sx q[0];
rz(2.0154542) q[0];
rz(-2.0938734) q[1];
sx q[1];
rz(-0.62647096) q[1];
sx q[1];
rz(1.610021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0275912) q[0];
sx q[0];
rz(-3.1182814) q[0];
sx q[0];
rz(-0.059904532) q[0];
rz(1.3381027) q[2];
sx q[2];
rz(-1.1232716) q[2];
sx q[2];
rz(-1.1648503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46343741) q[1];
sx q[1];
rz(-0.10129866) q[1];
sx q[1];
rz(0.8763635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6592386) q[3];
sx q[3];
rz(-2.1312661) q[3];
sx q[3];
rz(2.0129528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62022007) q[2];
sx q[2];
rz(-0.81338254) q[2];
sx q[2];
rz(2.9353471) q[2];
rz(-2.5824052) q[3];
sx q[3];
rz(-1.6181479) q[3];
sx q[3];
rz(1.1243813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6013721) q[0];
sx q[0];
rz(-1.843353) q[0];
sx q[0];
rz(-2.699615) q[0];
rz(-1.1454918) q[1];
sx q[1];
rz(-1.835123) q[1];
sx q[1];
rz(1.8225972) q[1];
rz(0.27362846) q[2];
sx q[2];
rz(-1.8406592) q[2];
sx q[2];
rz(-2.1218268) q[2];
rz(-0.72412422) q[3];
sx q[3];
rz(-1.4403314) q[3];
sx q[3];
rz(-2.452313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
