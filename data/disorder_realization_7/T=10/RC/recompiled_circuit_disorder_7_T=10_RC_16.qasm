OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(-0.20110826) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-2.1366871) q[0];
rz(0.5958545) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(-2.8516304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1605692) q[1];
sx q[1];
rz(-1.8796646) q[1];
sx q[1];
rz(-0.94998756) q[1];
rz(-pi) q[2];
rz(-0.17625531) q[3];
sx q[3];
rz(-1.3000543) q[3];
sx q[3];
rz(-0.83812974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.25519249) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(-2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24793808) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(3.0342297) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1194612) q[0];
sx q[0];
rz(-2.4298926) q[0];
sx q[0];
rz(2.2244781) q[0];
x q[1];
rz(2.7254843) q[2];
sx q[2];
rz(-2.5118718) q[2];
sx q[2];
rz(2.1925418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3423691) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(2.6745822) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3932088) q[3];
sx q[3];
rz(-2.7244096) q[3];
sx q[3];
rz(1.8821017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169516) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(2.4798933) q[2];
rz(0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-1.0823762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45535116) q[0];
sx q[0];
rz(-2.4903114) q[0];
sx q[0];
rz(-0.63182802) q[0];
x q[1];
rz(-1.1788549) q[2];
sx q[2];
rz(-2.760663) q[2];
sx q[2];
rz(0.71289635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0116926) q[1];
sx q[1];
rz(-1.7923647) q[1];
sx q[1];
rz(1.3264015) q[1];
rz(0.2228959) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(2.5877) q[2];
rz(-1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(0.035382263) q[0];
rz(0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-2.6534973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(-0.086829348) q[0];
rz(-pi) q[1];
rz(-2.343802) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(-2.4385902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0636053) q[1];
sx q[1];
rz(-1.6997153) q[1];
sx q[1];
rz(1.7267666) q[1];
rz(-pi) q[2];
rz(-1.0933236) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35486832) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(0.3702634) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(0.19651861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949961) q[0];
sx q[0];
rz(-1.8509522) q[0];
sx q[0];
rz(2.8269935) q[0];
x q[1];
rz(2.0318803) q[2];
sx q[2];
rz(-1.5994674) q[2];
sx q[2];
rz(1.2210786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4938441) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(0.26248787) q[1];
x q[2];
rz(0.27627857) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(-2.939558) q[0];
rz(-2.0687885) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(-1.3938168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545647) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(1.4902671) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9453085) q[2];
sx q[2];
rz(-1.6746582) q[2];
sx q[2];
rz(1.9749157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60336941) q[1];
sx q[1];
rz(-1.1829611) q[1];
sx q[1];
rz(-2.8549854) q[1];
rz(-pi) q[2];
rz(-2.6447547) q[3];
sx q[3];
rz(-0.79998575) q[3];
sx q[3];
rz(-0.62832384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(-2.2992772) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(-2.1202309) q[0];
rz(-2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(2.9313415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410941) q[0];
sx q[0];
rz(-2.9919618) q[0];
sx q[0];
rz(-1.11087) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9779151) q[2];
sx q[2];
rz(-1.5058766) q[2];
sx q[2];
rz(-2.373873) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0381895) q[1];
sx q[1];
rz(-1.8870755) q[1];
sx q[1];
rz(1.3693387) q[1];
x q[2];
rz(0.57742124) q[3];
sx q[3];
rz(-2.2038659) q[3];
sx q[3];
rz(-3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(2.3794877) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(-2.5543509) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2864838) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(2.8141862) q[0];
rz(-pi) q[1];
rz(3.061053) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(-1.5494407) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50423065) q[1];
sx q[1];
rz(-1.1334051) q[1];
sx q[1];
rz(2.7573542) q[1];
rz(-pi) q[2];
rz(-2.1920491) q[3];
sx q[3];
rz(-0.81406677) q[3];
sx q[3];
rz(0.66063389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(-2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-2.2424973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9102525) q[0];
sx q[0];
rz(-1.9950584) q[0];
sx q[0];
rz(-2.2146241) q[0];
x q[1];
rz(2.8092438) q[2];
sx q[2];
rz(-2.1045448) q[2];
sx q[2];
rz(-0.74408434) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2832665) q[1];
sx q[1];
rz(-0.35052931) q[1];
sx q[1];
rz(-1.9117029) q[1];
rz(-0.00021342834) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(-0.40094024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(-0.83834046) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(-0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(2.7446279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871646) q[0];
sx q[0];
rz(-1.0930601) q[0];
sx q[0];
rz(0.54425311) q[0];
rz(2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(3.0416833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1399999) q[1];
sx q[1];
rz(-1.6591751) q[1];
sx q[1];
rz(1.2338851) q[1];
rz(-pi) q[2];
rz(-0.40542094) q[3];
sx q[3];
rz(-2.0339031) q[3];
sx q[3];
rz(-2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(1.8718406) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(0.66432129) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(1.312064) q[3];
sx q[3];
rz(-1.4553634) q[3];
sx q[3];
rz(2.8298557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];