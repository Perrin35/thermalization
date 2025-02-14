OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57920116) q[0];
sx q[0];
rz(-0.76051036) q[0];
sx q[0];
rz(1.0150681) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(3.562869) q[1];
sx q[1];
rz(8.1864551) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27313156) q[0];
sx q[0];
rz(-2.1885311) q[0];
sx q[0];
rz(-2.4988089) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2767645) q[2];
sx q[2];
rz(-2.6061686) q[2];
sx q[2];
rz(3.0147417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4290966) q[1];
sx q[1];
rz(-2.0331906) q[1];
sx q[1];
rz(-2.7346947) q[1];
rz(-pi) q[2];
rz(-2.5776776) q[3];
sx q[3];
rz(-1.7132572) q[3];
sx q[3];
rz(-3.0092653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(2.2859196) q[3];
sx q[3];
rz(-2.296505) q[3];
sx q[3];
rz(1.8446911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6913476) q[0];
sx q[0];
rz(-1.0991199) q[0];
sx q[0];
rz(0.74478373) q[0];
rz(-1.5738457) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(-0.94211284) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0364931) q[0];
sx q[0];
rz(-0.78723037) q[0];
sx q[0];
rz(2.8544642) q[0];
x q[1];
rz(3.0622185) q[2];
sx q[2];
rz(-0.81981507) q[2];
sx q[2];
rz(-2.306983) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5444138) q[1];
sx q[1];
rz(-2.326818) q[1];
sx q[1];
rz(-0.18715231) q[1];
x q[2];
rz(-1.510624) q[3];
sx q[3];
rz(-2.1561557) q[3];
sx q[3];
rz(-0.28096052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2524903) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(1.4787632) q[3];
sx q[3];
rz(-2.3507599) q[3];
sx q[3];
rz(1.7553294) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304994) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(-1.4053364) q[0];
rz(-1.3966712) q[1];
sx q[1];
rz(-1.3537355) q[1];
sx q[1];
rz(0.90000802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4858682) q[0];
sx q[0];
rz(-1.2085755) q[0];
sx q[0];
rz(2.3535504) q[0];
rz(0.75863691) q[2];
sx q[2];
rz(-0.73497811) q[2];
sx q[2];
rz(1.0204362) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.450728) q[1];
sx q[1];
rz(-2.303586) q[1];
sx q[1];
rz(-1.8045252) q[1];
rz(-pi) q[2];
rz(2.4056959) q[3];
sx q[3];
rz(-2.3319339) q[3];
sx q[3];
rz(2.5084605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6806543) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(-2.1920965) q[2];
rz(-2.5203868) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(-1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42156521) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(0.048728745) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(2.8299832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2545812) q[0];
sx q[0];
rz(-1.100228) q[0];
sx q[0];
rz(-2.8356617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7976102) q[2];
sx q[2];
rz(-2.0728353) q[2];
sx q[2];
rz(-0.12884049) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8383465) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(1.1619199) q[1];
rz(2.0120088) q[3];
sx q[3];
rz(-1.3008683) q[3];
sx q[3];
rz(-2.2306311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.17778808) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(1.6507899) q[2];
rz(2.7866411) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(-1.8748913) q[0];
rz(-0.1162687) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(1.6273392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0851885) q[0];
sx q[0];
rz(-1.8833816) q[0];
sx q[0];
rz(-2.2645135) q[0];
rz(-pi) q[1];
rz(1.7909796) q[2];
sx q[2];
rz(-1.4823969) q[2];
sx q[2];
rz(-0.92313405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5942638) q[1];
sx q[1];
rz(-1.0272861) q[1];
sx q[1];
rz(-3.0819403) q[1];
x q[2];
rz(-0.84228869) q[3];
sx q[3];
rz(-1.5986048) q[3];
sx q[3];
rz(-2.9942346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9194455) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(1.1754645) q[2];
rz(-2.2373824) q[3];
sx q[3];
rz(-1.2491106) q[3];
sx q[3];
rz(0.97833943) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745875) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(2.7401155) q[0];
rz(-2.0145156) q[1];
sx q[1];
rz(-1.3104442) q[1];
sx q[1];
rz(2.3289767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881587) q[0];
sx q[0];
rz(-1.5749314) q[0];
sx q[0];
rz(-0.93664108) q[0];
x q[1];
rz(2.866686) q[2];
sx q[2];
rz(-1.3454352) q[2];
sx q[2];
rz(-1.7079086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8725295) q[1];
sx q[1];
rz(-1.7751964) q[1];
sx q[1];
rz(-2.5356631) q[1];
x q[2];
rz(-2.3964368) q[3];
sx q[3];
rz(-1.1155012) q[3];
sx q[3];
rz(-0.32768341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(-1.5852488) q[2];
rz(-0.19041348) q[3];
sx q[3];
rz(-2.3868581) q[3];
sx q[3];
rz(1.93369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.82454005) q[0];
sx q[0];
rz(-0.38924488) q[0];
sx q[0];
rz(-3.0188766) q[0];
rz(-1.1931194) q[1];
sx q[1];
rz(-2.2331608) q[1];
sx q[1];
rz(0.76748031) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90119367) q[0];
sx q[0];
rz(-1.6220105) q[0];
sx q[0];
rz(-2.2177494) q[0];
rz(-0.38184719) q[2];
sx q[2];
rz(-1.6497465) q[2];
sx q[2];
rz(0.71982924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9373405) q[1];
sx q[1];
rz(-0.91751999) q[1];
sx q[1];
rz(-0.055257992) q[1];
x q[2];
rz(-1.5662868) q[3];
sx q[3];
rz(-1.9770375) q[3];
sx q[3];
rz(0.81886983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(1.5519315) q[2];
rz(1.4895561) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(2.0261197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234696) q[0];
sx q[0];
rz(-2.7796845) q[0];
sx q[0];
rz(2.7662011) q[0];
rz(-2.2581532) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(-0.72867957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.323384) q[0];
sx q[0];
rz(-0.94037708) q[0];
sx q[0];
rz(0.66108836) q[0];
x q[1];
rz(2.905718) q[2];
sx q[2];
rz(-1.7348924) q[2];
sx q[2];
rz(2.1498888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.091944253) q[1];
sx q[1];
rz(-0.94251213) q[1];
sx q[1];
rz(-0.44558427) q[1];
rz(-pi) q[2];
rz(2.4861195) q[3];
sx q[3];
rz(-0.87067662) q[3];
sx q[3];
rz(-1.164013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6931307) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(-0.7473839) q[2];
rz(1.735431) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2324227) q[0];
sx q[0];
rz(-0.94038832) q[0];
sx q[0];
rz(0.7255834) q[0];
rz(-1.8046509) q[1];
sx q[1];
rz(-1.888211) q[1];
sx q[1];
rz(0.57428378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5142412) q[0];
sx q[0];
rz(-2.9459369) q[0];
sx q[0];
rz(0.95558863) q[0];
x q[1];
rz(2.1404187) q[2];
sx q[2];
rz(-1.9211848) q[2];
sx q[2];
rz(-2.1621665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7567819) q[1];
sx q[1];
rz(-0.73721209) q[1];
sx q[1];
rz(-1.739915) q[1];
rz(-2.5602407) q[3];
sx q[3];
rz(-0.87331088) q[3];
sx q[3];
rz(-2.4331828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9743222) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(0.66217011) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26226703) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(-1.7437438) q[0];
rz(-1.0700048) q[1];
sx q[1];
rz(-2.013423) q[1];
sx q[1];
rz(1.3319344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86517143) q[0];
sx q[0];
rz(-1.4201453) q[0];
sx q[0];
rz(0.030929203) q[0];
rz(-pi) q[1];
rz(-2.0831624) q[2];
sx q[2];
rz(-2.67423) q[2];
sx q[2];
rz(2.1455554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9402658) q[1];
sx q[1];
rz(-1.3405217) q[1];
sx q[1];
rz(0.23576945) q[1];
rz(2.2666107) q[3];
sx q[3];
rz(-2.0379645) q[3];
sx q[3];
rz(2.9206729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8269044) q[2];
sx q[2];
rz(-2.0411699) q[2];
sx q[2];
rz(0.17987128) q[2];
rz(-1.3966857) q[3];
sx q[3];
rz(-1.4254009) q[3];
sx q[3];
rz(2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(-0.63498598) q[2];
sx q[2];
rz(-1.9188835) q[2];
sx q[2];
rz(-1.9139482) q[2];
rz(0.64416364) q[3];
sx q[3];
rz(-0.5683036) q[3];
sx q[3];
rz(-0.34556942) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
