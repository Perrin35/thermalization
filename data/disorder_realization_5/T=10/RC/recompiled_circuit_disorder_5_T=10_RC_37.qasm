OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(2.6846057) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63859343) q[0];
sx q[0];
rz(-0.25829691) q[0];
sx q[0];
rz(0.55708416) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4957501) q[2];
sx q[2];
rz(-1.6610166) q[2];
sx q[2];
rz(0.98352945) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.16986019) q[1];
sx q[1];
rz(-1.5105376) q[1];
sx q[1];
rz(-1.7737284) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2184578) q[3];
sx q[3];
rz(-1.2520408) q[3];
sx q[3];
rz(1.2649328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0695641) q[0];
sx q[0];
rz(-1.1840128) q[0];
sx q[0];
rz(2.8974429) q[0];
rz(-pi) q[1];
rz(-1.7581802) q[2];
sx q[2];
rz(-2.6307081) q[2];
sx q[2];
rz(2.5201706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74233494) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(-1.4820815) q[1];
rz(-0.12985274) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(-3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(-2.6248698) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(1.3396938) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08683603) q[0];
sx q[0];
rz(-1.1491547) q[0];
sx q[0];
rz(-2.3479168) q[0];
x q[1];
rz(-0.3091829) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(3.0004629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9078319) q[1];
sx q[1];
rz(-1.7572174) q[1];
sx q[1];
rz(0.16252653) q[1];
rz(-1.9698824) q[3];
sx q[3];
rz(-0.92150021) q[3];
sx q[3];
rz(2.2513575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.7335256) q[2];
rz(0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(0.66657153) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4144856) q[0];
sx q[0];
rz(-1.3654728) q[0];
sx q[0];
rz(-1.3798825) q[0];
x q[1];
rz(2.4273657) q[2];
sx q[2];
rz(-1.3492609) q[2];
sx q[2];
rz(2.5941338) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92723144) q[1];
sx q[1];
rz(-0.9775369) q[1];
sx q[1];
rz(-0.67770047) q[1];
rz(1.7188844) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(1.577391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(2.3045325) q[0];
rz(-2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.9821092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5975208) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(-1.9835299) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7972838) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(3.0938873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3573208) q[1];
sx q[1];
rz(-1.1457774) q[1];
sx q[1];
rz(0.97370633) q[1];
rz(-pi) q[2];
rz(-0.67622185) q[3];
sx q[3];
rz(-2.3822504) q[3];
sx q[3];
rz(-0.29952213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5669252) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-3.040722) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2642919) q[0];
sx q[0];
rz(-1.722324) q[0];
sx q[0];
rz(2.507188) q[0];
x q[1];
rz(2.6542873) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(-3.0234408) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.069177376) q[1];
sx q[1];
rz(-0.25670708) q[1];
sx q[1];
rz(0.91255811) q[1];
rz(2.7363339) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(0.98291558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(2.4538453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15935414) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(2.9100145) q[0];
rz(0.22600941) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(1.4357391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5016985) q[1];
sx q[1];
rz(-1.5475029) q[1];
sx q[1];
rz(3.0122019) q[1];
rz(-pi) q[2];
rz(-1.0024768) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-2.2040099) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90889701) q[0];
sx q[0];
rz(-2.7017653) q[0];
sx q[0];
rz(0.16172414) q[0];
x q[1];
rz(1.7940117) q[2];
sx q[2];
rz(-2.7701021) q[2];
sx q[2];
rz(-1.5471293) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.291154) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(-0.50317851) q[1];
rz(2.6506181) q[3];
sx q[3];
rz(-2.6262865) q[3];
sx q[3];
rz(1.5332424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(-2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.8607128) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(-0.58832204) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(0.9334329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342455) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(-3.0349915) q[0];
rz(-0.8546631) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(1.2235299) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.583608) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(-0.058137356) q[1];
x q[2];
rz(0.65255717) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(-2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(2.771647) q[2];
rz(-0.89921078) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(-2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1856954) q[0];
sx q[0];
rz(-0.16793185) q[0];
sx q[0];
rz(-0.76783617) q[0];
rz(2.3046938) q[2];
sx q[2];
rz(-2.1914334) q[2];
sx q[2];
rz(1.2037954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2059584) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(2.2531855) q[1];
rz(3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(-1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.7741268) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(2.7783685) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(2.6828962) q[2];
sx q[2];
rz(-1.9445322) q[2];
sx q[2];
rz(-1.4945488) q[2];
rz(1.8300874) q[3];
sx q[3];
rz(-0.6060096) q[3];
sx q[3];
rz(0.93939645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
