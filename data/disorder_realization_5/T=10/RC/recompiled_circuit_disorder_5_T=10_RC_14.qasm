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
rz(-0.45698693) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(-1.8656123) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306643) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(-1.4320089) q[0];
x q[1];
rz(-0.090473526) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(-0.59404101) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0253804) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(-1.2799353) q[1];
rz(-2.071322) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(-0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(2.9764552) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(2.3388458) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0695641) q[0];
sx q[0];
rz(-1.1840128) q[0];
sx q[0];
rz(-0.2441498) q[0];
x q[1];
rz(-0.10404189) q[2];
sx q[2];
rz(-2.0718956) q[2];
sx q[2];
rz(0.83545557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3853559) q[1];
sx q[1];
rz(-1.5192351) q[1];
sx q[1];
rz(0.95183422) q[1];
rz(3.0117399) q[3];
sx q[3];
rz(-2.2359071) q[3];
sx q[3];
rz(-0.065404281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.8018988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1012672) q[0];
sx q[0];
rz(-0.87653941) q[0];
sx q[0];
rz(-0.56157748) q[0];
x q[1];
rz(1.9402903) q[2];
sx q[2];
rz(-0.72417799) q[2];
sx q[2];
rz(0.61831123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9078319) q[1];
sx q[1];
rz(-1.3843752) q[1];
sx q[1];
rz(-2.9790661) q[1];
rz(-pi) q[2];
rz(0.6891059) q[3];
sx q[3];
rz(-1.2561241) q[3];
sx q[3];
rz(0.93022197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4703935) q[2];
sx q[2];
rz(-0.931804) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-2.4750211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4857805) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(0.73894545) q[0];
rz(-pi) q[1];
rz(2.4273657) q[2];
sx q[2];
rz(-1.3492609) q[2];
sx q[2];
rz(-0.54745882) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8909292) q[1];
sx q[1];
rz(-2.2729985) q[1];
sx q[1];
rz(2.3198747) q[1];
x q[2];
rz(-1.4227082) q[3];
sx q[3];
rz(-1.423096) q[3];
sx q[3];
rz(-1.577391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(0.83706013) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.9821092) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5975208) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(1.9835299) q[0];
rz(-pi) q[1];
rz(-2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(1.3364524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3573208) q[1];
sx q[1];
rz(-1.1457774) q[1];
sx q[1];
rz(-0.97370633) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.106835) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(-2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(1.2716028) q[2];
rz(-2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(3.040722) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8773008) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(2.507188) q[0];
x q[1];
rz(2.6542873) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(-3.0234408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74360352) q[1];
sx q[1];
rz(-1.7730224) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50562596) q[3];
sx q[3];
rz(-2.6856832) q[3];
sx q[3];
rz(1.0491919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15935414) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(0.23157816) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2528) q[2];
sx q[2];
rz(-1.3557634) q[2];
sx q[2];
rz(-3.0766578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75377611) q[1];
sx q[1];
rz(-3.0101335) q[1];
sx q[1];
rz(0.17863518) q[1];
rz(-pi) q[2];
rz(2.1391159) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(2.3426567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(0.93758279) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(1.1212564) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(1.1987196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872902) q[0];
sx q[0];
rz(-1.1371005) q[0];
sx q[0];
rz(-1.646423) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.347581) q[2];
sx q[2];
rz(-2.7701021) q[2];
sx q[2];
rz(1.5944634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3938013) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(1.5138813) q[1];
rz(-pi) q[2];
rz(-0.49097455) q[3];
sx q[3];
rz(-2.6262865) q[3];
sx q[3];
rz(1.5332424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.2808799) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(2.5532706) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(-0.9334329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8342455) q[0];
sx q[0];
rz(-1.630162) q[0];
sx q[0];
rz(0.10660118) q[0];
x q[1];
rz(1.1912212) q[2];
sx q[2];
rz(-1.8827202) q[2];
sx q[2];
rz(-1.003007) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0837299) q[1];
sx q[1];
rz(-1.607556) q[1];
sx q[1];
rz(-2.4576393) q[1];
rz(-pi) q[2];
rz(-0.65255717) q[3];
sx q[3];
rz(-1.8592632) q[3];
sx q[3];
rz(-2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(2.771647) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(2.9577589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606124) q[0];
sx q[0];
rz(-1.6913337) q[0];
sx q[0];
rz(1.4535849) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76639558) q[2];
sx q[2];
rz(-2.1470214) q[2];
sx q[2];
rz(3.0255084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(-0.88840719) q[1];
rz(-pi) q[2];
x q[2];
rz(0.023601836) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(-2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(-0.20239057) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595173) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-0.45869641) q[2];
sx q[2];
rz(-1.9445322) q[2];
sx q[2];
rz(-1.4945488) q[2];
rz(-0.17584569) q[3];
sx q[3];
rz(-0.98777117) q[3];
sx q[3];
rz(-2.5143757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
