OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(4.7441626) q[0];
sx q[0];
rz(6.8466853) q[0];
sx q[0];
rz(9.8817649) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(-1.8656123) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63859343) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(2.5845085) q[0];
rz(-pi) q[1];
rz(1.4957501) q[2];
sx q[2];
rz(-1.4805761) q[2];
sx q[2];
rz(2.1580632) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.773355) q[1];
sx q[1];
rz(-3.0800746) q[1];
rz(-pi) q[2];
rz(-2.071322) q[3];
sx q[3];
rz(-2.4300017) q[3];
sx q[3];
rz(-2.4430032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(-0.16513744) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0720285) q[0];
sx q[0];
rz(-1.1840128) q[0];
sx q[0];
rz(-0.2441498) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0741834) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(2.3561321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75623679) q[1];
sx q[1];
rz(-1.6223575) q[1];
sx q[1];
rz(-2.1897584) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12985274) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(-3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.95214343) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.8018988) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2633318) q[0];
sx q[0];
rz(-0.86255951) q[0];
sx q[0];
rz(1.0017298) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88111934) q[2];
sx q[2];
rz(-1.8124049) q[2];
sx q[2];
rz(1.2348262) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(-1.3819441) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47311802) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(0.66657153) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19569451) q[0];
sx q[0];
rz(-1.3839405) q[0];
sx q[0];
rz(-0.20901434) q[0];
x q[1];
rz(0.33118403) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(0.77510288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2143612) q[1];
sx q[1];
rz(-0.9775369) q[1];
sx q[1];
rz(-0.67770047) q[1];
rz(2.3603667) q[3];
sx q[3];
rz(-0.20877148) q[3];
sx q[3];
rz(0.77199948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8551222) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(2.3045325) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.1594835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54407185) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(1.1580628) q[0];
rz(-pi) q[1];
rz(1.7972838) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(3.0938873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0817464) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(-0.50077011) q[1];
x q[2];
rz(-1.0347576) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7065457) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-0.10087068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8773008) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(2.507188) q[0];
rz(-pi) q[1];
rz(2.6249044) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(-1.237243) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85943356) q[1];
sx q[1];
rz(-1.4148501) q[1];
sx q[1];
rz(1.7755416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8039861) q[3];
sx q[3];
rz(-1.1753852) q[3];
sx q[3];
rz(2.6449441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15935414) q[0];
sx q[0];
rz(-1.9040477) q[0];
sx q[0];
rz(-0.23157816) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96104788) q[2];
sx q[2];
rz(-0.38182048) q[2];
sx q[2];
rz(2.080999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3878165) q[1];
sx q[1];
rz(-0.13145914) q[1];
sx q[1];
rz(-0.17863518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0024768) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(3.0623073) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(1.942873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51533651) q[0];
sx q[0];
rz(-1.6394098) q[0];
sx q[0];
rz(0.4347883) q[0];
rz(-pi) q[1];
rz(1.7940117) q[2];
sx q[2];
rz(-2.7701021) q[2];
sx q[2];
rz(-1.5471293) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3938013) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(1.5138813) q[1];
x q[2];
rz(-0.46320398) q[3];
sx q[3];
rz(-1.3363046) q[3];
sx q[3];
rz(0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(-0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3073472) q[0];
sx q[0];
rz(-1.630162) q[0];
sx q[0];
rz(-0.10660118) q[0];
x q[1];
rz(-2.80745) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(-2.6956218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(3.0834553) q[1];
x q[2];
rz(-1.2133573) q[3];
sx q[3];
rz(-0.94944984) q[3];
sx q[3];
rz(-2.3693717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(0.36994568) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(0.95329154) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(2.488234) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1856954) q[0];
sx q[0];
rz(-0.16793185) q[0];
sx q[0];
rz(2.3737565) q[0];
x q[1];
rz(0.76639558) q[2];
sx q[2];
rz(-0.99457127) q[2];
sx q[2];
rz(0.11608427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4348236) q[1];
sx q[1];
rz(-2.4425047) q[1];
sx q[1];
rz(-1.8301151) q[1];
rz(2.5582383) q[3];
sx q[3];
rz(-1.5577941) q[3];
sx q[3];
rz(2.0875723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.4595173) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(0.72485483) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(2.965747) q[3];
sx q[3];
rz(-0.98777117) q[3];
sx q[3];
rz(-2.5143757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
