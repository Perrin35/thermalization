OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39014434) q[0];
sx q[0];
rz(-1.7062618) q[0];
sx q[0];
rz(2.9209903) q[0];
rz(-2.4495507) q[2];
sx q[2];
rz(-3.0243052) q[2];
sx q[2];
rz(0.28809822) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16986019) q[1];
sx q[1];
rz(-1.6310551) q[1];
sx q[1];
rz(1.7737284) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2184578) q[3];
sx q[3];
rz(-1.2520408) q[3];
sx q[3];
rz(-1.2649328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052658) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(-2.8024407) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6532324) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(1.0351719) q[0];
x q[1];
rz(1.3834125) q[2];
sx q[2];
rz(-2.6307081) q[2];
sx q[2];
rz(-0.62142205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74233494) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-pi) q[2];
rz(-3.0117399) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(-0.065404281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(1.8018988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08683603) q[0];
sx q[0];
rz(-1.1491547) q[0];
sx q[0];
rz(-0.7936759) q[0];
rz(2.8324098) q[2];
sx q[2];
rz(-0.90484607) q[2];
sx q[2];
rz(0.14112976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2337607) q[1];
sx q[1];
rz(-1.3843752) q[1];
sx q[1];
rz(-2.9790661) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4524868) q[3];
sx q[3];
rz(-1.8854685) q[3];
sx q[3];
rz(0.93022197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(3.0696707) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.727107) q[0];
sx q[0];
rz(-1.7761199) q[0];
sx q[0];
rz(-1.3798825) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71422691) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(0.54745882) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(-0.82171792) q[1];
rz(-pi) q[2];
rz(0.14931071) q[3];
sx q[3];
rz(-1.7172604) q[3];
sx q[3];
rz(-0.028545054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-2.4812223) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(0.83706013) q[0];
rz(-2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.9821092) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54407185) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(-1.9835299) q[0];
x q[1];
rz(1.3443089) q[2];
sx q[2];
rz(-2.531257) q[2];
sx q[2];
rz(-0.047705334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0598462) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(2.6408225) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.106835) q[3];
sx q[3];
rz(-1.0039819) q[3];
sx q[3];
rz(-1.1353726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5669252) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-3.040722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4911762) q[0];
sx q[0];
rz(-0.6498148) q[0];
sx q[0];
rz(-0.25214904) q[0];
rz(-2.6542873) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(-0.11815182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3979891) q[1];
sx q[1];
rz(-1.7730224) q[1];
sx q[1];
rz(-0.15921758) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50562596) q[3];
sx q[3];
rz(-0.45590948) q[3];
sx q[3];
rz(-2.0924007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.5161783) q[2];
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
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(-1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(-0.68774736) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(1.9125008) q[0];
rz(-pi) q[1];
rz(1.2528) q[2];
sx q[2];
rz(-1.7858292) q[2];
sx q[2];
rz(3.0766578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.63989418) q[1];
sx q[1];
rz(-1.5940897) q[1];
sx q[1];
rz(-3.0122019) q[1];
rz(-1.5958105) q[3];
sx q[3];
rz(-0.56846148) q[3];
sx q[3];
rz(0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(1.1987196) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2326956) q[0];
sx q[0];
rz(-2.7017653) q[0];
sx q[0];
rz(0.16172414) q[0];
x q[1];
rz(-1.2077246) q[2];
sx q[2];
rz(-1.6512401) q[2];
sx q[2];
rz(0.23210873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7477914) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(-1.6277114) q[1];
rz(-0.49097455) q[3];
sx q[3];
rz(-0.51530616) q[3];
sx q[3];
rz(1.6083503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(2.5532706) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(0.9334329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844922) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(1.5110925) q[0];
rz(-pi) q[1];
rz(-2.2869296) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(1.9180627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0837299) q[1];
sx q[1];
rz(-1.5340367) q[1];
sx q[1];
rz(0.68395331) q[1];
rz(0.65255717) q[3];
sx q[3];
rz(-1.8592632) q[3];
sx q[3];
rz(-1.0126589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(-0.36994568) q[2];
rz(0.89921078) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(2.488234) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-0.18383372) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18098022) q[0];
sx q[0];
rz(-1.6913337) q[0];
sx q[0];
rz(-1.4535849) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75283639) q[2];
sx q[2];
rz(-2.2194153) q[2];
sx q[2];
rz(-0.93968771) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93563423) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(0.88840719) q[1];
x q[2];
rz(0.58335431) q[3];
sx q[3];
rz(-1.5577941) q[3];
sx q[3];
rz(-2.0875723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-2.9392021) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.9831052) q[2];
sx q[2];
rz(-1.9956797) q[2];
sx q[2];
rz(3.0394625) q[2];
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