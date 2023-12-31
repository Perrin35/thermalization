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
rz(1.2759804) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306643) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(1.7095837) q[0];
rz(-pi) q[1];
rz(3.0511191) q[2];
sx q[2];
rz(-1.6455368) q[2];
sx q[2];
rz(-2.5475516) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(3.0800746) q[1];
rz(-pi) q[2];
rz(1.0702707) q[3];
sx q[3];
rz(-2.4300017) q[3];
sx q[3];
rz(0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(2.9764552) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(2.3993649) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40753663) q[0];
sx q[0];
rz(-1.7965839) q[0];
sx q[0];
rz(-1.9682103) q[0];
rz(-pi) q[1];
rz(-1.7581802) q[2];
sx q[2];
rz(-0.51088453) q[2];
sx q[2];
rz(0.62142205) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3992577) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-pi) q[2];
rz(0.12985274) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95214343) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-0.28513518) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1012672) q[0];
sx q[0];
rz(-0.87653941) q[0];
sx q[0];
rz(-2.5800152) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2013024) q[2];
sx q[2];
rz(-2.4174147) q[2];
sx q[2];
rz(2.5232814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(1.3819441) q[1];
x q[2];
rz(-0.47311802) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(0.28120041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.7335256) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(0.66657153) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9458981) q[0];
sx q[0];
rz(-1.3839405) q[0];
sx q[0];
rz(-0.20901434) q[0];
x q[1];
rz(-0.71422691) q[2];
sx q[2];
rz(-1.3492609) q[2];
sx q[2];
rz(-0.54745882) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8909292) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(-0.82171792) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3603667) q[3];
sx q[3];
rz(-2.9328212) q[3];
sx q[3];
rz(-2.3695932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6289571) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.9821092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75393049) q[0];
sx q[0];
rz(-1.8844863) q[0];
sx q[0];
rz(-0.73648209) q[0];
rz(1.3443089) q[2];
sx q[2];
rz(-2.531257) q[2];
sx q[2];
rz(-0.047705334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9002478) q[1];
sx q[1];
rz(-2.4240139) q[1];
sx q[1];
rz(0.89300968) q[1];
x q[2];
rz(2.4653708) q[3];
sx q[3];
rz(-0.75934221) q[3];
sx q[3];
rz(0.29952213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7065457) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(1.8699899) q[2];
rz(1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3374702) q[0];
sx q[0];
rz(-2.1967948) q[0];
sx q[0];
rz(1.3834329) q[0];
rz(-pi) q[1];
rz(2.6249044) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(1.9043497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0724153) q[1];
sx q[1];
rz(-0.25670708) q[1];
sx q[1];
rz(-0.91255811) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40525873) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(-2.1586771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.293175) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(-1.2290918) q[0];
x q[1];
rz(2.1805448) q[2];
sx q[2];
rz(-0.38182048) q[2];
sx q[2];
rz(2.080999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93393275) q[1];
sx q[1];
rz(-1.7001517) q[1];
sx q[1];
rz(-1.594286) q[1];
rz(-0.01597605) q[3];
sx q[3];
rz(-2.139058) q[3];
sx q[3];
rz(-2.361134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.1987196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6262561) q[0];
sx q[0];
rz(-1.5021828) q[0];
sx q[0];
rz(-2.7068044) q[0];
rz(1.7940117) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(-1.5944634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63003507) q[1];
sx q[1];
rz(-0.50543284) q[1];
sx q[1];
rz(3.0384484) q[1];
rz(-2.6506181) q[3];
sx q[3];
rz(-0.51530616) q[3];
sx q[3];
rz(1.5332424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(-2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(0.9334329) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8844922) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(-1.5110925) q[0];
rz(-pi) q[1];
rz(-1.1912212) q[2];
sx q[2];
rz(-1.8827202) q[2];
sx q[2];
rz(-2.1385857) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0837299) q[1];
sx q[1];
rz(-1.607556) q[1];
sx q[1];
rz(-0.68395331) q[1];
rz(-pi) q[2];
rz(0.65255717) q[3];
sx q[3];
rz(-1.8592632) q[3];
sx q[3];
rz(2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.062722) q[2];
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
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(2.9577589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95589721) q[0];
sx q[0];
rz(-0.16793185) q[0];
sx q[0];
rz(2.3737565) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3887563) q[2];
sx q[2];
rz(-0.92217731) q[2];
sx q[2];
rz(2.2019049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(0.88840719) q[1];
x q[2];
rz(-0.58335431) q[3];
sx q[3];
rz(-1.5837985) q[3];
sx q[3];
rz(1.0540203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.7741268) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
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
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(2.6828962) q[2];
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
