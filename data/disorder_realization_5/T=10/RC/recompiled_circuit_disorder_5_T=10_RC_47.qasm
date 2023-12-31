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
rz(-1.8656123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39014434) q[0];
sx q[0];
rz(-1.7062618) q[0];
sx q[0];
rz(2.9209903) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4495507) q[2];
sx q[2];
rz(-3.0243052) q[2];
sx q[2];
rz(2.8534944) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1162122) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(-1.8616574) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2184578) q[3];
sx q[3];
rz(-1.8895518) q[3];
sx q[3];
rz(-1.8766599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(-1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(-2.3388458) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.734056) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(-1.1733823) q[0];
rz(-1.3834125) q[2];
sx q[2];
rz(-2.6307081) q[2];
sx q[2];
rz(-2.5201706) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2903366) q[1];
sx q[1];
rz(-0.95278059) q[1];
sx q[1];
rz(0.063277146) q[1];
rz(-pi) q[2];
rz(0.12985274) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(-0.60570335) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.8018988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08683603) q[0];
sx q[0];
rz(-1.9924379) q[0];
sx q[0];
rz(-0.7936759) q[0];
rz(-pi) q[1];
rz(2.8324098) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(-0.14112976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(-1.3819441) q[1];
x q[2];
rz(-0.6891059) q[3];
sx q[3];
rz(-1.8854685) q[3];
sx q[3];
rz(-2.2113707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4703935) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(2.9273422) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.727107) q[0];
sx q[0];
rz(-1.7761199) q[0];
sx q[0];
rz(1.3798825) q[0];
x q[1];
rz(2.4273657) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(0.54745882) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2143612) q[1];
sx q[1];
rz(-0.9775369) q[1];
sx q[1];
rz(-0.67770047) q[1];
rz(-2.3603667) q[3];
sx q[3];
rz(-2.9328212) q[3];
sx q[3];
rz(-2.3695932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(2.3045325) q[0];
rz(2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(-1.9821092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54407185) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(-1.9835299) q[0];
x q[1];
rz(1.7972838) q[2];
sx q[2];
rz(-2.531257) q[2];
sx q[2];
rz(0.047705334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3573208) q[1];
sx q[1];
rz(-1.9958152) q[1];
sx q[1];
rz(-0.97370633) q[1];
rz(-pi) q[2];
rz(-1.0347576) q[3];
sx q[3];
rz(-1.0039819) q[3];
sx q[3];
rz(-2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(1.2716028) q[2];
rz(-2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(-3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-3.040722) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2642919) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(0.63440462) q[0];
rz(2.6249044) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(1.9043497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.069177376) q[1];
sx q[1];
rz(-0.25670708) q[1];
sx q[1];
rz(0.91255811) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8039861) q[3];
sx q[3];
rz(-1.1753852) q[3];
sx q[3];
rz(-2.6449441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(1.2290918) q[0];
x q[1];
rz(-0.22600941) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(1.7058536) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63989418) q[1];
sx q[1];
rz(-1.5940897) q[1];
sx q[1];
rz(-0.1293907) q[1];
rz(-pi) q[2];
rz(-0.01597605) q[3];
sx q[3];
rz(-1.0025347) q[3];
sx q[3];
rz(-0.78045867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(0.93758279) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(-2.3576095) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(-1.942873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90889701) q[0];
sx q[0];
rz(-0.43982738) q[0];
sx q[0];
rz(-0.16172414) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7940117) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(-1.5471293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7477914) q[1];
sx q[1];
rz(-1.0683021) q[1];
sx q[1];
rz(-1.5138813) q[1];
x q[2];
rz(0.46320398) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(-0.9334329) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8342455) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(3.0349915) q[0];
rz(-pi) q[1];
rz(-2.2869296) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(-1.2235299) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(3.0834553) q[1];
x q[2];
rz(1.9282354) q[3];
sx q[3];
rz(-0.94944984) q[3];
sx q[3];
rz(-2.3693717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
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
rz(pi/2) q[3];
sx q[3];
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
rz(-1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-0.18383372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1856954) q[0];
sx q[0];
rz(-2.9736608) q[0];
sx q[0];
rz(-2.3737565) q[0];
rz(-2.3751971) q[2];
sx q[2];
rz(-0.99457127) q[2];
sx q[2];
rz(-3.0255084) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(-0.88840719) q[1];
x q[2];
rz(1.5863745) q[3];
sx q[3];
rz(-2.1540949) q[3];
sx q[3];
rz(2.616236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68207537) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(1.1584875) q[2];
sx q[2];
rz(-1.145913) q[2];
sx q[2];
rz(-0.10213012) q[2];
rz(-1.8300874) q[3];
sx q[3];
rz(-2.5355831) q[3];
sx q[3];
rz(-2.2021962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
