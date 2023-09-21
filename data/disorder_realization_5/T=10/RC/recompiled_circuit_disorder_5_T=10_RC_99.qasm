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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306643) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(1.7095837) q[0];
x q[1];
rz(1.4957501) q[2];
sx q[2];
rz(-1.6610166) q[2];
sx q[2];
rz(-2.1580632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.773355) q[1];
sx q[1];
rz(-3.0800746) q[1];
x q[2];
rz(-2.7492417) q[3];
sx q[3];
rz(-0.96071834) q[3];
sx q[3];
rz(3.0685134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(2.9764552) q[2];
rz(-1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.76666) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4883603) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(-1.0351719) q[0];
x q[1];
rz(-3.0375508) q[2];
sx q[2];
rz(-2.0718956) q[2];
sx q[2];
rz(-0.83545557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3992577) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-pi) q[2];
rz(1.734415) q[3];
sx q[3];
rz(-2.46582) q[3];
sx q[3];
rz(2.8676652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(-0.60570335) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(1.8018988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0403255) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(-2.5800152) q[0];
rz(-pi) q[1];
rz(0.3091829) q[2];
sx q[2];
rz(-0.90484607) q[2];
sx q[2];
rz(-0.14112976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(-1.7596485) q[1];
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
rz(-0.931804) q[2];
sx q[2];
rz(1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(-0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-2.4750211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65581215) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(0.73894545) q[0];
rz(-1.8604943) q[2];
sx q[2];
rz(-0.87756598) q[2];
sx q[2];
rz(1.9300269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2143612) q[1];
sx q[1];
rz(-2.1640558) q[1];
sx q[1];
rz(-2.4638922) q[1];
rz(-1.4227082) q[3];
sx q[3];
rz(-1.423096) q[3];
sx q[3];
rz(1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-2.58113) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(-1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5975208) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(1.1580628) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3443089) q[2];
sx q[2];
rz(-2.531257) q[2];
sx q[2];
rz(3.0938873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0817464) q[1];
sx q[1];
rz(-1.0330327) q[1];
sx q[1];
rz(0.50077011) q[1];
rz(-pi) q[2];
rz(2.106835) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(-1.1353726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(3.040722) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6504165) q[0];
sx q[0];
rz(-2.4917779) q[0];
sx q[0];
rz(2.8894436) q[0];
rz(1.8632554) q[2];
sx q[2];
rz(-1.1012226) q[2];
sx q[2];
rz(1.8243607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.069177376) q[1];
sx q[1];
rz(-2.8848856) q[1];
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
rz(1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-2.1784901) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15935414) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(-0.23157816) q[0];
rz(-pi) q[1];
rz(-0.96104788) q[2];
sx q[2];
rz(-2.7597722) q[2];
sx q[2];
rz(-2.080999) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93393275) q[1];
sx q[1];
rz(-1.4414409) q[1];
sx q[1];
rz(-1.5473066) q[1];
x q[2];
rz(0.01597605) q[3];
sx q[3];
rz(-1.0025347) q[3];
sx q[3];
rz(0.78045867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(0.93758279) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.2032578) q[1];
sx q[1];
rz(1.942873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543025) q[0];
sx q[0];
rz(-2.0044921) q[0];
sx q[0];
rz(1.646423) q[0];
rz(-pi) q[1];
x q[1];
rz(1.347581) q[2];
sx q[2];
rz(-2.7701021) q[2];
sx q[2];
rz(-1.5944634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.291154) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(-2.6384141) q[1];
rz(-pi) q[2];
rz(2.6783887) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(2.6686058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(-0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.566074) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(2.2081597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76970869) q[0];
sx q[0];
rz(-0.121962) q[0];
sx q[0];
rz(-2.6321649) q[0];
rz(-pi) q[1];
rz(1.9503715) q[2];
sx q[2];
rz(-1.2588725) q[2];
sx q[2];
rz(2.1385857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(-3.0834553) q[1];
rz(1.9282354) q[3];
sx q[3];
rz(-2.1921428) q[3];
sx q[3];
rz(2.3693717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(-2.771647) q[2];
rz(-0.89921078) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-2.9577589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659347) q[0];
sx q[0];
rz(-1.6871534) q[0];
sx q[0];
rz(3.0202306) q[0];
rz(-pi) q[1];
rz(-0.76639558) q[2];
sx q[2];
rz(-2.1470214) q[2];
sx q[2];
rz(-3.0255084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76845968) q[1];
sx q[1];
rz(-0.8994973) q[1];
sx q[1];
rz(0.21233227) q[1];
rz(-0.58335431) q[3];
sx q[3];
rz(-1.5837985) q[3];
sx q[3];
rz(-2.0875723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
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
rz(-2.1609859) q[3];
sx q[3];
rz(-1.7173613) q[3];
sx q[3];
rz(2.2955256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
