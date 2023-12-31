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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5029992) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(-2.5845085) q[0];
rz(1.6458426) q[2];
sx q[2];
rz(-1.6610166) q[2];
sx q[2];
rz(2.1580632) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0253804) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(1.8616574) q[1];
rz(2.2184578) q[3];
sx q[3];
rz(-1.8895518) q[3];
sx q[3];
rz(-1.8766599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(-1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(3.0965565) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6532324) q[0];
sx q[0];
rz(-0.45408861) q[0];
sx q[0];
rz(1.0351719) q[0];
rz(1.0674092) q[2];
sx q[2];
rz(-1.6620087) q[2];
sx q[2];
rz(0.78546055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2903366) q[1];
sx q[1];
rz(-2.1888121) q[1];
sx q[1];
rz(0.063277146) q[1];
rz(-2.2400168) q[3];
sx q[3];
rz(-1.4687317) q[3];
sx q[3];
rz(-1.7166184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(-0.60570335) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95214343) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-0.28513518) q[0];
rz(-2.6248698) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(-1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0547566) q[0];
sx q[0];
rz(-1.9924379) q[0];
sx q[0];
rz(-2.3479168) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88111934) q[2];
sx q[2];
rz(-1.8124049) q[2];
sx q[2];
rz(1.2348262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(1.7596485) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4524868) q[3];
sx q[3];
rz(-1.8854685) q[3];
sx q[3];
rz(-2.2113707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(1.7335256) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(2.9273422) q[0];
rz(-0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-2.4750211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9458981) q[0];
sx q[0];
rz(-1.7576522) q[0];
sx q[0];
rz(-0.20901434) q[0];
x q[1];
rz(-2.8104086) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(0.77510288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9208593) q[1];
sx q[1];
rz(-2.1174866) q[1];
sx q[1];
rz(0.85733719) q[1];
rz(2.3603667) q[3];
sx q[3];
rz(-2.9328212) q[3];
sx q[3];
rz(-0.77199948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(1.9541698) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-2.3045325) q[0];
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
rz(-2.5975208) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(1.1580628) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(1.3364524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0598462) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(-2.6408225) q[1];
rz(-pi) q[2];
rz(1.0347576) q[3];
sx q[3];
rz(-1.0039819) q[3];
sx q[3];
rz(2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7065457) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-0.10087068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6504165) q[0];
sx q[0];
rz(-2.4917779) q[0];
sx q[0];
rz(0.25214904) q[0];
rz(-2.6542873) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(3.0234408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.069177376) q[1];
sx q[1];
rz(-0.25670708) q[1];
sx q[1];
rz(2.2290345) q[1];
x q[2];
rz(-2.7363339) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(2.1586771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(0.96310258) q[0];
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
rz(-1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(-1.9125008) q[0];
rz(-pi) q[1];
rz(1.8887927) q[2];
sx q[2];
rz(-1.7858292) q[2];
sx q[2];
rz(-3.0766578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2076599) q[1];
sx q[1];
rz(-1.4414409) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(2.3576095) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.942873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90889701) q[0];
sx q[0];
rz(-0.43982738) q[0];
sx q[0];
rz(0.16172414) q[0];
rz(-0.086026831) q[2];
sx q[2];
rz(-1.2089529) q[2];
sx q[2];
rz(-1.8334243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3938013) q[1];
sx q[1];
rz(-1.0683021) q[1];
sx q[1];
rz(-1.5138813) q[1];
rz(-pi) q[2];
rz(-1.8317322) q[3];
sx q[3];
rz(-2.0203777) q[3];
sx q[3];
rz(-0.98228067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-0.58832204) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(0.9334329) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25710042) q[0];
sx q[0];
rz(-1.677209) q[0];
sx q[0];
rz(-1.5110925) q[0];
rz(-0.8546631) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(-1.9180627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(-0.058137356) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9282354) q[3];
sx q[3];
rz(-2.1921428) q[3];
sx q[3];
rz(0.77222094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(-2.771647) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(-0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(2.9577589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95589721) q[0];
sx q[0];
rz(-0.16793185) q[0];
sx q[0];
rz(-0.76783617) q[0];
x q[1];
rz(0.76639558) q[2];
sx q[2];
rz(-0.99457127) q[2];
sx q[2];
rz(0.11608427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7067691) q[1];
sx q[1];
rz(-2.4425047) q[1];
sx q[1];
rz(-1.3114776) q[1];
rz(-1.5552181) q[3];
sx q[3];
rz(-0.9874978) q[3];
sx q[3];
rz(-2.616236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-2.9392021) q[2];
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
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595173) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(2.7783685) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(-2.6828962) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(1.6470439) q[2];
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
