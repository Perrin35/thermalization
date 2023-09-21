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
rz(2.5198088) q[1];
sx q[1];
rz(3.8222651) q[1];
sx q[1];
rz(8.1487976) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63859343) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(-0.55708416) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6458426) q[2];
sx q[2];
rz(-1.4805761) q[2];
sx q[2];
rz(-0.98352945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7282655) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(3.0800746) q[1];
x q[2];
rz(0.39235093) q[3];
sx q[3];
rz(-2.1808743) q[3];
sx q[3];
rz(0.073079212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(2.9764552) q[2];
rz(-1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(-2.3388458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4883603) q[0];
sx q[0];
rz(-0.45408861) q[0];
sx q[0];
rz(-1.0351719) q[0];
x q[1];
rz(-1.3834125) q[2];
sx q[2];
rz(-0.51088453) q[2];
sx q[2];
rz(-0.62142205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85125605) q[1];
sx q[1];
rz(-2.1888121) q[1];
sx q[1];
rz(-0.063277146) q[1];
x q[2];
rz(-1.4071776) q[3];
sx q[3];
rz(-2.46582) q[3];
sx q[3];
rz(-0.27392745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7654483) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.7155898) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(-2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.3396938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1012672) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(0.56157748) q[0];
x q[1];
rz(-0.88111934) q[2];
sx q[2];
rz(-1.8124049) q[2];
sx q[2];
rz(1.2348262) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(1.7596485) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9698824) q[3];
sx q[3];
rz(-0.92150021) q[3];
sx q[3];
rz(-2.2513575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(1.408067) q[2];
rz(-0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-3.0696707) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-2.4750211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9458981) q[0];
sx q[0];
rz(-1.7576522) q[0];
sx q[0];
rz(0.20901434) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8604943) q[2];
sx q[2];
rz(-2.2640267) q[2];
sx q[2];
rz(1.9300269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22073332) q[1];
sx q[1];
rz(-1.024106) q[1];
sx q[1];
rz(-2.2842555) q[1];
rz(-0.14931071) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-0.028545054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(-2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.9821092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54407185) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(-1.1580628) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7972838) q[2];
sx q[2];
rz(-2.531257) q[2];
sx q[2];
rz(0.047705334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78427181) q[1];
sx q[1];
rz(-1.1457774) q[1];
sx q[1];
rz(-0.97370633) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67622185) q[3];
sx q[3];
rz(-2.3822504) q[3];
sx q[3];
rz(0.29952213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(0.10087068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4911762) q[0];
sx q[0];
rz(-2.4917779) q[0];
sx q[0];
rz(-2.8894436) q[0];
x q[1];
rz(-2.6542873) q[2];
sx q[2];
rz(-1.8308182) q[2];
sx q[2];
rz(0.11815182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74360352) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40525873) q[3];
sx q[3];
rz(-1.7856981) q[3];
sx q[3];
rz(0.98291558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-2.2163088) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46524099) q[0];
sx q[0];
rz(-2.7382593) q[0];
sx q[0];
rz(-2.1562735) q[0];
rz(0.96104788) q[2];
sx q[2];
rz(-0.38182048) q[2];
sx q[2];
rz(-2.080999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63989418) q[1];
sx q[1];
rz(-1.5940897) q[1];
sx q[1];
rz(-0.1293907) q[1];
x q[2];
rz(3.1256166) q[3];
sx q[3];
rz(-1.0025347) q[3];
sx q[3];
rz(2.361134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(0.93758279) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872902) q[0];
sx q[0];
rz(-2.0044921) q[0];
sx q[0];
rz(1.4951697) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(1.5471293) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3938013) q[1];
sx q[1];
rz(-1.0683021) q[1];
sx q[1];
rz(1.6277114) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8317322) q[3];
sx q[3];
rz(-1.121215) q[3];
sx q[3];
rz(-2.159312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(0.9334329) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.371884) q[0];
sx q[0];
rz(-3.0196307) q[0];
sx q[0];
rz(2.6321649) q[0];
rz(-0.33414267) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(2.6956218) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0578627) q[1];
sx q[1];
rz(-1.5340367) q[1];
sx q[1];
rz(-2.4576393) q[1];
x q[2];
rz(-2.4890355) q[3];
sx q[3];
rz(-1.8592632) q[3];
sx q[3];
rz(-1.0126589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(-0.36994568) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-2.1883011) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606124) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(-1.6880077) q[0];
rz(-pi) q[1];
rz(-0.75283639) q[2];
sx q[2];
rz(-0.92217731) q[2];
sx q[2];
rz(-2.2019049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76845968) q[1];
sx q[1];
rz(-2.2420954) q[1];
sx q[1];
rz(0.21233227) q[1];
x q[2];
rz(0.023601836) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-2.9392021) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(-1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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