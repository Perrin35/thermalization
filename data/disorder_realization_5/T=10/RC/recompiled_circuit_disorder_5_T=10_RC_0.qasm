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
rz(1.9306643) q[0];
sx q[0];
rz(-1.7893447) q[0];
sx q[0];
rz(-1.4320089) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6458426) q[2];
sx q[2];
rz(-1.4805761) q[2];
sx q[2];
rz(-0.98352945) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7282655) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(3.0800746) q[1];
x q[2];
rz(-0.92313487) q[3];
sx q[3];
rz(-1.8895518) q[3];
sx q[3];
rz(-1.8766599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(-0.80274686) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0720285) q[0];
sx q[0];
rz(-1.1840128) q[0];
sx q[0];
rz(0.2441498) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10404189) q[2];
sx q[2];
rz(-2.0718956) q[2];
sx q[2];
rz(-2.3061371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3992577) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90157585) q[3];
sx q[3];
rz(-1.672861) q[3];
sx q[3];
rz(-1.4249742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(-1.8018988) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547566) q[0];
sx q[0];
rz(-1.1491547) q[0];
sx q[0];
rz(0.7936759) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2013024) q[2];
sx q[2];
rz(-2.4174147) q[2];
sx q[2];
rz(2.5232814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9078319) q[1];
sx q[1];
rz(-1.3843752) q[1];
sx q[1];
rz(2.9790661) q[1];
x q[2];
rz(0.47311802) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(-0.28120041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4703935) q[2];
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
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-2.4750211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4857805) q[0];
sx q[0];
rz(-0.27944767) q[0];
sx q[0];
rz(-0.73894545) q[0];
x q[1];
rz(0.33118403) q[2];
sx q[2];
rz(-0.74197717) q[2];
sx q[2];
rz(-0.77510288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22073332) q[1];
sx q[1];
rz(-2.1174866) q[1];
sx q[1];
rz(2.2842555) q[1];
x q[2];
rz(1.4227082) q[3];
sx q[3];
rz(-1.423096) q[3];
sx q[3];
rz(1.577391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-2.4812223) q[2];
rz(1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-0.66851139) q[1];
sx q[1];
rz(1.9821092) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3876622) q[0];
sx q[0];
rz(-1.2571063) q[0];
sx q[0];
rz(0.73648209) q[0];
rz(-1.3443089) q[2];
sx q[2];
rz(-2.531257) q[2];
sx q[2];
rz(-3.0938873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0817464) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(2.6408225) q[1];
rz(0.67622185) q[3];
sx q[3];
rz(-2.3822504) q[3];
sx q[3];
rz(-2.8420705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(0.043958157) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(0.10087068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80412241) q[0];
sx q[0];
rz(-0.94479783) q[0];
sx q[0];
rz(-1.3834329) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6542873) q[2];
sx q[2];
rz(-1.8308182) q[2];
sx q[2];
rz(-3.0234408) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74360352) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-pi) q[2];
rz(-1.8039861) q[3];
sx q[3];
rz(-1.1753852) q[3];
sx q[3];
rz(-2.6449441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(-1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46524099) q[0];
sx q[0];
rz(-2.7382593) q[0];
sx q[0];
rz(2.1562735) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2528) q[2];
sx q[2];
rz(-1.7858292) q[2];
sx q[2];
rz(-3.0766578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5016985) q[1];
sx q[1];
rz(-1.5475029) q[1];
sx q[1];
rz(-3.0122019) q[1];
rz(-pi) q[2];
rz(1.0024768) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.942873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51533651) q[0];
sx q[0];
rz(-1.5021828) q[0];
sx q[0];
rz(-2.7068044) q[0];
rz(-pi) q[1];
rz(1.347581) q[2];
sx q[2];
rz(-2.7701021) q[2];
sx q[2];
rz(-1.5944634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.8504387) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(-2.6384141) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25710042) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(1.6305001) q[0];
x q[1];
rz(1.9503715) q[2];
sx q[2];
rz(-1.8827202) q[2];
sx q[2];
rz(1.003007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(0.058137356) q[1];
rz(0.45456072) q[3];
sx q[3];
rz(-2.436736) q[3];
sx q[3];
rz(-2.9398033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(2.771647) q[2];
rz(0.89921078) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-0.18383372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856954) q[0];
sx q[0];
rz(-2.9736608) q[0];
sx q[0];
rz(-0.76783617) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75283639) q[2];
sx q[2];
rz(-2.2194153) q[2];
sx q[2];
rz(2.2019049) q[2];
rz(-pi) q[3];
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
rz(-1.7365672) q[1];
sx q[1];
rz(2.2531855) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5582383) q[3];
sx q[3];
rz(-1.5577941) q[3];
sx q[3];
rz(1.0540203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-0.72485483) q[2];
sx q[2];
rz(-2.5584494) q[2];
sx q[2];
rz(0.7128788) q[2];
rz(2.1609859) q[3];
sx q[3];
rz(-1.4242314) q[3];
sx q[3];
rz(-0.84606708) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];