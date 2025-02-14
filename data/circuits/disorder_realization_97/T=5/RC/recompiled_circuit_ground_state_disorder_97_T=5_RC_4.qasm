OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.360541) q[0];
sx q[0];
rz(-0.6147576) q[0];
sx q[0];
rz(-1.0714666) q[0];
rz(-0.40845025) q[1];
sx q[1];
rz(-2.2041359) q[1];
sx q[1];
rz(1.6931005) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1030185) q[0];
sx q[0];
rz(-1.9644613) q[0];
sx q[0];
rz(2.066924) q[0];
rz(0.50183588) q[2];
sx q[2];
rz(-2.5060839) q[2];
sx q[2];
rz(-1.0226344) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2768184) q[1];
sx q[1];
rz(-1.4815104) q[1];
sx q[1];
rz(2.9922723) q[1];
rz(-2.2812649) q[3];
sx q[3];
rz(-1.9091879) q[3];
sx q[3];
rz(0.30358728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.19746) q[2];
sx q[2];
rz(-0.84386533) q[2];
sx q[2];
rz(-0.79305631) q[2];
rz(-2.872725) q[3];
sx q[3];
rz(-1.8157418) q[3];
sx q[3];
rz(2.1813006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307513) q[0];
sx q[0];
rz(-2.7590867) q[0];
sx q[0];
rz(2.261396) q[0];
rz(-0.63086069) q[1];
sx q[1];
rz(-0.81268251) q[1];
sx q[1];
rz(-2.0426483) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1504333) q[0];
sx q[0];
rz(-0.17788237) q[0];
sx q[0];
rz(-3.1129897) q[0];
rz(-pi) q[1];
rz(2.1889792) q[2];
sx q[2];
rz(-1.7064377) q[2];
sx q[2];
rz(-0.73276765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6827755) q[1];
sx q[1];
rz(-1.4084245) q[1];
sx q[1];
rz(2.008325) q[1];
rz(2.1707525) q[3];
sx q[3];
rz(-2.4106541) q[3];
sx q[3];
rz(0.39492861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7546996) q[2];
sx q[2];
rz(-1.6563481) q[2];
sx q[2];
rz(2.5999542) q[2];
rz(0.70332876) q[3];
sx q[3];
rz(-0.43916217) q[3];
sx q[3];
rz(-0.37731236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90848732) q[0];
sx q[0];
rz(-0.87924123) q[0];
sx q[0];
rz(3.044686) q[0];
rz(0.58755177) q[1];
sx q[1];
rz(-2.2511626) q[1];
sx q[1];
rz(0.60120916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58734632) q[0];
sx q[0];
rz(-2.4317435) q[0];
sx q[0];
rz(0.27940936) q[0];
x q[1];
rz(1.634609) q[2];
sx q[2];
rz(-1.4087965) q[2];
sx q[2];
rz(-2.8500003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2612348) q[1];
sx q[1];
rz(-2.6812892) q[1];
sx q[1];
rz(-3.1241547) q[1];
rz(-pi) q[2];
rz(-2.2279068) q[3];
sx q[3];
rz(-1.2877712) q[3];
sx q[3];
rz(-1.0237657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2931557) q[2];
sx q[2];
rz(-1.4508672) q[2];
sx q[2];
rz(2.4270774) q[2];
rz(1.4826639) q[3];
sx q[3];
rz(-0.27479333) q[3];
sx q[3];
rz(2.7627435) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42245427) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(-1.5041014) q[0];
rz(-1.0927041) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(-2.5561996) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4248237) q[0];
sx q[0];
rz(-0.94858525) q[0];
sx q[0];
rz(-1.2347925) q[0];
x q[1];
rz(1.4705067) q[2];
sx q[2];
rz(-3.0368251) q[2];
sx q[2];
rz(-0.47276326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3964557) q[1];
sx q[1];
rz(-1.0978699) q[1];
sx q[1];
rz(0.84348444) q[1];
x q[2];
rz(0.62111698) q[3];
sx q[3];
rz(-2.3986772) q[3];
sx q[3];
rz(-0.37534227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66712159) q[2];
sx q[2];
rz(-0.67255628) q[2];
sx q[2];
rz(-2.7453864) q[2];
rz(2.951156) q[3];
sx q[3];
rz(-2.2021144) q[3];
sx q[3];
rz(-2.6207391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025295479) q[0];
sx q[0];
rz(-0.46023661) q[0];
sx q[0];
rz(-0.39475557) q[0];
rz(1.0074298) q[1];
sx q[1];
rz(-1.9837244) q[1];
sx q[1];
rz(1.5347068) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43869492) q[0];
sx q[0];
rz(-0.83897018) q[0];
sx q[0];
rz(1.1870866) q[0];
rz(-pi) q[1];
rz(0.7163064) q[2];
sx q[2];
rz(-1.2898462) q[2];
sx q[2];
rz(1.0208875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6794656) q[1];
sx q[1];
rz(-1.7662342) q[1];
sx q[1];
rz(1.8848757) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0232361) q[3];
sx q[3];
rz(-1.1471738) q[3];
sx q[3];
rz(-2.771559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4452303) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(-0.17808476) q[2];
rz(0.94545025) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(0.43615714) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4431385) q[0];
sx q[0];
rz(-0.66265023) q[0];
sx q[0];
rz(-0.15289256) q[0];
rz(-1.0180417) q[1];
sx q[1];
rz(-2.1986304) q[1];
sx q[1];
rz(2.5315888) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65839889) q[0];
sx q[0];
rz(-0.81419277) q[0];
sx q[0];
rz(1.7479595) q[0];
rz(1.5581846) q[2];
sx q[2];
rz(-2.5676845) q[2];
sx q[2];
rz(-1.2645666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9023561) q[1];
sx q[1];
rz(-0.78241759) q[1];
sx q[1];
rz(-2.8066325) q[1];
rz(-pi) q[2];
rz(-0.097886622) q[3];
sx q[3];
rz(-1.6523419) q[3];
sx q[3];
rz(-1.7718601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6792004) q[2];
sx q[2];
rz(-1.2827164) q[2];
sx q[2];
rz(-2.882615) q[2];
rz(-2.9955043) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(-0.63905382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6530957) q[0];
sx q[0];
rz(-2.275832) q[0];
sx q[0];
rz(2.3897032) q[0];
rz(2.409626) q[1];
sx q[1];
rz(-1.7611793) q[1];
sx q[1];
rz(-0.52925777) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7649496) q[0];
sx q[0];
rz(-1.8762454) q[0];
sx q[0];
rz(-2.9719388) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0257607) q[2];
sx q[2];
rz(-1.6935529) q[2];
sx q[2];
rz(-1.794874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4715035) q[1];
sx q[1];
rz(-1.5836599) q[1];
sx q[1];
rz(1.6520581) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1629754) q[3];
sx q[3];
rz(-2.6088723) q[3];
sx q[3];
rz(-2.5902093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8818714) q[2];
sx q[2];
rz(-2.112381) q[2];
sx q[2];
rz(0.79505801) q[2];
rz(1.2231539) q[3];
sx q[3];
rz(-1.3431679) q[3];
sx q[3];
rz(-0.77643001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.935598) q[0];
sx q[0];
rz(-1.5482276) q[0];
sx q[0];
rz(0.11496168) q[0];
rz(0.089275442) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(0.1549391) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46613559) q[0];
sx q[0];
rz(-2.977006) q[0];
sx q[0];
rz(-2.229722) q[0];
x q[1];
rz(2.5903715) q[2];
sx q[2];
rz(-2.1586426) q[2];
sx q[2];
rz(1.2767222) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.47688) q[1];
sx q[1];
rz(-0.47213337) q[1];
sx q[1];
rz(-2.5499198) q[1];
rz(2.2959034) q[3];
sx q[3];
rz(-1.6250623) q[3];
sx q[3];
rz(0.76322633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(0.6883626) q[2];
rz(2.5734731) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(-2.4585371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883009) q[0];
sx q[0];
rz(-2.330307) q[0];
sx q[0];
rz(1.2055093) q[0];
rz(2.0506809) q[1];
sx q[1];
rz(-2.5542407) q[1];
sx q[1];
rz(1.5376512) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086370416) q[0];
sx q[0];
rz(-1.1366664) q[0];
sx q[0];
rz(-0.40121292) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1386991) q[2];
sx q[2];
rz(-1.5712103) q[2];
sx q[2];
rz(1.3430719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22828211) q[1];
sx q[1];
rz(-2.3582715) q[1];
sx q[1];
rz(1.4128774) q[1];
x q[2];
rz(-1.0991715) q[3];
sx q[3];
rz(-1.0256919) q[3];
sx q[3];
rz(-2.4555912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18498147) q[2];
sx q[2];
rz(-2.0272171) q[2];
sx q[2];
rz(-1.1448917) q[2];
rz(-1.6635118) q[3];
sx q[3];
rz(-2.3749115) q[3];
sx q[3];
rz(-2.0293106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803273) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(1.435745) q[0];
rz(-2.0044633) q[1];
sx q[1];
rz(-2.4500193) q[1];
sx q[1];
rz(0.93708509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1151871) q[0];
sx q[0];
rz(-1.6583024) q[0];
sx q[0];
rz(2.8390272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50855277) q[2];
sx q[2];
rz(-1.1932696) q[2];
sx q[2];
rz(-0.14103954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2236589) q[1];
sx q[1];
rz(-1.5451225) q[1];
sx q[1];
rz(-2.3707546) q[1];
x q[2];
rz(0.61599515) q[3];
sx q[3];
rz(-1.8503891) q[3];
sx q[3];
rz(-3.0699025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81637853) q[2];
sx q[2];
rz(-2.4983695) q[2];
sx q[2];
rz(-2.8312259) q[2];
rz(-0.54272932) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(-2.824596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19001374) q[0];
sx q[0];
rz(-1.0931451) q[0];
sx q[0];
rz(-2.6059294) q[0];
rz(-2.426892) q[1];
sx q[1];
rz(-1.2477881) q[1];
sx q[1];
rz(1.4664149) q[1];
rz(1.9541478) q[2];
sx q[2];
rz(-1.162983) q[2];
sx q[2];
rz(2.4833124) q[2];
rz(-0.067848005) q[3];
sx q[3];
rz(-2.4115457) q[3];
sx q[3];
rz(-1.5822165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
