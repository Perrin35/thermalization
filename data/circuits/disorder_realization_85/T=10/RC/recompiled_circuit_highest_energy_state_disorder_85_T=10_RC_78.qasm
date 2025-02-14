OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.149481) q[0];
sx q[0];
rz(-0.7465201) q[0];
sx q[0];
rz(2.574805) q[0];
rz(-1.8911288) q[1];
sx q[1];
rz(-1.1486147) q[1];
sx q[1];
rz(-2.221938) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56432589) q[0];
sx q[0];
rz(-2.3503605) q[0];
sx q[0];
rz(-1.7971695) q[0];
rz(-pi) q[1];
rz(-0.72534277) q[2];
sx q[2];
rz(-1.1769466) q[2];
sx q[2];
rz(0.79532901) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7417042) q[1];
sx q[1];
rz(-2.2117858) q[1];
sx q[1];
rz(-1.5128894) q[1];
rz(-pi) q[2];
rz(-2.9018974) q[3];
sx q[3];
rz(-2.5928232) q[3];
sx q[3];
rz(1.5769928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93996843) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(-1.0533818) q[2];
rz(-2.4273704) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(-1.8478954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6222222) q[0];
sx q[0];
rz(-0.70835963) q[0];
sx q[0];
rz(0.57193065) q[0];
rz(1.8427303) q[1];
sx q[1];
rz(-2.3426901) q[1];
sx q[1];
rz(0.78972185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30366722) q[0];
sx q[0];
rz(-1.2695489) q[0];
sx q[0];
rz(-0.10060723) q[0];
x q[1];
rz(-0.20767539) q[2];
sx q[2];
rz(-2.4768314) q[2];
sx q[2];
rz(1.1663811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29282031) q[1];
sx q[1];
rz(-2.1629984) q[1];
sx q[1];
rz(-1.9100389) q[1];
x q[2];
rz(-2.5511207) q[3];
sx q[3];
rz(-2.3749224) q[3];
sx q[3];
rz(1.1696512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(1.3624462) q[2];
rz(2.8507774) q[3];
sx q[3];
rz(-1.3664061) q[3];
sx q[3];
rz(-3.0349558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13880759) q[0];
sx q[0];
rz(-2.0464351) q[0];
sx q[0];
rz(-0.7005257) q[0];
rz(-2.6558212) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-2.9170759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.016708) q[0];
sx q[0];
rz(-1.9210235) q[0];
sx q[0];
rz(1.7279529) q[0];
rz(-pi) q[1];
rz(2.1491884) q[2];
sx q[2];
rz(-1.4881721) q[2];
sx q[2];
rz(-2.3985661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77188801) q[1];
sx q[1];
rz(-1.6460895) q[1];
sx q[1];
rz(-2.1545707) q[1];
rz(-pi) q[2];
x q[2];
rz(1.443601) q[3];
sx q[3];
rz(-2.4089097) q[3];
sx q[3];
rz(0.81058433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38480467) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(3.0745506) q[3];
sx q[3];
rz(-1.1578355) q[3];
sx q[3];
rz(0.66655794) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662358) q[0];
sx q[0];
rz(-2.5844564) q[0];
sx q[0];
rz(-0.019158451) q[0];
rz(-1.7823904) q[1];
sx q[1];
rz(-1.7117585) q[1];
sx q[1];
rz(0.0044936831) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.27994) q[0];
sx q[0];
rz(-2.3706382) q[0];
sx q[0];
rz(-2.3907135) q[0];
rz(-pi) q[1];
x q[1];
rz(1.121663) q[2];
sx q[2];
rz(-1.1156848) q[2];
sx q[2];
rz(-3.1246076) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5623206) q[1];
sx q[1];
rz(-1.4021789) q[1];
sx q[1];
rz(-1.9344485) q[1];
rz(2.717797) q[3];
sx q[3];
rz(-2.9188699) q[3];
sx q[3];
rz(0.84171695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6898474) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(1.6881662) q[2];
rz(1.2545741) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(2.5793251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8101863) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(-1.9796665) q[0];
rz(-2.3792073) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(-0.42957482) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31072894) q[0];
sx q[0];
rz(-2.7318044) q[0];
sx q[0];
rz(1.662246) q[0];
x q[1];
rz(7/(8*pi)) q[2];
sx q[2];
rz(-2.3633011) q[2];
sx q[2];
rz(-0.24519224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67744327) q[1];
sx q[1];
rz(-0.92532571) q[1];
sx q[1];
rz(-2.7149625) q[1];
rz(-pi) q[2];
rz(-1.1384835) q[3];
sx q[3];
rz(-2.4412182) q[3];
sx q[3];
rz(2.5770503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21542159) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(1.311709) q[2];
rz(-2.6808776) q[3];
sx q[3];
rz(-2.1381133) q[3];
sx q[3];
rz(-3.0084012) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(-2.0352236) q[0];
rz(3.0270992) q[1];
sx q[1];
rz(-2.0888927) q[1];
sx q[1];
rz(3.0879424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.69218) q[0];
sx q[0];
rz(-2.8859038) q[0];
sx q[0];
rz(1.5480255) q[0];
rz(-pi) q[1];
rz(1.6752376) q[2];
sx q[2];
rz(-0.76600961) q[2];
sx q[2];
rz(3.0867689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.304639) q[1];
sx q[1];
rz(-1.0598247) q[1];
sx q[1];
rz(-3.0797019) q[1];
rz(-pi) q[2];
rz(-2.1838004) q[3];
sx q[3];
rz(-1.0292064) q[3];
sx q[3];
rz(-2.0939746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(2.6344521) q[2];
rz(-1.3745314) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(-1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5116665) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(0.041286904) q[0];
rz(-0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(2.1521177) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7558003) q[0];
sx q[0];
rz(-0.72052252) q[0];
sx q[0];
rz(0.81554834) q[0];
x q[1];
rz(2.6648952) q[2];
sx q[2];
rz(-2.7455491) q[2];
sx q[2];
rz(0.11644289) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75162269) q[1];
sx q[1];
rz(-2.1643359) q[1];
sx q[1];
rz(-2.4291527) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0416077) q[3];
sx q[3];
rz(-2.4306477) q[3];
sx q[3];
rz(-2.8266852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5993293) q[2];
sx q[2];
rz(-1.0487391) q[2];
sx q[2];
rz(-2.2824724) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(3.0288127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82666731) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(-2.9097606) q[0];
rz(-0.58206093) q[1];
sx q[1];
rz(-1.8128017) q[1];
sx q[1];
rz(1.9225165) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9052637) q[0];
sx q[0];
rz(-1.4538611) q[0];
sx q[0];
rz(-2.5811483) q[0];
rz(-0.35301669) q[2];
sx q[2];
rz(-1.5215877) q[2];
sx q[2];
rz(2.4491058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2197709) q[1];
sx q[1];
rz(-1.3857462) q[1];
sx q[1];
rz(-2.9875523) q[1];
rz(-pi) q[2];
rz(-1.0295415) q[3];
sx q[3];
rz(-2.2502021) q[3];
sx q[3];
rz(-2.816659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1845188) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(2.3547724) q[2];
rz(1.6400853) q[3];
sx q[3];
rz(-2.7101176) q[3];
sx q[3];
rz(1.7155581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805434) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(2.2013262) q[0];
rz(-2.6967948) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(2.99517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6165926) q[0];
sx q[0];
rz(-1.6967275) q[0];
sx q[0];
rz(-0.30091826) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8690879) q[2];
sx q[2];
rz(-1.8159869) q[2];
sx q[2];
rz(-1.5281072) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89288974) q[1];
sx q[1];
rz(-1.8161402) q[1];
sx q[1];
rz(0.69458436) q[1];
rz(-pi) q[2];
rz(1.4692471) q[3];
sx q[3];
rz(-0.78208215) q[3];
sx q[3];
rz(2.8080432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0545097) q[2];
sx q[2];
rz(-1.6528218) q[2];
sx q[2];
rz(-0.83756891) q[2];
rz(2.2086823) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(-2.3605409) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26302108) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(1.7735057) q[0];
rz(0.076016501) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(-2.0215633) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.66066) q[0];
sx q[0];
rz(-1.2613861) q[0];
sx q[0];
rz(0.27746986) q[0];
rz(-pi) q[1];
rz(-2.4330288) q[2];
sx q[2];
rz(-1.8052793) q[2];
sx q[2];
rz(1.3402779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6516797) q[1];
sx q[1];
rz(-1.8615684) q[1];
sx q[1];
rz(-0.099822961) q[1];
x q[2];
rz(1.7342274) q[3];
sx q[3];
rz(-1.6728969) q[3];
sx q[3];
rz(2.7354376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7398305) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(-2.8945727) q[2];
rz(-2.5714827) q[3];
sx q[3];
rz(-0.48738042) q[3];
sx q[3];
rz(-2.0232239) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0433255) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(-0.098516057) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(1.1578887) q[2];
sx q[2];
rz(-1.558424) q[2];
sx q[2];
rz(1.3267714) q[2];
rz(0.5876585) q[3];
sx q[3];
rz(-2.6937204) q[3];
sx q[3];
rz(3.1093521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
