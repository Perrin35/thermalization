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
rz(1.0032049) q[0];
sx q[0];
rz(-1.2004852) q[0];
sx q[0];
rz(1.0263654) q[0];
rz(-2.8431471) q[1];
sx q[1];
rz(-1.5084074) q[1];
sx q[1];
rz(2.1389979) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63051414) q[0];
sx q[0];
rz(-1.3670232) q[0];
sx q[0];
rz(0.023292631) q[0];
rz(3.0709362) q[2];
sx q[2];
rz(-2.4774654) q[2];
sx q[2];
rz(2.9717367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1977495) q[1];
sx q[1];
rz(-1.7043866) q[1];
sx q[1];
rz(-2.751256) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75768023) q[3];
sx q[3];
rz(-0.79059383) q[3];
sx q[3];
rz(-2.9942715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5528494) q[2];
sx q[2];
rz(-1.0007977) q[2];
sx q[2];
rz(3.106485) q[2];
rz(1.2037753) q[3];
sx q[3];
rz(-1.3222008) q[3];
sx q[3];
rz(-2.8573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069000706) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(-2.1373855) q[0];
rz(-3.0191811) q[1];
sx q[1];
rz(-1.6740084) q[1];
sx q[1];
rz(-2.4845128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96913494) q[0];
sx q[0];
rz(-1.9577049) q[0];
sx q[0];
rz(1.659214) q[0];
x q[1];
rz(2.9802965) q[2];
sx q[2];
rz(-1.2968924) q[2];
sx q[2];
rz(2.1333926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2942336) q[1];
sx q[1];
rz(-2.2072133) q[1];
sx q[1];
rz(2.5381303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3260001) q[3];
sx q[3];
rz(-2.3395633) q[3];
sx q[3];
rz(2.0662465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5848026) q[2];
sx q[2];
rz(-0.57791296) q[2];
sx q[2];
rz(2.54971) q[2];
rz(1.4764192) q[3];
sx q[3];
rz(-1.6074601) q[3];
sx q[3];
rz(2.2679451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1378491) q[0];
sx q[0];
rz(-2.9326404) q[0];
sx q[0];
rz(1.9269706) q[0];
rz(-0.95642033) q[1];
sx q[1];
rz(-1.950187) q[1];
sx q[1];
rz(1.9699684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75367613) q[0];
sx q[0];
rz(-1.1769341) q[0];
sx q[0];
rz(2.0772165) q[0];
x q[1];
rz(2.2713216) q[2];
sx q[2];
rz(-1.4129801) q[2];
sx q[2];
rz(1.6788788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70832789) q[1];
sx q[1];
rz(-1.2339795) q[1];
sx q[1];
rz(2.5418806) q[1];
x q[2];
rz(-0.089853386) q[3];
sx q[3];
rz(-1.615534) q[3];
sx q[3];
rz(-0.60715946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20022915) q[2];
sx q[2];
rz(-2.0912632) q[2];
sx q[2];
rz(-1.3776779) q[2];
rz(-2.7348943) q[3];
sx q[3];
rz(-2.4923057) q[3];
sx q[3];
rz(-0.46050921) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7895301) q[0];
sx q[0];
rz(-1.2593513) q[0];
sx q[0];
rz(1.2404741) q[0];
rz(-0.71818304) q[1];
sx q[1];
rz(-1.8156464) q[1];
sx q[1];
rz(0.2389508) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093042308) q[0];
sx q[0];
rz(-1.4497546) q[0];
sx q[0];
rz(2.8516475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3117606) q[2];
sx q[2];
rz(-0.90832635) q[2];
sx q[2];
rz(-1.5892513) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0480369) q[1];
sx q[1];
rz(-0.29680064) q[1];
sx q[1];
rz(-1.0637299) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4364868) q[3];
sx q[3];
rz(-2.1013936) q[3];
sx q[3];
rz(2.1710896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51603985) q[2];
sx q[2];
rz(-1.7480787) q[2];
sx q[2];
rz(1.7986521) q[2];
rz(-2.2602153) q[3];
sx q[3];
rz(-0.16888976) q[3];
sx q[3];
rz(-2.3605997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.4137022) q[0];
sx q[0];
rz(-2.8975633) q[0];
sx q[0];
rz(-1.6572886) q[0];
rz(0.65385747) q[1];
sx q[1];
rz(-1.6203208) q[1];
sx q[1];
rz(0.13793764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1178095) q[0];
sx q[0];
rz(-1.9236757) q[0];
sx q[0];
rz(2.9830052) q[0];
x q[1];
rz(2.8213812) q[2];
sx q[2];
rz(-0.67812014) q[2];
sx q[2];
rz(-1.4411639) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5522462) q[1];
sx q[1];
rz(-1.668743) q[1];
sx q[1];
rz(-1.1393273) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.82965) q[3];
sx q[3];
rz(-1.5167243) q[3];
sx q[3];
rz(1.3075706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5518034) q[2];
sx q[2];
rz(-1.9686331) q[2];
sx q[2];
rz(0.087184437) q[2];
rz(0.11463556) q[3];
sx q[3];
rz(-0.81511027) q[3];
sx q[3];
rz(2.7242928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7755985) q[0];
sx q[0];
rz(-1.6021148) q[0];
sx q[0];
rz(0.48496801) q[0];
rz(2.1980749) q[1];
sx q[1];
rz(-2.3285995) q[1];
sx q[1];
rz(1.9224723) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4480369) q[0];
sx q[0];
rz(-0.72759923) q[0];
sx q[0];
rz(-1.4240525) q[0];
rz(1.2487683) q[2];
sx q[2];
rz(-1.5552943) q[2];
sx q[2];
rz(0.44336927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.74138) q[1];
sx q[1];
rz(-1.4948436) q[1];
sx q[1];
rz(1.6981359) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0756453) q[3];
sx q[3];
rz(-2.2643914) q[3];
sx q[3];
rz(-1.949273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.24011) q[2];
sx q[2];
rz(-1.0604475) q[2];
sx q[2];
rz(2.2002952) q[2];
rz(-2.8955722) q[3];
sx q[3];
rz(-1.6451719) q[3];
sx q[3];
rz(1.3601607) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098002382) q[0];
sx q[0];
rz(-1.101838) q[0];
sx q[0];
rz(-1.164042) q[0];
rz(-1.3580492) q[1];
sx q[1];
rz(-2.2750504) q[1];
sx q[1];
rz(-1.7485626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37993615) q[0];
sx q[0];
rz(-1.6448061) q[0];
sx q[0];
rz(-1.9475027) q[0];
x q[1];
rz(-0.54945182) q[2];
sx q[2];
rz(-1.3260815) q[2];
sx q[2];
rz(-0.050087226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5389312) q[1];
sx q[1];
rz(-1.4653413) q[1];
sx q[1];
rz(-0.47172539) q[1];
rz(-pi) q[2];
rz(-2.249516) q[3];
sx q[3];
rz(-1.7463272) q[3];
sx q[3];
rz(2.9316559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67036575) q[2];
sx q[2];
rz(-1.4281851) q[2];
sx q[2];
rz(-0.29895374) q[2];
rz(-2.7361338) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(0.56431842) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4295171) q[0];
sx q[0];
rz(-2.8686664) q[0];
sx q[0];
rz(-0.78936973) q[0];
rz(2.7366267) q[1];
sx q[1];
rz(-1.3507495) q[1];
sx q[1];
rz(-2.6466323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4433561) q[0];
sx q[0];
rz(-2.3670417) q[0];
sx q[0];
rz(-2.2775035) q[0];
rz(-pi) q[1];
rz(-1.1210006) q[2];
sx q[2];
rz(-2.5230222) q[2];
sx q[2];
rz(0.44396675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5304184) q[1];
sx q[1];
rz(-1.3278759) q[1];
sx q[1];
rz(-2.3526885) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1405419) q[3];
sx q[3];
rz(-2.3214139) q[3];
sx q[3];
rz(-2.6626183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1627545) q[2];
sx q[2];
rz(-1.2252204) q[2];
sx q[2];
rz(2.8863353) q[2];
rz(-0.034959547) q[3];
sx q[3];
rz(-1.6182263) q[3];
sx q[3];
rz(2.8953654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73664767) q[0];
sx q[0];
rz(-2.078779) q[0];
sx q[0];
rz(2.431562) q[0];
rz(1.0011477) q[1];
sx q[1];
rz(-0.89718693) q[1];
sx q[1];
rz(0.62090105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7707405) q[0];
sx q[0];
rz(-1.9319879) q[0];
sx q[0];
rz(0.06451635) q[0];
x q[1];
rz(-2.484637) q[2];
sx q[2];
rz(-1.959942) q[2];
sx q[2];
rz(-0.78788131) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8709594) q[1];
sx q[1];
rz(-0.3644202) q[1];
sx q[1];
rz(-0.30521133) q[1];
x q[2];
rz(2.8594703) q[3];
sx q[3];
rz(-2.3446313) q[3];
sx q[3];
rz(-2.7180501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3322477) q[2];
sx q[2];
rz(-1.5231909) q[2];
sx q[2];
rz(-2.8524032) q[2];
rz(2.0293763) q[3];
sx q[3];
rz(-2.8465392) q[3];
sx q[3];
rz(-1.0203863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6638829) q[0];
sx q[0];
rz(-1.7998671) q[0];
sx q[0];
rz(0.92673242) q[0];
rz(-1.1070739) q[1];
sx q[1];
rz(-1.6278382) q[1];
sx q[1];
rz(0.56799299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999789) q[0];
sx q[0];
rz(-1.2547412) q[0];
sx q[0];
rz(-2.5074682) q[0];
x q[1];
rz(1.307042) q[2];
sx q[2];
rz(-1.0737952) q[2];
sx q[2];
rz(2.618034) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24577877) q[1];
sx q[1];
rz(-1.2061283) q[1];
sx q[1];
rz(2.2321537) q[1];
rz(-pi) q[2];
rz(-2.0906914) q[3];
sx q[3];
rz(-2.1445159) q[3];
sx q[3];
rz(0.75999505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.052281436) q[2];
sx q[2];
rz(-2.3312882) q[2];
sx q[2];
rz(1.4386162) q[2];
rz(0.29664052) q[3];
sx q[3];
rz(-0.34160015) q[3];
sx q[3];
rz(2.6945485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.457837) q[0];
sx q[0];
rz(-2.044027) q[0];
sx q[0];
rz(-2.5720163) q[0];
rz(-2.8029022) q[1];
sx q[1];
rz(-1.5298005) q[1];
sx q[1];
rz(-1.6385967) q[1];
rz(1.7768301) q[2];
sx q[2];
rz(-1.3246957) q[2];
sx q[2];
rz(-0.93304721) q[2];
rz(1.8453311) q[3];
sx q[3];
rz(-1.7332776) q[3];
sx q[3];
rz(0.83040614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
