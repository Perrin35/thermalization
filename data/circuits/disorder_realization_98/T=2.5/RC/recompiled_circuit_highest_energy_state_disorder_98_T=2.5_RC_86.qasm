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
rz(3.0124445) q[0];
sx q[0];
rz(-1.4718066) q[0];
sx q[0];
rz(2.1762525) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(-1.3774011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7962153) q[0];
sx q[0];
rz(-0.85897972) q[0];
sx q[0];
rz(0.9839566) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1093311) q[2];
sx q[2];
rz(-0.92515495) q[2];
sx q[2];
rz(-1.0021068) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7143615) q[1];
sx q[1];
rz(-2.9032482) q[1];
sx q[1];
rz(2.4803216) q[1];
x q[2];
rz(-3.1226452) q[3];
sx q[3];
rz(-2.3403882) q[3];
sx q[3];
rz(-2.612118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2237902) q[2];
sx q[2];
rz(-1.6518355) q[2];
sx q[2];
rz(1.6415143) q[2];
rz(-0.7817868) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(-2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016945275) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(2.1710904) q[0];
rz(-1.3264725) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(-2.2296947) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2925948) q[0];
sx q[0];
rz(-2.1540897) q[0];
sx q[0];
rz(-1.769577) q[0];
rz(-pi) q[1];
rz(-1.3746337) q[2];
sx q[2];
rz(-1.7387783) q[2];
sx q[2];
rz(0.84174918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32102958) q[1];
sx q[1];
rz(-2.2707175) q[1];
sx q[1];
rz(-3.0882443) q[1];
x q[2];
rz(-1.6025869) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(-1.8546752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.2328925) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(-2.1523037) q[2];
rz(2.7214637) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1437538) q[0];
sx q[0];
rz(-1.1622575) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(0.85820091) q[1];
sx q[1];
rz(-2.61519) q[1];
sx q[1];
rz(2.9720378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1829878) q[0];
sx q[0];
rz(-2.9979994) q[0];
sx q[0];
rz(-1.5426226) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1716258) q[2];
sx q[2];
rz(-1.2720275) q[2];
sx q[2];
rz(2.5577161) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21903164) q[1];
sx q[1];
rz(-1.2777849) q[1];
sx q[1];
rz(-2.0907563) q[1];
rz(-2.9870728) q[3];
sx q[3];
rz(-1.9294881) q[3];
sx q[3];
rz(1.4160938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9618591) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(-2.7355984) q[2];
rz(-0.77754846) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(-0.8146666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4985713) q[0];
sx q[0];
rz(-1.2914456) q[0];
sx q[0];
rz(-0.92798573) q[0];
rz(-3.0827177) q[1];
sx q[1];
rz(-1.4480271) q[1];
sx q[1];
rz(1.4307129) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3979523) q[0];
sx q[0];
rz(-2.7317606) q[0];
sx q[0];
rz(1.3751956) q[0];
x q[1];
rz(2.7785382) q[2];
sx q[2];
rz(-1.3842693) q[2];
sx q[2];
rz(-0.26939738) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7652755) q[1];
sx q[1];
rz(-0.54794725) q[1];
sx q[1];
rz(-1.1330695) q[1];
x q[2];
rz(-0.9417504) q[3];
sx q[3];
rz(-1.3247196) q[3];
sx q[3];
rz(1.6501282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4579939) q[2];
sx q[2];
rz(-1.4983838) q[2];
sx q[2];
rz(1.9963473) q[2];
rz(-2.2942719) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97750807) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(0.384828) q[0];
rz(-0.79967868) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(-1.5481366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.528397) q[0];
sx q[0];
rz(-1.3269182) q[0];
sx q[0];
rz(2.8923678) q[0];
x q[1];
rz(-1.9673011) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(-2.3421974) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.3272737) q[1];
sx q[1];
rz(-0.47879915) q[1];
sx q[1];
rz(-2.4079066) q[1];
x q[2];
rz(1.0209072) q[3];
sx q[3];
rz(-1.8916777) q[3];
sx q[3];
rz(-2.4060017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7657713) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(3.0613464) q[2];
rz(-0.56096983) q[3];
sx q[3];
rz(-1.6601945) q[3];
sx q[3];
rz(2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-2.4764562) q[0];
sx q[0];
rz(1.8810062) q[0];
rz(0.27443019) q[1];
sx q[1];
rz(-1.6703037) q[1];
sx q[1];
rz(2.4226277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711766) q[0];
sx q[0];
rz(-2.377802) q[0];
sx q[0];
rz(1.0538573) q[0];
rz(-pi) q[1];
rz(2.3173213) q[2];
sx q[2];
rz(-1.5635412) q[2];
sx q[2];
rz(1.5985009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2430035) q[1];
sx q[1];
rz(-0.92783992) q[1];
sx q[1];
rz(-1.7220294) q[1];
x q[2];
rz(-0.218504) q[3];
sx q[3];
rz(-1.2077304) q[3];
sx q[3];
rz(-0.6930815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6545973) q[2];
sx q[2];
rz(-1.2664814) q[2];
sx q[2];
rz(-0.60301644) q[2];
rz(-1.4740137) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6239768) q[0];
sx q[0];
rz(-1.6083953) q[0];
sx q[0];
rz(-2.9543167) q[0];
rz(0.97224832) q[1];
sx q[1];
rz(-2.9863803) q[1];
sx q[1];
rz(0.42047277) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3492476) q[0];
sx q[0];
rz(-2.2195243) q[0];
sx q[0];
rz(-1.0386085) q[0];
rz(2.114061) q[2];
sx q[2];
rz(-1.2422556) q[2];
sx q[2];
rz(-1.3628886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6217611) q[1];
sx q[1];
rz(-2.6606307) q[1];
sx q[1];
rz(0.93641187) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99655788) q[3];
sx q[3];
rz(-1.3769994) q[3];
sx q[3];
rz(-2.7791948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(-2.8908758) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.3919132) q[3];
sx q[3];
rz(-2.5417476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.8641758) q[0];
sx q[0];
rz(-2.5499948) q[0];
sx q[0];
rz(0.94171062) q[0];
rz(-0.038453728) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(-0.11631913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021892) q[0];
sx q[0];
rz(-0.84745126) q[0];
sx q[0];
rz(0.96203104) q[0];
rz(-pi) q[1];
rz(2.2519926) q[2];
sx q[2];
rz(-1.4228914) q[2];
sx q[2];
rz(3.0513632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47919905) q[1];
sx q[1];
rz(-0.31649193) q[1];
sx q[1];
rz(-2.5280158) q[1];
rz(-pi) q[2];
rz(0.38920684) q[3];
sx q[3];
rz(-1.4688244) q[3];
sx q[3];
rz(1.9736279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2945127) q[2];
sx q[2];
rz(-0.81081644) q[2];
sx q[2];
rz(-1.026356) q[2];
rz(-1.9319084) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(0.9616372) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66861361) q[0];
sx q[0];
rz(-2.0675779) q[0];
sx q[0];
rz(-2.7630254) q[0];
rz(-1.1096795) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(0.027677061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6768764) q[0];
sx q[0];
rz(-2.4188571) q[0];
sx q[0];
rz(0.56761543) q[0];
rz(-pi) q[1];
rz(-0.1849298) q[2];
sx q[2];
rz(-1.7835296) q[2];
sx q[2];
rz(-0.30480584) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3125449) q[1];
sx q[1];
rz(-1.9282189) q[1];
sx q[1];
rz(0.19695671) q[1];
x q[2];
rz(0.94407336) q[3];
sx q[3];
rz(-0.51261307) q[3];
sx q[3];
rz(2.1533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3045584) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(2.7216116) q[2];
rz(-2.3000681) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(-0.37874547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(2.8564659) q[0];
rz(0.20052234) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(2.3060422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0203945) q[0];
sx q[0];
rz(-2.2403952) q[0];
sx q[0];
rz(-2.5479937) q[0];
x q[1];
rz(-3.1329536) q[2];
sx q[2];
rz(-0.64707478) q[2];
sx q[2];
rz(1.4467913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3491464) q[1];
sx q[1];
rz(-2.0172152) q[1];
sx q[1];
rz(-1.940215) q[1];
rz(-1.6099168) q[3];
sx q[3];
rz(-1.8062544) q[3];
sx q[3];
rz(-0.92604107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7150813) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(0.69768989) q[2];
rz(2.6155124) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(-1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1220916) q[0];
sx q[0];
rz(-0.87845907) q[0];
sx q[0];
rz(-0.39977951) q[0];
rz(-1.9922235) q[1];
sx q[1];
rz(-2.414357) q[1];
sx q[1];
rz(2.3946708) q[1];
rz(2.2483038) q[2];
sx q[2];
rz(-1.8468241) q[2];
sx q[2];
rz(2.1459864) q[2];
rz(-0.4646085) q[3];
sx q[3];
rz(-0.59904848) q[3];
sx q[3];
rz(-0.11943331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
