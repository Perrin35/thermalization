OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(1.8811037) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7605654) q[0];
sx q[0];
rz(-2.0187223) q[0];
sx q[0];
rz(-2.7678124) q[0];
rz(-pi) q[1];
x q[1];
rz(1.107723) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(0.83480922) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1957789) q[1];
sx q[1];
rz(-1.073277) q[1];
sx q[1];
rz(-1.0832018) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5039623) q[3];
sx q[3];
rz(-0.59083592) q[3];
sx q[3];
rz(1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(0.16201924) q[2];
rz(-2.2062733) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.686036) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8287795) q[0];
sx q[0];
rz(-1.7377186) q[0];
sx q[0];
rz(-2.1524327) q[0];
rz(-pi) q[1];
rz(-2.935264) q[2];
sx q[2];
rz(-2.7596139) q[2];
sx q[2];
rz(-1.0155592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93745366) q[1];
sx q[1];
rz(-2.0935645) q[1];
sx q[1];
rz(-0.015603113) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0807642) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(-3.1128186) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5867509) q[0];
sx q[0];
rz(-1.8790073) q[0];
sx q[0];
rz(2.5462333) q[0];
rz(-2.795479) q[2];
sx q[2];
rz(-2.3556404) q[2];
sx q[2];
rz(-1.2197989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6781569) q[1];
sx q[1];
rz(-1.0532182) q[1];
sx q[1];
rz(-2.8606158) q[1];
rz(-pi) q[2];
rz(1.7964732) q[3];
sx q[3];
rz(-0.25697069) q[3];
sx q[3];
rz(-0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0353161) q[0];
sx q[0];
rz(-1.9395394) q[0];
sx q[0];
rz(0.62555255) q[0];
rz(-pi) q[1];
rz(1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(2.0563682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3829141) q[1];
sx q[1];
rz(-0.74712336) q[1];
sx q[1];
rz(-1.9488571) q[1];
x q[2];
rz(-0.13682271) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(-2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.6500641) q[0];
rz(-0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46734738) q[0];
sx q[0];
rz(-1.6773946) q[0];
sx q[0];
rz(2.9664413) q[0];
x q[1];
rz(3.132658) q[2];
sx q[2];
rz(-1.8736471) q[2];
sx q[2];
rz(-2.6645899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3757513) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(-2.8326616) q[1];
rz(-1.7080073) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(-1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(-1.4250925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8558559) q[0];
sx q[0];
rz(-2.5897313) q[0];
sx q[0];
rz(-2.7049271) q[0];
rz(-1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(-2.595682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1849991) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(1.1207629) q[1];
rz(-2.1720042) q[3];
sx q[3];
rz(-1.117327) q[3];
sx q[3];
rz(0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99188995) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(-2.4292612) q[0];
x q[1];
rz(2.2981811) q[2];
sx q[2];
rz(-1.9551829) q[2];
sx q[2];
rz(-0.20993983) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8824749) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(-2.4005753) q[1];
rz(-pi) q[2];
rz(-2.5542198) q[3];
sx q[3];
rz(-1.8837187) q[3];
sx q[3];
rz(1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3372779) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-0.24469963) q[2];
rz(3.0120567) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.8364505) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96852036) q[0];
sx q[0];
rz(-0.99181306) q[0];
sx q[0];
rz(1.9645683) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9403946) q[2];
sx q[2];
rz(-2.2410789) q[2];
sx q[2];
rz(3.0073462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95272428) q[1];
sx q[1];
rz(-2.3844068) q[1];
sx q[1];
rz(2.5785179) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6373789) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(2.7293527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(1.0772609) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7548783) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(0.3219147) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(2.4386491) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84530572) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(1.9592497) q[0];
rz(-pi) q[1];
x q[1];
rz(0.039673294) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(2.2850125) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.326509) q[1];
sx q[1];
rz(-1.0352967) q[1];
sx q[1];
rz(2.9498847) q[1];
x q[2];
rz(2.1861107) q[3];
sx q[3];
rz(-1.2947047) q[3];
sx q[3];
rz(-0.51481065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9194473) q[0];
sx q[0];
rz(-1.3627909) q[0];
sx q[0];
rz(2.749445) q[0];
rz(-pi) q[1];
rz(1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(0.38244837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.66099) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(-1.8250699) q[1];
x q[2];
rz(-2.3234899) q[3];
sx q[3];
rz(-2.5460498) q[3];
sx q[3];
rz(-0.18225741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(0.11463595) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(-1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464012) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-0.8016349) q[2];
sx q[2];
rz(-2.2710706) q[2];
sx q[2];
rz(2.0588277) q[2];
rz(-3.0631089) q[3];
sx q[3];
rz(-2.2208636) q[3];
sx q[3];
rz(0.29490864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
