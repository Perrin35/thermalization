OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.779939) q[0];
sx q[0];
rz(-0.10804636) q[0];
sx q[0];
rz(1.894423) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(-1.127004) q[1];
sx q[1];
rz(-0.74581528) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818258) q[0];
sx q[0];
rz(-0.39224658) q[0];
sx q[0];
rz(2.0779209) q[0];
rz(-pi) q[1];
rz(3.132944) q[2];
sx q[2];
rz(-1.7315355) q[2];
sx q[2];
rz(-0.5102821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2850029) q[1];
sx q[1];
rz(-1.7816646) q[1];
sx q[1];
rz(-0.1182571) q[1];
x q[2];
rz(-1.4036101) q[3];
sx q[3];
rz(-1.0568398) q[3];
sx q[3];
rz(-1.592976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1625533) q[2];
sx q[2];
rz(-0.59430846) q[2];
sx q[2];
rz(1.7621367) q[2];
rz(2.4123794) q[3];
sx q[3];
rz(-2.3766434) q[3];
sx q[3];
rz(-0.929207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58296975) q[0];
sx q[0];
rz(-2.0780777) q[0];
sx q[0];
rz(-0.24120086) q[0];
rz(-2.8855715) q[1];
sx q[1];
rz(-1.0216917) q[1];
sx q[1];
rz(2.3604438) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9140919) q[0];
sx q[0];
rz(-2.1315398) q[0];
sx q[0];
rz(-2.1888211) q[0];
rz(-pi) q[1];
rz(0.56480572) q[2];
sx q[2];
rz(-1.6132659) q[2];
sx q[2];
rz(-3.0926585) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1684188) q[1];
sx q[1];
rz(-2.2464464) q[1];
sx q[1];
rz(2.3526117) q[1];
rz(-1.2595176) q[3];
sx q[3];
rz(-2.9972509) q[3];
sx q[3];
rz(2.8022769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8476734) q[2];
sx q[2];
rz(-1.5896229) q[2];
sx q[2];
rz(1.5576564) q[2];
rz(3.1368351) q[3];
sx q[3];
rz(-0.59499756) q[3];
sx q[3];
rz(0.34576542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3874338) q[0];
sx q[0];
rz(-1.4492946) q[0];
sx q[0];
rz(-0.19101983) q[0];
rz(-0.35107958) q[1];
sx q[1];
rz(-2.7103238) q[1];
sx q[1];
rz(-1.4494337) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80181771) q[0];
sx q[0];
rz(-1.2317941) q[0];
sx q[0];
rz(1.5168241) q[0];
rz(-pi) q[1];
rz(0.44523737) q[2];
sx q[2];
rz(-1.9539333) q[2];
sx q[2];
rz(-0.68510011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3684769) q[1];
sx q[1];
rz(-2.0252899) q[1];
sx q[1];
rz(1.2746432) q[1];
rz(-pi) q[2];
rz(-2.5496629) q[3];
sx q[3];
rz(-1.8576277) q[3];
sx q[3];
rz(-1.5805707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2991952) q[2];
sx q[2];
rz(-1.7123875) q[2];
sx q[2];
rz(3.0136285) q[2];
rz(-1.9757102) q[3];
sx q[3];
rz(-2.2539625) q[3];
sx q[3];
rz(-0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.8562451) q[0];
sx q[0];
rz(-2.0576394) q[0];
sx q[0];
rz(1.6612843) q[0];
rz(0.081143204) q[1];
sx q[1];
rz(-2.0442043) q[1];
sx q[1];
rz(2.2611484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9362088) q[0];
sx q[0];
rz(-3.1160389) q[0];
sx q[0];
rz(-2.8841302) q[0];
rz(-0.95194503) q[2];
sx q[2];
rz(-1.610872) q[2];
sx q[2];
rz(0.71499707) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8740011) q[1];
sx q[1];
rz(-2.1028215) q[1];
sx q[1];
rz(-1.4417955) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.038136884) q[3];
sx q[3];
rz(-2.5835369) q[3];
sx q[3];
rz(-1.7779153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5389898) q[2];
sx q[2];
rz(-2.4675641) q[2];
sx q[2];
rz(-1.4269786) q[2];
rz(-1.9849298) q[3];
sx q[3];
rz(-1.4256698) q[3];
sx q[3];
rz(0.15779933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.91931525) q[0];
sx q[0];
rz(-2.3265525) q[0];
sx q[0];
rz(-0.94006938) q[0];
rz(2.6433511) q[1];
sx q[1];
rz(-0.91612852) q[1];
sx q[1];
rz(-1.1246276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8064708) q[0];
sx q[0];
rz(-2.2735519) q[0];
sx q[0];
rz(-0.051856144) q[0];
x q[1];
rz(1.1246936) q[2];
sx q[2];
rz(-0.44844018) q[2];
sx q[2];
rz(2.5288127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4044721) q[1];
sx q[1];
rz(-1.1386765) q[1];
sx q[1];
rz(0.67373709) q[1];
x q[2];
rz(3.1203007) q[3];
sx q[3];
rz(-0.7856889) q[3];
sx q[3];
rz(1.1710492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4336832) q[2];
sx q[2];
rz(-0.83432546) q[2];
sx q[2];
rz(-1.679812) q[2];
rz(2.8953569) q[3];
sx q[3];
rz(-2.2610531) q[3];
sx q[3];
rz(-2.06854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6980625) q[0];
sx q[0];
rz(-1.227523) q[0];
sx q[0];
rz(-0.45502934) q[0];
rz(-2.9235234) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(0.47854447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6821176) q[0];
sx q[0];
rz(-2.2322901) q[0];
sx q[0];
rz(-1.4465998) q[0];
rz(2.6784513) q[2];
sx q[2];
rz(-0.69370334) q[2];
sx q[2];
rz(2.5460473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2854954) q[1];
sx q[1];
rz(-1.2427274) q[1];
sx q[1];
rz(-1.0873919) q[1];
x q[2];
rz(0.6932186) q[3];
sx q[3];
rz(-0.81855921) q[3];
sx q[3];
rz(-2.4002176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1856508) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(-0.68598023) q[2];
rz(-1.8020804) q[3];
sx q[3];
rz(-1.1687665) q[3];
sx q[3];
rz(1.7495988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5550391) q[0];
sx q[0];
rz(-1.7551442) q[0];
sx q[0];
rz(-2.6626124) q[0];
rz(2.2827177) q[1];
sx q[1];
rz(-0.54194599) q[1];
sx q[1];
rz(-3.0368793) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59706748) q[0];
sx q[0];
rz(-1.9550208) q[0];
sx q[0];
rz(1.7756419) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9373879) q[2];
sx q[2];
rz(-1.1127923) q[2];
sx q[2];
rz(1.4187494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4762209) q[1];
sx q[1];
rz(-1.7867861) q[1];
sx q[1];
rz(-0.98276999) q[1];
rz(1.0927622) q[3];
sx q[3];
rz(-1.6217188) q[3];
sx q[3];
rz(-1.4560521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0844448) q[2];
sx q[2];
rz(-1.8859325) q[2];
sx q[2];
rz(-0.88195938) q[2];
rz(0.070934892) q[3];
sx q[3];
rz(-1.8788012) q[3];
sx q[3];
rz(2.1766263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98081812) q[0];
sx q[0];
rz(-0.47176281) q[0];
sx q[0];
rz(1.8881352) q[0];
rz(-2.4313633) q[1];
sx q[1];
rz(-1.1980779) q[1];
sx q[1];
rz(-2.0517147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3364351) q[0];
sx q[0];
rz(-1.1602279) q[0];
sx q[0];
rz(-0.76469501) q[0];
rz(-pi) q[1];
rz(-2.5023191) q[2];
sx q[2];
rz(-2.6007923) q[2];
sx q[2];
rz(1.6317473) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1068258) q[1];
sx q[1];
rz(-2.2112101) q[1];
sx q[1];
rz(-0.57042112) q[1];
x q[2];
rz(1.4185216) q[3];
sx q[3];
rz(-0.99968761) q[3];
sx q[3];
rz(-0.53779569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7465884) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(-1.0224226) q[2];
rz(-2.5452781) q[3];
sx q[3];
rz(-2.3583581) q[3];
sx q[3];
rz(2.7885126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8280867) q[0];
sx q[0];
rz(-0.49740288) q[0];
sx q[0];
rz(1.585438) q[0];
rz(-1.6515139) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(0.99686399) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6108902) q[0];
sx q[0];
rz(-0.35323745) q[0];
sx q[0];
rz(-0.20521407) q[0];
x q[1];
rz(-1.385594) q[2];
sx q[2];
rz(-1.5723937) q[2];
sx q[2];
rz(1.314437) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0744508) q[1];
sx q[1];
rz(-2.0623226) q[1];
sx q[1];
rz(2.3105826) q[1];
rz(-0.50052754) q[3];
sx q[3];
rz(-2.1033035) q[3];
sx q[3];
rz(0.73126572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31522125) q[2];
sx q[2];
rz(-0.73306495) q[2];
sx q[2];
rz(-2.9054387) q[2];
rz(-0.30596966) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(-1.6570305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9295101) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(0.65810743) q[0];
rz(1.2492389) q[1];
sx q[1];
rz(-2.55195) q[1];
sx q[1];
rz(-2.6499937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015305681) q[0];
sx q[0];
rz(-2.0184061) q[0];
sx q[0];
rz(-2.5372895) q[0];
rz(-pi) q[1];
rz(-1.3479606) q[2];
sx q[2];
rz(-2.1249944) q[2];
sx q[2];
rz(1.8775307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2053615) q[1];
sx q[1];
rz(-2.865701) q[1];
sx q[1];
rz(3.107168) q[1];
rz(3.0594143) q[3];
sx q[3];
rz(-1.0834457) q[3];
sx q[3];
rz(-0.10457071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3601941) q[2];
sx q[2];
rz(-0.34803826) q[2];
sx q[2];
rz(1.6602824) q[2];
rz(0.71634746) q[3];
sx q[3];
rz(-2.4951388) q[3];
sx q[3];
rz(-2.9334478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3502055) q[0];
sx q[0];
rz(-1.5443784) q[0];
sx q[0];
rz(-2.7271893) q[0];
rz(0.95175891) q[1];
sx q[1];
rz(-1.8487683) q[1];
sx q[1];
rz(1.9602736) q[1];
rz(-0.69618113) q[2];
sx q[2];
rz(-1.0715108) q[2];
sx q[2];
rz(0.61658695) q[2];
rz(-2.2496102) q[3];
sx q[3];
rz(-1.7041698) q[3];
sx q[3];
rz(1.1052122) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
