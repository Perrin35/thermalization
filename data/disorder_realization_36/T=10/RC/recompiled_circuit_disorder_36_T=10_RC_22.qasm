OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429457) q[0];
sx q[0];
rz(-0.67910128) q[0];
sx q[0];
rz(-2.8773017) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3876786) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-0.66893259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9211728) q[1];
sx q[1];
rz(-2.6487659) q[1];
sx q[1];
rz(-2.1652031) q[1];
rz(-pi) q[2];
rz(-2.4926315) q[3];
sx q[3];
rz(-2.0994086) q[3];
sx q[3];
rz(1.1543857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(-0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-0.84567436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613522) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(-2.4387226) q[0];
rz(1.7010062) q[2];
sx q[2];
rz(-1.7034334) q[2];
sx q[2];
rz(0.33078937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0963124) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(2.6394096) q[1];
rz(0.50176974) q[3];
sx q[3];
rz(-2.1309149) q[3];
sx q[3];
rz(-0.027651699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(1.7885615) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-0.95124328) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7032996) q[0];
sx q[0];
rz(-0.46143954) q[0];
sx q[0];
rz(1.5745387) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7227206) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(0.34130794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.31049) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(0.28880854) q[1];
x q[2];
rz(-0.17351405) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(-0.80843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0125668) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(-2.7584372) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(3.0932328) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.4455459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955129) q[0];
sx q[0];
rz(-2.4582986) q[0];
sx q[0];
rz(2.8542095) q[0];
rz(1.6608597) q[2];
sx q[2];
rz(-1.7863331) q[2];
sx q[2];
rz(0.71722523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1452892) q[1];
sx q[1];
rz(-2.9128296) q[1];
sx q[1];
rz(0.28782515) q[1];
x q[2];
rz(2.0914145) q[3];
sx q[3];
rz(-1.319066) q[3];
sx q[3];
rz(2.7057735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0754898) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4860977) q[0];
sx q[0];
rz(-2.4276519) q[0];
sx q[0];
rz(-1.5833202) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48798497) q[2];
sx q[2];
rz(-0.93163604) q[2];
sx q[2];
rz(-1.1346863) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8728719) q[1];
sx q[1];
rz(-2.1495021) q[1];
sx q[1];
rz(-0.89923664) q[1];
rz(-pi) q[2];
rz(1.0075931) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(1.0790881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(3.0153826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431625) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(-0.14847319) q[0];
rz(-pi) q[1];
rz(-0.77504471) q[2];
sx q[2];
rz(-1.2942874) q[2];
sx q[2];
rz(-2.0691878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5050161) q[1];
sx q[1];
rz(-2.2762183) q[1];
sx q[1];
rz(-0.24800639) q[1];
rz(1.24228) q[3];
sx q[3];
rz(-1.3887172) q[3];
sx q[3];
rz(1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-0.77511707) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59259748) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-2.5820406) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643499) q[0];
sx q[0];
rz(-0.61315216) q[0];
sx q[0];
rz(2.9186547) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6227116) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.3494929) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2110062) q[1];
sx q[1];
rz(-1.6582656) q[1];
sx q[1];
rz(-1.7006111) q[1];
rz(2.5114602) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(0.4160479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7155834) q[0];
sx q[0];
rz(-1.4953488) q[0];
sx q[0];
rz(1.5619318) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34198728) q[2];
sx q[2];
rz(-2.7139398) q[2];
sx q[2];
rz(2.4004186) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.65054446) q[1];
sx q[1];
rz(-0.77495134) q[1];
sx q[1];
rz(-0.52853711) q[1];
x q[2];
rz(-1.3685162) q[3];
sx q[3];
rz(-0.61198046) q[3];
sx q[3];
rz(-2.4806541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-1.139572) q[2];
rz(-1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(1.5244012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7302007) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(-1.9090396) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47713251) q[2];
sx q[2];
rz(-1.7413365) q[2];
sx q[2];
rz(2.0810623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2251687) q[1];
sx q[1];
rz(-1.4125707) q[1];
sx q[1];
rz(-1.2374864) q[1];
rz(-pi) q[2];
rz(1.8839621) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(-2.9477011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(2.8881853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5288552) q[0];
sx q[0];
rz(-0.40898541) q[0];
sx q[0];
rz(1.8324052) q[0];
rz(0.16742736) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(-1.6850922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(2.5245689) q[1];
rz(0.6110544) q[3];
sx q[3];
rz(-1.2130034) q[3];
sx q[3];
rz(-0.65294453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59166756) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(0.7111711) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(1.8297557) q[2];
sx q[2];
rz(-1.2987518) q[2];
sx q[2];
rz(-2.1546298) q[2];
rz(-2.3861804) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
