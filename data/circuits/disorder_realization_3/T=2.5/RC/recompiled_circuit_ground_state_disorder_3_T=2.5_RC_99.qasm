OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2663015) q[0];
sx q[0];
rz(-0.40691352) q[0];
sx q[0];
rz(0.20242515) q[0];
rz(-1.3099194) q[1];
sx q[1];
rz(-1.1453495) q[1];
sx q[1];
rz(-2.1529012) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017254596) q[0];
sx q[0];
rz(-1.7600804) q[0];
sx q[0];
rz(2.3293077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0568809) q[2];
sx q[2];
rz(-2.9836453) q[2];
sx q[2];
rz(2.081209) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6472345) q[1];
sx q[1];
rz(-2.1796759) q[1];
sx q[1];
rz(-0.27642864) q[1];
rz(-pi) q[2];
rz(-1.3019788) q[3];
sx q[3];
rz(-0.43987396) q[3];
sx q[3];
rz(-1.1209821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6053091) q[2];
sx q[2];
rz(-0.80159694) q[2];
sx q[2];
rz(-1.0514222) q[2];
rz(1.6384404) q[3];
sx q[3];
rz(-1.4132615) q[3];
sx q[3];
rz(-0.2352636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94238344) q[0];
sx q[0];
rz(-2.1114045) q[0];
sx q[0];
rz(2.8675766) q[0];
rz(0.38385299) q[1];
sx q[1];
rz(-0.63261384) q[1];
sx q[1];
rz(1.8208108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21652554) q[0];
sx q[0];
rz(-1.7293879) q[0];
sx q[0];
rz(3.0761883) q[0];
x q[1];
rz(1.8019559) q[2];
sx q[2];
rz(-2.0589239) q[2];
sx q[2];
rz(1.7423992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7410424) q[1];
sx q[1];
rz(-0.84851754) q[1];
sx q[1];
rz(0.063401179) q[1];
x q[2];
rz(-2.3521496) q[3];
sx q[3];
rz(-0.94067162) q[3];
sx q[3];
rz(-0.9073782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5968093) q[2];
sx q[2];
rz(-1.7388672) q[2];
sx q[2];
rz(1.7421494) q[2];
rz(-0.60747373) q[3];
sx q[3];
rz(-1.6739269) q[3];
sx q[3];
rz(-0.39052159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6121599) q[0];
sx q[0];
rz(-2.2859892) q[0];
sx q[0];
rz(2.8223619) q[0];
rz(-3.0901129) q[1];
sx q[1];
rz(-1.9620644) q[1];
sx q[1];
rz(-1.0565588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.257911) q[0];
sx q[0];
rz(-1.5192832) q[0];
sx q[0];
rz(-1.7345588) q[0];
rz(-0.048398971) q[2];
sx q[2];
rz(-0.86586398) q[2];
sx q[2];
rz(2.8466109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3650546) q[1];
sx q[1];
rz(-0.94887304) q[1];
sx q[1];
rz(-2.9261241) q[1];
x q[2];
rz(-1.9740417) q[3];
sx q[3];
rz(-1.4856899) q[3];
sx q[3];
rz(0.81604609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9487379) q[2];
sx q[2];
rz(-1.7098018) q[2];
sx q[2];
rz(-1.3445492) q[2];
rz(2.2236842) q[3];
sx q[3];
rz(-2.5036006) q[3];
sx q[3];
rz(1.2558254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767839) q[0];
sx q[0];
rz(-0.5905686) q[0];
sx q[0];
rz(1.6326686) q[0];
rz(-0.50679874) q[1];
sx q[1];
rz(-1.9281887) q[1];
sx q[1];
rz(-2.7510344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4669663) q[0];
sx q[0];
rz(-1.5780562) q[0];
sx q[0];
rz(1.5750336) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7616169) q[2];
sx q[2];
rz(-1.3936498) q[2];
sx q[2];
rz(-2.4823657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4949172) q[1];
sx q[1];
rz(-0.61958379) q[1];
sx q[1];
rz(-1.5219206) q[1];
rz(-pi) q[2];
rz(1.6612442) q[3];
sx q[3];
rz(-1.7285943) q[3];
sx q[3];
rz(-3.0441847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8916696) q[2];
sx q[2];
rz(-2.2107783) q[2];
sx q[2];
rz(0.8963975) q[2];
rz(0.91626382) q[3];
sx q[3];
rz(-0.93583411) q[3];
sx q[3];
rz(0.96074444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5498098) q[0];
sx q[0];
rz(-0.35583219) q[0];
sx q[0];
rz(-0.47019666) q[0];
rz(1.7377986) q[1];
sx q[1];
rz(-0.70039582) q[1];
sx q[1];
rz(-1.2276924) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.550023) q[0];
sx q[0];
rz(-1.8409022) q[0];
sx q[0];
rz(3.1349934) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.50493) q[2];
sx q[2];
rz(-2.4877254) q[2];
sx q[2];
rz(0.9034397) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44547651) q[1];
sx q[1];
rz(-1.4951147) q[1];
sx q[1];
rz(-2.2161525) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6032734) q[3];
sx q[3];
rz(-1.3019058) q[3];
sx q[3];
rz(-1.5017832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97051835) q[2];
sx q[2];
rz(-2.4288364) q[2];
sx q[2];
rz(-1.3842545) q[2];
rz(0.61128831) q[3];
sx q[3];
rz(-1.7620112) q[3];
sx q[3];
rz(2.7454564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9983845) q[0];
sx q[0];
rz(-1.0118326) q[0];
sx q[0];
rz(-0.55214733) q[0];
rz(0.72969189) q[1];
sx q[1];
rz(-2.153986) q[1];
sx q[1];
rz(-2.139835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4491548) q[0];
sx q[0];
rz(-1.4375411) q[0];
sx q[0];
rz(-0.081234531) q[0];
x q[1];
rz(2.6438786) q[2];
sx q[2];
rz(-0.98955284) q[2];
sx q[2];
rz(1.1377942) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.525108) q[1];
sx q[1];
rz(-1.387038) q[1];
sx q[1];
rz(-2.3701369) q[1];
x q[2];
rz(2.4625307) q[3];
sx q[3];
rz(-1.9894532) q[3];
sx q[3];
rz(-2.1931894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6843162) q[2];
sx q[2];
rz(-1.3490973) q[2];
sx q[2];
rz(-1.3433733) q[2];
rz(2.4172879) q[3];
sx q[3];
rz(-1.2035921) q[3];
sx q[3];
rz(-0.30266416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90172076) q[0];
sx q[0];
rz(-1.3322823) q[0];
sx q[0];
rz(-2.4161762) q[0];
rz(1.8431009) q[1];
sx q[1];
rz(-0.45227554) q[1];
sx q[1];
rz(-0.27870146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.168047) q[0];
sx q[0];
rz(-2.5802045) q[0];
sx q[0];
rz(0.28837236) q[0];
rz(2.2594035) q[2];
sx q[2];
rz(-0.49748245) q[2];
sx q[2];
rz(-0.51940268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6141521) q[1];
sx q[1];
rz(-0.88178712) q[1];
sx q[1];
rz(-2.7428738) q[1];
rz(-0.82741957) q[3];
sx q[3];
rz(-1.9636834) q[3];
sx q[3];
rz(-3.06638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38410386) q[2];
sx q[2];
rz(-2.8748547) q[2];
sx q[2];
rz(-2.1534446) q[2];
rz(3.0366376) q[3];
sx q[3];
rz(-1.8433439) q[3];
sx q[3];
rz(-2.9102563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4526378) q[0];
sx q[0];
rz(-2.3404558) q[0];
sx q[0];
rz(-2.5222006) q[0];
rz(3.1267005) q[1];
sx q[1];
rz(-0.89034021) q[1];
sx q[1];
rz(-2.4598918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8312163) q[0];
sx q[0];
rz(-1.8583603) q[0];
sx q[0];
rz(-1.0108666) q[0];
rz(-pi) q[1];
rz(0.080198535) q[2];
sx q[2];
rz(-0.90723443) q[2];
sx q[2];
rz(2.4303183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3368083) q[1];
sx q[1];
rz(-2.7899335) q[1];
sx q[1];
rz(3.1218887) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7813788) q[3];
sx q[3];
rz(-0.74880785) q[3];
sx q[3];
rz(1.9702828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90998489) q[2];
sx q[2];
rz(-0.80314988) q[2];
sx q[2];
rz(2.8019359) q[2];
rz(-0.61399442) q[3];
sx q[3];
rz(-1.3795779) q[3];
sx q[3];
rz(1.5617255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851819) q[0];
sx q[0];
rz(-2.7934533) q[0];
sx q[0];
rz(-0.85439318) q[0];
rz(-2.5482381) q[1];
sx q[1];
rz(-1.0320458) q[1];
sx q[1];
rz(2.8462483) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6969374) q[0];
sx q[0];
rz(-1.5155223) q[0];
sx q[0];
rz(-1.6037206) q[0];
x q[1];
rz(-0.02252553) q[2];
sx q[2];
rz(-1.9884544) q[2];
sx q[2];
rz(0.45271046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.01646811) q[1];
sx q[1];
rz(-1.09519) q[1];
sx q[1];
rz(2.3549902) q[1];
rz(-pi) q[2];
rz(-0.11225742) q[3];
sx q[3];
rz(-1.8786544) q[3];
sx q[3];
rz(-0.69428919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.16193834) q[2];
sx q[2];
rz(-1.4583959) q[2];
sx q[2];
rz(2.4597607) q[2];
rz(-1.2006867) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(-0.43584263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0308663) q[0];
sx q[0];
rz(-1.2799355) q[0];
sx q[0];
rz(2.6357292) q[0];
rz(-2.1234296) q[1];
sx q[1];
rz(-1.4351427) q[1];
sx q[1];
rz(2.2811269) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31193908) q[0];
sx q[0];
rz(-1.6180493) q[0];
sx q[0];
rz(-2.465399) q[0];
rz(0.14218743) q[2];
sx q[2];
rz(-1.8609797) q[2];
sx q[2];
rz(0.99925257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3597534) q[1];
sx q[1];
rz(-1.5379496) q[1];
sx q[1];
rz(-0.24928045) q[1];
x q[2];
rz(-0.56475909) q[3];
sx q[3];
rz(-1.0339206) q[3];
sx q[3];
rz(-0.91954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0818417) q[2];
sx q[2];
rz(-1.1885208) q[2];
sx q[2];
rz(-0.58319485) q[2];
rz(3.1405295) q[3];
sx q[3];
rz(-2.21057) q[3];
sx q[3];
rz(2.6482436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6625593) q[0];
sx q[0];
rz(-1.4190577) q[0];
sx q[0];
rz(0.90202913) q[0];
rz(-2.5262911) q[1];
sx q[1];
rz(-1.6533783) q[1];
sx q[1];
rz(1.3396214) q[1];
rz(2.7131296) q[2];
sx q[2];
rz(-0.56461538) q[2];
sx q[2];
rz(2.3587894) q[2];
rz(2.9499182) q[3];
sx q[3];
rz(-1.0036461) q[3];
sx q[3];
rz(-0.20958513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
