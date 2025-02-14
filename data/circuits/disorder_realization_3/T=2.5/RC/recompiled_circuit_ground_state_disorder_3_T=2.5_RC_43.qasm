OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8752911) q[0];
sx q[0];
rz(3.5485062) q[0];
sx q[0];
rz(9.2223528) q[0];
rz(1.8316733) q[1];
sx q[1];
rz(4.2869422) q[1];
sx q[1];
rz(11.577679) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017254596) q[0];
sx q[0];
rz(-1.3815123) q[0];
sx q[0];
rz(0.81228492) q[0];
rz(1.4329756) q[2];
sx q[2];
rz(-1.4933961) q[2];
sx q[2];
rz(3.1397506) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0571756) q[1];
sx q[1];
rz(-1.345008) q[1];
sx q[1];
rz(2.1979519) q[1];
x q[2];
rz(1.1448496) q[3];
sx q[3];
rz(-1.4574582) q[3];
sx q[3];
rz(2.9360611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6053091) q[2];
sx q[2];
rz(-0.80159694) q[2];
sx q[2];
rz(-2.0901704) q[2];
rz(1.5031523) q[3];
sx q[3];
rz(-1.4132615) q[3];
sx q[3];
rz(0.2352636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94238344) q[0];
sx q[0];
rz(-1.0301882) q[0];
sx q[0];
rz(-2.8675766) q[0];
rz(-0.38385299) q[1];
sx q[1];
rz(-2.5089788) q[1];
sx q[1];
rz(-1.3207818) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1766168) q[0];
sx q[0];
rz(-2.9701485) q[0];
sx q[0];
rz(-1.9587396) q[0];
x q[1];
rz(-0.4993778) q[2];
sx q[2];
rz(-1.3670397) q[2];
sx q[2];
rz(2.8600542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49628851) q[1];
sx q[1];
rz(-2.4170365) q[1];
sx q[1];
rz(1.4990119) q[1];
rz(2.3521496) q[3];
sx q[3];
rz(-2.200921) q[3];
sx q[3];
rz(-0.9073782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5447834) q[2];
sx q[2];
rz(-1.4027255) q[2];
sx q[2];
rz(-1.7421494) q[2];
rz(2.5341189) q[3];
sx q[3];
rz(-1.6739269) q[3];
sx q[3];
rz(2.7510711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6121599) q[0];
sx q[0];
rz(-2.2859892) q[0];
sx q[0];
rz(2.8223619) q[0];
rz(-0.051479738) q[1];
sx q[1];
rz(-1.9620644) q[1];
sx q[1];
rz(-2.0850339) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.257911) q[0];
sx q[0];
rz(-1.6223094) q[0];
sx q[0];
rz(-1.7345588) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86528565) q[2];
sx q[2];
rz(-1.6076536) q[2];
sx q[2];
rz(1.3071905) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33256862) q[1];
sx q[1];
rz(-1.7454552) q[1];
sx q[1];
rz(2.2038384) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3567188) q[3];
sx q[3];
rz(-0.41164648) q[3];
sx q[3];
rz(0.95142196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9487379) q[2];
sx q[2];
rz(-1.4317908) q[2];
sx q[2];
rz(-1.3445492) q[2];
rz(-2.2236842) q[3];
sx q[3];
rz(-2.5036006) q[3];
sx q[3];
rz(-1.2558254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767839) q[0];
sx q[0];
rz(-2.5510241) q[0];
sx q[0];
rz(1.508924) q[0];
rz(-2.6347939) q[1];
sx q[1];
rz(-1.2134039) q[1];
sx q[1];
rz(0.39055821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9386425) q[0];
sx q[0];
rz(-0.0084059518) q[0];
sx q[0];
rz(-0.52830835) q[0];
rz(0.37997577) q[2];
sx q[2];
rz(-1.3936498) q[2];
sx q[2];
rz(-2.4823657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4349065) q[1];
sx q[1];
rz(-0.95206407) q[1];
sx q[1];
rz(-3.1067585) q[1];
rz(-2.6253013) q[3];
sx q[3];
rz(-2.9598979) q[3];
sx q[3];
rz(2.5207375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8916696) q[2];
sx q[2];
rz(-2.2107783) q[2];
sx q[2];
rz(2.2451952) q[2];
rz(2.2253288) q[3];
sx q[3];
rz(-2.2057585) q[3];
sx q[3];
rz(-2.1808482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59178281) q[0];
sx q[0];
rz(-2.7857605) q[0];
sx q[0];
rz(-2.671396) q[0];
rz(-1.7377986) q[1];
sx q[1];
rz(-0.70039582) q[1];
sx q[1];
rz(1.2276924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0225343) q[0];
sx q[0];
rz(-1.5644363) q[0];
sx q[0];
rz(-1.8409078) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1433021) q[2];
sx q[2];
rz(-1.059747) q[2];
sx q[2];
rz(1.6533899) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0684546) q[1];
sx q[1];
rz(-0.92759354) q[1];
sx q[1];
rz(3.0469608) q[1];
rz(-pi) q[2];
rz(0.53831921) q[3];
sx q[3];
rz(-1.3019058) q[3];
sx q[3];
rz(-1.6398095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1710743) q[2];
sx q[2];
rz(-0.71275622) q[2];
sx q[2];
rz(-1.7573382) q[2];
rz(2.5303043) q[3];
sx q[3];
rz(-1.7620112) q[3];
sx q[3];
rz(0.39613625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9983845) q[0];
sx q[0];
rz(-2.12976) q[0];
sx q[0];
rz(-0.55214733) q[0];
rz(-0.72969189) q[1];
sx q[1];
rz(-2.153986) q[1];
sx q[1];
rz(2.139835) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14267966) q[0];
sx q[0];
rz(-0.15593869) q[0];
sx q[0];
rz(-2.1151311) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2127953) q[2];
sx q[2];
rz(-1.1603519) q[2];
sx q[2];
rz(0.14308077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3730091) q[1];
sx q[1];
rz(-0.78861744) q[1];
sx q[1];
rz(-2.8810701) q[1];
rz(-pi) q[2];
rz(0.67906191) q[3];
sx q[3];
rz(-1.1521395) q[3];
sx q[3];
rz(-2.1931894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6843162) q[2];
sx q[2];
rz(-1.3490973) q[2];
sx q[2];
rz(1.7982193) q[2];
rz(2.4172879) q[3];
sx q[3];
rz(-1.9380006) q[3];
sx q[3];
rz(-2.8389285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90172076) q[0];
sx q[0];
rz(-1.8093103) q[0];
sx q[0];
rz(2.4161762) q[0];
rz(1.2984917) q[1];
sx q[1];
rz(-2.6893171) q[1];
sx q[1];
rz(2.8628912) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9735457) q[0];
sx q[0];
rz(-2.5802045) q[0];
sx q[0];
rz(2.8532203) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33228525) q[2];
sx q[2];
rz(-1.1934308) q[2];
sx q[2];
rz(1.8695631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0289474) q[1];
sx q[1];
rz(-2.3621846) q[1];
sx q[1];
rz(-2.0112627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82741957) q[3];
sx q[3];
rz(-1.1779092) q[3];
sx q[3];
rz(3.06638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38410386) q[2];
sx q[2];
rz(-2.8748547) q[2];
sx q[2];
rz(-2.1534446) q[2];
rz(-3.0366376) q[3];
sx q[3];
rz(-1.2982488) q[3];
sx q[3];
rz(0.23133639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.4526378) q[0];
sx q[0];
rz(-0.80113688) q[0];
sx q[0];
rz(0.6193921) q[0];
rz(3.1267005) q[1];
sx q[1];
rz(-2.2512524) q[1];
sx q[1];
rz(-0.68170086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68530689) q[0];
sx q[0];
rz(-2.5192252) q[0];
sx q[0];
rz(-2.078889) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.080198535) q[2];
sx q[2];
rz(-2.2343582) q[2];
sx q[2];
rz(2.4303183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.926103) q[1];
sx q[1];
rz(-1.5640096) q[1];
sx q[1];
rz(0.35159638) q[1];
rz(-0.83311773) q[3];
sx q[3];
rz(-1.4280115) q[3];
sx q[3];
rz(0.24417434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2316078) q[2];
sx q[2];
rz(-2.3384428) q[2];
sx q[2];
rz(-0.33965674) q[2];
rz(-0.61399442) q[3];
sx q[3];
rz(-1.7620148) q[3];
sx q[3];
rz(-1.5617255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9851819) q[0];
sx q[0];
rz(-0.34813938) q[0];
sx q[0];
rz(-0.85439318) q[0];
rz(2.5482381) q[1];
sx q[1];
rz(-2.1095468) q[1];
sx q[1];
rz(-0.29534435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12796062) q[0];
sx q[0];
rz(-1.5379224) q[0];
sx q[0];
rz(0.055303975) q[0];
x q[1];
rz(3.1190671) q[2];
sx q[2];
rz(-1.1531382) q[2];
sx q[2];
rz(-0.45271046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.015739) q[1];
sx q[1];
rz(-0.89198128) q[1];
sx q[1];
rz(0.62894459) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.909524) q[3];
sx q[3];
rz(-0.32707387) q[3];
sx q[3];
rz(-1.0504521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16193834) q[2];
sx q[2];
rz(-1.6831968) q[2];
sx q[2];
rz(0.68183199) q[2];
rz(1.9409059) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(-0.43584263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1107263) q[0];
sx q[0];
rz(-1.2799355) q[0];
sx q[0];
rz(2.6357292) q[0];
rz(-2.1234296) q[1];
sx q[1];
rz(-1.4351427) q[1];
sx q[1];
rz(2.2811269) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2967401) q[0];
sx q[0];
rz(-2.2460947) q[0];
sx q[0];
rz(1.5102415) q[0];
rz(-pi) q[1];
rz(-2.0138836) q[2];
sx q[2];
rz(-2.8193316) q[2];
sx q[2];
rz(1.4631504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0588433) q[1];
sx q[1];
rz(-2.8902021) q[1];
sx q[1];
rz(0.13240929) q[1];
rz(-pi) q[2];
rz(2.303185) q[3];
sx q[3];
rz(-2.3831019) q[3];
sx q[3];
rz(-3.1137147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0818417) q[2];
sx q[2];
rz(-1.9530719) q[2];
sx q[2];
rz(-2.5583978) q[2];
rz(-3.1405295) q[3];
sx q[3];
rz(-2.21057) q[3];
sx q[3];
rz(0.49334905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4790333) q[0];
sx q[0];
rz(-1.722535) q[0];
sx q[0];
rz(-2.2395635) q[0];
rz(0.61530151) q[1];
sx q[1];
rz(-1.6533783) q[1];
sx q[1];
rz(1.3396214) q[1];
rz(1.3134708) q[2];
sx q[2];
rz(-2.0791292) q[2];
sx q[2];
rz(2.8544478) q[2];
rz(-1.8614122) q[3];
sx q[3];
rz(-2.5463085) q[3];
sx q[3];
rz(0.13704722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
