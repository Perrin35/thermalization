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
rz(-2.7346791) q[0];
sx q[0];
rz(-0.20242515) q[0];
rz(1.8316733) q[1];
sx q[1];
rz(4.2869422) q[1];
sx q[1];
rz(11.577679) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3770954) q[0];
sx q[0];
rz(-0.82907739) q[0];
sx q[0];
rz(-2.8835456) q[0];
rz(2.0847118) q[2];
sx q[2];
rz(-2.9836453) q[2];
sx q[2];
rz(1.0603836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4943581) q[1];
sx q[1];
rz(-2.1796759) q[1];
sx q[1];
rz(0.27642864) q[1];
rz(-pi) q[2];
rz(3.0172431) q[3];
sx q[3];
rz(-1.1477587) q[3];
sx q[3];
rz(-1.4165341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53628355) q[2];
sx q[2];
rz(-0.80159694) q[2];
sx q[2];
rz(1.0514222) q[2];
rz(1.6384404) q[3];
sx q[3];
rz(-1.7283311) q[3];
sx q[3];
rz(-2.9063291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94238344) q[0];
sx q[0];
rz(-2.1114045) q[0];
sx q[0];
rz(-0.27401608) q[0];
rz(0.38385299) q[1];
sx q[1];
rz(-0.63261384) q[1];
sx q[1];
rz(1.8208108) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21652554) q[0];
sx q[0];
rz(-1.4122047) q[0];
sx q[0];
rz(0.065404371) q[0];
rz(-pi) q[1];
rz(1.3396367) q[2];
sx q[2];
rz(-2.0589239) q[2];
sx q[2];
rz(1.3991935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6453041) q[1];
sx q[1];
rz(-2.4170365) q[1];
sx q[1];
rz(-1.4990119) q[1];
rz(-0.78944309) q[3];
sx q[3];
rz(-0.94067162) q[3];
sx q[3];
rz(-2.2342145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5968093) q[2];
sx q[2];
rz(-1.4027255) q[2];
sx q[2];
rz(1.7421494) q[2];
rz(-2.5341189) q[3];
sx q[3];
rz(-1.6739269) q[3];
sx q[3];
rz(-2.7510711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6121599) q[0];
sx q[0];
rz(-0.85560346) q[0];
sx q[0];
rz(-2.8223619) q[0];
rz(0.051479738) q[1];
sx q[1];
rz(-1.1795283) q[1];
sx q[1];
rz(1.0565588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38504544) q[0];
sx q[0];
rz(-0.17160417) q[0];
sx q[0];
rz(-1.8770939) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0931937) q[2];
sx q[2];
rz(-2.2757287) q[2];
sx q[2];
rz(2.8466109) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.776538) q[1];
sx q[1];
rz(-0.94887304) q[1];
sx q[1];
rz(-2.9261241) q[1];
rz(-pi) q[2];
rz(1.7848739) q[3];
sx q[3];
rz(-0.41164648) q[3];
sx q[3];
rz(0.95142196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1928548) q[2];
sx q[2];
rz(-1.7098018) q[2];
sx q[2];
rz(1.7970435) q[2];
rz(0.91790849) q[3];
sx q[3];
rz(-2.5036006) q[3];
sx q[3];
rz(-1.2558254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767839) q[0];
sx q[0];
rz(-0.5905686) q[0];
sx q[0];
rz(1.508924) q[0];
rz(2.6347939) q[1];
sx q[1];
rz(-1.2134039) q[1];
sx q[1];
rz(-0.39055821) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2029501) q[0];
sx q[0];
rz(-3.1331867) q[0];
sx q[0];
rz(2.6132843) q[0];
rz(-2.7616169) q[2];
sx q[2];
rz(-1.3936498) q[2];
sx q[2];
rz(0.65922696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11568053) q[1];
sx q[1];
rz(-1.5424219) q[1];
sx q[1];
rz(2.1898153) q[1];
rz(0.15843491) q[3];
sx q[3];
rz(-1.6601175) q[3];
sx q[3];
rz(-1.6824555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8916696) q[2];
sx q[2];
rz(-0.93081433) q[2];
sx q[2];
rz(0.8963975) q[2];
rz(2.2253288) q[3];
sx q[3];
rz(-0.93583411) q[3];
sx q[3];
rz(2.1808482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59178281) q[0];
sx q[0];
rz(-2.7857605) q[0];
sx q[0];
rz(-0.47019666) q[0];
rz(1.4037941) q[1];
sx q[1];
rz(-0.70039582) q[1];
sx q[1];
rz(-1.9139003) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5915697) q[0];
sx q[0];
rz(-1.3006905) q[0];
sx q[0];
rz(3.1349934) q[0];
rz(-pi) q[1];
x q[1];
rz(2.589354) q[2];
sx q[2];
rz(-1.9408014) q[2];
sx q[2];
rz(0.13665119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44547651) q[1];
sx q[1];
rz(-1.646478) q[1];
sx q[1];
rz(-0.92544011) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8813666) q[3];
sx q[3];
rz(-1.0538007) q[3];
sx q[3];
rz(-0.22634144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97051835) q[2];
sx q[2];
rz(-0.71275622) q[2];
sx q[2];
rz(-1.7573382) q[2];
rz(-0.61128831) q[3];
sx q[3];
rz(-1.3795815) q[3];
sx q[3];
rz(2.7454564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1432081) q[0];
sx q[0];
rz(-2.12976) q[0];
sx q[0];
rz(0.55214733) q[0];
rz(0.72969189) q[1];
sx q[1];
rz(-0.98760664) q[1];
sx q[1];
rz(-1.0017576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2524179) q[0];
sx q[0];
rz(-1.4902835) q[0];
sx q[0];
rz(-1.4371055) q[0];
rz(-0.4977141) q[2];
sx q[2];
rz(-0.98955284) q[2];
sx q[2];
rz(1.1377942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3730091) q[1];
sx q[1];
rz(-0.78861744) q[1];
sx q[1];
rz(-2.8810701) q[1];
rz(-0.61638636) q[3];
sx q[3];
rz(-0.77993837) q[3];
sx q[3];
rz(0.15579789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.45727649) q[2];
sx q[2];
rz(-1.3490973) q[2];
sx q[2];
rz(1.7982193) q[2];
rz(-2.4172879) q[3];
sx q[3];
rz(-1.2035921) q[3];
sx q[3];
rz(-2.8389285) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2398719) q[0];
sx q[0];
rz(-1.3322823) q[0];
sx q[0];
rz(-0.72541642) q[0];
rz(1.8431009) q[1];
sx q[1];
rz(-2.6893171) q[1];
sx q[1];
rz(0.27870146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50510079) q[0];
sx q[0];
rz(-1.0351702) q[0];
sx q[0];
rz(-1.3938175) q[0];
rz(1.1737661) q[2];
sx q[2];
rz(-1.8788998) q[2];
sx q[2];
rz(-2.7163519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52744059) q[1];
sx q[1];
rz(-0.88178712) q[1];
sx q[1];
rz(2.7428738) q[1];
x q[2];
rz(0.82741957) q[3];
sx q[3];
rz(-1.9636834) q[3];
sx q[3];
rz(3.06638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7574888) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68895483) q[0];
sx q[0];
rz(-2.3404558) q[0];
sx q[0];
rz(-2.5222006) q[0];
rz(-0.014892189) q[1];
sx q[1];
rz(-0.89034021) q[1];
sx q[1];
rz(0.68170086) q[1];
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
x q[1];
rz(-1.4686844) q[2];
sx q[2];
rz(-2.4739304) q[2];
sx q[2];
rz(-0.84102902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8047844) q[1];
sx q[1];
rz(-2.7899335) q[1];
sx q[1];
rz(0.019703948) q[1];
rz(-pi) q[2];
rz(-0.83311773) q[3];
sx q[3];
rz(-1.7135812) q[3];
sx q[3];
rz(-0.24417434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90998489) q[2];
sx q[2];
rz(-2.3384428) q[2];
sx q[2];
rz(-2.8019359) q[2];
rz(0.61399442) q[3];
sx q[3];
rz(-1.3795779) q[3];
sx q[3];
rz(-1.5617255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851819) q[0];
sx q[0];
rz(-2.7934533) q[0];
sx q[0];
rz(-2.2871995) q[0];
rz(-0.59335452) q[1];
sx q[1];
rz(-2.1095468) q[1];
sx q[1];
rz(-0.29534435) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4446553) q[0];
sx q[0];
rz(-1.5155223) q[0];
sx q[0];
rz(1.6037206) q[0];
rz(1.9885485) q[2];
sx q[2];
rz(-1.5502073) q[2];
sx q[2];
rz(1.127224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.01646811) q[1];
sx q[1];
rz(-1.09519) q[1];
sx q[1];
rz(-0.7866025) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2611102) q[3];
sx q[3];
rz(-1.4638374) q[3];
sx q[3];
rz(2.2992319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9796543) q[2];
sx q[2];
rz(-1.4583959) q[2];
sx q[2];
rz(-2.4597607) q[2];
rz(-1.9409059) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(0.43584263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1107263) q[0];
sx q[0];
rz(-1.8616572) q[0];
sx q[0];
rz(-2.6357292) q[0];
rz(-2.1234296) q[1];
sx q[1];
rz(-1.4351427) q[1];
sx q[1];
rz(-0.86046576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8448526) q[0];
sx q[0];
rz(-2.2460947) q[0];
sx q[0];
rz(1.6313511) q[0];
rz(-pi) q[1];
x q[1];
rz(1.127709) q[2];
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
rz(-1.7818393) q[1];
sx q[1];
rz(-1.5379496) q[1];
sx q[1];
rz(0.24928045) q[1];
rz(2.5768336) q[3];
sx q[3];
rz(-1.0339206) q[3];
sx q[3];
rz(2.2220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0597509) q[2];
sx q[2];
rz(-1.9530719) q[2];
sx q[2];
rz(2.5583978) q[2];
rz(-0.0010631337) q[3];
sx q[3];
rz(-0.93102264) q[3];
sx q[3];
rz(0.49334905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6625593) q[0];
sx q[0];
rz(-1.722535) q[0];
sx q[0];
rz(-2.2395635) q[0];
rz(2.5262911) q[1];
sx q[1];
rz(-1.4882144) q[1];
sx q[1];
rz(-1.8019713) q[1];
rz(-0.42846305) q[2];
sx q[2];
rz(-0.56461538) q[2];
sx q[2];
rz(2.3587894) q[2];
rz(-2.9499182) q[3];
sx q[3];
rz(-2.1379466) q[3];
sx q[3];
rz(2.9320075) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
