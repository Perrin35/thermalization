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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(1.992835) q[0];
rz(-2.5454638) q[1];
sx q[1];
rz(-2.7355173) q[1];
sx q[1];
rz(-1.6869071) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5243149) q[0];
sx q[0];
rz(-1.9348659) q[0];
sx q[0];
rz(2.2979256) q[0];
rz(-pi) q[1];
rz(-3.0266255) q[2];
sx q[2];
rz(-0.76412725) q[2];
sx q[2];
rz(-0.72944356) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0815989) q[1];
sx q[1];
rz(-0.76548701) q[1];
sx q[1];
rz(-0.19123921) q[1];
rz(-3.0166854) q[3];
sx q[3];
rz(-2.4797399) q[3];
sx q[3];
rz(3.007909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9601606) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(-1.206548) q[2];
rz(-1.9801961) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(1.5709706) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777765) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(-2.1695082) q[0];
rz(-0.26845911) q[1];
sx q[1];
rz(-1.4413036) q[1];
sx q[1];
rz(-0.45480248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80538023) q[0];
sx q[0];
rz(-0.32355645) q[0];
sx q[0];
rz(-1.9791043) q[0];
rz(1.201638) q[2];
sx q[2];
rz(-2.7281986) q[2];
sx q[2];
rz(-1.8231376) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6356335) q[1];
sx q[1];
rz(-1.1065496) q[1];
sx q[1];
rz(1.5134769) q[1];
rz(2.11449) q[3];
sx q[3];
rz(-2.0683401) q[3];
sx q[3];
rz(3.080201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0316169) q[2];
sx q[2];
rz(-0.72828186) q[2];
sx q[2];
rz(-2.1419683) q[2];
rz(0.085974606) q[3];
sx q[3];
rz(-1.5558259) q[3];
sx q[3];
rz(-0.011367817) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407532) q[0];
sx q[0];
rz(-1.6735621) q[0];
sx q[0];
rz(2.9174347) q[0];
rz(-1.5466746) q[1];
sx q[1];
rz(-1.5432576) q[1];
sx q[1];
rz(1.3067783) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6431893) q[0];
sx q[0];
rz(-2.1191594) q[0];
sx q[0];
rz(1.4275309) q[0];
rz(2.3036912) q[2];
sx q[2];
rz(-0.97992491) q[2];
sx q[2];
rz(1.6333579) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0780454) q[1];
sx q[1];
rz(-1.1544466) q[1];
sx q[1];
rz(-1.1148808) q[1];
rz(2.8257224) q[3];
sx q[3];
rz(-2.4474553) q[3];
sx q[3];
rz(0.71198502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33009067) q[2];
sx q[2];
rz(-1.5218488) q[2];
sx q[2];
rz(2.1211993) q[2];
rz(-1.329782) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(-1.5248732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40600768) q[0];
sx q[0];
rz(-2.0841053) q[0];
sx q[0];
rz(-1.1757346) q[0];
rz(-0.27578393) q[1];
sx q[1];
rz(-0.83668721) q[1];
sx q[1];
rz(-0.63794678) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9383134) q[0];
sx q[0];
rz(-1.5814476) q[0];
sx q[0];
rz(3.108884) q[0];
rz(-0.3023382) q[2];
sx q[2];
rz(-2.4983642) q[2];
sx q[2];
rz(1.1272205) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1012816) q[1];
sx q[1];
rz(-2.4263546) q[1];
sx q[1];
rz(-0.50846993) q[1];
x q[2];
rz(0.27880554) q[3];
sx q[3];
rz(-2.8099217) q[3];
sx q[3];
rz(1.3607894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.63567579) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(1.6944616) q[2];
rz(0.21008374) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(-0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90784812) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(2.0727169) q[0];
rz(1.2999889) q[1];
sx q[1];
rz(-1.735382) q[1];
sx q[1];
rz(1.9357505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.962709) q[0];
sx q[0];
rz(-2.7307352) q[0];
sx q[0];
rz(-1.370048) q[0];
x q[1];
rz(-0.0095179518) q[2];
sx q[2];
rz(-2.2235907) q[2];
sx q[2];
rz(-1.643484) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1002378) q[1];
sx q[1];
rz(-2.2861028) q[1];
sx q[1];
rz(-3.0042786) q[1];
rz(-pi) q[2];
rz(0.48128328) q[3];
sx q[3];
rz(-1.4247155) q[3];
sx q[3];
rz(-1.8446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24870366) q[2];
sx q[2];
rz(-2.7619669) q[2];
sx q[2];
rz(-3.0730263) q[2];
rz(-0.93962234) q[3];
sx q[3];
rz(-1.7328123) q[3];
sx q[3];
rz(1.2287963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762961) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(-0.5214386) q[0];
rz(1.5154845) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(-2.5159786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2170048) q[0];
sx q[0];
rz(-0.5408322) q[0];
sx q[0];
rz(-1.5103673) q[0];
rz(-pi) q[1];
rz(-0.31452532) q[2];
sx q[2];
rz(-0.39681602) q[2];
sx q[2];
rz(-1.0971958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.292899) q[1];
sx q[1];
rz(-1.1176511) q[1];
sx q[1];
rz(2.921656) q[1];
rz(-pi) q[2];
rz(1.5099105) q[3];
sx q[3];
rz(-1.1495483) q[3];
sx q[3];
rz(2.6256068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1336512) q[2];
sx q[2];
rz(-2.5456754) q[2];
sx q[2];
rz(3.0360743) q[2];
rz(0.96406913) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(0.9849557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.069020011) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(-1.7286812) q[0];
rz(-2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(-0.55317318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39554292) q[0];
sx q[0];
rz(-2.2043214) q[0];
sx q[0];
rz(2.6468011) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4804391) q[2];
sx q[2];
rz(-1.7117662) q[2];
sx q[2];
rz(-0.80328926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87482086) q[1];
sx q[1];
rz(-1.8541341) q[1];
sx q[1];
rz(-1.6127519) q[1];
rz(-1.8314701) q[3];
sx q[3];
rz(-1.7672774) q[3];
sx q[3];
rz(0.5082265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46099123) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(-2.7294532) q[2];
rz(-2.4814217) q[3];
sx q[3];
rz(-0.5414525) q[3];
sx q[3];
rz(-2.6485543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93765813) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(2.9397021) q[0];
rz(-2.0837325) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(0.26022628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79151151) q[0];
sx q[0];
rz(-2.4790194) q[0];
sx q[0];
rz(2.4988079) q[0];
x q[1];
rz(1.0551532) q[2];
sx q[2];
rz(-1.9571575) q[2];
sx q[2];
rz(-0.48395115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59108666) q[1];
sx q[1];
rz(-2.6896853) q[1];
sx q[1];
rz(-2.4585548) q[1];
x q[2];
rz(1.4020355) q[3];
sx q[3];
rz(-0.91723727) q[3];
sx q[3];
rz(-0.26168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1883833) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(-0.58237135) q[2];
rz(-2.6992056) q[3];
sx q[3];
rz(-1.235639) q[3];
sx q[3];
rz(0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9198832) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(-1.7145994) q[0];
rz(1.7859979) q[1];
sx q[1];
rz(-1.6537063) q[1];
sx q[1];
rz(1.4367163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6528931) q[0];
sx q[0];
rz(-0.32036361) q[0];
sx q[0];
rz(1.5071177) q[0];
rz(-1.8709917) q[2];
sx q[2];
rz(-1.1087024) q[2];
sx q[2];
rz(1.0754881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.88594) q[1];
sx q[1];
rz(-2.8999834) q[1];
sx q[1];
rz(-1.5109946) q[1];
rz(-pi) q[2];
rz(-1.837265) q[3];
sx q[3];
rz(-2.5992166) q[3];
sx q[3];
rz(1.2790542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9068678) q[2];
sx q[2];
rz(-1.2805254) q[2];
sx q[2];
rz(1.5391763) q[2];
rz(-2.5336044) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(0.68979818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749087) q[0];
sx q[0];
rz(-2.4803949) q[0];
sx q[0];
rz(-2.3181424) q[0];
rz(0.89556328) q[1];
sx q[1];
rz(-1.7915553) q[1];
sx q[1];
rz(2.7604738) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6247647) q[0];
sx q[0];
rz(-2.4104558) q[0];
sx q[0];
rz(-0.69964377) q[0];
x q[1];
rz(-2.7398749) q[2];
sx q[2];
rz(-1.7251996) q[2];
sx q[2];
rz(-0.32333514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.95914074) q[1];
sx q[1];
rz(-0.63196665) q[1];
sx q[1];
rz(-2.0672634) q[1];
x q[2];
rz(-0.70638871) q[3];
sx q[3];
rz(-2.3800142) q[3];
sx q[3];
rz(2.0743845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33619189) q[2];
sx q[2];
rz(-0.41173428) q[2];
sx q[2];
rz(-2.1545048) q[2];
rz(-1.6030715) q[3];
sx q[3];
rz(-2.5542732) q[3];
sx q[3];
rz(0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22454746) q[0];
sx q[0];
rz(-1.2397091) q[0];
sx q[0];
rz(1.5201257) q[0];
rz(2.4317901) q[1];
sx q[1];
rz(-0.55955049) q[1];
sx q[1];
rz(1.4581663) q[1];
rz(-0.29322704) q[2];
sx q[2];
rz(-0.37920375) q[2];
sx q[2];
rz(-1.3398021) q[2];
rz(0.62623528) q[3];
sx q[3];
rz(-0.84624419) q[3];
sx q[3];
rz(-0.14583896) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
