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
rz(0.69831508) q[0];
sx q[0];
rz(2.7628216) q[0];
sx q[0];
rz(8.2674352) q[0];
rz(10.209822) q[1];
sx q[1];
rz(0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049126712) q[0];
sx q[0];
rz(-1.8688339) q[0];
sx q[0];
rz(-2.1597693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73424642) q[2];
sx q[2];
rz(-2.1484005) q[2];
sx q[2];
rz(-0.89370382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68702811) q[1];
sx q[1];
rz(-0.44378456) q[1];
sx q[1];
rz(-0.084974809) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.111638) q[3];
sx q[3];
rz(-2.1308793) q[3];
sx q[3];
rz(-0.64489472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.779125) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(-2.989952) q[2];
rz(-2.9996297) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66459429) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(2.3240996) q[0];
rz(-0.11257653) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(2.2419825) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0257638) q[0];
sx q[0];
rz(-1.2983822) q[0];
sx q[0];
rz(1.8194356) q[0];
rz(-pi) q[1];
rz(2.4075872) q[2];
sx q[2];
rz(-0.81799928) q[2];
sx q[2];
rz(-0.98246511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1020722) q[1];
sx q[1];
rz(-1.5145434) q[1];
sx q[1];
rz(-2.392676) q[1];
rz(-pi) q[2];
rz(0.95590653) q[3];
sx q[3];
rz(-2.3645176) q[3];
sx q[3];
rz(-1.9333378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4704935) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(1.0279083) q[2];
rz(2.4043064) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(0.024287311) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1360433) q[0];
sx q[0];
rz(-0.5558973) q[0];
sx q[0];
rz(1.1935724) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(1.6568291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.68736) q[0];
sx q[0];
rz(-3.0890565) q[0];
sx q[0];
rz(2.0665069) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3593) q[2];
sx q[2];
rz(-1.4471608) q[2];
sx q[2];
rz(2.5619626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.130467) q[1];
sx q[1];
rz(-0.36917205) q[1];
sx q[1];
rz(0.38291957) q[1];
x q[2];
rz(1.3562702) q[3];
sx q[3];
rz(-2.7709922) q[3];
sx q[3];
rz(1.724898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16126157) q[2];
sx q[2];
rz(-1.223246) q[2];
sx q[2];
rz(2.4397591) q[2];
rz(-3.0724604) q[3];
sx q[3];
rz(-0.88874236) q[3];
sx q[3];
rz(0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0858916) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(-2.5162146) q[0];
rz(-1.0944132) q[1];
sx q[1];
rz(-1.5233327) q[1];
sx q[1];
rz(-1.3166924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65526456) q[0];
sx q[0];
rz(-1.0287971) q[0];
sx q[0];
rz(-0.16714759) q[0];
rz(0.28202348) q[2];
sx q[2];
rz(-2.8057753) q[2];
sx q[2];
rz(-1.5329602) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.229828) q[1];
sx q[1];
rz(-0.87512866) q[1];
sx q[1];
rz(1.3352446) q[1];
rz(0.98424498) q[3];
sx q[3];
rz(-2.2797425) q[3];
sx q[3];
rz(1.1701442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40654287) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(-3.1114846) q[2];
rz(0.41664577) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(-3.0863975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77867126) q[0];
sx q[0];
rz(-1.2782949) q[0];
sx q[0];
rz(1.9238506) q[0];
rz(2.9171004) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(2.7405558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825077) q[0];
sx q[0];
rz(-1.9422724) q[0];
sx q[0];
rz(-2.0187733) q[0];
rz(-pi) q[1];
rz(-1.2000167) q[2];
sx q[2];
rz(-0.94964281) q[2];
sx q[2];
rz(-2.7948684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2485473) q[1];
sx q[1];
rz(-1.5233938) q[1];
sx q[1];
rz(2.9045312) q[1];
rz(-pi) q[2];
rz(2.5878536) q[3];
sx q[3];
rz(-2.9512292) q[3];
sx q[3];
rz(-2.7680264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84247056) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(2.6341338) q[2];
rz(-0.92042813) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(-0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64666635) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(-2.1221509) q[0];
rz(1.9693718) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(1.2219465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2766314) q[0];
sx q[0];
rz(-2.2233133) q[0];
sx q[0];
rz(-0.44256532) q[0];
rz(2.0856306) q[2];
sx q[2];
rz(-1.5400572) q[2];
sx q[2];
rz(-1.6248425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6060562) q[1];
sx q[1];
rz(-2.5185985) q[1];
sx q[1];
rz(0.48320233) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29767146) q[3];
sx q[3];
rz(-2.1250013) q[3];
sx q[3];
rz(-0.25458529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6601861) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(-1.0931724) q[2];
rz(0.53705755) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(-0.66948906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401232) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(1.4599266) q[0];
rz(1.5096674) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(-3.0622838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96937856) q[0];
sx q[0];
rz(-1.8076118) q[0];
sx q[0];
rz(0.083061465) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4334636) q[2];
sx q[2];
rz(-1.3690565) q[2];
sx q[2];
rz(-1.84795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39187688) q[1];
sx q[1];
rz(-1.2971216) q[1];
sx q[1];
rz(2.7065212) q[1];
rz(2.4238911) q[3];
sx q[3];
rz(-1.6738179) q[3];
sx q[3];
rz(-0.75005248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0038393) q[2];
sx q[2];
rz(-1.1649818) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(-1.8462935) q[3];
sx q[3];
rz(-0.92419878) q[3];
sx q[3];
rz(-0.41954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9846648) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(0.28537634) q[0];
rz(3.0889619) q[1];
sx q[1];
rz(-1.4080518) q[1];
sx q[1];
rz(-0.1836798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7516458) q[0];
sx q[0];
rz(-3.0245298) q[0];
sx q[0];
rz(-1.879891) q[0];
rz(-pi) q[1];
rz(1.4168315) q[2];
sx q[2];
rz(-1.8530669) q[2];
sx q[2];
rz(-2.8090614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.86183) q[1];
sx q[1];
rz(-2.1643736) q[1];
sx q[1];
rz(1.4788126) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7643572) q[3];
sx q[3];
rz(-0.68359112) q[3];
sx q[3];
rz(2.205276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8949184) q[2];
sx q[2];
rz(-1.4072714) q[2];
sx q[2];
rz(-1.0651917) q[2];
rz(1.4554321) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(2.5559032) q[3];
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
rz(-0.10373779) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(1.6023585) q[0];
rz(2.1922951) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.4303713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1827138) q[0];
sx q[0];
rz(-2.4005167) q[0];
sx q[0];
rz(0.7345906) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2995641) q[2];
sx q[2];
rz(-1.6269738) q[2];
sx q[2];
rz(-2.2420355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0956968) q[1];
sx q[1];
rz(-1.8664845) q[1];
sx q[1];
rz(-1.6017833) q[1];
x q[2];
rz(1.0046047) q[3];
sx q[3];
rz(-0.49612415) q[3];
sx q[3];
rz(-1.8593018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1324233) q[2];
sx q[2];
rz(-1.8370266) q[2];
sx q[2];
rz(3.0150748) q[2];
rz(0.33454076) q[3];
sx q[3];
rz(-2.8821475) q[3];
sx q[3];
rz(-2.2693966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3009406) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(-0.51710039) q[0];
rz(0.05489796) q[1];
sx q[1];
rz(-1.5093191) q[1];
sx q[1];
rz(3.0518234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9835998) q[0];
sx q[0];
rz(-1.9741892) q[0];
sx q[0];
rz(2.1703815) q[0];
x q[1];
rz(-3.0788306) q[2];
sx q[2];
rz(-2.1307935) q[2];
sx q[2];
rz(-0.43207174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.43347142) q[1];
sx q[1];
rz(-0.29954391) q[1];
sx q[1];
rz(-0.56510651) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9338984) q[3];
sx q[3];
rz(-1.3570367) q[3];
sx q[3];
rz(0.0068706415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6101997) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(0.72719491) q[2];
rz(0.7555035) q[3];
sx q[3];
rz(-0.38755363) q[3];
sx q[3];
rz(2.9288647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.94201921) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(-1.3923116) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(-3.0346995) q[2];
sx q[2];
rz(-2.2724367) q[2];
sx q[2];
rz(2.5913749) q[2];
rz(-0.94115067) q[3];
sx q[3];
rz(-2.0654021) q[3];
sx q[3];
rz(0.13765814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
