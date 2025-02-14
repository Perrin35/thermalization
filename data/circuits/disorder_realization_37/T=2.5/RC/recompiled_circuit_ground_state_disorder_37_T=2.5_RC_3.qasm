OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3545544) q[0];
sx q[0];
rz(-2.3784502) q[0];
sx q[0];
rz(-0.95844498) q[0];
rz(-2.7148442) q[1];
sx q[1];
rz(3.5659748) q[1];
sx q[1];
rz(10.397484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63494191) q[0];
sx q[0];
rz(-1.6641064) q[0];
sx q[0];
rz(-0.20185457) q[0];
rz(1.6390267) q[2];
sx q[2];
rz(-1.3989324) q[2];
sx q[2];
rz(-1.1504988) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2610361) q[1];
sx q[1];
rz(-2.2734005) q[1];
sx q[1];
rz(-1.9126585) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28222059) q[3];
sx q[3];
rz(-1.422554) q[3];
sx q[3];
rz(2.6765598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.080231754) q[2];
sx q[2];
rz(-1.9539359) q[2];
sx q[2];
rz(0.80002552) q[2];
rz(3.0112265) q[3];
sx q[3];
rz(-2.3733807) q[3];
sx q[3];
rz(-0.31726328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3404563) q[0];
sx q[0];
rz(-1.2202593) q[0];
sx q[0];
rz(0.98943797) q[0];
rz(-0.89775741) q[1];
sx q[1];
rz(-1.0549301) q[1];
sx q[1];
rz(0.78136939) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1290478) q[0];
sx q[0];
rz(-0.54586403) q[0];
sx q[0];
rz(-1.9485628) q[0];
rz(-pi) q[1];
rz(-1.4825253) q[2];
sx q[2];
rz(-1.8039743) q[2];
sx q[2];
rz(-0.424343) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20618379) q[1];
sx q[1];
rz(-2.0061135) q[1];
sx q[1];
rz(2.101442) q[1];
x q[2];
rz(0.3831634) q[3];
sx q[3];
rz(-1.2943315) q[3];
sx q[3];
rz(-3.059946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.117131) q[2];
sx q[2];
rz(-1.7639561) q[2];
sx q[2];
rz(-0.77829877) q[2];
rz(2.4513169) q[3];
sx q[3];
rz(-2.2199151) q[3];
sx q[3];
rz(1.5604304) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86163259) q[0];
sx q[0];
rz(-2.3065688) q[0];
sx q[0];
rz(-0.41197187) q[0];
rz(2.1906134) q[1];
sx q[1];
rz(-2.488766) q[1];
sx q[1];
rz(-1.9297809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6854343) q[0];
sx q[0];
rz(-0.051101772) q[0];
sx q[0];
rz(0.85275485) q[0];
rz(-1.1583382) q[2];
sx q[2];
rz(-0.24246267) q[2];
sx q[2];
rz(-1.860577) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6784994) q[1];
sx q[1];
rz(-2.2804567) q[1];
sx q[1];
rz(-1.6072431) q[1];
rz(-pi) q[2];
rz(1.7318299) q[3];
sx q[3];
rz(-1.2098243) q[3];
sx q[3];
rz(-1.4396813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1872306) q[2];
sx q[2];
rz(-2.6095641) q[2];
sx q[2];
rz(-1.6383891) q[2];
rz(-3.1014118) q[3];
sx q[3];
rz(-0.92622042) q[3];
sx q[3];
rz(2.9068376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9784341) q[0];
sx q[0];
rz(-1.5143159) q[0];
sx q[0];
rz(0.12983313) q[0];
rz(1.8561329) q[1];
sx q[1];
rz(-1.1760271) q[1];
sx q[1];
rz(2.1066378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0823815) q[0];
sx q[0];
rz(-0.9810582) q[0];
sx q[0];
rz(2.608631) q[0];
x q[1];
rz(-0.1088684) q[2];
sx q[2];
rz(-0.90990674) q[2];
sx q[2];
rz(2.3488597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1521406) q[1];
sx q[1];
rz(-1.3229645) q[1];
sx q[1];
rz(1.4825443) q[1];
x q[2];
rz(-2.2620151) q[3];
sx q[3];
rz(-1.1926418) q[3];
sx q[3];
rz(-2.27431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7776362) q[2];
sx q[2];
rz(-1.7549606) q[2];
sx q[2];
rz(-0.094495471) q[2];
rz(0.42094055) q[3];
sx q[3];
rz(-1.4237483) q[3];
sx q[3];
rz(-1.7054935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41973758) q[0];
sx q[0];
rz(-2.5696745) q[0];
sx q[0];
rz(-3.1255334) q[0];
rz(2.2968966) q[1];
sx q[1];
rz(-0.41929308) q[1];
sx q[1];
rz(-2.2339581) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3428872) q[0];
sx q[0];
rz(-1.3488975) q[0];
sx q[0];
rz(0.7733487) q[0];
rz(0.24036562) q[2];
sx q[2];
rz(-1.9800926) q[2];
sx q[2];
rz(1.586692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6542977) q[1];
sx q[1];
rz(-2.0072486) q[1];
sx q[1];
rz(0.36466332) q[1];
rz(-3.0673712) q[3];
sx q[3];
rz(-0.34498131) q[3];
sx q[3];
rz(0.94946194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5661261) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(-3.0779085) q[2];
rz(-0.56695402) q[3];
sx q[3];
rz(-0.83438116) q[3];
sx q[3];
rz(-0.59757346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8276234) q[0];
sx q[0];
rz(-2.9997928) q[0];
sx q[0];
rz(1.9837448) q[0];
rz(-1.4843548) q[1];
sx q[1];
rz(-1.6969705) q[1];
sx q[1];
rz(-2.8725913) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7662748) q[0];
sx q[0];
rz(-0.14785375) q[0];
sx q[0];
rz(0.30884858) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9018752) q[2];
sx q[2];
rz(-0.28713687) q[2];
sx q[2];
rz(1.3763217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8527166) q[1];
sx q[1];
rz(-1.7501131) q[1];
sx q[1];
rz(-0.10097558) q[1];
x q[2];
rz(1.7888467) q[3];
sx q[3];
rz(-1.9761128) q[3];
sx q[3];
rz(2.1637555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3645662) q[2];
sx q[2];
rz(-0.51743162) q[2];
sx q[2];
rz(0.88097921) q[2];
rz(0.12428728) q[3];
sx q[3];
rz(-1.4275987) q[3];
sx q[3];
rz(-2.4833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86150259) q[0];
sx q[0];
rz(-0.26743356) q[0];
sx q[0];
rz(2.4524443) q[0];
rz(0.46734494) q[1];
sx q[1];
rz(-2.5655949) q[1];
sx q[1];
rz(-1.0980094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2240967) q[0];
sx q[0];
rz(-2.9584329) q[0];
sx q[0];
rz(-0.44395019) q[0];
x q[1];
rz(-2.5937366) q[2];
sx q[2];
rz(-1.8286053) q[2];
sx q[2];
rz(-1.1950891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5913691) q[1];
sx q[1];
rz(-0.5412296) q[1];
sx q[1];
rz(3.079097) q[1];
rz(-pi) q[2];
rz(-2.1995898) q[3];
sx q[3];
rz(-1.9900454) q[3];
sx q[3];
rz(-0.85070953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.565862) q[2];
sx q[2];
rz(-2.2634759) q[2];
sx q[2];
rz(0.6305387) q[2];
rz(-0.60728836) q[3];
sx q[3];
rz(-2.2346965) q[3];
sx q[3];
rz(-1.0926532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31834114) q[0];
sx q[0];
rz(-2.9956151) q[0];
sx q[0];
rz(-1.3463705) q[0];
rz(1.775555) q[1];
sx q[1];
rz(-2.4925158) q[1];
sx q[1];
rz(-2.5128561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36489428) q[0];
sx q[0];
rz(-1.3140251) q[0];
sx q[0];
rz(-2.5904694) q[0];
rz(-pi) q[1];
rz(0.5127617) q[2];
sx q[2];
rz(-2.1760941) q[2];
sx q[2];
rz(0.61601854) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6639316) q[1];
sx q[1];
rz(-2.0545425) q[1];
sx q[1];
rz(-2.5238634) q[1];
rz(1.2090861) q[3];
sx q[3];
rz(-0.21245281) q[3];
sx q[3];
rz(-0.3750876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2748572) q[2];
sx q[2];
rz(-1.3254415) q[2];
sx q[2];
rz(-0.87230116) q[2];
rz(0.21271475) q[3];
sx q[3];
rz(-1.745696) q[3];
sx q[3];
rz(2.5581636) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9048731) q[0];
sx q[0];
rz(-1.5308335) q[0];
sx q[0];
rz(-1.2637631) q[0];
rz(-1.0172552) q[1];
sx q[1];
rz(-0.89824289) q[1];
sx q[1];
rz(0.57473007) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86704937) q[0];
sx q[0];
rz(-0.5349656) q[0];
sx q[0];
rz(-0.127129) q[0];
rz(2.3129914) q[2];
sx q[2];
rz(-1.7485837) q[2];
sx q[2];
rz(2.1091903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1803429) q[1];
sx q[1];
rz(-1.3836134) q[1];
sx q[1];
rz(-2.768152) q[1];
rz(2.0323864) q[3];
sx q[3];
rz(-2.6127041) q[3];
sx q[3];
rz(-2.703107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42402521) q[2];
sx q[2];
rz(-1.8342476) q[2];
sx q[2];
rz(3.1094816) q[2];
rz(1.1561681) q[3];
sx q[3];
rz(-0.38438946) q[3];
sx q[3];
rz(1.2110075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0403989) q[0];
sx q[0];
rz(-0.54553425) q[0];
sx q[0];
rz(-1.1241166) q[0];
rz(2.595937) q[1];
sx q[1];
rz(-0.35174313) q[1];
sx q[1];
rz(-2.5901332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41697793) q[0];
sx q[0];
rz(-1.7584086) q[0];
sx q[0];
rz(-3.0155988) q[0];
rz(-pi) q[1];
rz(1.7573171) q[2];
sx q[2];
rz(-2.8628446) q[2];
sx q[2];
rz(-2.2997625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.76186) q[1];
sx q[1];
rz(-1.1071536) q[1];
sx q[1];
rz(0.54370466) q[1];
rz(2.8352883) q[3];
sx q[3];
rz(-1.4147537) q[3];
sx q[3];
rz(-0.60221218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5896899) q[2];
sx q[2];
rz(-2.6305113) q[2];
sx q[2];
rz(-1.4116633) q[2];
rz(-0.74665135) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(0.97486973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31155561) q[0];
sx q[0];
rz(-1.5811601) q[0];
sx q[0];
rz(1.5720221) q[0];
rz(2.742782) q[1];
sx q[1];
rz(-1.8067982) q[1];
sx q[1];
rz(-1.6511818) q[1];
rz(-1.082303) q[2];
sx q[2];
rz(-1.563896) q[2];
sx q[2];
rz(0.36120216) q[2];
rz(-1.87181) q[3];
sx q[3];
rz(-2.3475937) q[3];
sx q[3];
rz(0.54318843) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
