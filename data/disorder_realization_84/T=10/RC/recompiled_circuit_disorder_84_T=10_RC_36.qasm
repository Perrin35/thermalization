OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(-2.3875561) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(1.1608646) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0275118) q[0];
sx q[0];
rz(-2.049276) q[0];
sx q[0];
rz(1.5374684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8765175) q[2];
sx q[2];
rz(-0.49760488) q[2];
sx q[2];
rz(1.9203609) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7293284) q[1];
sx q[1];
rz(-2.4630479) q[1];
sx q[1];
rz(-0.88339427) q[1];
x q[2];
rz(-1.4046895) q[3];
sx q[3];
rz(-1.3940092) q[3];
sx q[3];
rz(0.12648957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(2.6317821) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.3382834) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1424944) q[0];
sx q[0];
rz(-0.72744838) q[0];
sx q[0];
rz(-2.8508458) q[0];
rz(-0.83805214) q[2];
sx q[2];
rz(-1.8997955) q[2];
sx q[2];
rz(2.7555562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6515793) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(-0.3117805) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9017341) q[3];
sx q[3];
rz(-2.5407102) q[3];
sx q[3];
rz(2.2440653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(2.1792049) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798379) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-2.8082074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282368) q[0];
sx q[0];
rz(-2.2524999) q[0];
sx q[0];
rz(2.8744254) q[0];
rz(-pi) q[1];
rz(-1.1836428) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(1.7533592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15409878) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(3.0497453) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7942579) q[3];
sx q[3];
rz(-0.94933214) q[3];
sx q[3];
rz(-2.7623451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(-2.9004167) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(-0.10087092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872831) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(0.56418443) q[0];
rz(-2.5634649) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(2.633458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3343737) q[0];
sx q[0];
rz(-0.59893543) q[0];
sx q[0];
rz(-2.3301624) q[0];
rz(-1.5925243) q[2];
sx q[2];
rz(-0.93291514) q[2];
sx q[2];
rz(-0.69597543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6259039) q[1];
sx q[1];
rz(-1.3090033) q[1];
sx q[1];
rz(-0.7598676) q[1];
rz(-0.18211256) q[3];
sx q[3];
rz(-1.84387) q[3];
sx q[3];
rz(-1.9338716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-2.508146) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(0.19609837) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(-2.0702147) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8349985) q[0];
sx q[0];
rz(-0.95933611) q[0];
sx q[0];
rz(-2.6131265) q[0];
rz(-1.1159665) q[2];
sx q[2];
rz(-0.78274667) q[2];
sx q[2];
rz(-0.46403971) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(0.96697076) q[1];
rz(2.3134872) q[3];
sx q[3];
rz(-0.58907408) q[3];
sx q[3];
rz(2.7554054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(0.98199797) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.8765607) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(0.21662724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075561698) q[0];
sx q[0];
rz(-0.78853411) q[0];
sx q[0];
rz(-1.9476452) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40041133) q[2];
sx q[2];
rz(-2.4715804) q[2];
sx q[2];
rz(-0.77202051) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10303282) q[1];
sx q[1];
rz(-1.235504) q[1];
sx q[1];
rz(-0.0086987728) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2413007) q[3];
sx q[3];
rz(-1.7219464) q[3];
sx q[3];
rz(-2.8063262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(-1.023863) q[2];
rz(2.2653545) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(-0.52893692) q[0];
rz(1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(2.0524009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3750317) q[0];
sx q[0];
rz(-1.8683109) q[0];
sx q[0];
rz(-3.1356698) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94605298) q[2];
sx q[2];
rz(-0.81753987) q[2];
sx q[2];
rz(-2.3600876) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.293922) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(-3.1084204) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.180021) q[3];
sx q[3];
rz(-0.15768356) q[3];
sx q[3];
rz(-2.3458583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(1.4922173) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(-2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(-0.067226974) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(-2.9464088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0474437) q[0];
sx q[0];
rz(-2.7820245) q[0];
sx q[0];
rz(0.6028428) q[0];
rz(0.40114258) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(-2.0502979) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7797797) q[1];
sx q[1];
rz(-2.7098512) q[1];
sx q[1];
rz(-1.4713431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94594749) q[3];
sx q[3];
rz(-0.56560707) q[3];
sx q[3];
rz(0.92029508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(-2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72682196) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(-3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(-2.7862766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32896341) q[0];
sx q[0];
rz(-1.2041429) q[0];
sx q[0];
rz(3.0366412) q[0];
rz(-2.0673429) q[2];
sx q[2];
rz(-1.3585919) q[2];
sx q[2];
rz(0.70763904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4081817) q[1];
sx q[1];
rz(-1.2473885) q[1];
sx q[1];
rz(-0.059493382) q[1];
rz(-1.7570417) q[3];
sx q[3];
rz(-2.1631141) q[3];
sx q[3];
rz(-2.4663962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(3.1353531) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.86826098) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(1.2982752) q[0];
rz(0.4459933) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(2.840852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7556831) q[0];
sx q[0];
rz(-1.3433546) q[0];
sx q[0];
rz(-1.0230416) q[0];
rz(-1.1624596) q[2];
sx q[2];
rz(-1.0620585) q[2];
sx q[2];
rz(0.12847729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86665323) q[1];
sx q[1];
rz(-1.5870952) q[1];
sx q[1];
rz(-0.29124041) q[1];
x q[2];
rz(-2.7281076) q[3];
sx q[3];
rz(-0.95260145) q[3];
sx q[3];
rz(1.9181044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62844244) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-1.1395617) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(-1.4032455) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(-1.5724814) q[2];
sx q[2];
rz(-1.5148074) q[2];
sx q[2];
rz(-2.2956216) q[2];
rz(-1.2107522) q[3];
sx q[3];
rz(-2.4891709) q[3];
sx q[3];
rz(-1.8579033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];