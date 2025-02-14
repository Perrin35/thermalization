OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(2.6074183) q[0];
sx q[0];
rz(8.3745126) q[0];
rz(-1.2644816) q[1];
sx q[1];
rz(-2.2883132) q[1];
sx q[1];
rz(2.5929911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48196822) q[0];
sx q[0];
rz(-0.42185703) q[0];
sx q[0];
rz(-2.3530988) q[0];
rz(-pi) q[1];
rz(1.5919551) q[2];
sx q[2];
rz(-0.97351626) q[2];
sx q[2];
rz(2.6279272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2254991) q[1];
sx q[1];
rz(-1.0584944) q[1];
sx q[1];
rz(0.47173791) q[1];
rz(2.2150346) q[3];
sx q[3];
rz(-0.38739714) q[3];
sx q[3];
rz(0.54033632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7301664) q[2];
sx q[2];
rz(-2.7226518) q[2];
sx q[2];
rz(-2.1298998) q[2];
rz(2.2569979) q[3];
sx q[3];
rz(-1.1583637) q[3];
sx q[3];
rz(-1.2456892) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0678134) q[0];
sx q[0];
rz(-0.99869204) q[0];
sx q[0];
rz(2.3538537) q[0];
rz(-2.1630321) q[1];
sx q[1];
rz(-1.9906882) q[1];
sx q[1];
rz(1.2501134) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.015991) q[0];
sx q[0];
rz(-2.0569394) q[0];
sx q[0];
rz(-3.019633) q[0];
rz(-pi) q[1];
rz(-2.8626094) q[2];
sx q[2];
rz(-1.9595385) q[2];
sx q[2];
rz(1.1065567) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82455222) q[1];
sx q[1];
rz(-2.3352156) q[1];
sx q[1];
rz(0.033407465) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9216209) q[3];
sx q[3];
rz(-1.7860054) q[3];
sx q[3];
rz(0.53088354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0222187) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(3.019943) q[2];
rz(2.4308448) q[3];
sx q[3];
rz(-2.9080279) q[3];
sx q[3];
rz(2.5392883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1043333) q[0];
sx q[0];
rz(-0.72838825) q[0];
sx q[0];
rz(-0.65993586) q[0];
rz(2.624699) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(-1.0924115) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7832062) q[0];
sx q[0];
rz(-0.45166812) q[0];
sx q[0];
rz(-0.93675254) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2683892) q[2];
sx q[2];
rz(-1.3371468) q[2];
sx q[2];
rz(0.2327118) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9398492) q[1];
sx q[1];
rz(-2.1437316) q[1];
sx q[1];
rz(0.64676379) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64869439) q[3];
sx q[3];
rz(-1.3645384) q[3];
sx q[3];
rz(1.0252531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2163781) q[2];
sx q[2];
rz(-2.7984012) q[2];
sx q[2];
rz(1.7197616) q[2];
rz(-2.6575994) q[3];
sx q[3];
rz(-1.4839987) q[3];
sx q[3];
rz(1.7234195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5928818) q[0];
sx q[0];
rz(-0.084267862) q[0];
sx q[0];
rz(1.8362057) q[0];
rz(3.1343754) q[1];
sx q[1];
rz(-2.9199298) q[1];
sx q[1];
rz(-0.73297393) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6306046) q[0];
sx q[0];
rz(-2.9633491) q[0];
sx q[0];
rz(2.3018738) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8470331) q[2];
sx q[2];
rz(-1.6798875) q[2];
sx q[2];
rz(2.8677169) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8165255) q[1];
sx q[1];
rz(-0.96615929) q[1];
sx q[1];
rz(2.6139392) q[1];
rz(0.72317883) q[3];
sx q[3];
rz(-1.6537602) q[3];
sx q[3];
rz(-1.7570329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2923773) q[2];
sx q[2];
rz(-1.7376309) q[2];
sx q[2];
rz(0.89548573) q[2];
rz(2.605947) q[3];
sx q[3];
rz(-1.5629385) q[3];
sx q[3];
rz(-3.079788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8931005) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(-0.50791159) q[0];
rz(1.6912564) q[1];
sx q[1];
rz(-2.129887) q[1];
sx q[1];
rz(1.7844261) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0944509) q[0];
sx q[0];
rz(-0.12221065) q[0];
sx q[0];
rz(-2.727174) q[0];
rz(-pi) q[1];
rz(-1.1512027) q[2];
sx q[2];
rz(-1.1418248) q[2];
sx q[2];
rz(-1.2871413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0101945) q[1];
sx q[1];
rz(-1.4631738) q[1];
sx q[1];
rz(-1.1788998) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97220631) q[3];
sx q[3];
rz(-1.8892885) q[3];
sx q[3];
rz(1.1725669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42673972) q[2];
sx q[2];
rz(-1.667495) q[2];
sx q[2];
rz(1.1023785) q[2];
rz(-2.0057996) q[3];
sx q[3];
rz(-1.9250684) q[3];
sx q[3];
rz(-2.3701325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8847467) q[0];
sx q[0];
rz(-2.1717635) q[0];
sx q[0];
rz(2.2127175) q[0];
rz(0.38201395) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(-1.4487723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3253052) q[0];
sx q[0];
rz(-1.0950118) q[0];
sx q[0];
rz(1.2289117) q[0];
rz(-pi) q[1];
rz(0.48488481) q[2];
sx q[2];
rz(-1.6248091) q[2];
sx q[2];
rz(1.206813) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.12940059) q[1];
sx q[1];
rz(-1.4621648) q[1];
sx q[1];
rz(1.9700178) q[1];
rz(1.9453064) q[3];
sx q[3];
rz(-1.4099965) q[3];
sx q[3];
rz(-0.2674777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53943071) q[2];
sx q[2];
rz(-2.7950037) q[2];
sx q[2];
rz(2.3740785) q[2];
rz(-1.8481988) q[3];
sx q[3];
rz(-0.69197217) q[3];
sx q[3];
rz(-0.20492157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2061737) q[0];
sx q[0];
rz(-0.17429166) q[0];
sx q[0];
rz(-2.8906004) q[0];
rz(2.9639066) q[1];
sx q[1];
rz(-1.4439986) q[1];
sx q[1];
rz(0.65690717) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85992438) q[0];
sx q[0];
rz(-0.31552464) q[0];
sx q[0];
rz(2.2918227) q[0];
x q[1];
rz(2.9365402) q[2];
sx q[2];
rz(-0.66407138) q[2];
sx q[2];
rz(2.1560706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6823947) q[1];
sx q[1];
rz(-0.49811813) q[1];
sx q[1];
rz(2.47704) q[1];
rz(-pi) q[2];
rz(-2.4961996) q[3];
sx q[3];
rz(-1.925996) q[3];
sx q[3];
rz(-1.3076289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3333007) q[2];
sx q[2];
rz(-2.5623645) q[2];
sx q[2];
rz(-0.91326886) q[2];
rz(-0.67780668) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(-1.0031797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0388357) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(-2.5586149) q[0];
rz(-1.673117) q[1];
sx q[1];
rz(-0.99248326) q[1];
sx q[1];
rz(-0.34117064) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.720168) q[0];
sx q[0];
rz(-2.9377794) q[0];
sx q[0];
rz(1.782062) q[0];
rz(-pi) q[1];
rz(-1.2874574) q[2];
sx q[2];
rz(-1.8036799) q[2];
sx q[2];
rz(2.8849059) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.063733) q[1];
sx q[1];
rz(-1.9394411) q[1];
sx q[1];
rz(-0.39409448) q[1];
rz(2.0604793) q[3];
sx q[3];
rz(-2.4016909) q[3];
sx q[3];
rz(2.5748753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1950281) q[2];
sx q[2];
rz(-0.82411689) q[2];
sx q[2];
rz(-0.80671802) q[2];
rz(-0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(2.257982) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44613999) q[0];
sx q[0];
rz(-2.0841632) q[0];
sx q[0];
rz(-2.9651508) q[0];
rz(2.8314619) q[1];
sx q[1];
rz(-1.6519494) q[1];
sx q[1];
rz(-2.2055221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061358) q[0];
sx q[0];
rz(-1.8671037) q[0];
sx q[0];
rz(1.4151876) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1175198) q[2];
sx q[2];
rz(-1.0220811) q[2];
sx q[2];
rz(-2.0224151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99451274) q[1];
sx q[1];
rz(-1.6638866) q[1];
sx q[1];
rz(-1.9398111) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8027169) q[3];
sx q[3];
rz(-2.6565564) q[3];
sx q[3];
rz(0.95279653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0940242) q[2];
sx q[2];
rz(-1.3022283) q[2];
sx q[2];
rz(-3.1330718) q[2];
rz(-0.73355567) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(3.0146397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215288) q[0];
sx q[0];
rz(-0.41291741) q[0];
sx q[0];
rz(2.371696) q[0];
rz(0.13433111) q[1];
sx q[1];
rz(-0.7901935) q[1];
sx q[1];
rz(-1.92164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3111356) q[0];
sx q[0];
rz(-2.756098) q[0];
sx q[0];
rz(0.66118447) q[0];
rz(1.8656857) q[2];
sx q[2];
rz(-0.98433009) q[2];
sx q[2];
rz(-0.49582729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7999477) q[1];
sx q[1];
rz(-1.1468256) q[1];
sx q[1];
rz(1.9862224) q[1];
rz(-pi) q[2];
rz(0.071000428) q[3];
sx q[3];
rz(-1.7216847) q[3];
sx q[3];
rz(-1.7689266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.28025815) q[2];
sx q[2];
rz(-0.65797776) q[2];
sx q[2];
rz(0.78557837) q[2];
rz(2.806459) q[3];
sx q[3];
rz(-0.98888713) q[3];
sx q[3];
rz(-0.82211632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.817374) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(2.4556976) q[1];
sx q[1];
rz(-2.5527725) q[1];
sx q[1];
rz(1.7935161) q[1];
rz(-1.5625207) q[2];
sx q[2];
rz(-0.58115056) q[2];
sx q[2];
rz(-1.5528408) q[2];
rz(3.104628) q[3];
sx q[3];
rz(-1.8633458) q[3];
sx q[3];
rz(2.9316791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
