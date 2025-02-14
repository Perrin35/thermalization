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
rz(0.74547493) q[0];
sx q[0];
rz(2.5479654) q[0];
sx q[0];
rz(9.0802703) q[0];
rz(-1.966882) q[1];
sx q[1];
rz(-0.36829683) q[1];
sx q[1];
rz(-0.60454291) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6942351) q[0];
sx q[0];
rz(-1.6295054) q[0];
sx q[0];
rz(0.82253463) q[0];
x q[1];
rz(-1.749332) q[2];
sx q[2];
rz(-1.5024619) q[2];
sx q[2];
rz(0.19388895) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.39554) q[1];
sx q[1];
rz(-2.722858) q[1];
sx q[1];
rz(0.0050425649) q[1];
rz(0.56960241) q[3];
sx q[3];
rz(-2.0853373) q[3];
sx q[3];
rz(-1.2698184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.979226) q[2];
sx q[2];
rz(-1.484885) q[2];
sx q[2];
rz(1.7171198) q[2];
rz(-0.10719565) q[3];
sx q[3];
rz(-2.9469979) q[3];
sx q[3];
rz(-2.181982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2694038) q[0];
sx q[0];
rz(-0.93282455) q[0];
sx q[0];
rz(-1.3099571) q[0];
rz(-0.45920363) q[1];
sx q[1];
rz(-1.2745067) q[1];
sx q[1];
rz(2.8244663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8175791) q[0];
sx q[0];
rz(-1.5234103) q[0];
sx q[0];
rz(3.0035278) q[0];
x q[1];
rz(-0.049811157) q[2];
sx q[2];
rz(-1.5731817) q[2];
sx q[2];
rz(1.1394237) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0786653) q[1];
sx q[1];
rz(-1.0243653) q[1];
sx q[1];
rz(1.648047) q[1];
rz(-0.22038118) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(3.0241682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7121048) q[2];
sx q[2];
rz(-1.8450582) q[2];
sx q[2];
rz(-0.8075766) q[2];
rz(-0.56337041) q[3];
sx q[3];
rz(-2.6475776) q[3];
sx q[3];
rz(1.0529244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.5092369) q[0];
sx q[0];
rz(-2.95166) q[0];
sx q[0];
rz(0.83455363) q[0];
rz(-0.50476152) q[1];
sx q[1];
rz(-2.7752462) q[1];
sx q[1];
rz(1.4588446) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6989649) q[0];
sx q[0];
rz(-2.1219398) q[0];
sx q[0];
rz(1.8056662) q[0];
x q[1];
rz(-0.73503942) q[2];
sx q[2];
rz(-0.95434084) q[2];
sx q[2];
rz(0.40555813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8275833) q[1];
sx q[1];
rz(-1.4928668) q[1];
sx q[1];
rz(0.72300903) q[1];
rz(0.97455131) q[3];
sx q[3];
rz(-2.2828492) q[3];
sx q[3];
rz(-2.2012868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0316281) q[2];
sx q[2];
rz(-1.6692903) q[2];
sx q[2];
rz(0.36299452) q[2];
rz(1.1686769) q[3];
sx q[3];
rz(-0.76255885) q[3];
sx q[3];
rz(0.76538435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4204243) q[0];
sx q[0];
rz(-0.52614251) q[0];
sx q[0];
rz(-0.59071937) q[0];
rz(-2.4144454) q[1];
sx q[1];
rz(-2.7666481) q[1];
sx q[1];
rz(2.4730543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0787857) q[0];
sx q[0];
rz(-1.6123839) q[0];
sx q[0];
rz(-1.6522626) q[0];
x q[1];
rz(0.48769571) q[2];
sx q[2];
rz(-0.82789153) q[2];
sx q[2];
rz(2.6726892) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82610735) q[1];
sx q[1];
rz(-2.25137) q[1];
sx q[1];
rz(-0.95113753) q[1];
rz(0.7410243) q[3];
sx q[3];
rz(-1.1644568) q[3];
sx q[3];
rz(-0.63593731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10355243) q[2];
sx q[2];
rz(-1.870564) q[2];
sx q[2];
rz(-1.5070149) q[2];
rz(-0.42190894) q[3];
sx q[3];
rz(-0.76015893) q[3];
sx q[3];
rz(2.4222477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5108532) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(-1.0605633) q[0];
rz(-0.063449055) q[1];
sx q[1];
rz(-0.63065204) q[1];
sx q[1];
rz(-2.5877156) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83186524) q[0];
sx q[0];
rz(-0.99156443) q[0];
sx q[0];
rz(-2.6088891) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80762005) q[2];
sx q[2];
rz(-2.3381026) q[2];
sx q[2];
rz(2.9832911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3556695) q[1];
sx q[1];
rz(-0.3057978) q[1];
sx q[1];
rz(3.0088916) q[1];
rz(-1.2451671) q[3];
sx q[3];
rz(-1.6998597) q[3];
sx q[3];
rz(1.8092524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8744897) q[2];
sx q[2];
rz(-2.5449982) q[2];
sx q[2];
rz(2.5388517) q[2];
rz(0.069325773) q[3];
sx q[3];
rz(-1.5963138) q[3];
sx q[3];
rz(-2.9749405) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9423187) q[0];
sx q[0];
rz(-0.60180226) q[0];
sx q[0];
rz(-2.6875575) q[0];
rz(1.8858689) q[1];
sx q[1];
rz(-2.1818826) q[1];
sx q[1];
rz(1.4383291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5944477) q[0];
sx q[0];
rz(-1.9701574) q[0];
sx q[0];
rz(-0.85712503) q[0];
x q[1];
rz(-2.8642162) q[2];
sx q[2];
rz(-1.8581333) q[2];
sx q[2];
rz(0.6544906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0918947) q[1];
sx q[1];
rz(-2.0573074) q[1];
sx q[1];
rz(2.2177433) q[1];
x q[2];
rz(1.0694396) q[3];
sx q[3];
rz(-0.92235288) q[3];
sx q[3];
rz(-1.4251777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85256514) q[2];
sx q[2];
rz(-1.8708517) q[2];
sx q[2];
rz(1.3555869) q[2];
rz(-1.2191314) q[3];
sx q[3];
rz(-1.6617323) q[3];
sx q[3];
rz(1.5752972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8236302) q[0];
sx q[0];
rz(-0.83331236) q[0];
sx q[0];
rz(1.747636) q[0];
rz(0.45669237) q[1];
sx q[1];
rz(-0.87867457) q[1];
sx q[1];
rz(-1.7535694) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7308424) q[0];
sx q[0];
rz(-1.1880945) q[0];
sx q[0];
rz(2.9992473) q[0];
rz(2.9406383) q[2];
sx q[2];
rz(-1.4575053) q[2];
sx q[2];
rz(-2.410881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2051831) q[1];
sx q[1];
rz(-0.35297063) q[1];
sx q[1];
rz(-0.38683968) q[1];
x q[2];
rz(-0.61856602) q[3];
sx q[3];
rz(-1.1453218) q[3];
sx q[3];
rz(-2.8989603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8689416) q[2];
sx q[2];
rz(-2.36918) q[2];
sx q[2];
rz(-0.19115494) q[2];
rz(1.8027421) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(-0.52201456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.987748) q[0];
sx q[0];
rz(-0.67180434) q[0];
sx q[0];
rz(1.9195358) q[0];
rz(2.1568495) q[1];
sx q[1];
rz(-2.1538815) q[1];
sx q[1];
rz(-0.15880671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75750763) q[0];
sx q[0];
rz(-2.473978) q[0];
sx q[0];
rz(0.13170858) q[0];
rz(-pi) q[1];
rz(-0.42694636) q[2];
sx q[2];
rz(-1.5265577) q[2];
sx q[2];
rz(2.6298863) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.040222283) q[1];
sx q[1];
rz(-1.985038) q[1];
sx q[1];
rz(2.657402) q[1];
x q[2];
rz(2.5582223) q[3];
sx q[3];
rz(-2.8900263) q[3];
sx q[3];
rz(2.5473197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1211991) q[2];
sx q[2];
rz(-2.8195916) q[2];
sx q[2];
rz(2.2484692) q[2];
rz(2.761306) q[3];
sx q[3];
rz(-1.9392574) q[3];
sx q[3];
rz(-1.7072385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0685843) q[0];
sx q[0];
rz(-2.6717654) q[0];
sx q[0];
rz(-1.6597066) q[0];
rz(2.1549554) q[1];
sx q[1];
rz(-1.4886798) q[1];
sx q[1];
rz(-0.557244) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065345848) q[0];
sx q[0];
rz(-1.6827415) q[0];
sx q[0];
rz(1.5165138) q[0];
rz(1.9300013) q[2];
sx q[2];
rz(-1.85382) q[2];
sx q[2];
rz(0.67365269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4367974) q[1];
sx q[1];
rz(-2.3782502) q[1];
sx q[1];
rz(-0.82398681) q[1];
rz(1.581963) q[3];
sx q[3];
rz(-1.5538994) q[3];
sx q[3];
rz(-0.12728413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8318994) q[2];
sx q[2];
rz(-2.5688186) q[2];
sx q[2];
rz(-1.0581623) q[2];
rz(1.7582827) q[3];
sx q[3];
rz(-2.1027095) q[3];
sx q[3];
rz(0.96341187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6365373) q[0];
sx q[0];
rz(-1.1543244) q[0];
sx q[0];
rz(-0.43168798) q[0];
rz(-0.70026669) q[1];
sx q[1];
rz(-1.9077178) q[1];
sx q[1];
rz(-0.69469992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2852531) q[0];
sx q[0];
rz(-2.4972557) q[0];
sx q[0];
rz(-0.9327234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0236597) q[2];
sx q[2];
rz(-1.6647571) q[2];
sx q[2];
rz(0.31459034) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52076653) q[1];
sx q[1];
rz(-1.0559901) q[1];
sx q[1];
rz(0.52052814) q[1];
rz(-pi) q[2];
rz(1.5273109) q[3];
sx q[3];
rz(-1.3924034) q[3];
sx q[3];
rz(-0.94387335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6557287) q[2];
sx q[2];
rz(-2.8317917) q[2];
sx q[2];
rz(-0.41717213) q[2];
rz(-0.74472767) q[3];
sx q[3];
rz(-1.3436907) q[3];
sx q[3];
rz(0.97851306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3485296) q[0];
sx q[0];
rz(-1.2611669) q[0];
sx q[0];
rz(1.4766759) q[0];
rz(-1.9229802) q[1];
sx q[1];
rz(-1.1398133) q[1];
sx q[1];
rz(-2.555991) q[1];
rz(-0.55773363) q[2];
sx q[2];
rz(-2.1699868) q[2];
sx q[2];
rz(1.5716519) q[2];
rz(-0.95016913) q[3];
sx q[3];
rz(-1.8425377) q[3];
sx q[3];
rz(0.80811926) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
