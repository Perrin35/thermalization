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
rz(3.0577793) q[0];
sx q[0];
rz(2.7901791) q[0];
sx q[0];
rz(8.3793381) q[0];
rz(0.84402973) q[1];
sx q[1];
rz(-2.5764155) q[1];
sx q[1];
rz(1.3619818) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365252) q[0];
sx q[0];
rz(-2.1224667) q[0];
sx q[0];
rz(1.9776634) q[0];
x q[1];
rz(0.092081618) q[2];
sx q[2];
rz(-1.144334) q[2];
sx q[2];
rz(1.2464166) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1782383) q[1];
sx q[1];
rz(-0.46666086) q[1];
sx q[1];
rz(2.6299824) q[1];
x q[2];
rz(0.33884873) q[3];
sx q[3];
rz(-2.04755) q[3];
sx q[3];
rz(-2.2173405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43209806) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(0.70322767) q[2];
rz(0.75913366) q[3];
sx q[3];
rz(-2.1742564) q[3];
sx q[3];
rz(-2.4659992) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264939) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(2.3469927) q[0];
rz(-0.82851797) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(0.4812831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7672509) q[0];
sx q[0];
rz(-1.4515948) q[0];
sx q[0];
rz(1.5425372) q[0];
x q[1];
rz(2.3929446) q[2];
sx q[2];
rz(-1.949913) q[2];
sx q[2];
rz(-2.9826255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0121721) q[1];
sx q[1];
rz(-0.24298619) q[1];
sx q[1];
rz(1.245973) q[1];
rz(-2.4886143) q[3];
sx q[3];
rz(-1.1140454) q[3];
sx q[3];
rz(1.2596318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1543697) q[2];
sx q[2];
rz(-0.96174806) q[2];
sx q[2];
rz(0.50755429) q[2];
rz(2.4842723) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(1.9067732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364612) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(0.24612799) q[0];
rz(0.71324619) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(0.86457843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893892) q[0];
sx q[0];
rz(-1.8789904) q[0];
sx q[0];
rz(-2.9257923) q[0];
rz(-pi) q[1];
rz(2.2717996) q[2];
sx q[2];
rz(-0.73852362) q[2];
sx q[2];
rz(3.125762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1684226) q[1];
sx q[1];
rz(-1.4900786) q[1];
sx q[1];
rz(-2.0328841) q[1];
rz(-pi) q[2];
rz(-0.78737463) q[3];
sx q[3];
rz(-1.4630166) q[3];
sx q[3];
rz(-0.4103578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41512179) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(-1.2064639) q[2];
rz(-3.0237899) q[3];
sx q[3];
rz(-2.466187) q[3];
sx q[3];
rz(3.0142768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29663157) q[0];
sx q[0];
rz(-1.4609818) q[0];
sx q[0];
rz(-0.46517459) q[0];
rz(-0.90452114) q[1];
sx q[1];
rz(-1.0135087) q[1];
sx q[1];
rz(1.4314338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2003882) q[0];
sx q[0];
rz(-0.53576195) q[0];
sx q[0];
rz(2.2493811) q[0];
x q[1];
rz(-1.0776005) q[2];
sx q[2];
rz(-1.3921129) q[2];
sx q[2];
rz(-2.5692232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91421724) q[1];
sx q[1];
rz(-0.7694811) q[1];
sx q[1];
rz(-2.6282267) q[1];
x q[2];
rz(2.7560985) q[3];
sx q[3];
rz(-0.76498308) q[3];
sx q[3];
rz(-1.0631595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5005834) q[2];
sx q[2];
rz(-2.9312129) q[2];
sx q[2];
rz(1.496199) q[2];
rz(3.1032041) q[3];
sx q[3];
rz(-1.8206785) q[3];
sx q[3];
rz(-2.3827609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043623) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(-1.5414365) q[0];
rz(1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(-1.5024332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1609057) q[0];
sx q[0];
rz(-2.0269927) q[0];
sx q[0];
rz(-2.6260757) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2404332) q[2];
sx q[2];
rz(-0.56686646) q[2];
sx q[2];
rz(2.4971003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0509661) q[1];
sx q[1];
rz(-0.6399261) q[1];
sx q[1];
rz(-0.76211318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5556031) q[3];
sx q[3];
rz(-3.0137339) q[3];
sx q[3];
rz(-3.1154261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1372823) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(0.38055554) q[2];
rz(-0.54008326) q[3];
sx q[3];
rz(-1.8916847) q[3];
sx q[3];
rz(1.9758457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49749097) q[0];
sx q[0];
rz(-0.044450132) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(-0.050845536) q[1];
sx q[1];
rz(-2.6652002) q[1];
sx q[1];
rz(1.4643889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8576524) q[0];
sx q[0];
rz(-0.63125718) q[0];
sx q[0];
rz(2.3842978) q[0];
x q[1];
rz(0.91422407) q[2];
sx q[2];
rz(-0.79812917) q[2];
sx q[2];
rz(1.4887336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6932043) q[1];
sx q[1];
rz(-0.12388661) q[1];
sx q[1];
rz(-2.3883567) q[1];
x q[2];
rz(2.4568632) q[3];
sx q[3];
rz(-1.4103508) q[3];
sx q[3];
rz(-1.993865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9836318) q[2];
sx q[2];
rz(-0.76405683) q[2];
sx q[2];
rz(2.2264886) q[2];
rz(0.30465952) q[3];
sx q[3];
rz(-2.4200078) q[3];
sx q[3];
rz(-1.6662395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.3604928) q[0];
sx q[0];
rz(-0.89291328) q[0];
sx q[0];
rz(-1.1156981) q[0];
rz(0.64104331) q[1];
sx q[1];
rz(-1.2572925) q[1];
sx q[1];
rz(1.8152274) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4531739) q[0];
sx q[0];
rz(-2.6368615) q[0];
sx q[0];
rz(2.8674528) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37441476) q[2];
sx q[2];
rz(-1.4157989) q[2];
sx q[2];
rz(0.58954504) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39001993) q[1];
sx q[1];
rz(-2.2604001) q[1];
sx q[1];
rz(-0.71038664) q[1];
rz(-pi) q[2];
rz(2.2828045) q[3];
sx q[3];
rz(-2.1534216) q[3];
sx q[3];
rz(1.7395541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5850087) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(1.2674241) q[2];
rz(-0.061138717) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(0.88501969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0201482) q[0];
sx q[0];
rz(-2.3763438) q[0];
sx q[0];
rz(2.6159317) q[0];
rz(1.4875745) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(-2.9590327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199969) q[0];
sx q[0];
rz(-0.47843009) q[0];
sx q[0];
rz(2.6493376) q[0];
rz(2.0261954) q[2];
sx q[2];
rz(-1.2898852) q[2];
sx q[2];
rz(-1.8473491) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8570188) q[1];
sx q[1];
rz(-2.1942687) q[1];
sx q[1];
rz(-0.2116043) q[1];
rz(-2.8556999) q[3];
sx q[3];
rz(-2.0276311) q[3];
sx q[3];
rz(-1.4057807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21060264) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(3.0328499) q[2];
rz(0.10793081) q[3];
sx q[3];
rz(-2.7863672) q[3];
sx q[3];
rz(-1.7598553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0906319) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(-2.4884124) q[0];
rz(1.2921035) q[1];
sx q[1];
rz(-0.2457681) q[1];
sx q[1];
rz(0.076315708) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6600403) q[0];
sx q[0];
rz(-1.3306434) q[0];
sx q[0];
rz(1.2368519) q[0];
rz(0.033360783) q[2];
sx q[2];
rz(-1.4028869) q[2];
sx q[2];
rz(1.0238992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2321978) q[1];
sx q[1];
rz(-1.7095209) q[1];
sx q[1];
rz(-2.3826284) q[1];
rz(1.1006484) q[3];
sx q[3];
rz(-0.68504928) q[3];
sx q[3];
rz(0.73891008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2448347) q[2];
sx q[2];
rz(-1.1847757) q[2];
sx q[2];
rz(0.99311382) q[2];
rz(-1.3541597) q[3];
sx q[3];
rz(-0.5391776) q[3];
sx q[3];
rz(-1.7206934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83574522) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(2.848023) q[0];
rz(-2.1096032) q[1];
sx q[1];
rz(-0.76621619) q[1];
sx q[1];
rz(-2.8242677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.246884) q[0];
sx q[0];
rz(-0.63011677) q[0];
sx q[0];
rz(2.9219841) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55462448) q[2];
sx q[2];
rz(-1.3828424) q[2];
sx q[2];
rz(0.30890572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.369243) q[1];
sx q[1];
rz(-1.1398106) q[1];
sx q[1];
rz(-1.1915156) q[1];
x q[2];
rz(-2.844048) q[3];
sx q[3];
rz(-1.625522) q[3];
sx q[3];
rz(-1.7136337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62852922) q[2];
sx q[2];
rz(-0.3231914) q[2];
sx q[2];
rz(-0.31488669) q[2];
rz(-0.033673938) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(1.7033738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2122129) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(0.20044151) q[1];
sx q[1];
rz(-2.0574175) q[1];
sx q[1];
rz(2.5331694) q[1];
rz(2.5949316) q[2];
sx q[2];
rz(-0.056686747) q[2];
sx q[2];
rz(-2.0591339) q[2];
rz(2.759886) q[3];
sx q[3];
rz(-2.5504938) q[3];
sx q[3];
rz(2.9989178) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
