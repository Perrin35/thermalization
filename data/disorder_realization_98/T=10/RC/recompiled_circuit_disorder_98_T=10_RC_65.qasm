OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(-0.9057522) q[0];
x q[1];
rz(2.3157273) q[2];
sx q[2];
rz(-2.4561433) q[2];
sx q[2];
rz(0.65537383) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61923164) q[1];
sx q[1];
rz(-1.3904966) q[1];
sx q[1];
rz(-0.038832263) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2655067) q[3];
sx q[3];
rz(-2.5698834) q[3];
sx q[3];
rz(-1.7892464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95603847) q[0];
sx q[0];
rz(-1.599405) q[0];
sx q[0];
rz(-1.5831328) q[0];
rz(-0.89276887) q[2];
sx q[2];
rz(-0.36913482) q[2];
sx q[2];
rz(-2.6004651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7420885) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(-0.29466596) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2964301) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(-2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0911672) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(-1.3495548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916546) q[0];
sx q[0];
rz(-2.2313801) q[0];
sx q[0];
rz(-1.8750989) q[0];
x q[1];
rz(-1.6367958) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(-0.50022349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.459356) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(2.6911246) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1317741) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(-0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(-1.7416471) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(3.1052123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.13658) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(1.8965917) q[2];
sx q[2];
rz(-2.1278283) q[2];
sx q[2];
rz(1.4404802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.034678) q[1];
sx q[1];
rz(-2.2639096) q[1];
sx q[1];
rz(0.17130674) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4724457) q[3];
sx q[3];
rz(-1.8042759) q[3];
sx q[3];
rz(2.1360306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(2.5454583) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(-2.9751076) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(0.23434815) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4011824) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(0.58437225) q[0];
rz(1.7469823) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(1.2622152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.036451642) q[1];
sx q[1];
rz(-0.76117939) q[1];
sx q[1];
rz(2.9245604) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3393199) q[3];
sx q[3];
rz(-2.8095062) q[3];
sx q[3];
rz(-2.092923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.1436499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1724388) q[0];
sx q[0];
rz(-1.6108496) q[0];
sx q[0];
rz(0.37102951) q[0];
x q[1];
rz(1.8311062) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(1.7973763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5864582) q[1];
sx q[1];
rz(-2.0124334) q[1];
sx q[1];
rz(-0.099979062) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46623047) q[3];
sx q[3];
rz(-2.5737408) q[3];
sx q[3];
rz(-0.0044435244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75366655) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(-0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(-1.2333262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9702643) q[0];
sx q[0];
rz(-1.5597222) q[0];
sx q[0];
rz(1.1867255) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2592505) q[2];
sx q[2];
rz(-2.2922278) q[2];
sx q[2];
rz(1.3426069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8334956) q[1];
sx q[1];
rz(-1.4596246) q[1];
sx q[1];
rz(-1.5566467) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60621467) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.258237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90342605) q[0];
sx q[0];
rz(-1.705372) q[0];
sx q[0];
rz(-2.0057136) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4969205) q[2];
sx q[2];
rz(-2.1053227) q[2];
sx q[2];
rz(2.1246186) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.313169) q[1];
sx q[1];
rz(-1.0272044) q[1];
sx q[1];
rz(-1.4491175) q[1];
x q[2];
rz(-0.85556742) q[3];
sx q[3];
rz(-1.9868402) q[3];
sx q[3];
rz(2.1895529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(0.02773157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6513034) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(-0.62255967) q[0];
rz(-0.8717732) q[2];
sx q[2];
rz(-1.36424) q[2];
sx q[2];
rz(-0.86824647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2003277) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(1.4501249) q[1];
rz(0.025891993) q[3];
sx q[3];
rz(-1.6373487) q[3];
sx q[3];
rz(0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.998418) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(2.6877158) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2430902) q[0];
sx q[0];
rz(-1.9542964) q[0];
sx q[0];
rz(-0.48454185) q[0];
x q[1];
rz(1.0742513) q[2];
sx q[2];
rz(-1.0814582) q[2];
sx q[2];
rz(-0.41000965) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.232302) q[1];
sx q[1];
rz(-1.5564939) q[1];
sx q[1];
rz(-0.69625744) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24210838) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(-2.9374591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(0.47484067) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(0.67129927) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(2.0541035) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
