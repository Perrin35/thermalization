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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(2.8454236) q[1];
sx q[1];
rz(12.256395) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5910019) q[0];
sx q[0];
rz(-2.7059116) q[0];
sx q[0];
rz(-2.8753619) q[0];
rz(-3.0344738) q[2];
sx q[2];
rz(-2.9351882) q[2];
sx q[2];
rz(-2.602488) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9604608) q[1];
sx q[1];
rz(-0.38552654) q[1];
sx q[1];
rz(1.0268282) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83155379) q[3];
sx q[3];
rz(-1.420701) q[3];
sx q[3];
rz(-2.7287366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4787204) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(0.09566801) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(1.4314502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-1.2410759) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(-0.61479968) q[0];
rz(2.2531033) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(-0.49957553) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.963012) q[0];
sx q[0];
rz(-2.7546429) q[0];
sx q[0];
rz(2.1301756) q[0];
x q[1];
rz(-0.52421661) q[2];
sx q[2];
rz(-1.7078064) q[2];
sx q[2];
rz(1.7643339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.412498) q[1];
sx q[1];
rz(-2.4654268) q[1];
sx q[1];
rz(1.2926284) q[1];
rz(-2.0604443) q[3];
sx q[3];
rz(-0.64125618) q[3];
sx q[3];
rz(-1.412751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8932314) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(-0.020817967) q[2];
rz(1.2402395) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(-1.9564691) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5027387) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(0.33367208) q[0];
rz(0.73792136) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(0.12942448) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0282818) q[0];
sx q[0];
rz(-0.89963642) q[0];
sx q[0];
rz(-2.9379803) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(-2.3455623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8379535) q[1];
sx q[1];
rz(-0.15175125) q[1];
sx q[1];
rz(-2.2084153) q[1];
x q[2];
rz(2.0468726) q[3];
sx q[3];
rz(-2.2942846) q[3];
sx q[3];
rz(-1.0915966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0176598) q[2];
sx q[2];
rz(-1.5256226) q[2];
sx q[2];
rz(1.7714436) q[2];
rz(-0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(-0.082775041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.133404) q[0];
sx q[0];
rz(-1.9318102) q[0];
sx q[0];
rz(0.56513894) q[0];
rz(-2.7124229) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(-2.3133004) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.013731) q[0];
sx q[0];
rz(-2.3963701) q[0];
sx q[0];
rz(3.0237314) q[0];
rz(2.8548334) q[2];
sx q[2];
rz(-2.0398524) q[2];
sx q[2];
rz(-3.0931851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.284621) q[1];
sx q[1];
rz(-1.0880252) q[1];
sx q[1];
rz(1.2699782) q[1];
x q[2];
rz(-2.8935562) q[3];
sx q[3];
rz(-2.0741077) q[3];
sx q[3];
rz(-2.4865884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68507489) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(-1.431541) q[2];
rz(-2.2020014) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(3.140894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058896) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(-1.3294719) q[0];
rz(1.9895408) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(-0.00024814127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084293289) q[0];
sx q[0];
rz(-1.7547742) q[0];
sx q[0];
rz(1.8198245) q[0];
rz(-pi) q[1];
rz(1.0365965) q[2];
sx q[2];
rz(-1.4956258) q[2];
sx q[2];
rz(2.6435564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2622288) q[1];
sx q[1];
rz(-1.2254224) q[1];
sx q[1];
rz(-2.097258) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73152929) q[3];
sx q[3];
rz(-1.6247182) q[3];
sx q[3];
rz(-1.1889072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0033215) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(3.1357583) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(2.0148923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936845) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(-1.6538612) q[0];
rz(-3.1112025) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(2.1991275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93118942) q[0];
sx q[0];
rz(-1.382834) q[0];
sx q[0];
rz(0.24954777) q[0];
rz(-pi) q[1];
rz(-2.4320388) q[2];
sx q[2];
rz(-0.38564577) q[2];
sx q[2];
rz(2.8288159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2561431) q[1];
sx q[1];
rz(-1.1032618) q[1];
sx q[1];
rz(2.8199151) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0464155) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(1.3568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5421062) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(-1.8360809) q[2];
rz(-2.5944338) q[3];
sx q[3];
rz(-1.2055509) q[3];
sx q[3];
rz(-0.80593306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18043537) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(-0.85711342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094497546) q[0];
sx q[0];
rz(-2.408354) q[0];
sx q[0];
rz(-0.49219699) q[0];
rz(-3.0761511) q[2];
sx q[2];
rz(-0.64787728) q[2];
sx q[2];
rz(0.53109785) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3639314) q[1];
sx q[1];
rz(-2.6865733) q[1];
sx q[1];
rz(2.6520686) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1825652) q[3];
sx q[3];
rz(-0.60743466) q[3];
sx q[3];
rz(-1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40019217) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(-2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(1.2790595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.361146) q[0];
sx q[0];
rz(-3.0028711) q[0];
sx q[0];
rz(-3.1106023) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(-1.8642289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5073858) q[0];
sx q[0];
rz(-1.9045826) q[0];
sx q[0];
rz(0.35315634) q[0];
rz(-1.1436815) q[2];
sx q[2];
rz(-2.6027205) q[2];
sx q[2];
rz(2.8123735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9121418) q[1];
sx q[1];
rz(-1.0849285) q[1];
sx q[1];
rz(0.082991675) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2856917) q[3];
sx q[3];
rz(-1.5114956) q[3];
sx q[3];
rz(1.3467195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.10451) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(2.6541397) q[2];
rz(-2.4449352) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-0.35097861) q[0];
rz(-2.9674496) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(1.7399656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49920666) q[0];
sx q[0];
rz(-2.8642352) q[0];
sx q[0];
rz(2.6794898) q[0];
rz(2.6675111) q[2];
sx q[2];
rz(-2.3326121) q[2];
sx q[2];
rz(2.7681729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.735504) q[1];
sx q[1];
rz(-1.5417409) q[1];
sx q[1];
rz(1.7543704) q[1];
rz(-pi) q[2];
rz(2.5085658) q[3];
sx q[3];
rz(-2.2928975) q[3];
sx q[3];
rz(-0.71114572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(2.5049211) q[2];
rz(2.0311671) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(-2.6231664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7745895) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(-0.038657945) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(1.0640594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23180873) q[0];
sx q[0];
rz(-1.5752821) q[0];
sx q[0];
rz(-1.5741328) q[0];
rz(-pi) q[1];
rz(0.27711192) q[2];
sx q[2];
rz(-1.849035) q[2];
sx q[2];
rz(-2.3269175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0088996) q[1];
sx q[1];
rz(-1.510396) q[1];
sx q[1];
rz(1.6113043) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6056772) q[3];
sx q[3];
rz(-2.0300755) q[3];
sx q[3];
rz(1.9495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(0.074020298) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(-2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1114125) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(0.61182712) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(1.189497) q[2];
sx q[2];
rz(-1.0972037) q[2];
sx q[2];
rz(0.9158132) q[2];
rz(1.0985804) q[3];
sx q[3];
rz(-0.87971148) q[3];
sx q[3];
rz(1.9426027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
