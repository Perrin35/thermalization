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
rz(-1.9266204) q[0];
sx q[0];
rz(-1.5705234) q[0];
sx q[0];
rz(-2.7745752) q[0];
rz(3.1103599) q[1];
sx q[1];
rz(0.90489689) q[1];
sx q[1];
rz(8.7483258) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1679992) q[0];
sx q[0];
rz(-1.8449835) q[0];
sx q[0];
rz(1.6674158) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1284087) q[2];
sx q[2];
rz(-0.73650515) q[2];
sx q[2];
rz(-2.1188096) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45647168) q[1];
sx q[1];
rz(-1.1053797) q[1];
sx q[1];
rz(0.19486787) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6630456) q[3];
sx q[3];
rz(-1.4395096) q[3];
sx q[3];
rz(-1.4660975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.61554945) q[2];
sx q[2];
rz(-2.575826) q[2];
sx q[2];
rz(-0.88741285) q[2];
rz(0.081485661) q[3];
sx q[3];
rz(-1.2469651) q[3];
sx q[3];
rz(-0.87765774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3990729) q[0];
sx q[0];
rz(-2.6953473) q[0];
sx q[0];
rz(-1.2267858) q[0];
rz(-2.1555105) q[1];
sx q[1];
rz(-1.6796651) q[1];
sx q[1];
rz(2.2373534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2806633) q[0];
sx q[0];
rz(-0.21337803) q[0];
sx q[0];
rz(2.5819874) q[0];
x q[1];
rz(2.9071525) q[2];
sx q[2];
rz(-0.16142217) q[2];
sx q[2];
rz(-1.8983253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0457669) q[1];
sx q[1];
rz(-0.44039044) q[1];
sx q[1];
rz(0.38118036) q[1];
rz(-pi) q[2];
rz(-2.9504602) q[3];
sx q[3];
rz(-2.2588552) q[3];
sx q[3];
rz(0.21746527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51332372) q[2];
sx q[2];
rz(-2.8474035) q[2];
sx q[2];
rz(2.7163556) q[2];
rz(-2.6325295) q[3];
sx q[3];
rz(-2.013701) q[3];
sx q[3];
rz(1.7019255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2929448) q[0];
sx q[0];
rz(-0.078462891) q[0];
sx q[0];
rz(1.5516094) q[0];
rz(-1.4750922) q[1];
sx q[1];
rz(-1.0868797) q[1];
sx q[1];
rz(2.1045254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3939115) q[0];
sx q[0];
rz(-2.5153749) q[0];
sx q[0];
rz(-1.3514723) q[0];
rz(-pi) q[1];
rz(0.53735781) q[2];
sx q[2];
rz(-1.4721057) q[2];
sx q[2];
rz(1.7627929) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.711896) q[1];
sx q[1];
rz(-2.8968004) q[1];
sx q[1];
rz(2.1866287) q[1];
rz(-0.45944233) q[3];
sx q[3];
rz(-1.407335) q[3];
sx q[3];
rz(1.3889988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7324149) q[2];
sx q[2];
rz(-2.7517509) q[2];
sx q[2];
rz(1.2624435) q[2];
rz(-3.0564195) q[3];
sx q[3];
rz(-2.6933935) q[3];
sx q[3];
rz(-1.1220773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20268102) q[0];
sx q[0];
rz(-2.1124463) q[0];
sx q[0];
rz(-0.19126782) q[0];
rz(2.250504) q[1];
sx q[1];
rz(-2.8137408) q[1];
sx q[1];
rz(1.0506312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3252661) q[0];
sx q[0];
rz(-2.3065595) q[0];
sx q[0];
rz(0.28006552) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2810535) q[2];
sx q[2];
rz(-1.8122753) q[2];
sx q[2];
rz(0.84404404) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6022608) q[1];
sx q[1];
rz(-2.6289399) q[1];
sx q[1];
rz(-2.4668478) q[1];
x q[2];
rz(2.1414143) q[3];
sx q[3];
rz(-1.4550536) q[3];
sx q[3];
rz(-1.1463349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8656859) q[2];
sx q[2];
rz(-2.9624717) q[2];
sx q[2];
rz(0.56037819) q[2];
rz(-0.066702453) q[3];
sx q[3];
rz(-1.3542465) q[3];
sx q[3];
rz(2.1393447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10750833) q[0];
sx q[0];
rz(-1.1488687) q[0];
sx q[0];
rz(0.68503553) q[0];
rz(-0.3512474) q[1];
sx q[1];
rz(-2.5076187) q[1];
sx q[1];
rz(1.8018855) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.811243) q[0];
sx q[0];
rz(-1.0817391) q[0];
sx q[0];
rz(1.5633718) q[0];
x q[1];
rz(-0.5184104) q[2];
sx q[2];
rz(-0.74096402) q[2];
sx q[2];
rz(0.10185845) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1692002) q[1];
sx q[1];
rz(-2.9771617) q[1];
sx q[1];
rz(1.753767) q[1];
rz(2.4248554) q[3];
sx q[3];
rz(-0.61541688) q[3];
sx q[3];
rz(-0.36726609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54905218) q[2];
sx q[2];
rz(-1.668674) q[2];
sx q[2];
rz(-0.215691) q[2];
rz(0.14821626) q[3];
sx q[3];
rz(-0.18311466) q[3];
sx q[3];
rz(0.7557925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.6052979) q[0];
sx q[0];
rz(-1.4301825) q[0];
sx q[0];
rz(-2.9624665) q[0];
rz(2.1419549) q[1];
sx q[1];
rz(-2.3536286) q[1];
sx q[1];
rz(-0.16440186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59578204) q[0];
sx q[0];
rz(-0.82909697) q[0];
sx q[0];
rz(-1.8442283) q[0];
x q[1];
rz(0.51785179) q[2];
sx q[2];
rz(-1.4538063) q[2];
sx q[2];
rz(-1.3530664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9008847) q[1];
sx q[1];
rz(-1.5635033) q[1];
sx q[1];
rz(0.65967719) q[1];
rz(-pi) q[2];
rz(2.7751924) q[3];
sx q[3];
rz(-1.783964) q[3];
sx q[3];
rz(2.873796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0187692) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(1.4336047) q[2];
rz(0.99099365) q[3];
sx q[3];
rz(-1.7510479) q[3];
sx q[3];
rz(1.8979134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3864415) q[0];
sx q[0];
rz(-0.069510892) q[0];
sx q[0];
rz(-2.884927) q[0];
rz(3.0811884) q[1];
sx q[1];
rz(-0.81788617) q[1];
sx q[1];
rz(-2.7108257) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11837932) q[0];
sx q[0];
rz(-2.1617365) q[0];
sx q[0];
rz(0.33765461) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97051986) q[2];
sx q[2];
rz(-0.97504967) q[2];
sx q[2];
rz(-0.76825324) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6182729) q[1];
sx q[1];
rz(-2.225481) q[1];
sx q[1];
rz(-2.1140845) q[1];
rz(-2.2604094) q[3];
sx q[3];
rz(-1.6828711) q[3];
sx q[3];
rz(-1.4209233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7595235) q[2];
sx q[2];
rz(-0.52983317) q[2];
sx q[2];
rz(2.1257373) q[2];
rz(1.2782512) q[3];
sx q[3];
rz(-1.2334373) q[3];
sx q[3];
rz(-2.5281233) q[3];
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
rz(2.627219) q[0];
sx q[0];
rz(-0.52424279) q[0];
sx q[0];
rz(0.71994495) q[0];
rz(-0.68079692) q[1];
sx q[1];
rz(-0.38377181) q[1];
sx q[1];
rz(-2.0489571) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2075296) q[0];
sx q[0];
rz(-2.060256) q[0];
sx q[0];
rz(1.5402769) q[0];
rz(-1.2038379) q[2];
sx q[2];
rz(-0.99571823) q[2];
sx q[2];
rz(-0.42189156) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5515308) q[1];
sx q[1];
rz(-2.2023337) q[1];
sx q[1];
rz(2.2582891) q[1];
rz(-pi) q[2];
rz(-0.44050782) q[3];
sx q[3];
rz(-0.28127174) q[3];
sx q[3];
rz(-0.59556304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0138268) q[2];
sx q[2];
rz(-2.3769145) q[2];
sx q[2];
rz(-0.71316767) q[2];
rz(1.7662883) q[3];
sx q[3];
rz(-1.9153374) q[3];
sx q[3];
rz(1.3885952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16667287) q[0];
sx q[0];
rz(-2.9627934) q[0];
sx q[0];
rz(1.7145702) q[0];
rz(2.7339281) q[1];
sx q[1];
rz(-1.0382321) q[1];
sx q[1];
rz(-2.8020249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7106066) q[0];
sx q[0];
rz(-2.2453968) q[0];
sx q[0];
rz(-0.33682025) q[0];
rz(-1.6194862) q[2];
sx q[2];
rz(-1.9867989) q[2];
sx q[2];
rz(2.4456519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.133278) q[1];
sx q[1];
rz(-2.5247658) q[1];
sx q[1];
rz(0.99442579) q[1];
rz(-pi) q[2];
rz(-2.5224116) q[3];
sx q[3];
rz(-1.4805438) q[3];
sx q[3];
rz(-2.7803286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0582383) q[2];
sx q[2];
rz(-2.3149172) q[2];
sx q[2];
rz(-2.6969686) q[2];
rz(-0.63117635) q[3];
sx q[3];
rz(-1.3381713) q[3];
sx q[3];
rz(-1.9806503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6889965) q[0];
sx q[0];
rz(-1.3225553) q[0];
sx q[0];
rz(-0.026206503) q[0];
rz(1.7474489) q[1];
sx q[1];
rz(-1.1528287) q[1];
sx q[1];
rz(1.1411512) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885006) q[0];
sx q[0];
rz(-1.1682751) q[0];
sx q[0];
rz(2.5293468) q[0];
x q[1];
rz(3.0442002) q[2];
sx q[2];
rz(-0.44008128) q[2];
sx q[2];
rz(-0.68005622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.929229) q[1];
sx q[1];
rz(-2.432669) q[1];
sx q[1];
rz(2.2422748) q[1];
x q[2];
rz(-2.908292) q[3];
sx q[3];
rz(-0.58252347) q[3];
sx q[3];
rz(1.7482116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.074284) q[2];
sx q[2];
rz(-1.1003541) q[2];
sx q[2];
rz(-1.9791774) q[2];
rz(-0.31636604) q[3];
sx q[3];
rz(-1.1917944) q[3];
sx q[3];
rz(-1.4828064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3027073) q[0];
sx q[0];
rz(-0.067262983) q[0];
sx q[0];
rz(2.5433232) q[0];
rz(-0.0089946714) q[1];
sx q[1];
rz(-2.3458377) q[1];
sx q[1];
rz(-1.5942106) q[1];
rz(-0.66214421) q[2];
sx q[2];
rz(-1.6877996) q[2];
sx q[2];
rz(0.43588426) q[2];
rz(0.056445382) q[3];
sx q[3];
rz(-0.79462449) q[3];
sx q[3];
rz(0.36804646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
