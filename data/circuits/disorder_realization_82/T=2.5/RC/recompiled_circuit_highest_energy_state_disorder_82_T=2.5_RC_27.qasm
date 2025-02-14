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
rz(-0.43871969) q[0];
sx q[0];
rz(3.7511711) q[0];
sx q[0];
rz(10.114976) q[0];
rz(-1.8742427) q[1];
sx q[1];
rz(-2.6522418) q[1];
sx q[1];
rz(0.34440053) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716838) q[0];
sx q[0];
rz(-0.2580041) q[0];
sx q[0];
rz(-2.1841316) q[0];
x q[1];
rz(-1.3252668) q[2];
sx q[2];
rz(-1.1628764) q[2];
sx q[2];
rz(-1.7657042) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8616383) q[1];
sx q[1];
rz(-0.32740232) q[1];
sx q[1];
rz(-0.91307098) q[1];
rz(-2.9696483) q[3];
sx q[3];
rz(-2.056582) q[3];
sx q[3];
rz(2.7051089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.043618) q[2];
sx q[2];
rz(-0.93150413) q[2];
sx q[2];
rz(1.0574868) q[2];
rz(0.99772325) q[3];
sx q[3];
rz(-1.3910553) q[3];
sx q[3];
rz(1.1725496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4443896) q[0];
sx q[0];
rz(-0.81962219) q[0];
sx q[0];
rz(2.3890553) q[0];
rz(-1.2731816) q[1];
sx q[1];
rz(-1.1538785) q[1];
sx q[1];
rz(0.85535991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0146926) q[0];
sx q[0];
rz(-1.3570667) q[0];
sx q[0];
rz(-3.041211) q[0];
x q[1];
rz(0.94903058) q[2];
sx q[2];
rz(-1.3675642) q[2];
sx q[2];
rz(1.0946314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29578698) q[1];
sx q[1];
rz(-1.9545022) q[1];
sx q[1];
rz(-1.3944666) q[1];
rz(0.19892502) q[3];
sx q[3];
rz(-1.595605) q[3];
sx q[3];
rz(2.5682784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7684043) q[2];
sx q[2];
rz(-1.6104001) q[2];
sx q[2];
rz(-1.1676403) q[2];
rz(3.043637) q[3];
sx q[3];
rz(-2.8601213) q[3];
sx q[3];
rz(1.9907985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34404594) q[0];
sx q[0];
rz(-0.56586376) q[0];
sx q[0];
rz(1.6746445) q[0];
rz(0.037874669) q[1];
sx q[1];
rz(-1.2118309) q[1];
sx q[1];
rz(2.0272592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5112202) q[0];
sx q[0];
rz(-1.1603705) q[0];
sx q[0];
rz(-3.0188796) q[0];
rz(0.75598209) q[2];
sx q[2];
rz(-0.70394433) q[2];
sx q[2];
rz(-2.9005425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.072664) q[1];
sx q[1];
rz(-2.0004099) q[1];
sx q[1];
rz(-2.0868029) q[1];
x q[2];
rz(-0.30348482) q[3];
sx q[3];
rz(-2.1305829) q[3];
sx q[3];
rz(-1.8523077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36180878) q[2];
sx q[2];
rz(-0.59541687) q[2];
sx q[2];
rz(-2.6105866) q[2];
rz(-0.94046721) q[3];
sx q[3];
rz(-2.2898424) q[3];
sx q[3];
rz(-1.4550335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.72857147) q[0];
sx q[0];
rz(-0.41550264) q[0];
sx q[0];
rz(0.79237932) q[0];
rz(3.0089695) q[1];
sx q[1];
rz(-0.68999973) q[1];
sx q[1];
rz(2.2161868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26811582) q[0];
sx q[0];
rz(-0.24304427) q[0];
sx q[0];
rz(0.8732218) q[0];
rz(-0.058387832) q[2];
sx q[2];
rz(-2.0217253) q[2];
sx q[2];
rz(2.1356867) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2823688) q[1];
sx q[1];
rz(-1.7710905) q[1];
sx q[1];
rz(-0.97334169) q[1];
rz(-1.2424281) q[3];
sx q[3];
rz(-2.066062) q[3];
sx q[3];
rz(-1.424448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7856019) q[2];
sx q[2];
rz(-1.8428558) q[2];
sx q[2];
rz(-1.402727) q[2];
rz(1.2518903) q[3];
sx q[3];
rz(-1.8391049) q[3];
sx q[3];
rz(-0.29786626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70016015) q[0];
sx q[0];
rz(-0.97235632) q[0];
sx q[0];
rz(-0.98883072) q[0];
rz(-1.5160457) q[1];
sx q[1];
rz(-1.4715618) q[1];
sx q[1];
rz(0.54738799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1158041) q[0];
sx q[0];
rz(-2.0594547) q[0];
sx q[0];
rz(-1.5343094) q[0];
x q[1];
rz(-1.4896738) q[2];
sx q[2];
rz(-1.2077959) q[2];
sx q[2];
rz(0.93349071) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0802287) q[1];
sx q[1];
rz(-1.1669817) q[1];
sx q[1];
rz(2.3263127) q[1];
rz(-2.8655254) q[3];
sx q[3];
rz(-1.4763586) q[3];
sx q[3];
rz(0.51968473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0543412) q[2];
sx q[2];
rz(-2.6332899) q[2];
sx q[2];
rz(0.15092078) q[2];
rz(-1.5058676) q[3];
sx q[3];
rz(-1.3003636) q[3];
sx q[3];
rz(1.4668363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6557789) q[0];
sx q[0];
rz(-1.1968311) q[0];
sx q[0];
rz(2.6659513) q[0];
rz(0.42426839) q[1];
sx q[1];
rz(-1.2356267) q[1];
sx q[1];
rz(-1.3729399) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17775336) q[0];
sx q[0];
rz(-0.31597695) q[0];
sx q[0];
rz(-2.55992) q[0];
rz(-pi) q[1];
rz(0.67694725) q[2];
sx q[2];
rz(-0.71093762) q[2];
sx q[2];
rz(2.9589911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9109505) q[1];
sx q[1];
rz(-1.9021209) q[1];
sx q[1];
rz(2.1403007) q[1];
rz(-pi) q[2];
rz(3.025029) q[3];
sx q[3];
rz(-2.0102083) q[3];
sx q[3];
rz(1.1573462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0515392) q[2];
sx q[2];
rz(-0.97344437) q[2];
sx q[2];
rz(0.092183979) q[2];
rz(-1.9696382) q[3];
sx q[3];
rz(-0.53297526) q[3];
sx q[3];
rz(2.2251825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5343269) q[0];
sx q[0];
rz(-2.3054275) q[0];
sx q[0];
rz(-2.109206) q[0];
rz(0.1296002) q[1];
sx q[1];
rz(-1.7584636) q[1];
sx q[1];
rz(-1.3547156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808896) q[0];
sx q[0];
rz(-0.96183813) q[0];
sx q[0];
rz(-0.97228284) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9117891) q[2];
sx q[2];
rz(-1.908506) q[2];
sx q[2];
rz(2.4203398) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2191085) q[1];
sx q[1];
rz(-2.0754083) q[1];
sx q[1];
rz(-0.64448661) q[1];
rz(-pi) q[2];
x q[2];
rz(0.078523648) q[3];
sx q[3];
rz(-2.0451945) q[3];
sx q[3];
rz(2.7767088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20595343) q[2];
sx q[2];
rz(-0.26198584) q[2];
sx q[2];
rz(-1.4865173) q[2];
rz(0.67241159) q[3];
sx q[3];
rz(-2.0629864) q[3];
sx q[3];
rz(3.0730754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2023778) q[0];
sx q[0];
rz(-0.12513932) q[0];
sx q[0];
rz(-2.8676721) q[0];
rz(-0.1637474) q[1];
sx q[1];
rz(-1.8293019) q[1];
sx q[1];
rz(2.7408677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3794601) q[0];
sx q[0];
rz(-0.5665938) q[0];
sx q[0];
rz(0.48948009) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2873307) q[2];
sx q[2];
rz(-0.41850433) q[2];
sx q[2];
rz(2.037627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0060746) q[1];
sx q[1];
rz(-1.589314) q[1];
sx q[1];
rz(0.89864267) q[1];
x q[2];
rz(-1.1464045) q[3];
sx q[3];
rz(-1.2514858) q[3];
sx q[3];
rz(-0.13892787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2694232) q[2];
sx q[2];
rz(-2.9672406) q[2];
sx q[2];
rz(1.4738458) q[2];
rz(3.040124) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(-1.7469223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87429109) q[0];
sx q[0];
rz(-0.75031459) q[0];
sx q[0];
rz(-0.0016203298) q[0];
rz(-2.5770309) q[1];
sx q[1];
rz(-1.3357342) q[1];
sx q[1];
rz(0.50216215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9410132) q[0];
sx q[0];
rz(-2.1920125) q[0];
sx q[0];
rz(-0.70141478) q[0];
x q[1];
rz(-2.5259168) q[2];
sx q[2];
rz(-0.64794316) q[2];
sx q[2];
rz(0.16251646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6356023) q[1];
sx q[1];
rz(-2.5924304) q[1];
sx q[1];
rz(-2.7628187) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21162866) q[3];
sx q[3];
rz(-2.0787079) q[3];
sx q[3];
rz(0.25903364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.042171176) q[2];
sx q[2];
rz(-1.9440938) q[2];
sx q[2];
rz(1.8966804) q[2];
rz(0.20243195) q[3];
sx q[3];
rz(-1.3916241) q[3];
sx q[3];
rz(-0.99948731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56364432) q[0];
sx q[0];
rz(-0.36643323) q[0];
sx q[0];
rz(-1.4165437) q[0];
rz(-1.1997403) q[1];
sx q[1];
rz(-1.1734633) q[1];
sx q[1];
rz(-2.8660668) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4352132) q[0];
sx q[0];
rz(-1.2093028) q[0];
sx q[0];
rz(2.778591) q[0];
rz(-0.24015719) q[2];
sx q[2];
rz(-0.36061812) q[2];
sx q[2];
rz(0.92133969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78026315) q[1];
sx q[1];
rz(-2.0320914) q[1];
sx q[1];
rz(0.023246846) q[1];
rz(-1.6943135) q[3];
sx q[3];
rz(-0.97573167) q[3];
sx q[3];
rz(1.7473011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3755017) q[2];
sx q[2];
rz(-1.4644863) q[2];
sx q[2];
rz(-0.38140934) q[2];
rz(2.8219847) q[3];
sx q[3];
rz(-0.8173129) q[3];
sx q[3];
rz(-2.4735425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1565336) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(1.001724) q[1];
sx q[1];
rz(-1.4717419) q[1];
sx q[1];
rz(-0.91632661) q[1];
rz(-2.9849595) q[2];
sx q[2];
rz(-2.4703783) q[2];
sx q[2];
rz(-3.0117161) q[2];
rz(-2.9974964) q[3];
sx q[3];
rz(-2.477705) q[3];
sx q[3];
rz(-2.4841819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
