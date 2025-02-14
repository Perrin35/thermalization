OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7118536) q[0];
sx q[0];
rz(-1.1450333) q[0];
sx q[0];
rz(1.0166919) q[0];
rz(2.297205) q[1];
sx q[1];
rz(-2.3250186) q[1];
sx q[1];
rz(-0.56632298) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071958) q[0];
sx q[0];
rz(-1.473794) q[0];
sx q[0];
rz(0.96235458) q[0];
rz(-pi) q[1];
rz(1.7687321) q[2];
sx q[2];
rz(-1.6467484) q[2];
sx q[2];
rz(1.4966485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8493718) q[1];
sx q[1];
rz(-1.9927036) q[1];
sx q[1];
rz(1.4825691) q[1];
rz(-pi) q[2];
rz(-0.58523607) q[3];
sx q[3];
rz(-1.397831) q[3];
sx q[3];
rz(1.1542785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5113968) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(1.1279747) q[2];
rz(1.6389716) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(2.716841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2750435) q[0];
sx q[0];
rz(-0.4568704) q[0];
sx q[0];
rz(-3.1006815) q[0];
rz(2.5919137) q[1];
sx q[1];
rz(-2.7626541) q[1];
sx q[1];
rz(3.0612225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5588829) q[0];
sx q[0];
rz(-2.2930995) q[0];
sx q[0];
rz(-1.5651907) q[0];
rz(2.2533312) q[2];
sx q[2];
rz(-2.3739034) q[2];
sx q[2];
rz(1.7861799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.064172) q[1];
sx q[1];
rz(-1.5158093) q[1];
sx q[1];
rz(2.0559146) q[1];
x q[2];
rz(-2.8714058) q[3];
sx q[3];
rz(-2.1572497) q[3];
sx q[3];
rz(1.0932066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4139159) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(0.28398871) q[2];
rz(-1.6207638) q[3];
sx q[3];
rz(-1.5512543) q[3];
sx q[3];
rz(2.3782597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336998) q[0];
sx q[0];
rz(-2.9871873) q[0];
sx q[0];
rz(-0.41686091) q[0];
rz(2.2082224) q[1];
sx q[1];
rz(-0.88637543) q[1];
sx q[1];
rz(1.0881759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8899208) q[0];
sx q[0];
rz(-1.8150738) q[0];
sx q[0];
rz(-1.0142465) q[0];
rz(-pi) q[1];
rz(1.5303262) q[2];
sx q[2];
rz(-1.010561) q[2];
sx q[2];
rz(-0.76228415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58055604) q[1];
sx q[1];
rz(-1.2966178) q[1];
sx q[1];
rz(-2.8304173) q[1];
x q[2];
rz(1.635664) q[3];
sx q[3];
rz(-1.2838351) q[3];
sx q[3];
rz(-0.6395517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1209968) q[2];
sx q[2];
rz(-2.4211297) q[2];
sx q[2];
rz(-0.32568112) q[2];
rz(-2.7459512) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(2.4237848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2242103) q[0];
sx q[0];
rz(-2.2069187) q[0];
sx q[0];
rz(3.000946) q[0];
rz(2.7463101) q[1];
sx q[1];
rz(-1.5970767) q[1];
sx q[1];
rz(-0.43221727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062586322) q[0];
sx q[0];
rz(-0.64029988) q[0];
sx q[0];
rz(2.2878134) q[0];
rz(-pi) q[1];
rz(3.0555058) q[2];
sx q[2];
rz(-2.0755092) q[2];
sx q[2];
rz(1.9534257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2905457) q[1];
sx q[1];
rz(-1.6509027) q[1];
sx q[1];
rz(2.8677031) q[1];
x q[2];
rz(2.3234576) q[3];
sx q[3];
rz(-0.77518089) q[3];
sx q[3];
rz(1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29493368) q[2];
sx q[2];
rz(-1.1120956) q[2];
sx q[2];
rz(-2.6726932) q[2];
rz(2.9851959) q[3];
sx q[3];
rz(-0.96894914) q[3];
sx q[3];
rz(-1.2676574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0676607) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(-2.3611948) q[0];
rz(-2.228915) q[1];
sx q[1];
rz(-0.78881216) q[1];
sx q[1];
rz(-1.3891634) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71043832) q[0];
sx q[0];
rz(-1.9214993) q[0];
sx q[0];
rz(1.2855269) q[0];
x q[1];
rz(-0.19540968) q[2];
sx q[2];
rz(-2.5511046) q[2];
sx q[2];
rz(0.54412847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.011574419) q[1];
sx q[1];
rz(-0.60638529) q[1];
sx q[1];
rz(0.19488402) q[1];
rz(-pi) q[2];
rz(-1.0205998) q[3];
sx q[3];
rz(-0.76348272) q[3];
sx q[3];
rz(-0.17624763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3461561) q[2];
sx q[2];
rz(-0.6554335) q[2];
sx q[2];
rz(2.9620841) q[2];
rz(-1.6620212) q[3];
sx q[3];
rz(-2.447465) q[3];
sx q[3];
rz(-3.1365385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96106225) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(1.6777212) q[0];
rz(3.1278817) q[1];
sx q[1];
rz(-1.6746215) q[1];
sx q[1];
rz(-0.94508583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4714519) q[0];
sx q[0];
rz(-1.9197122) q[0];
sx q[0];
rz(0.80100153) q[0];
rz(-0.31107549) q[2];
sx q[2];
rz(-1.7698423) q[2];
sx q[2];
rz(-2.8761326) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3774776) q[1];
sx q[1];
rz(-2.4903653) q[1];
sx q[1];
rz(0.73173827) q[1];
rz(-2.6078647) q[3];
sx q[3];
rz(-2.0954359) q[3];
sx q[3];
rz(0.49688767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1549687) q[2];
sx q[2];
rz(-2.0782317) q[2];
sx q[2];
rz(2.092579) q[2];
rz(-1.5341885) q[3];
sx q[3];
rz(-2.1395855) q[3];
sx q[3];
rz(-1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55117115) q[0];
sx q[0];
rz(-2.6541002) q[0];
sx q[0];
rz(-2.3792939) q[0];
rz(-0.47362348) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(0.90829888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0001091) q[0];
sx q[0];
rz(-2.3628652) q[0];
sx q[0];
rz(1.5634006) q[0];
x q[1];
rz(1.3495096) q[2];
sx q[2];
rz(-1.7402116) q[2];
sx q[2];
rz(2.7986023) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6544853) q[1];
sx q[1];
rz(-1.5853549) q[1];
sx q[1];
rz(1.3338939) q[1];
rz(-1.5987464) q[3];
sx q[3];
rz(-1.5582435) q[3];
sx q[3];
rz(2.8421228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0198387) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(-2.1799977) q[2];
rz(-2.614295) q[3];
sx q[3];
rz(-1.3798102) q[3];
sx q[3];
rz(2.5808047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2335662) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(-2.7954234) q[0];
rz(1.8442122) q[1];
sx q[1];
rz(-0.66645122) q[1];
sx q[1];
rz(0.60874879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93701836) q[0];
sx q[0];
rz(-1.4868344) q[0];
sx q[0];
rz(0.89421009) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7815963) q[2];
sx q[2];
rz(-1.8328875) q[2];
sx q[2];
rz(2.0767625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37571463) q[1];
sx q[1];
rz(-0.27395136) q[1];
sx q[1];
rz(0.25744827) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5777588) q[3];
sx q[3];
rz(-0.878351) q[3];
sx q[3];
rz(0.13217029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79494563) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-2.4897599) q[2];
rz(2.4018304) q[3];
sx q[3];
rz(-1.0659822) q[3];
sx q[3];
rz(-2.2447926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5935434) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-0.59762534) q[0];
rz(-1.301544) q[1];
sx q[1];
rz(-1.5664682) q[1];
sx q[1];
rz(0.32599932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6582002) q[0];
sx q[0];
rz(-2.3478386) q[0];
sx q[0];
rz(-1.9824337) q[0];
x q[1];
rz(-2.8140961) q[2];
sx q[2];
rz(-1.3262889) q[2];
sx q[2];
rz(1.2933019) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44931901) q[1];
sx q[1];
rz(-0.92781583) q[1];
sx q[1];
rz(-2.8161956) q[1];
x q[2];
rz(-0.80504412) q[3];
sx q[3];
rz(-0.16318233) q[3];
sx q[3];
rz(1.2393348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43805435) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(-0.34029141) q[2];
rz(1.5002286) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(-2.5435508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.930645) q[0];
sx q[0];
rz(-2.0238545) q[0];
sx q[0];
rz(-0.043638226) q[0];
rz(0.71815193) q[1];
sx q[1];
rz(-1.8225881) q[1];
sx q[1];
rz(-1.4290379) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315009) q[0];
sx q[0];
rz(-2.8866308) q[0];
sx q[0];
rz(2.405317) q[0];
rz(2.9731644) q[2];
sx q[2];
rz(-0.68466869) q[2];
sx q[2];
rz(-1.9569966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70245269) q[1];
sx q[1];
rz(-1.0826775) q[1];
sx q[1];
rz(-1.6678651) q[1];
rz(-pi) q[2];
rz(0.36504443) q[3];
sx q[3];
rz(-1.924558) q[3];
sx q[3];
rz(2.7937074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4622197) q[2];
sx q[2];
rz(-2.2023109) q[2];
sx q[2];
rz(-2.3636554) q[2];
rz(-1.0434307) q[3];
sx q[3];
rz(-1.611004) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935681) q[0];
sx q[0];
rz(-2.1401736) q[0];
sx q[0];
rz(-2.370136) q[0];
rz(1.363516) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(-2.099299) q[2];
sx q[2];
rz(-2.2215371) q[2];
sx q[2];
rz(2.110832) q[2];
rz(-2.2424768) q[3];
sx q[3];
rz(-1.1976783) q[3];
sx q[3];
rz(-1.4730361) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
