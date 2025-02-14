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
rz(-2.5858606) q[0];
sx q[0];
rz(-1.2795804) q[0];
sx q[0];
rz(-2.8137141) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(-1.0110923) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7231434) q[0];
sx q[0];
rz(-1.7307502) q[0];
sx q[0];
rz(2.5470995) q[0];
rz(-pi) q[1];
rz(0.39530547) q[2];
sx q[2];
rz(-0.43839851) q[2];
sx q[2];
rz(2.663118) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7465621) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(1.0814352) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1527083) q[3];
sx q[3];
rz(-1.7888513) q[3];
sx q[3];
rz(-2.5421028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8622387) q[2];
sx q[2];
rz(-0.78247672) q[2];
sx q[2];
rz(1.5834825) q[2];
rz(-2.8067348) q[3];
sx q[3];
rz(-1.0840651) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8849477) q[0];
sx q[0];
rz(-0.56459752) q[0];
sx q[0];
rz(-2.8096492) q[0];
rz(0.360082) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(2.8538381) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.677414) q[0];
sx q[0];
rz(-0.3218284) q[0];
sx q[0];
rz(2.4610956) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2505202) q[2];
sx q[2];
rz(-0.8647635) q[2];
sx q[2];
rz(1.6635739) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8043878) q[1];
sx q[1];
rz(-1.241695) q[1];
sx q[1];
rz(1.9371402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5310181) q[3];
sx q[3];
rz(-0.27974162) q[3];
sx q[3];
rz(0.66204643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.687279) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(-1.7506556) q[2];
rz(-2.2823997) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(-0.19101645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39621064) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(-0.23042738) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(0.28883019) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5342425) q[0];
sx q[0];
rz(-1.6301883) q[0];
sx q[0];
rz(-1.1618105) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14634653) q[2];
sx q[2];
rz(-0.47292559) q[2];
sx q[2];
rz(2.2951916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7728665) q[1];
sx q[1];
rz(-1.9802367) q[1];
sx q[1];
rz(3.0837497) q[1];
rz(-pi) q[2];
rz(-0.90694859) q[3];
sx q[3];
rz(-2.8606599) q[3];
sx q[3];
rz(-2.8317766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1032224) q[2];
sx q[2];
rz(-0.73870814) q[2];
sx q[2];
rz(-0.040741097) q[2];
rz(-3.0401958) q[3];
sx q[3];
rz(-1.8056168) q[3];
sx q[3];
rz(2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7539702) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(-2.907584) q[0];
rz(2.9435844) q[1];
sx q[1];
rz(-2.104069) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7957243) q[0];
sx q[0];
rz(-2.0540753) q[0];
sx q[0];
rz(0.048075284) q[0];
x q[1];
rz(-1.4794691) q[2];
sx q[2];
rz(-1.0589561) q[2];
sx q[2];
rz(-0.78081607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23468221) q[1];
sx q[1];
rz(-2.249472) q[1];
sx q[1];
rz(-0.80295103) q[1];
rz(-pi) q[2];
rz(-1.9601106) q[3];
sx q[3];
rz(-1.6004171) q[3];
sx q[3];
rz(2.9148341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6733751) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(-0.77486983) q[2];
rz(1.8153927) q[3];
sx q[3];
rz(-1.1893136) q[3];
sx q[3];
rz(0.51138043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020017) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(-0.032489754) q[0];
rz(-1.912311) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(3.0585739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1047528) q[0];
sx q[0];
rz(-0.51114619) q[0];
sx q[0];
rz(-2.6698861) q[0];
rz(-pi) q[1];
rz(1.4013702) q[2];
sx q[2];
rz(-1.3063141) q[2];
sx q[2];
rz(1.5096017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8578945) q[1];
sx q[1];
rz(-0.91111983) q[1];
sx q[1];
rz(-2.1723351) q[1];
rz(-pi) q[2];
rz(2.2076603) q[3];
sx q[3];
rz(-1.5583002) q[3];
sx q[3];
rz(-0.10904551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3518389) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(0.16072533) q[2];
rz(-1.1927346) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(2.3737657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47650325) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(-2.8506668) q[0];
rz(-1.7581958) q[1];
sx q[1];
rz(-1.8427126) q[1];
sx q[1];
rz(-2.6417522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4668149) q[0];
sx q[0];
rz(-1.5272015) q[0];
sx q[0];
rz(1.528426) q[0];
rz(-pi) q[1];
rz(-0.96237125) q[2];
sx q[2];
rz(-2.314724) q[2];
sx q[2];
rz(2.8679071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.19491296) q[1];
sx q[1];
rz(-0.99338594) q[1];
sx q[1];
rz(-2.1140079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1736702) q[3];
sx q[3];
rz(-0.80652666) q[3];
sx q[3];
rz(-1.6953118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9342039) q[2];
sx q[2];
rz(-0.61177212) q[2];
sx q[2];
rz(-0.70154166) q[2];
rz(-0.97405854) q[3];
sx q[3];
rz(-0.8166703) q[3];
sx q[3];
rz(0.90132236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262254) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(-0.65959626) q[0];
rz(1.9024128) q[1];
sx q[1];
rz(-1.0737123) q[1];
sx q[1];
rz(-1.0741796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0461615) q[0];
sx q[0];
rz(-3.1218596) q[0];
sx q[0];
rz(2.31953) q[0];
rz(-pi) q[1];
rz(-1.1848106) q[2];
sx q[2];
rz(-1.9705271) q[2];
sx q[2];
rz(-1.0077493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85107899) q[1];
sx q[1];
rz(-2.004262) q[1];
sx q[1];
rz(-2.9032533) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3775565) q[3];
sx q[3];
rz(-1.2437399) q[3];
sx q[3];
rz(0.59248176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85019511) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(-0.58471739) q[2];
rz(2.0753453) q[3];
sx q[3];
rz(-1.2277579) q[3];
sx q[3];
rz(0.40759459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(-0.48072746) q[0];
rz(-1.9302543) q[1];
sx q[1];
rz(-1.5296661) q[1];
sx q[1];
rz(-2.5239677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4803333) q[0];
sx q[0];
rz(-2.8164589) q[0];
sx q[0];
rz(-1.5839769) q[0];
x q[1];
rz(-1.0477871) q[2];
sx q[2];
rz(-2.1126502) q[2];
sx q[2];
rz(0.43719342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94719515) q[1];
sx q[1];
rz(-2.7734967) q[1];
sx q[1];
rz(-0.93666623) q[1];
rz(-0.17980735) q[3];
sx q[3];
rz(-2.9463065) q[3];
sx q[3];
rz(2.2168468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1787313) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(-2.2507131) q[2];
rz(2.1122872) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(1.5503957) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7903098) q[0];
sx q[0];
rz(-0.93515486) q[0];
sx q[0];
rz(1.9679605) q[0];
rz(1.952518) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(0.48114166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136776) q[0];
sx q[0];
rz(-0.60773931) q[0];
sx q[0];
rz(0.39791664) q[0];
rz(-2.9815707) q[2];
sx q[2];
rz(-2.5759187) q[2];
sx q[2];
rz(-2.6765347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18080841) q[1];
sx q[1];
rz(-2.167302) q[1];
sx q[1];
rz(0.1779314) q[1];
x q[2];
rz(-0.44409306) q[3];
sx q[3];
rz(-2.0578579) q[3];
sx q[3];
rz(-0.57898486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9159307) q[2];
sx q[2];
rz(-1.91232) q[2];
sx q[2];
rz(-1.6602328) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.1888844) q[3];
sx q[3];
rz(1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34761053) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(-3.0718497) q[0];
rz(-0.28930411) q[1];
sx q[1];
rz(-1.2846839) q[1];
sx q[1];
rz(-1.0522254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69221321) q[0];
sx q[0];
rz(-1.894763) q[0];
sx q[0];
rz(2.762336) q[0];
rz(-pi) q[1];
rz(1.7830816) q[2];
sx q[2];
rz(-0.67321482) q[2];
sx q[2];
rz(2.5517983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1919162) q[1];
sx q[1];
rz(-1.1947617) q[1];
sx q[1];
rz(-0.38825775) q[1];
rz(-pi) q[2];
rz(0.42843786) q[3];
sx q[3];
rz(-1.9588545) q[3];
sx q[3];
rz(2.4840499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59794402) q[2];
sx q[2];
rz(-1.1755377) q[2];
sx q[2];
rz(-3.0808466) q[2];
rz(-2.8525823) q[3];
sx q[3];
rz(-1.776639) q[3];
sx q[3];
rz(1.8665159) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.65171) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(1.0501077) q[1];
sx q[1];
rz(-1.0625912) q[1];
sx q[1];
rz(-1.2407632) q[1];
rz(1.1012668) q[2];
sx q[2];
rz(-1.4360089) q[2];
sx q[2];
rz(0.77965005) q[2];
rz(-1.0084739) q[3];
sx q[3];
rz(-0.49839603) q[3];
sx q[3];
rz(-2.8965542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
