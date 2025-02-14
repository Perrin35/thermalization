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
rz(-3.1336194) q[0];
sx q[0];
rz(-2.70533) q[0];
sx q[0];
rz(-2.1313957) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(7.0893256) q[1];
sx q[1];
rz(10.448699) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14903325) q[0];
sx q[0];
rz(-1.5776538) q[0];
sx q[0];
rz(1.0829145) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3918565) q[2];
sx q[2];
rz(-1.7855682) q[2];
sx q[2];
rz(-1.880876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0144452) q[1];
sx q[1];
rz(-1.5513485) q[1];
sx q[1];
rz(1.4102742) q[1];
x q[2];
rz(-1.4825304) q[3];
sx q[3];
rz(-1.3219943) q[3];
sx q[3];
rz(3.1052447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94850928) q[2];
sx q[2];
rz(-0.32863363) q[2];
sx q[2];
rz(2.9060717) q[2];
rz(2.113302) q[3];
sx q[3];
rz(-1.2582658) q[3];
sx q[3];
rz(1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4368206) q[0];
sx q[0];
rz(-1.7828159) q[0];
sx q[0];
rz(-2.219668) q[0];
rz(-0.1164662) q[1];
sx q[1];
rz(-1.7200229) q[1];
sx q[1];
rz(-2.2540653) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015892643) q[0];
sx q[0];
rz(-1.3845446) q[0];
sx q[0];
rz(3.0086579) q[0];
x q[1];
rz(1.6294125) q[2];
sx q[2];
rz(-1.2208856) q[2];
sx q[2];
rz(2.0417803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16819363) q[1];
sx q[1];
rz(-1.0032986) q[1];
sx q[1];
rz(-1.7048111) q[1];
x q[2];
rz(-1.4027897) q[3];
sx q[3];
rz(-0.6989494) q[3];
sx q[3];
rz(-2.2022171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.55887115) q[2];
sx q[2];
rz(-1.0251986) q[2];
sx q[2];
rz(-2.9026046) q[2];
rz(0.72238266) q[3];
sx q[3];
rz(-0.84010774) q[3];
sx q[3];
rz(-1.9766138) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31579414) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(-0.78829366) q[0];
rz(-1.7514508) q[1];
sx q[1];
rz(-2.7056521) q[1];
sx q[1];
rz(1.4424666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1628796) q[0];
sx q[0];
rz(-1.7121234) q[0];
sx q[0];
rz(-3.0706329) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8462538) q[2];
sx q[2];
rz(-1.2111693) q[2];
sx q[2];
rz(-2.4725898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.07787598) q[1];
sx q[1];
rz(-2.2635679) q[1];
sx q[1];
rz(-0.31097842) q[1];
x q[2];
rz(1.888186) q[3];
sx q[3];
rz(-2.6167777) q[3];
sx q[3];
rz(-1.2167236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73551377) q[2];
sx q[2];
rz(-1.9903851) q[2];
sx q[2];
rz(2.8847412) q[2];
rz(-0.86366051) q[3];
sx q[3];
rz(-2.05859) q[3];
sx q[3];
rz(-2.5679722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0470806) q[0];
sx q[0];
rz(-0.032945078) q[0];
sx q[0];
rz(1.5091913) q[0];
rz(0.31653658) q[1];
sx q[1];
rz(-1.5455952) q[1];
sx q[1];
rz(1.0155771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1986562) q[0];
sx q[0];
rz(-1.1186677) q[0];
sx q[0];
rz(-0.52152918) q[0];
x q[1];
rz(-0.79358856) q[2];
sx q[2];
rz(-1.7670791) q[2];
sx q[2];
rz(1.0283141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3104183) q[1];
sx q[1];
rz(-1.1929346) q[1];
sx q[1];
rz(2.8997594) q[1];
x q[2];
rz(-2.9731941) q[3];
sx q[3];
rz(-2.3108498) q[3];
sx q[3];
rz(-0.5357045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8429026) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(-1.451937) q[2];
rz(1.2339833) q[3];
sx q[3];
rz(-1.0943639) q[3];
sx q[3];
rz(-2.2598677) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7538309) q[0];
sx q[0];
rz(-1.3085288) q[0];
sx q[0];
rz(-2.0552788) q[0];
rz(-1.7408675) q[1];
sx q[1];
rz(-0.81175214) q[1];
sx q[1];
rz(-3.1382255) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2732179) q[0];
sx q[0];
rz(-1.5766607) q[0];
sx q[0];
rz(-0.59451367) q[0];
rz(1.0900159) q[2];
sx q[2];
rz(-0.64489105) q[2];
sx q[2];
rz(2.1457248) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7351384) q[1];
sx q[1];
rz(-2.6310096) q[1];
sx q[1];
rz(3.0264454) q[1];
rz(-pi) q[2];
rz(-3.1037381) q[3];
sx q[3];
rz(-2.541171) q[3];
sx q[3];
rz(-0.63943938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9747666) q[2];
sx q[2];
rz(-0.70009118) q[2];
sx q[2];
rz(2.8066446) q[2];
rz(-0.29609984) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(-0.12224841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.629338) q[0];
sx q[0];
rz(-2.6423995) q[0];
sx q[0];
rz(-2.7235624) q[0];
rz(0.66608518) q[1];
sx q[1];
rz(-1.2434375) q[1];
sx q[1];
rz(-2.8935249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9739415) q[0];
sx q[0];
rz(-0.60567666) q[0];
sx q[0];
rz(-0.96346897) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0866805) q[2];
sx q[2];
rz(-2.0198727) q[2];
sx q[2];
rz(-2.1841846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0937178) q[1];
sx q[1];
rz(-1.95922) q[1];
sx q[1];
rz(-2.1234496) q[1];
rz(-pi) q[2];
rz(0.25603981) q[3];
sx q[3];
rz(-2.1845792) q[3];
sx q[3];
rz(-2.7308794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4416113) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(-1.1629533) q[2];
rz(-2.5476088) q[3];
sx q[3];
rz(-1.5409527) q[3];
sx q[3];
rz(0.99015132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.537259) q[0];
sx q[0];
rz(-2.8626677) q[0];
sx q[0];
rz(0.064924031) q[0];
rz(-0.36554947) q[1];
sx q[1];
rz(-1.7100916) q[1];
sx q[1];
rz(2.4117267) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9655617) q[0];
sx q[0];
rz(-2.4906127) q[0];
sx q[0];
rz(-1.3887547) q[0];
rz(-pi) q[1];
rz(-1.7634298) q[2];
sx q[2];
rz(-2.3245077) q[2];
sx q[2];
rz(2.0157004) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.262558) q[1];
sx q[1];
rz(-1.015519) q[1];
sx q[1];
rz(0.11206514) q[1];
x q[2];
rz(-0.49298005) q[3];
sx q[3];
rz(-2.3928309) q[3];
sx q[3];
rz(1.4593339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87021747) q[2];
sx q[2];
rz(-1.368618) q[2];
sx q[2];
rz(2.0591056) q[2];
rz(-1.4879976) q[3];
sx q[3];
rz(-2.7797785) q[3];
sx q[3];
rz(-2.7910119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38700405) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(-0.3311232) q[0];
rz(-1.6638727) q[1];
sx q[1];
rz(-0.96550566) q[1];
sx q[1];
rz(2.8211735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5950111) q[0];
sx q[0];
rz(-3.0394331) q[0];
sx q[0];
rz(2.8861396) q[0];
x q[1];
rz(-1.4244979) q[2];
sx q[2];
rz(-1.4182036) q[2];
sx q[2];
rz(2.1426147) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.85727188) q[1];
sx q[1];
rz(-1.9969988) q[1];
sx q[1];
rz(-1.5537474) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5005169) q[3];
sx q[3];
rz(-0.68551862) q[3];
sx q[3];
rz(-2.7932515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42935818) q[2];
sx q[2];
rz(-1.4171436) q[2];
sx q[2];
rz(2.7799535) q[2];
rz(0.9137736) q[3];
sx q[3];
rz(-0.29699609) q[3];
sx q[3];
rz(-2.698212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906032) q[0];
sx q[0];
rz(-0.78762233) q[0];
sx q[0];
rz(-0.11601624) q[0];
rz(1.8138255) q[1];
sx q[1];
rz(-1.1632183) q[1];
sx q[1];
rz(1.2001002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3108771) q[0];
sx q[0];
rz(-0.19234622) q[0];
sx q[0];
rz(-3.0430692) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3692672) q[2];
sx q[2];
rz(-0.69567615) q[2];
sx q[2];
rz(-1.8303378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.42182339) q[1];
sx q[1];
rz(-1.0092111) q[1];
sx q[1];
rz(1.5535068) q[1];
x q[2];
rz(-1.8103792) q[3];
sx q[3];
rz(-1.7322043) q[3];
sx q[3];
rz(-0.39612285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3001331) q[2];
sx q[2];
rz(-1.7068784) q[2];
sx q[2];
rz(-0.33162281) q[2];
rz(2.8578109) q[3];
sx q[3];
rz(-1.1016568) q[3];
sx q[3];
rz(-2.7200429) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3109741) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(0.16657883) q[0];
rz(2.2919948) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(-0.073624484) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4723785) q[0];
sx q[0];
rz(-1.3166691) q[0];
sx q[0];
rz(0.73331244) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1207091) q[2];
sx q[2];
rz(-1.1556018) q[2];
sx q[2];
rz(0.5116432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8371997) q[1];
sx q[1];
rz(-1.0719258) q[1];
sx q[1];
rz(-1.9552597) q[1];
x q[2];
rz(1.0898548) q[3];
sx q[3];
rz(-1.0006128) q[3];
sx q[3];
rz(-0.4057623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9061884) q[2];
sx q[2];
rz(-0.81934682) q[2];
sx q[2];
rz(-2.4833615) q[2];
rz(1.8002347) q[3];
sx q[3];
rz(-2.1737289) q[3];
sx q[3];
rz(2.2906176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3157208) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(0.52147621) q[1];
sx q[1];
rz(-1.5745402) q[1];
sx q[1];
rz(1.5706617) q[1];
rz(-2.6084788) q[2];
sx q[2];
rz(-2.2775912) q[2];
sx q[2];
rz(-2.9698042) q[2];
rz(-2.9110649) q[3];
sx q[3];
rz(-2.5067825) q[3];
sx q[3];
rz(1.4092177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
