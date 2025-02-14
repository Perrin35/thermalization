OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(-2.8328083) q[0];
sx q[0];
rz(0.2398332) q[0];
rz(-0.56107768) q[1];
sx q[1];
rz(-2.3993888) q[1];
sx q[1];
rz(-1.5129369) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1334907) q[0];
sx q[0];
rz(-1.2043796) q[0];
sx q[0];
rz(-0.13704637) q[0];
x q[1];
rz(-1.0885029) q[2];
sx q[2];
rz(-0.52327195) q[2];
sx q[2];
rz(-0.54973093) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0061683361) q[1];
sx q[1];
rz(-1.0841771) q[1];
sx q[1];
rz(-1.4044589) q[1];
rz(-1.739108) q[3];
sx q[3];
rz(-1.6279847) q[3];
sx q[3];
rz(0.38874278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1538887) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(2.0858916) q[2];
rz(2.740247) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(-1.6373985) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317247) q[0];
sx q[0];
rz(-0.26917502) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(2.3954605) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(3.0768118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822183) q[0];
sx q[0];
rz(-1.1752593) q[0];
sx q[0];
rz(-2.6526124) q[0];
x q[1];
rz(0.20938907) q[2];
sx q[2];
rz(-2.4509666) q[2];
sx q[2];
rz(-2.2356981) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83103115) q[1];
sx q[1];
rz(-1.4906724) q[1];
sx q[1];
rz(-1.4270272) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.719789) q[3];
sx q[3];
rz(-0.21315609) q[3];
sx q[3];
rz(-0.085863559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88227162) q[2];
sx q[2];
rz(-2.6210625) q[2];
sx q[2];
rz(3.1003013) q[2];
rz(1.1193554) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(2.0004499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(0.69931716) q[0];
rz(-2.518867) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(0.25845382) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42180443) q[0];
sx q[0];
rz(-2.5076619) q[0];
sx q[0];
rz(0.020672043) q[0];
rz(-2.223143) q[2];
sx q[2];
rz(-1.47965) q[2];
sx q[2];
rz(0.80207588) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6221348) q[1];
sx q[1];
rz(-1.5909275) q[1];
sx q[1];
rz(2.6926671) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1408125) q[3];
sx q[3];
rz(-1.6039492) q[3];
sx q[3];
rz(1.6044817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93566018) q[2];
sx q[2];
rz(-1.6814597) q[2];
sx q[2];
rz(0.70651954) q[2];
rz(-2.5126854) q[3];
sx q[3];
rz(-0.8258515) q[3];
sx q[3];
rz(-1.1227054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(-0.97417796) q[0];
rz(-1.4840508) q[1];
sx q[1];
rz(-2.0482792) q[1];
sx q[1];
rz(1.0008224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580622) q[0];
sx q[0];
rz(-0.62072004) q[0];
sx q[0];
rz(-1.3002943) q[0];
x q[1];
rz(-0.77032178) q[2];
sx q[2];
rz(-0.97409596) q[2];
sx q[2];
rz(1.53231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3127733) q[1];
sx q[1];
rz(-2.3760894) q[1];
sx q[1];
rz(-0.70154066) q[1];
rz(0.14814143) q[3];
sx q[3];
rz(-2.8175948) q[3];
sx q[3];
rz(-1.9213541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68690825) q[2];
sx q[2];
rz(-1.800622) q[2];
sx q[2];
rz(-2.3243813) q[2];
rz(2.5456083) q[3];
sx q[3];
rz(-1.417336) q[3];
sx q[3];
rz(-1.7433085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3572094) q[0];
sx q[0];
rz(-1.4056118) q[0];
sx q[0];
rz(2.0696409) q[0];
rz(-1.1373854) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(2.9683108) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3629775) q[0];
sx q[0];
rz(-2.5234875) q[0];
sx q[0];
rz(-0.72160665) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7543147) q[2];
sx q[2];
rz(-1.1114632) q[2];
sx q[2];
rz(2.0044553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2318423) q[1];
sx q[1];
rz(-2.3810003) q[1];
sx q[1];
rz(-1.4666345) q[1];
rz(2.9786879) q[3];
sx q[3];
rz(-1.4789733) q[3];
sx q[3];
rz(0.11549982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6413573) q[2];
sx q[2];
rz(-2.6884029) q[2];
sx q[2];
rz(-2.267061) q[2];
rz(0.66679653) q[3];
sx q[3];
rz(-1.6801445) q[3];
sx q[3];
rz(-0.82505208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2905529) q[0];
sx q[0];
rz(-2.9930826) q[0];
sx q[0];
rz(-0.17459757) q[0];
rz(-0.54706508) q[1];
sx q[1];
rz(-2.6262296) q[1];
sx q[1];
rz(-1.2778767) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7449943) q[0];
sx q[0];
rz(-1.6327259) q[0];
sx q[0];
rz(1.8577736) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0401947) q[2];
sx q[2];
rz(-2.2661327) q[2];
sx q[2];
rz(2.4027195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6656832) q[1];
sx q[1];
rz(-1.2194677) q[1];
sx q[1];
rz(2.372588) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2383591) q[3];
sx q[3];
rz(-2.0732911) q[3];
sx q[3];
rz(-0.63695217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3829019) q[2];
sx q[2];
rz(-0.26472696) q[2];
sx q[2];
rz(-0.36941377) q[2];
rz(0.82768011) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740042) q[0];
sx q[0];
rz(-0.91667691) q[0];
sx q[0];
rz(0.23707238) q[0];
rz(-1.7489307) q[1];
sx q[1];
rz(-1.1940414) q[1];
sx q[1];
rz(-1.0038092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2938376) q[0];
sx q[0];
rz(-1.7373573) q[0];
sx q[0];
rz(2.7012347) q[0];
x q[1];
rz(-2.4280274) q[2];
sx q[2];
rz(-1.2110365) q[2];
sx q[2];
rz(-0.52116115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29664111) q[1];
sx q[1];
rz(-1.9012723) q[1];
sx q[1];
rz(0.56092324) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4455261) q[3];
sx q[3];
rz(-2.4579774) q[3];
sx q[3];
rz(-0.044134951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5281333) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(2.6252739) q[2];
rz(2.3033219) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(-2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8539921) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(-2.2273492) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(-2.2343238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4205243) q[0];
sx q[0];
rz(-0.64890175) q[0];
sx q[0];
rz(-1.6008928) q[0];
rz(-2.5261717) q[2];
sx q[2];
rz(-1.713101) q[2];
sx q[2];
rz(-0.9535383) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0338194) q[1];
sx q[1];
rz(-0.93074742) q[1];
sx q[1];
rz(2.1442851) q[1];
x q[2];
rz(-2.8738638) q[3];
sx q[3];
rz(-0.98891034) q[3];
sx q[3];
rz(1.8393593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3153136) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(-0.86177525) q[2];
rz(1.2365384) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(1.2110565) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(1.030141) q[0];
rz(-2.8252699) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(2.8533459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41053) q[0];
sx q[0];
rz(-1.9707435) q[0];
sx q[0];
rz(1.4131141) q[0];
rz(2.4924303) q[2];
sx q[2];
rz(-1.9948975) q[2];
sx q[2];
rz(2.8555388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2001674) q[1];
sx q[1];
rz(-1.1516478) q[1];
sx q[1];
rz(2.2302385) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4699798) q[3];
sx q[3];
rz(-1.3928431) q[3];
sx q[3];
rz(-2.2554397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2890702) q[2];
sx q[2];
rz(-1.728629) q[2];
sx q[2];
rz(0.22656013) q[2];
rz(1.4551) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(-0.69211012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15014547) q[0];
sx q[0];
rz(-1.2757855) q[0];
sx q[0];
rz(-0.62514296) q[0];
rz(-1.9236247) q[1];
sx q[1];
rz(-2.4324799) q[1];
sx q[1];
rz(-1.9295173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0481794) q[0];
sx q[0];
rz(-2.163889) q[0];
sx q[0];
rz(2.1593447) q[0];
rz(-pi) q[1];
rz(1.8292959) q[2];
sx q[2];
rz(-2.2257559) q[2];
sx q[2];
rz(1.0884681) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0733302) q[1];
sx q[1];
rz(-1.1302599) q[1];
sx q[1];
rz(3.0777626) q[1];
rz(-pi) q[2];
rz(2.3732568) q[3];
sx q[3];
rz(-2.7339122) q[3];
sx q[3];
rz(2.3967495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7071699) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(2.0909069) q[2];
rz(0.55189842) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(1.2845854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(2.3552409) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(-1.9342061) q[2];
sx q[2];
rz(-1.6061693) q[2];
sx q[2];
rz(2.1314878) q[2];
rz(-2.4554853) q[3];
sx q[3];
rz(-1.5024363) q[3];
sx q[3];
rz(-0.92492044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
