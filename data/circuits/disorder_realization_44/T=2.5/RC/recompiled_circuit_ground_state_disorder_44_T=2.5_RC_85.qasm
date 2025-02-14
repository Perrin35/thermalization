OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39785102) q[0];
sx q[0];
rz(-2.1704817) q[0];
sx q[0];
rz(0.71075332) q[0];
rz(-2.1739668) q[1];
sx q[1];
rz(-0.014048014) q[1];
sx q[1];
rz(2.3486121) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4913038) q[0];
sx q[0];
rz(-0.86435917) q[0];
sx q[0];
rz(2.6992348) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2797727) q[2];
sx q[2];
rz(-3.0918192) q[2];
sx q[2];
rz(-0.39779824) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1268908) q[1];
sx q[1];
rz(-1.1387991) q[1];
sx q[1];
rz(2.4149681) q[1];
x q[2];
rz(2.3550911) q[3];
sx q[3];
rz(-0.77148998) q[3];
sx q[3];
rz(2.2685693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.27796081) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(-2.5521736) q[2];
rz(-0.74728084) q[3];
sx q[3];
rz(-0.092844754) q[3];
sx q[3];
rz(-2.5417627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.079161949) q[0];
sx q[0];
rz(-0.23167647) q[0];
sx q[0];
rz(-0.9011426) q[0];
rz(-0.98183739) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(2.1493886) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81047786) q[0];
sx q[0];
rz(-1.7077291) q[0];
sx q[0];
rz(2.1934319) q[0];
rz(-pi) q[1];
rz(-1.3937226) q[2];
sx q[2];
rz(-2.3858527) q[2];
sx q[2];
rz(-2.9301639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7957669) q[1];
sx q[1];
rz(-0.92928934) q[1];
sx q[1];
rz(-0.87916763) q[1];
rz(2.2295932) q[3];
sx q[3];
rz(-0.9992632) q[3];
sx q[3];
rz(0.4812355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.046752669) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(3.1165822) q[2];
rz(-0.95995861) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(-0.068232603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(2.4175194) q[0];
rz(3.030576) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(3.1287126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2483978) q[0];
sx q[0];
rz(-1.7832558) q[0];
sx q[0];
rz(-0.080470632) q[0];
rz(-2.3518042) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(-2.546026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4531969) q[1];
sx q[1];
rz(-1.4047785) q[1];
sx q[1];
rz(2.9012381) q[1];
rz(-pi) q[2];
rz(0.37498388) q[3];
sx q[3];
rz(-0.91911784) q[3];
sx q[3];
rz(0.76256547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24590242) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(2.8172909) q[2];
rz(0.4970099) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(-2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683891) q[0];
sx q[0];
rz(-0.16981801) q[0];
sx q[0];
rz(-2.66535) q[0];
rz(-0.4854804) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(0.21864299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9354585) q[0];
sx q[0];
rz(-1.1674321) q[0];
sx q[0];
rz(-2.8950188) q[0];
rz(-2.4179732) q[2];
sx q[2];
rz(-1.5800522) q[2];
sx q[2];
rz(1.1030359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8805595) q[1];
sx q[1];
rz(-1.1487242) q[1];
sx q[1];
rz(2.5690394) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8273638) q[3];
sx q[3];
rz(-1.7933729) q[3];
sx q[3];
rz(-1.4861388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.3525047) q[2];
sx q[2];
rz(-2.3829491) q[2];
sx q[2];
rz(2.5701806) q[2];
rz(2.2032951) q[3];
sx q[3];
rz(-2.5550227) q[3];
sx q[3];
rz(-0.17294426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5649696) q[0];
sx q[0];
rz(-0.40410703) q[0];
sx q[0];
rz(0.47082666) q[0];
rz(0.91104031) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(-0.43112531) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3456988) q[0];
sx q[0];
rz(-0.78471334) q[0];
sx q[0];
rz(-0.037952947) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2324824) q[2];
sx q[2];
rz(-2.4801284) q[2];
sx q[2];
rz(-2.4544883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4404802) q[1];
sx q[1];
rz(-1.639469) q[1];
sx q[1];
rz(-0.027570034) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60089941) q[3];
sx q[3];
rz(-0.24293262) q[3];
sx q[3];
rz(-0.082279131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33991459) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-0.63569331) q[3];
sx q[3];
rz(0.46975964) q[3];
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
rz(-0.85207283) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(-0.12292718) q[0];
rz(0.3854824) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(0.72365671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.052304) q[0];
sx q[0];
rz(-1.5255685) q[0];
sx q[0];
rz(1.7448698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0704109) q[2];
sx q[2];
rz(-0.37795174) q[2];
sx q[2];
rz(-0.17411451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8116709) q[1];
sx q[1];
rz(-0.99024665) q[1];
sx q[1];
rz(-2.4075721) q[1];
rz(2.0615929) q[3];
sx q[3];
rz(-2.7088968) q[3];
sx q[3];
rz(-3.1263292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75189292) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(-0.67328084) q[2];
rz(-2.8781387) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(-2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.70169705) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(2.6252966) q[0];
rz(-2.3856178) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(-0.36852536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719292) q[0];
sx q[0];
rz(-2.2145529) q[0];
sx q[0];
rz(-0.66515775) q[0];
rz(-pi) q[1];
rz(2.4368144) q[2];
sx q[2];
rz(-1.5411589) q[2];
sx q[2];
rz(0.88479048) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8766899) q[1];
sx q[1];
rz(-1.4646834) q[1];
sx q[1];
rz(1.4011032) q[1];
rz(-pi) q[2];
rz(-2.2880408) q[3];
sx q[3];
rz(-1.9738758) q[3];
sx q[3];
rz(-2.1385101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4964909) q[2];
sx q[2];
rz(-0.056702159) q[2];
sx q[2];
rz(2.6594095) q[2];
rz(-0.21969806) q[3];
sx q[3];
rz(-0.69742656) q[3];
sx q[3];
rz(0.68035948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41736233) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(2.505488) q[0];
rz(-0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(-2.2669534) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8281404) q[0];
sx q[0];
rz(-2.9813672) q[0];
sx q[0];
rz(-1.942722) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6804439) q[2];
sx q[2];
rz(-1.6911518) q[2];
sx q[2];
rz(0.14794825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2541209) q[1];
sx q[1];
rz(-1.9560818) q[1];
sx q[1];
rz(-1.1347358) q[1];
rz(-pi) q[2];
rz(-2.465807) q[3];
sx q[3];
rz(-0.99832557) q[3];
sx q[3];
rz(1.0990717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26802289) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(-2.2250309) q[2];
rz(-2.364184) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(-2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1189608) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(2.90888) q[0];
rz(2.4932056) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(0.56786215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647301) q[0];
sx q[0];
rz(-1.1608539) q[0];
sx q[0];
rz(-0.38500824) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90654208) q[2];
sx q[2];
rz(-0.34053206) q[2];
sx q[2];
rz(-0.9685002) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9514002) q[1];
sx q[1];
rz(-1.4245598) q[1];
sx q[1];
rz(1.3657354) q[1];
rz(-pi) q[2];
rz(0.73091032) q[3];
sx q[3];
rz(-0.95940351) q[3];
sx q[3];
rz(-0.17043336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(-0.5522716) q[2];
rz(-2.8596089) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(3.0388487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.929739) q[0];
sx q[0];
rz(-0.064082853) q[0];
sx q[0];
rz(3.0098359) q[0];
rz(2.6051104) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(0.90824711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2836766) q[0];
sx q[0];
rz(-0.87661298) q[0];
sx q[0];
rz(3.1026918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0591878) q[2];
sx q[2];
rz(-0.51967144) q[2];
sx q[2];
rz(-2.1888417) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7365807) q[1];
sx q[1];
rz(-0.93376505) q[1];
sx q[1];
rz(-0.65959658) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.578701) q[3];
sx q[3];
rz(-0.22592446) q[3];
sx q[3];
rz(0.28858063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50916719) q[2];
sx q[2];
rz(-2.2278892) q[2];
sx q[2];
rz(-0.2779648) q[2];
rz(-2.5748504) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(-2.22866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.3315898) q[0];
sx q[0];
rz(-1.7362052) q[0];
sx q[0];
rz(2.195634) q[0];
rz(2.6592061) q[1];
sx q[1];
rz(-1.9079897) q[1];
sx q[1];
rz(2.447396) q[1];
rz(0.98050465) q[2];
sx q[2];
rz(-1.7750778) q[2];
sx q[2];
rz(1.9951174) q[2];
rz(0.031470555) q[3];
sx q[3];
rz(-2.0768055) q[3];
sx q[3];
rz(0.16184645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
