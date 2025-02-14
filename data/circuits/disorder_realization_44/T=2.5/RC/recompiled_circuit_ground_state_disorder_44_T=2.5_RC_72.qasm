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
rz(0.96762586) q[1];
sx q[1];
rz(3.1556407) q[1];
sx q[1];
rz(10.217759) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4913038) q[0];
sx q[0];
rz(-2.2772335) q[0];
sx q[0];
rz(-0.44235787) q[0];
x q[1];
rz(-1.6184801) q[2];
sx q[2];
rz(-1.5565201) q[2];
sx q[2];
rz(-0.88231495) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1268908) q[1];
sx q[1];
rz(-2.0027936) q[1];
sx q[1];
rz(0.72662455) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7865016) q[3];
sx q[3];
rz(-0.77148998) q[3];
sx q[3];
rz(2.2685693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27796081) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(2.5521736) q[2];
rz(-0.74728084) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(2.5417627) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0624307) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(-2.2404501) q[0];
rz(-0.98183739) q[1];
sx q[1];
rz(-0.58126175) q[1];
sx q[1];
rz(0.99220401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3311148) q[0];
sx q[0];
rz(-1.7077291) q[0];
sx q[0];
rz(2.1934319) q[0];
rz(0.82291863) q[2];
sx q[2];
rz(-1.6919005) q[2];
sx q[2];
rz(-1.9117282) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8507311) q[1];
sx q[1];
rz(-0.90576101) q[1];
sx q[1];
rz(2.4348201) q[1];
rz(0.76073356) q[3];
sx q[3];
rz(-2.2983716) q[3];
sx q[3];
rz(2.6618877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-0.53043008) q[2];
sx q[2];
rz(3.1165822) q[2];
rz(2.181634) q[3];
sx q[3];
rz(-0.046006087) q[3];
sx q[3];
rz(-3.0733601) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18441021) q[0];
sx q[0];
rz(-3.130641) q[0];
sx q[0];
rz(0.72407323) q[0];
rz(3.030576) q[1];
sx q[1];
rz(-0.60581517) q[1];
sx q[1];
rz(0.012880005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89319481) q[0];
sx q[0];
rz(-1.7832558) q[0];
sx q[0];
rz(-0.080470632) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7628646) q[2];
sx q[2];
rz(-1.3837866) q[2];
sx q[2];
rz(-2.929304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4531969) q[1];
sx q[1];
rz(-1.4047785) q[1];
sx q[1];
rz(0.24035458) q[1];
rz(-pi) q[2];
rz(1.1231842) q[3];
sx q[3];
rz(-0.73799282) q[3];
sx q[3];
rz(0.18692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8956902) q[2];
sx q[2];
rz(-2.4195713) q[2];
sx q[2];
rz(-0.32430172) q[2];
rz(0.4970099) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(-2.9089109) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67320353) q[0];
sx q[0];
rz(-2.9717746) q[0];
sx q[0];
rz(2.66535) q[0];
rz(0.4854804) q[1];
sx q[1];
rz(-2.6464033) q[1];
sx q[1];
rz(0.21864299) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9354585) q[0];
sx q[0];
rz(-1.9741606) q[0];
sx q[0];
rz(-2.8950188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5584458) q[2];
sx q[2];
rz(-0.84721476) q[2];
sx q[2];
rz(-2.6820094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0899892) q[1];
sx q[1];
rz(-2.0877503) q[1];
sx q[1];
rz(-2.0614784) q[1];
rz(-2.5095482) q[3];
sx q[3];
rz(-0.38292745) q[3];
sx q[3];
rz(-2.6296089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.789088) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(2.5701806) q[2];
rz(0.93829751) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(-0.17294426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57662302) q[0];
sx q[0];
rz(-0.40410703) q[0];
sx q[0];
rz(2.670766) q[0];
rz(-2.2305523) q[1];
sx q[1];
rz(-2.7016579) q[1];
sx q[1];
rz(-2.7104673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8017641) q[0];
sx q[0];
rz(-1.5976115) q[0];
sx q[0];
rz(0.78435314) q[0];
rz(-pi) q[1];
rz(-2.2041915) q[2];
sx q[2];
rz(-1.3654815) q[2];
sx q[2];
rz(-0.61287731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31897983) q[1];
sx q[1];
rz(-0.07399223) q[1];
sx q[1];
rz(1.1896108) q[1];
rz(-pi) q[2];
rz(-0.2016368) q[3];
sx q[3];
rz(-1.4343702) q[3];
sx q[3];
rz(-1.0659983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33991459) q[2];
sx q[2];
rz(-3.1368003) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(-0.46975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.85207283) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(3.0186655) q[0];
rz(0.3854824) q[1];
sx q[1];
rz(-0.79817927) q[1];
sx q[1];
rz(2.4179359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0892886) q[0];
sx q[0];
rz(-1.6160242) q[0];
sx q[0];
rz(1.7448698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0704109) q[2];
sx q[2];
rz(-2.7636409) q[2];
sx q[2];
rz(-2.9674781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9229615) q[1];
sx q[1];
rz(-2.1653163) q[1];
sx q[1];
rz(2.2943952) q[1];
rz(2.0615929) q[3];
sx q[3];
rz(-2.7088968) q[3];
sx q[3];
rz(0.015263488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75189292) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(0.67328084) q[2];
rz(0.26345396) q[3];
sx q[3];
rz(-0.47502381) q[3];
sx q[3];
rz(2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4398956) q[0];
sx q[0];
rz(-0.79653996) q[0];
sx q[0];
rz(2.6252966) q[0];
rz(2.3856178) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(0.36852536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4004423) q[0];
sx q[0];
rz(-2.0870805) q[0];
sx q[0];
rz(-2.3325066) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0958648) q[2];
sx q[2];
rz(-2.4362982) q[2];
sx q[2];
rz(-0.7208342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3893551) q[1];
sx q[1];
rz(-2.9417243) q[1];
sx q[1];
rz(-1.0081069) q[1];
rz(-0.99530275) q[3];
sx q[3];
rz(-0.80484521) q[3];
sx q[3];
rz(-0.14508776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4964909) q[2];
sx q[2];
rz(-0.056702159) q[2];
sx q[2];
rz(2.6594095) q[2];
rz(-0.21969806) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(-0.68035948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.41736233) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(-2.505488) q[0];
rz(0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(2.2669534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8281404) q[0];
sx q[0];
rz(-0.16022542) q[0];
sx q[0];
rz(-1.1988706) q[0];
x q[1];
rz(-0.73546715) q[2];
sx q[2];
rz(-0.16263419) q[2];
sx q[2];
rz(-0.59413183) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85669163) q[1];
sx q[1];
rz(-1.1686348) q[1];
sx q[1];
rz(-0.42070893) q[1];
x q[2];
rz(2.3412744) q[3];
sx q[3];
rz(-2.2860675) q[3];
sx q[3];
rz(0.12249882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8735698) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(-0.91656172) q[2];
rz(-2.364184) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(-2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022631835) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(-0.23271261) q[0];
rz(2.4932056) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(-2.5737305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647301) q[0];
sx q[0];
rz(-1.9807388) q[0];
sx q[0];
rz(-2.7565844) q[0];
x q[1];
rz(0.90654208) q[2];
sx q[2];
rz(-0.34053206) q[2];
sx q[2];
rz(0.9685002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4109013) q[1];
sx q[1];
rz(-1.7736378) q[1];
sx q[1];
rz(0.14932015) q[1];
rz(-pi) q[2];
rz(-0.81553163) q[3];
sx q[3];
rz(-0.99247265) q[3];
sx q[3];
rz(-0.92507833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(2.589321) q[2];
rz(2.8596089) q[3];
sx q[3];
rz(-2.7115287) q[3];
sx q[3];
rz(3.0388487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2118537) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(-0.1317568) q[0];
rz(2.6051104) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-2.2333455) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344438) q[0];
sx q[0];
rz(-2.4465009) q[0];
sx q[0];
rz(-1.6174843) q[0];
rz(-pi) q[1];
rz(2.0386253) q[2];
sx q[2];
rz(-1.3356294) q[2];
sx q[2];
rz(0.18593341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6301292) q[1];
sx q[1];
rz(-0.88246934) q[1];
sx q[1];
rz(-0.87911112) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1397758) q[3];
sx q[3];
rz(-1.344879) q[3];
sx q[3];
rz(2.8449013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6324255) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(-2.8636279) q[2];
rz(-0.5667423) q[3];
sx q[3];
rz(-0.20213474) q[3];
sx q[3];
rz(0.91293269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81000281) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(0.48238659) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-0.24438582) q[2];
sx q[2];
rz(-0.99437154) q[2];
sx q[2];
rz(-2.852358) q[2];
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
