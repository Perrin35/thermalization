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
rz(2.9444115) q[0];
sx q[0];
rz(-1.0721167) q[0];
sx q[0];
rz(1.5443235) q[0];
rz(0.3577258) q[1];
sx q[1];
rz(4.975714) q[1];
sx q[1];
rz(9.8006021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.887325) q[0];
sx q[0];
rz(-1.642881) q[0];
sx q[0];
rz(0.76728009) q[0];
rz(-pi) q[1];
rz(-0.4514427) q[2];
sx q[2];
rz(-0.57381781) q[2];
sx q[2];
rz(-0.30893477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42106146) q[1];
sx q[1];
rz(-2.1481247) q[1];
sx q[1];
rz(-2.4666822) q[1];
rz(-pi) q[2];
rz(2.014767) q[3];
sx q[3];
rz(-0.68279632) q[3];
sx q[3];
rz(1.7449774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69349217) q[2];
sx q[2];
rz(-0.38303146) q[2];
sx q[2];
rz(-0.4105655) q[2];
rz(-0.52452123) q[3];
sx q[3];
rz(-0.21206847) q[3];
sx q[3];
rz(-1.1493523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18325026) q[0];
sx q[0];
rz(-2.9999833) q[0];
sx q[0];
rz(0.70567065) q[0];
rz(-1.2880098) q[1];
sx q[1];
rz(-0.57669222) q[1];
sx q[1];
rz(-2.0075331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9767446) q[0];
sx q[0];
rz(-2.3253092) q[0];
sx q[0];
rz(2.1378986) q[0];
rz(-pi) q[1];
rz(-2.6666849) q[2];
sx q[2];
rz(-0.38198745) q[2];
sx q[2];
rz(-1.8595921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4422851) q[1];
sx q[1];
rz(-2.2622427) q[1];
sx q[1];
rz(2.0705018) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5763578) q[3];
sx q[3];
rz(-1.7218523) q[3];
sx q[3];
rz(2.8961669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80295339) q[2];
sx q[2];
rz(-2.6889375) q[2];
sx q[2];
rz(-2.6111641) q[2];
rz(0.28702921) q[3];
sx q[3];
rz(-2.65286) q[3];
sx q[3];
rz(-0.73634017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257783) q[0];
sx q[0];
rz(-1.1017355) q[0];
sx q[0];
rz(-1.0544448) q[0];
rz(-1.6619445) q[1];
sx q[1];
rz(-0.92250657) q[1];
sx q[1];
rz(1.066801) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118059) q[0];
sx q[0];
rz(-1.3830782) q[0];
sx q[0];
rz(0.037506692) q[0];
rz(-1.4786903) q[2];
sx q[2];
rz(-2.2450441) q[2];
sx q[2];
rz(2.9363652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2640564) q[1];
sx q[1];
rz(-2.0320368) q[1];
sx q[1];
rz(0.25564927) q[1];
rz(-2.2852632) q[3];
sx q[3];
rz(-1.3420055) q[3];
sx q[3];
rz(-1.3135873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5532784) q[2];
sx q[2];
rz(-1.30043) q[2];
sx q[2];
rz(2.039457) q[2];
rz(-2.6831902) q[3];
sx q[3];
rz(-1.6130684) q[3];
sx q[3];
rz(-2.5100191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.691953) q[0];
sx q[0];
rz(-2.3821558) q[0];
sx q[0];
rz(-2.6547883) q[0];
rz(1.5961157) q[1];
sx q[1];
rz(-0.38914746) q[1];
sx q[1];
rz(-0.62852377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768054) q[0];
sx q[0];
rz(-2.724455) q[0];
sx q[0];
rz(-1.9030626) q[0];
x q[1];
rz(2.9942368) q[2];
sx q[2];
rz(-2.1389845) q[2];
sx q[2];
rz(-0.94765608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13884096) q[1];
sx q[1];
rz(-1.758444) q[1];
sx q[1];
rz(1.1739338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1194589) q[3];
sx q[3];
rz(-1.6319773) q[3];
sx q[3];
rz(-0.010329846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5635166) q[2];
sx q[2];
rz(-3.0264644) q[2];
sx q[2];
rz(-2.3797177) q[2];
rz(-0.21841194) q[3];
sx q[3];
rz(-0.29774791) q[3];
sx q[3];
rz(-2.0930321) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8258284) q[0];
sx q[0];
rz(-0.35218069) q[0];
sx q[0];
rz(0.6045565) q[0];
rz(1.8479337) q[1];
sx q[1];
rz(-2.3217521) q[1];
sx q[1];
rz(2.7972319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647756) q[0];
sx q[0];
rz(-1.6533543) q[0];
sx q[0];
rz(0.020391012) q[0];
rz(-pi) q[1];
rz(2.5141469) q[2];
sx q[2];
rz(-1.2973669) q[2];
sx q[2];
rz(-0.16579096) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3185259) q[1];
sx q[1];
rz(-1.112655) q[1];
sx q[1];
rz(1.0847102) q[1];
x q[2];
rz(-2.2456393) q[3];
sx q[3];
rz(-1.1955441) q[3];
sx q[3];
rz(-2.8651875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1206426) q[2];
sx q[2];
rz(-0.63434333) q[2];
sx q[2];
rz(2.9316588) q[2];
rz(2.0338992) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(0.2618739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27860677) q[0];
sx q[0];
rz(-0.59026533) q[0];
sx q[0];
rz(-0.069409542) q[0];
rz(0.033992652) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(-0.45190826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2632217) q[0];
sx q[0];
rz(-0.42455205) q[0];
sx q[0];
rz(0.39841899) q[0];
rz(-2.185427) q[2];
sx q[2];
rz(-0.98041222) q[2];
sx q[2];
rz(-2.0942795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44575366) q[1];
sx q[1];
rz(-1.5924982) q[1];
sx q[1];
rz(0.0071956102) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3159368) q[3];
sx q[3];
rz(-1.3075124) q[3];
sx q[3];
rz(1.7039434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41428337) q[2];
sx q[2];
rz(-2.3392129) q[2];
sx q[2];
rz(1.1451954) q[2];
rz(-2.1833873) q[3];
sx q[3];
rz(-1.2088135) q[3];
sx q[3];
rz(-0.27829471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3713643) q[0];
sx q[0];
rz(-2.6537277) q[0];
sx q[0];
rz(2.2801953) q[0];
rz(-2.1791012) q[1];
sx q[1];
rz(-0.95314127) q[1];
sx q[1];
rz(0.15693396) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10805865) q[0];
sx q[0];
rz(-2.3367051) q[0];
sx q[0];
rz(2.0867086) q[0];
rz(-pi) q[1];
rz(-1.8489776) q[2];
sx q[2];
rz(-0.37019324) q[2];
sx q[2];
rz(0.37281677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3536071) q[1];
sx q[1];
rz(-2.3448886) q[1];
sx q[1];
rz(-2.4322047) q[1];
rz(-pi) q[2];
rz(1.7926927) q[3];
sx q[3];
rz(-2.923344) q[3];
sx q[3];
rz(-1.9980007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12338766) q[2];
sx q[2];
rz(-0.79424477) q[2];
sx q[2];
rz(-2.0920853) q[2];
rz(0.47979313) q[3];
sx q[3];
rz(-0.19194651) q[3];
sx q[3];
rz(-2.979579) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70742339) q[0];
sx q[0];
rz(-0.24288414) q[0];
sx q[0];
rz(-2.6450787) q[0];
rz(1.4855344) q[1];
sx q[1];
rz(-2.2746706) q[1];
sx q[1];
rz(0.022580126) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2537947) q[0];
sx q[0];
rz(-2.8631475) q[0];
sx q[0];
rz(-0.87775567) q[0];
x q[1];
rz(-3.0737799) q[2];
sx q[2];
rz(-1.4275506) q[2];
sx q[2];
rz(-0.80729173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1068017) q[1];
sx q[1];
rz(-0.80134922) q[1];
sx q[1];
rz(-1.31557) q[1];
x q[2];
rz(-2.1053505) q[3];
sx q[3];
rz(-1.4009985) q[3];
sx q[3];
rz(-1.9056729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5217343) q[2];
sx q[2];
rz(-2.0010184) q[2];
sx q[2];
rz(0.015259585) q[2];
rz(-0.36733019) q[3];
sx q[3];
rz(-2.6945249) q[3];
sx q[3];
rz(-0.36146155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5527282) q[0];
sx q[0];
rz(-0.22607729) q[0];
sx q[0];
rz(3.0638301) q[0];
rz(-0.11847682) q[1];
sx q[1];
rz(-2.8006554) q[1];
sx q[1];
rz(2.4854787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7897624) q[0];
sx q[0];
rz(-1.8740204) q[0];
sx q[0];
rz(-3.0089761) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.838573) q[2];
sx q[2];
rz(-1.736044) q[2];
sx q[2];
rz(-2.6337773) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39086884) q[1];
sx q[1];
rz(-1.0926529) q[1];
sx q[1];
rz(2.6873846) q[1];
rz(-pi) q[2];
rz(0.55226294) q[3];
sx q[3];
rz(-1.1202462) q[3];
sx q[3];
rz(1.0737863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6180827) q[2];
sx q[2];
rz(-1.221035) q[2];
sx q[2];
rz(2.6005884) q[2];
rz(2.6909761) q[3];
sx q[3];
rz(-1.485598) q[3];
sx q[3];
rz(0.97644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4377874) q[0];
sx q[0];
rz(-1.3492275) q[0];
sx q[0];
rz(3.034814) q[0];
rz(-1.5497426) q[1];
sx q[1];
rz(-2.6390862) q[1];
sx q[1];
rz(-1.3776616) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1049809) q[0];
sx q[0];
rz(-2.689765) q[0];
sx q[0];
rz(1.4696448) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5306931) q[2];
sx q[2];
rz(-2.7654117) q[2];
sx q[2];
rz(3.0665325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71860524) q[1];
sx q[1];
rz(-1.3235301) q[1];
sx q[1];
rz(0.24691564) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.315134) q[3];
sx q[3];
rz(-2.4240806) q[3];
sx q[3];
rz(-0.96271587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0489073) q[2];
sx q[2];
rz(-1.4212298) q[2];
sx q[2];
rz(2.8468724) q[2];
rz(-0.60485351) q[3];
sx q[3];
rz(-1.9113144) q[3];
sx q[3];
rz(-3.0321583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(0.43815628) q[0];
sx q[0];
rz(-1.3180757) q[0];
sx q[0];
rz(-0.99676589) q[0];
rz(0.017771487) q[1];
sx q[1];
rz(-1.8780864) q[1];
sx q[1];
rz(-1.3395739) q[1];
rz(2.3370635) q[2];
sx q[2];
rz(-1.1297516) q[2];
sx q[2];
rz(-0.17937112) q[2];
rz(-0.65405398) q[3];
sx q[3];
rz(-0.97915836) q[3];
sx q[3];
rz(-0.98042515) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
