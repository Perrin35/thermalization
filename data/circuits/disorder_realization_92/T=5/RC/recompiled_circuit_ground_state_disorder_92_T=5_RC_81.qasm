OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.3761223) q[0];
sx q[0];
rz(2.0630615) q[0];
sx q[0];
rz(10.800092) q[0];
rz(3.5612192) q[1];
sx q[1];
rz(1.761275) q[1];
sx q[1];
rz(7.0786369) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6559443) q[0];
sx q[0];
rz(-0.38165452) q[0];
sx q[0];
rz(-0.21393805) q[0];
rz(-pi) q[1];
rz(-2.8654557) q[2];
sx q[2];
rz(-2.727446) q[2];
sx q[2];
rz(0.93246952) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7955543) q[1];
sx q[1];
rz(-0.70540326) q[1];
sx q[1];
rz(2.8349174) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0555435) q[3];
sx q[3];
rz(-1.8294191) q[3];
sx q[3];
rz(-1.8515967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2835283) q[2];
sx q[2];
rz(-2.8499446) q[2];
sx q[2];
rz(-2.7533599) q[2];
rz(-2.9325716) q[3];
sx q[3];
rz(-1.4744604) q[3];
sx q[3];
rz(2.2504263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.46880576) q[0];
sx q[0];
rz(-0.46872941) q[0];
sx q[0];
rz(-0.76954532) q[0];
rz(-1.0308713) q[1];
sx q[1];
rz(-2.2831235) q[1];
sx q[1];
rz(1.7147725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1772004) q[0];
sx q[0];
rz(-1.3528723) q[0];
sx q[0];
rz(-2.833018) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6278817) q[2];
sx q[2];
rz(-0.89623129) q[2];
sx q[2];
rz(-3.0115173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0116358) q[1];
sx q[1];
rz(-2.4677688) q[1];
sx q[1];
rz(-0.41566276) q[1];
rz(-pi) q[2];
rz(-0.053228037) q[3];
sx q[3];
rz(-0.8570519) q[3];
sx q[3];
rz(0.66859879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1398805) q[2];
sx q[2];
rz(-0.44846815) q[2];
sx q[2];
rz(-1.9072388) q[2];
rz(-0.52302304) q[3];
sx q[3];
rz(-2.0165899) q[3];
sx q[3];
rz(-2.4856868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457542) q[0];
sx q[0];
rz(-3.1390751) q[0];
sx q[0];
rz(-2.5939831) q[0];
rz(0.86946636) q[1];
sx q[1];
rz(-2.165386) q[1];
sx q[1];
rz(-1.4617823) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47785178) q[0];
sx q[0];
rz(-1.2292851) q[0];
sx q[0];
rz(2.6966901) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8274948) q[2];
sx q[2];
rz(-2.3163598) q[2];
sx q[2];
rz(-2.1113393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1717277) q[1];
sx q[1];
rz(-1.1369924) q[1];
sx q[1];
rz(-2.8126905) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8601772) q[3];
sx q[3];
rz(-2.0915935) q[3];
sx q[3];
rz(-1.2788683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72472858) q[2];
sx q[2];
rz(-0.29837307) q[2];
sx q[2];
rz(2.5073063) q[2];
rz(2.0554845) q[3];
sx q[3];
rz(-2.6606798) q[3];
sx q[3];
rz(1.1336795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5241765) q[0];
sx q[0];
rz(-0.15707459) q[0];
sx q[0];
rz(-1.5217391) q[0];
rz(-2.1583083) q[1];
sx q[1];
rz(-1.1055999) q[1];
sx q[1];
rz(-0.081667893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33591341) q[0];
sx q[0];
rz(-1.6821852) q[0];
sx q[0];
rz(2.1301444) q[0];
rz(-pi) q[1];
rz(-1.2377023) q[2];
sx q[2];
rz(-1.5577847) q[2];
sx q[2];
rz(-1.5934554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6165627) q[1];
sx q[1];
rz(-1.9989415) q[1];
sx q[1];
rz(-1.031428) q[1];
rz(-pi) q[2];
rz(-1.0809298) q[3];
sx q[3];
rz(-1.5233763) q[3];
sx q[3];
rz(-1.0256378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4982768) q[2];
sx q[2];
rz(-1.6643107) q[2];
sx q[2];
rz(2.2225883) q[2];
rz(2.8194341) q[3];
sx q[3];
rz(-1.6291658) q[3];
sx q[3];
rz(-0.34644103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0568986) q[0];
sx q[0];
rz(-1.1321122) q[0];
sx q[0];
rz(-1.3380916) q[0];
rz(2.7390506) q[1];
sx q[1];
rz(-1.2502547) q[1];
sx q[1];
rz(-1.5571627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65619233) q[0];
sx q[0];
rz(-1.8372522) q[0];
sx q[0];
rz(-2.0977667) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3473559) q[2];
sx q[2];
rz(-1.734724) q[2];
sx q[2];
rz(-1.9203223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30774035) q[1];
sx q[1];
rz(-0.55451143) q[1];
sx q[1];
rz(-0.28077287) q[1];
x q[2];
rz(-1.1013056) q[3];
sx q[3];
rz(-1.3302505) q[3];
sx q[3];
rz(-0.32215531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0255787) q[2];
sx q[2];
rz(-1.2776926) q[2];
sx q[2];
rz(0.73520994) q[2];
rz(2.0871346) q[3];
sx q[3];
rz(-2.1078883) q[3];
sx q[3];
rz(1.916877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16578199) q[0];
sx q[0];
rz(-2.7175792) q[0];
sx q[0];
rz(-1.8884416) q[0];
rz(-1.5660628) q[1];
sx q[1];
rz(-1.1352481) q[1];
sx q[1];
rz(0.53322405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5767007) q[0];
sx q[0];
rz(-0.78540388) q[0];
sx q[0];
rz(1.9068973) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27006939) q[2];
sx q[2];
rz(-0.52723072) q[2];
sx q[2];
rz(0.6747077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8044843) q[1];
sx q[1];
rz(-1.7586367) q[1];
sx q[1];
rz(-0.2262745) q[1];
rz(-pi) q[2];
rz(-0.76119555) q[3];
sx q[3];
rz(-2.4223339) q[3];
sx q[3];
rz(2.066589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6553361) q[2];
sx q[2];
rz(-2.76261) q[2];
sx q[2];
rz(0.8302702) q[2];
rz(1.8298979) q[3];
sx q[3];
rz(-1.0887086) q[3];
sx q[3];
rz(2.8970498) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.526392) q[0];
sx q[0];
rz(-1.6561693) q[0];
sx q[0];
rz(-2.3080589) q[0];
rz(2.1444881) q[1];
sx q[1];
rz(-1.3778967) q[1];
sx q[1];
rz(1.5514143) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9587904) q[0];
sx q[0];
rz(-2.0062464) q[0];
sx q[0];
rz(1.4884401) q[0];
rz(-2.1623108) q[2];
sx q[2];
rz(-1.8663317) q[2];
sx q[2];
rz(-2.599567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60234705) q[1];
sx q[1];
rz(-1.6825469) q[1];
sx q[1];
rz(1.4503327) q[1];
rz(-pi) q[2];
rz(-0.40682934) q[3];
sx q[3];
rz(-2.5689044) q[3];
sx q[3];
rz(2.0175479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3414063) q[2];
sx q[2];
rz(-2.1087346) q[2];
sx q[2];
rz(1.7436854) q[2];
rz(-0.37311113) q[3];
sx q[3];
rz(-2.2171376) q[3];
sx q[3];
rz(2.1227409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7031192) q[0];
sx q[0];
rz(-1.6699474) q[0];
sx q[0];
rz(-2.8224509) q[0];
rz(-1.1646264) q[1];
sx q[1];
rz(-1.171037) q[1];
sx q[1];
rz(-3.107531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8933786) q[0];
sx q[0];
rz(-2.3015018) q[0];
sx q[0];
rz(2.7284751) q[0];
rz(-pi) q[1];
rz(1.2581052) q[2];
sx q[2];
rz(-2.3159928) q[2];
sx q[2];
rz(-1.6751877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0553375) q[1];
sx q[1];
rz(-3.1172522) q[1];
sx q[1];
rz(-1.1175977) q[1];
rz(1.930792) q[3];
sx q[3];
rz(-1.0673095) q[3];
sx q[3];
rz(-2.9140123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0762735) q[2];
sx q[2];
rz(-1.3955782) q[2];
sx q[2];
rz(0.070153959) q[2];
rz(-1.8326727) q[3];
sx q[3];
rz(-2.2891243) q[3];
sx q[3];
rz(-3.0064228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604621) q[0];
sx q[0];
rz(-1.141893) q[0];
sx q[0];
rz(-2.6666226) q[0];
rz(-1.9604856) q[1];
sx q[1];
rz(-0.89473748) q[1];
sx q[1];
rz(-2.2244577) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.006893039) q[0];
sx q[0];
rz(-2.1297161) q[0];
sx q[0];
rz(1.8447645) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.195329) q[2];
sx q[2];
rz(-0.72403833) q[2];
sx q[2];
rz(-0.53508102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0168382) q[1];
sx q[1];
rz(-1.6528582) q[1];
sx q[1];
rz(1.7793307) q[1];
x q[2];
rz(0.59983045) q[3];
sx q[3];
rz(-1.8919865) q[3];
sx q[3];
rz(2.7293492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8721623) q[2];
sx q[2];
rz(-1.7717125) q[2];
sx q[2];
rz(-1.5838464) q[2];
rz(-2.9861084) q[3];
sx q[3];
rz(-1.0289861) q[3];
sx q[3];
rz(2.1492929) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601722) q[0];
sx q[0];
rz(-1.6654797) q[0];
sx q[0];
rz(0.36556622) q[0];
rz(-1.9407326) q[1];
sx q[1];
rz(-1.5301842) q[1];
sx q[1];
rz(2.0948804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7378254) q[0];
sx q[0];
rz(-2.1008472) q[0];
sx q[0];
rz(-0.60204317) q[0];
rz(-pi) q[1];
rz(1.0253941) q[2];
sx q[2];
rz(-1.1869245) q[2];
sx q[2];
rz(-2.0929558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2859058) q[1];
sx q[1];
rz(-2.1033204) q[1];
sx q[1];
rz(-2.3175815) q[1];
rz(-pi) q[2];
rz(-0.58308122) q[3];
sx q[3];
rz(-1.1016535) q[3];
sx q[3];
rz(-0.52419477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4552348) q[2];
sx q[2];
rz(-1.1150259) q[2];
sx q[2];
rz(0.70351797) q[2];
rz(-0.55758682) q[3];
sx q[3];
rz(-0.78123868) q[3];
sx q[3];
rz(-0.47344661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2932417) q[0];
sx q[0];
rz(-1.0307172) q[0];
sx q[0];
rz(0.95681554) q[0];
rz(-1.035607) q[1];
sx q[1];
rz(-2.5685812) q[1];
sx q[1];
rz(-2.0249637) q[1];
rz(1.0445486) q[2];
sx q[2];
rz(-2.1696354) q[2];
sx q[2];
rz(-1.3641691) q[2];
rz(-0.21743786) q[3];
sx q[3];
rz(-1.6119381) q[3];
sx q[3];
rz(1.987206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
