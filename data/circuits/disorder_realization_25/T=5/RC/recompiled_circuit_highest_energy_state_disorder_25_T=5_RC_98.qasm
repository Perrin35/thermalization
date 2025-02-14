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
rz(-0.93936062) q[0];
sx q[0];
rz(-2.0225749) q[0];
sx q[0];
rz(0.033666704) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(-1.6330947) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93532186) q[0];
sx q[0];
rz(-1.540872) q[0];
sx q[0];
rz(1.6851835) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.066285) q[2];
sx q[2];
rz(-1.8641553) q[2];
sx q[2];
rz(-1.8541898) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7693387) q[1];
sx q[1];
rz(-0.70972356) q[1];
sx q[1];
rz(0.22270568) q[1];
rz(-pi) q[2];
rz(0.81013362) q[3];
sx q[3];
rz(-1.0768693) q[3];
sx q[3];
rz(-2.877748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9713251) q[2];
sx q[2];
rz(-0.89189947) q[2];
sx q[2];
rz(2.3049221) q[2];
rz(-2.3352046) q[3];
sx q[3];
rz(-0.92525768) q[3];
sx q[3];
rz(0.18536082) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065501325) q[0];
sx q[0];
rz(-1.0668904) q[0];
sx q[0];
rz(-1.9245032) q[0];
rz(0.20054664) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(-2.5977871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018262176) q[0];
sx q[0];
rz(-1.8866492) q[0];
sx q[0];
rz(1.3297775) q[0];
rz(-pi) q[1];
rz(-0.020140127) q[2];
sx q[2];
rz(-2.5465362) q[2];
sx q[2];
rz(0.12660566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3873991) q[1];
sx q[1];
rz(-0.48739823) q[1];
sx q[1];
rz(0.91895603) q[1];
x q[2];
rz(-3.1065349) q[3];
sx q[3];
rz(-0.89588273) q[3];
sx q[3];
rz(-1.9652776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66889888) q[2];
sx q[2];
rz(-2.3987179) q[2];
sx q[2];
rz(-1.5354068) q[2];
rz(1.499048) q[3];
sx q[3];
rz(-2.5583772) q[3];
sx q[3];
rz(-1.8892052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1010308) q[0];
sx q[0];
rz(-2.5767548) q[0];
sx q[0];
rz(0.0095796883) q[0];
rz(1.0293695) q[1];
sx q[1];
rz(-1.5681489) q[1];
sx q[1];
rz(1.3363438) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1270008) q[0];
sx q[0];
rz(-0.83084241) q[0];
sx q[0];
rz(-0.71944351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0687683) q[2];
sx q[2];
rz(-1.8548551) q[2];
sx q[2];
rz(-1.9845923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6273986) q[1];
sx q[1];
rz(-1.3954867) q[1];
sx q[1];
rz(1.2199376) q[1];
rz(-pi) q[2];
rz(-1.1747178) q[3];
sx q[3];
rz(-1.4142591) q[3];
sx q[3];
rz(2.1780632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0336527) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(1.8281724) q[2];
rz(0.039947346) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(-1.1279172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0082598) q[0];
sx q[0];
rz(-1.8513716) q[0];
sx q[0];
rz(2.7449352) q[0];
rz(-0.89504129) q[1];
sx q[1];
rz(-1.7170693) q[1];
sx q[1];
rz(2.3779714) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0741054) q[0];
sx q[0];
rz(-0.25831902) q[0];
sx q[0];
rz(-1.8909569) q[0];
x q[1];
rz(3.0020048) q[2];
sx q[2];
rz(-3.0198823) q[2];
sx q[2];
rz(1.9892529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29429193) q[1];
sx q[1];
rz(-2.9025893) q[1];
sx q[1];
rz(1.4579178) q[1];
rz(2.850789) q[3];
sx q[3];
rz(-1.606308) q[3];
sx q[3];
rz(-2.5574656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7274196) q[2];
sx q[2];
rz(-1.6215934) q[2];
sx q[2];
rz(2.011389) q[2];
rz(1.1388418) q[3];
sx q[3];
rz(-1.1915221) q[3];
sx q[3];
rz(1.1677008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.8324757) q[0];
sx q[0];
rz(-0.55823767) q[0];
sx q[0];
rz(-0.48156893) q[0];
rz(-2.8140977) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(-2.17735) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5835579) q[0];
sx q[0];
rz(-2.6917771) q[0];
sx q[0];
rz(-2.0925453) q[0];
x q[1];
rz(1.8580074) q[2];
sx q[2];
rz(-2.5406844) q[2];
sx q[2];
rz(-1.7696966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9157) q[1];
sx q[1];
rz(-1.6927162) q[1];
sx q[1];
rz(-1.3242766) q[1];
x q[2];
rz(1.3664512) q[3];
sx q[3];
rz(-1.1875523) q[3];
sx q[3];
rz(0.97571532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1714736) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(1.6469275) q[2];
rz(1.0129048) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(2.6265898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9377015) q[0];
sx q[0];
rz(-1.0045811) q[0];
sx q[0];
rz(-2.0830233) q[0];
rz(2.2189498) q[1];
sx q[1];
rz(-2.1244996) q[1];
sx q[1];
rz(-0.10291084) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.218821) q[0];
sx q[0];
rz(-2.0334683) q[0];
sx q[0];
rz(0.86313049) q[0];
x q[1];
rz(0.40252588) q[2];
sx q[2];
rz(-2.5136097) q[2];
sx q[2];
rz(-1.8875853) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7988777) q[1];
sx q[1];
rz(-2.5763576) q[1];
sx q[1];
rz(1.2223627) q[1];
x q[2];
rz(1.5636121) q[3];
sx q[3];
rz(-0.71643351) q[3];
sx q[3];
rz(2.7897657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70204488) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(2.8257545) q[2];
rz(-0.99509197) q[3];
sx q[3];
rz(-1.2732048) q[3];
sx q[3];
rz(-0.72527138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(2.134326) q[0];
sx q[0];
rz(-2.8747929) q[0];
sx q[0];
rz(2.5079978) q[0];
rz(-2.7569356) q[1];
sx q[1];
rz(-1.1793914) q[1];
sx q[1];
rz(0.42919174) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4576526) q[0];
sx q[0];
rz(-0.54163092) q[0];
sx q[0];
rz(-1.0628048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6056988) q[2];
sx q[2];
rz(-2.2560824) q[2];
sx q[2];
rz(-1.8165782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0269484) q[1];
sx q[1];
rz(-2.1430795) q[1];
sx q[1];
rz(1.9932822) q[1];
rz(-pi) q[2];
rz(-2.1294534) q[3];
sx q[3];
rz(-1.9724761) q[3];
sx q[3];
rz(0.60810773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0523494) q[2];
sx q[2];
rz(-2.3640335) q[2];
sx q[2];
rz(2.8746129) q[2];
rz(0.071768196) q[3];
sx q[3];
rz(-2.1238582) q[3];
sx q[3];
rz(1.4192386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778022) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(1.4005533) q[0];
rz(1.9969223) q[1];
sx q[1];
rz(-2.0796227) q[1];
sx q[1];
rz(-0.58715075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42204866) q[0];
sx q[0];
rz(-1.9912155) q[0];
sx q[0];
rz(-0.69164135) q[0];
rz(-pi) q[1];
rz(-1.4540423) q[2];
sx q[2];
rz(-1.1551305) q[2];
sx q[2];
rz(2.1374201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59881951) q[1];
sx q[1];
rz(-2.2428998) q[1];
sx q[1];
rz(-1.2818976) q[1];
rz(-pi) q[2];
rz(-2.8244157) q[3];
sx q[3];
rz(-0.97517386) q[3];
sx q[3];
rz(0.45444876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2862386) q[2];
sx q[2];
rz(-1.9619532) q[2];
sx q[2];
rz(-0.41909763) q[2];
rz(1.6176443) q[3];
sx q[3];
rz(-0.91846275) q[3];
sx q[3];
rz(-1.5036478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54290259) q[0];
sx q[0];
rz(-1.221523) q[0];
sx q[0];
rz(-0.30938095) q[0];
rz(-1.7912553) q[1];
sx q[1];
rz(-0.58809909) q[1];
sx q[1];
rz(2.5877171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5039936) q[0];
sx q[0];
rz(-1.1553815) q[0];
sx q[0];
rz(1.2211772) q[0];
x q[1];
rz(-1.925681) q[2];
sx q[2];
rz(-2.7352374) q[2];
sx q[2];
rz(-2.2333437) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83912056) q[1];
sx q[1];
rz(-2.1697147) q[1];
sx q[1];
rz(0.37521404) q[1];
rz(-0.36333824) q[3];
sx q[3];
rz(-0.73550341) q[3];
sx q[3];
rz(0.21061646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89173335) q[2];
sx q[2];
rz(-1.5313818) q[2];
sx q[2];
rz(-1.7943133) q[2];
rz(0.76752082) q[3];
sx q[3];
rz(-1.9556421) q[3];
sx q[3];
rz(-0.84848136) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42221853) q[0];
sx q[0];
rz(-0.30000559) q[0];
sx q[0];
rz(-1.2147709) q[0];
rz(-2.1396554) q[1];
sx q[1];
rz(-1.6286214) q[1];
sx q[1];
rz(-0.92438662) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533646) q[0];
sx q[0];
rz(-1.9583252) q[0];
sx q[0];
rz(0.94382443) q[0];
rz(-0.027410313) q[2];
sx q[2];
rz(-2.1152585) q[2];
sx q[2];
rz(-1.4278864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1781285) q[1];
sx q[1];
rz(-0.22451065) q[1];
sx q[1];
rz(3.0369395) q[1];
rz(-0.59238418) q[3];
sx q[3];
rz(-2.1796556) q[3];
sx q[3];
rz(0.99623535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34710106) q[2];
sx q[2];
rz(-2.2091986) q[2];
sx q[2];
rz(-2.7744228) q[2];
rz(-1.9764887) q[3];
sx q[3];
rz(-2.1737183) q[3];
sx q[3];
rz(0.86618209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544871) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(0.26609303) q[1];
sx q[1];
rz(-0.91231822) q[1];
sx q[1];
rz(-1.4516861) q[1];
rz(0.75917203) q[2];
sx q[2];
rz(-2.7552032) q[2];
sx q[2];
rz(-0.0061782171) q[2];
rz(-1.2721619) q[3];
sx q[3];
rz(-2.5187026) q[3];
sx q[3];
rz(-1.8395195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
