OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17570198) q[0];
sx q[0];
rz(5.0134563) q[0];
sx q[0];
rz(11.391517) q[0];
rz(2.0570316) q[1];
sx q[1];
rz(-1.6697872) q[1];
sx q[1];
rz(-0.56301277) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969471) q[0];
sx q[0];
rz(-1.5176231) q[0];
sx q[0];
rz(-0.19617686) q[0];
x q[1];
rz(2.1839574) q[2];
sx q[2];
rz(-2.4102825) q[2];
sx q[2];
rz(-0.72145185) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9809499) q[1];
sx q[1];
rz(-1.8997571) q[1];
sx q[1];
rz(-1.4018707) q[1];
rz(-2.7607542) q[3];
sx q[3];
rz(-1.9054753) q[3];
sx q[3];
rz(0.0045485529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6072443) q[2];
sx q[2];
rz(-1.6515942) q[2];
sx q[2];
rz(1.6687757) q[2];
rz(-1.8913174) q[3];
sx q[3];
rz(-1.5273124) q[3];
sx q[3];
rz(-0.97243029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8984574) q[0];
sx q[0];
rz(-3.0861096) q[0];
sx q[0];
rz(2.6836416) q[0];
rz(3.0398439) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(2.8153458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3889623) q[0];
sx q[0];
rz(-1.9174308) q[0];
sx q[0];
rz(-1.1244357) q[0];
rz(-pi) q[1];
rz(-0.63133755) q[2];
sx q[2];
rz(-2.7777618) q[2];
sx q[2];
rz(-1.199786) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40405289) q[1];
sx q[1];
rz(-1.0657601) q[1];
sx q[1];
rz(-0.86143273) q[1];
rz(-pi) q[2];
rz(0.34290727) q[3];
sx q[3];
rz(-2.3686633) q[3];
sx q[3];
rz(2.8445304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.132167) q[2];
sx q[2];
rz(-0.99401179) q[2];
sx q[2];
rz(-1.2635292) q[2];
rz(0.020091232) q[3];
sx q[3];
rz(-2.3720436) q[3];
sx q[3];
rz(1.2300389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024071368) q[0];
sx q[0];
rz(-0.063952359) q[0];
sx q[0];
rz(-2.5653895) q[0];
rz(-1.4441215) q[1];
sx q[1];
rz(-1.2497808) q[1];
sx q[1];
rz(-1.8796657) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1133973) q[0];
sx q[0];
rz(-1.5388267) q[0];
sx q[0];
rz(-1.8799923) q[0];
x q[1];
rz(0.52719614) q[2];
sx q[2];
rz(-1.1644852) q[2];
sx q[2];
rz(0.57002178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0351959) q[1];
sx q[1];
rz(-0.96265618) q[1];
sx q[1];
rz(1.2713096) q[1];
x q[2];
rz(-1.0141054) q[3];
sx q[3];
rz(-2.2051349) q[3];
sx q[3];
rz(1.2403099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.55804092) q[2];
sx q[2];
rz(-2.5463107) q[2];
sx q[2];
rz(2.9884647) q[2];
rz(-0.88614744) q[3];
sx q[3];
rz(-1.359442) q[3];
sx q[3];
rz(-1.2019633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2578289) q[0];
sx q[0];
rz(-1.9572636) q[0];
sx q[0];
rz(-0.15876874) q[0];
rz(2.1067045) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(0.75611702) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0949815) q[0];
sx q[0];
rz(-1.621052) q[0];
sx q[0];
rz(-3.0836225) q[0];
rz(-1.6768084) q[2];
sx q[2];
rz(-0.25463018) q[2];
sx q[2];
rz(2.1531596) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9489574) q[1];
sx q[1];
rz(-1.3657943) q[1];
sx q[1];
rz(0.24688772) q[1];
rz(-2.9063734) q[3];
sx q[3];
rz(-1.7718214) q[3];
sx q[3];
rz(0.26409322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5990344) q[2];
sx q[2];
rz(-2.2381146) q[2];
sx q[2];
rz(0.8503882) q[2];
rz(1.6709857) q[3];
sx q[3];
rz(-1.8149523) q[3];
sx q[3];
rz(0.82408041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4849512) q[0];
sx q[0];
rz(-1.4018207) q[0];
sx q[0];
rz(-2.9630419) q[0];
rz(1.2483596) q[1];
sx q[1];
rz(-1.5310042) q[1];
sx q[1];
rz(0.53370968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2634791) q[0];
sx q[0];
rz(-1.7537259) q[0];
sx q[0];
rz(-0.35533743) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87119331) q[2];
sx q[2];
rz(-1.9661926) q[2];
sx q[2];
rz(1.4104539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7099784) q[1];
sx q[1];
rz(-1.6233161) q[1];
sx q[1];
rz(0.90265982) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0751791) q[3];
sx q[3];
rz(-1.4850067) q[3];
sx q[3];
rz(-0.64297167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3137714) q[2];
sx q[2];
rz(-1.8905819) q[2];
sx q[2];
rz(2.0884464) q[2];
rz(-2.961212) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(-0.36261305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5798222) q[0];
sx q[0];
rz(-2.1403911) q[0];
sx q[0];
rz(-1.2286105) q[0];
rz(0.46663943) q[1];
sx q[1];
rz(-1.9231611) q[1];
sx q[1];
rz(-0.12900464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59959847) q[0];
sx q[0];
rz(-0.16762531) q[0];
sx q[0];
rz(-0.92131281) q[0];
rz(-pi) q[1];
rz(0.1618218) q[2];
sx q[2];
rz(-2.2941049) q[2];
sx q[2];
rz(-0.79297334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.946917) q[1];
sx q[1];
rz(-1.2763775) q[1];
sx q[1];
rz(0.63398449) q[1];
rz(-pi) q[2];
rz(0.13997649) q[3];
sx q[3];
rz(-1.8365321) q[3];
sx q[3];
rz(-1.1418726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0064319) q[2];
sx q[2];
rz(-2.5401523) q[2];
sx q[2];
rz(2.637291) q[2];
rz(2.048061) q[3];
sx q[3];
rz(-1.8176327) q[3];
sx q[3];
rz(2.1849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4920376) q[0];
sx q[0];
rz(-1.858359) q[0];
sx q[0];
rz(-1.5851703) q[0];
rz(-2.5632312) q[1];
sx q[1];
rz(-1.2588986) q[1];
sx q[1];
rz(0.10231054) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033275398) q[0];
sx q[0];
rz(-0.96848291) q[0];
sx q[0];
rz(2.6541416) q[0];
rz(1.341171) q[2];
sx q[2];
rz(-1.183325) q[2];
sx q[2];
rz(1.4619399) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.426541) q[1];
sx q[1];
rz(-0.3449966) q[1];
sx q[1];
rz(-2.662602) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75179799) q[3];
sx q[3];
rz(-2.2357142) q[3];
sx q[3];
rz(-0.81937688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9461225) q[2];
sx q[2];
rz(-1.8778233) q[2];
sx q[2];
rz(-2.5275687) q[2];
rz(-0.92769235) q[3];
sx q[3];
rz(-1.5427019) q[3];
sx q[3];
rz(2.4641002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55352655) q[0];
sx q[0];
rz(-2.8480242) q[0];
sx q[0];
rz(-1.9422096) q[0];
rz(-1.9623494) q[1];
sx q[1];
rz(-2.113138) q[1];
sx q[1];
rz(-2.1882122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3971553) q[0];
sx q[0];
rz(-2.1040475) q[0];
sx q[0];
rz(-2.5313733) q[0];
rz(-pi) q[1];
x q[1];
rz(1.032845) q[2];
sx q[2];
rz(-1.6643833) q[2];
sx q[2];
rz(0.46582281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7443707) q[1];
sx q[1];
rz(-0.83852856) q[1];
sx q[1];
rz(2.9930755) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9913919) q[3];
sx q[3];
rz(-2.6552192) q[3];
sx q[3];
rz(1.2260162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74932468) q[2];
sx q[2];
rz(-0.12568036) q[2];
sx q[2];
rz(2.9607062) q[2];
rz(-2.0285897) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(2.7614312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0893294) q[0];
sx q[0];
rz(-0.90173975) q[0];
sx q[0];
rz(-2.3774636) q[0];
rz(2.1185421) q[1];
sx q[1];
rz(-1.5307531) q[1];
sx q[1];
rz(-1.6463564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9272629) q[0];
sx q[0];
rz(-1.8678878) q[0];
sx q[0];
rz(-1.731575) q[0];
rz(-pi) q[1];
rz(1.1029215) q[2];
sx q[2];
rz(-2.4870424) q[2];
sx q[2];
rz(-0.86328116) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0937522) q[1];
sx q[1];
rz(-1.6421659) q[1];
sx q[1];
rz(-2.755295) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71919433) q[3];
sx q[3];
rz(-0.56108863) q[3];
sx q[3];
rz(2.6230375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19227795) q[2];
sx q[2];
rz(-1.01769) q[2];
sx q[2];
rz(-0.44753543) q[2];
rz(-3.0053075) q[3];
sx q[3];
rz(-0.49301967) q[3];
sx q[3];
rz(-2.5435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3560781) q[0];
sx q[0];
rz(-1.0250174) q[0];
sx q[0];
rz(-0.4726952) q[0];
rz(0.39974943) q[1];
sx q[1];
rz(-0.70416299) q[1];
sx q[1];
rz(-2.3230816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8560573) q[0];
sx q[0];
rz(-1.7489079) q[0];
sx q[0];
rz(0.31383236) q[0];
x q[1];
rz(-2.632896) q[2];
sx q[2];
rz(-2.3983068) q[2];
sx q[2];
rz(2.7270779) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.012072969) q[1];
sx q[1];
rz(-0.32271472) q[1];
sx q[1];
rz(1.5606776) q[1];
x q[2];
rz(2.7943939) q[3];
sx q[3];
rz(-1.7337017) q[3];
sx q[3];
rz(-0.81219514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58273903) q[2];
sx q[2];
rz(-1.0820729) q[2];
sx q[2];
rz(1.8633128) q[2];
rz(-1.1996972) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(-2.8619158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3661135) q[0];
sx q[0];
rz(-1.4237325) q[0];
sx q[0];
rz(1.7970418) q[0];
rz(-1.7300425) q[1];
sx q[1];
rz(-2.330214) q[1];
sx q[1];
rz(-0.65705962) q[1];
rz(-1.5700339) q[2];
sx q[2];
rz(-1.7117725) q[2];
sx q[2];
rz(2.6374809) q[2];
rz(-2.2349002) q[3];
sx q[3];
rz(-1.2994874) q[3];
sx q[3];
rz(-0.52728925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
