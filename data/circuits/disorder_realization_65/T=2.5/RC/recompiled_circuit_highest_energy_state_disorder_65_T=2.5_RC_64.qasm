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
rz(0.12198099) q[0];
sx q[0];
rz(3.0966336) q[0];
sx q[0];
rz(9.3293204) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(-0.065923318) q[1];
sx q[1];
rz(-2.5479868) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9531692) q[0];
sx q[0];
rz(-2.0521518) q[0];
sx q[0];
rz(1.166422) q[0];
x q[1];
rz(2.9903052) q[2];
sx q[2];
rz(-1.4522284) q[2];
sx q[2];
rz(-1.2158082) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2954457) q[1];
sx q[1];
rz(-2.4762003) q[1];
sx q[1];
rz(0.34122463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.024385) q[3];
sx q[3];
rz(-2.5684169) q[3];
sx q[3];
rz(0.53861618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.060071271) q[2];
sx q[2];
rz(-1.2520484) q[2];
sx q[2];
rz(-2.579465) q[2];
rz(0.33956042) q[3];
sx q[3];
rz(-1.5226676) q[3];
sx q[3];
rz(2.0558527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.958441) q[0];
sx q[0];
rz(-2.9957132) q[0];
sx q[0];
rz(-1.9802144) q[0];
rz(-2.9005652) q[1];
sx q[1];
rz(-2.3190505) q[1];
sx q[1];
rz(-2.8384812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7890772) q[0];
sx q[0];
rz(-0.20931986) q[0];
sx q[0];
rz(1.2417488) q[0];
x q[1];
rz(1.4463521) q[2];
sx q[2];
rz(-1.4258988) q[2];
sx q[2];
rz(3.0595487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.0063340291) q[1];
sx q[1];
rz(-0.63992441) q[1];
sx q[1];
rz(2.695822) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0996898) q[3];
sx q[3];
rz(-2.0031824) q[3];
sx q[3];
rz(1.3392836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3031528) q[2];
sx q[2];
rz(-1.2829605) q[2];
sx q[2];
rz(-2.8008833) q[2];
rz(0.43863145) q[3];
sx q[3];
rz(-2.3352968) q[3];
sx q[3];
rz(0.058569245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230147) q[0];
sx q[0];
rz(-0.61841643) q[0];
sx q[0];
rz(-0.83786905) q[0];
rz(-2.2184929) q[1];
sx q[1];
rz(-1.9337312) q[1];
sx q[1];
rz(-0.34742483) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4283963) q[0];
sx q[0];
rz(-1.5335463) q[0];
sx q[0];
rz(-1.6646882) q[0];
x q[1];
rz(-1.4342257) q[2];
sx q[2];
rz(-0.25871181) q[2];
sx q[2];
rz(-1.2070036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76207896) q[1];
sx q[1];
rz(-2.3752932) q[1];
sx q[1];
rz(2.422782) q[1];
rz(-0.33073552) q[3];
sx q[3];
rz(-1.8285255) q[3];
sx q[3];
rz(2.7455519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.122494) q[2];
sx q[2];
rz(-3.0235897) q[2];
sx q[2];
rz(-0.72354358) q[2];
rz(-1.8540234) q[3];
sx q[3];
rz(-2.5836594) q[3];
sx q[3];
rz(-3.0961228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7751854) q[0];
sx q[0];
rz(-2.7025096) q[0];
sx q[0];
rz(1.084569) q[0];
rz(-1.9547801) q[1];
sx q[1];
rz(-1.2326406) q[1];
sx q[1];
rz(-1.4694227) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6220734) q[0];
sx q[0];
rz(-1.5876829) q[0];
sx q[0];
rz(0.030120912) q[0];
rz(-1.182556) q[2];
sx q[2];
rz(-1.2080384) q[2];
sx q[2];
rz(-1.4826536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.96257229) q[1];
sx q[1];
rz(-2.031759) q[1];
sx q[1];
rz(-2.712972) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3923548) q[3];
sx q[3];
rz(-0.47296528) q[3];
sx q[3];
rz(-1.3441031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0161418) q[2];
sx q[2];
rz(-2.2405388) q[2];
sx q[2];
rz(2.2060642) q[2];
rz(2.9901796) q[3];
sx q[3];
rz(-0.97984034) q[3];
sx q[3];
rz(0.51101959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9147515) q[0];
sx q[0];
rz(-1.5578569) q[0];
sx q[0];
rz(2.3874808) q[0];
rz(-2.5304645) q[1];
sx q[1];
rz(-1.1437623) q[1];
sx q[1];
rz(-2.4653844) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0128764) q[0];
sx q[0];
rz(-2.7749847) q[0];
sx q[0];
rz(1.104273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1242179) q[2];
sx q[2];
rz(-2.5141735) q[2];
sx q[2];
rz(-1.2469893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83482997) q[1];
sx q[1];
rz(-1.7422297) q[1];
sx q[1];
rz(-1.7333561) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0178278) q[3];
sx q[3];
rz(-2.5522557) q[3];
sx q[3];
rz(-0.86245727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1375229) q[2];
sx q[2];
rz(-2.4293032) q[2];
sx q[2];
rz(-2.3884657) q[2];
rz(-2.5774041) q[3];
sx q[3];
rz(-2.0337532) q[3];
sx q[3];
rz(-2.3851725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9516528) q[0];
sx q[0];
rz(-0.48007444) q[0];
sx q[0];
rz(2.9204364) q[0];
rz(-3.0033374) q[1];
sx q[1];
rz(-1.6860551) q[1];
sx q[1];
rz(0.55603212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4512999) q[0];
sx q[0];
rz(-2.8273819) q[0];
sx q[0];
rz(-0.67361535) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37555917) q[2];
sx q[2];
rz(-1.9347768) q[2];
sx q[2];
rz(-3.08422) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1631524) q[1];
sx q[1];
rz(-1.5500871) q[1];
sx q[1];
rz(0.15343517) q[1];
rz(-pi) q[2];
rz(-2.0222072) q[3];
sx q[3];
rz(-1.4638293) q[3];
sx q[3];
rz(-0.2052923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14744645) q[2];
sx q[2];
rz(-2.3799956) q[2];
sx q[2];
rz(0.93696326) q[2];
rz(0.41038904) q[3];
sx q[3];
rz(-0.97987163) q[3];
sx q[3];
rz(0.65569896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3963102) q[0];
sx q[0];
rz(-0.87438011) q[0];
sx q[0];
rz(-2.1790047) q[0];
rz(-2.4707322) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(2.3765391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6214692) q[0];
sx q[0];
rz(-1.4111821) q[0];
sx q[0];
rz(-1.9748715) q[0];
rz(-pi) q[1];
rz(0.069326055) q[2];
sx q[2];
rz(-1.6717807) q[2];
sx q[2];
rz(0.19420964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9468967) q[1];
sx q[1];
rz(-1.0651154) q[1];
sx q[1];
rz(-0.45241385) q[1];
x q[2];
rz(1.8840065) q[3];
sx q[3];
rz(-2.7566285) q[3];
sx q[3];
rz(2.0244618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.700231) q[2];
sx q[2];
rz(-0.92741489) q[2];
sx q[2];
rz(2.3564763) q[2];
rz(0.49078068) q[3];
sx q[3];
rz(-1.9793341) q[3];
sx q[3];
rz(-0.47500113) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9908776) q[0];
sx q[0];
rz(-0.10902037) q[0];
sx q[0];
rz(-2.9955067) q[0];
rz(2.8699005) q[1];
sx q[1];
rz(-2.3482595) q[1];
sx q[1];
rz(-0.21154107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95446996) q[0];
sx q[0];
rz(-0.7899219) q[0];
sx q[0];
rz(-0.23440897) q[0];
rz(1.2810938) q[2];
sx q[2];
rz(-2.7549676) q[2];
sx q[2];
rz(0.36687106) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22979888) q[1];
sx q[1];
rz(-1.8066563) q[1];
sx q[1];
rz(1.6423507) q[1];
rz(-3.0307496) q[3];
sx q[3];
rz(-2.1680084) q[3];
sx q[3];
rz(0.66652121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3286256) q[2];
sx q[2];
rz(-0.93572891) q[2];
sx q[2];
rz(-0.82175559) q[2];
rz(-1.8030608) q[3];
sx q[3];
rz(-1.0889564) q[3];
sx q[3];
rz(-0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61447918) q[0];
sx q[0];
rz(-2.4452657) q[0];
sx q[0];
rz(1.3592199) q[0];
rz(2.1837153) q[1];
sx q[1];
rz(-0.33708894) q[1];
sx q[1];
rz(0.27761308) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.255186) q[0];
sx q[0];
rz(-1.0346812) q[0];
sx q[0];
rz(-0.84250959) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2587373) q[2];
sx q[2];
rz(-0.81969374) q[2];
sx q[2];
rz(1.8783542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.902144) q[1];
sx q[1];
rz(-2.7056597) q[1];
sx q[1];
rz(2.8365342) q[1];
rz(-pi) q[2];
rz(1.7429535) q[3];
sx q[3];
rz(-1.1973698) q[3];
sx q[3];
rz(0.40478727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17310384) q[2];
sx q[2];
rz(-2.5471881) q[2];
sx q[2];
rz(-1.8572726) q[2];
rz(-2.3641018) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(-0.29977453) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1182627) q[0];
sx q[0];
rz(-1.6602004) q[0];
sx q[0];
rz(-0.0059286038) q[0];
rz(1.3395576) q[1];
sx q[1];
rz(-1.0362933) q[1];
sx q[1];
rz(-0.48066995) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5222134) q[0];
sx q[0];
rz(-1.5850768) q[0];
sx q[0];
rz(-0.016187384) q[0];
rz(-pi) q[1];
rz(-0.92222442) q[2];
sx q[2];
rz(-1.2184452) q[2];
sx q[2];
rz(2.1395526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25998589) q[1];
sx q[1];
rz(-0.98129683) q[1];
sx q[1];
rz(2.0719253) q[1];
rz(-pi) q[2];
rz(-3.0044286) q[3];
sx q[3];
rz(-1.3068294) q[3];
sx q[3];
rz(-0.19318737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19114755) q[2];
sx q[2];
rz(-2.5823249) q[2];
sx q[2];
rz(-0.17678235) q[2];
rz(2.2117129) q[3];
sx q[3];
rz(-1.9527438) q[3];
sx q[3];
rz(-1.0298347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4491691) q[0];
sx q[0];
rz(-1.7208736) q[0];
sx q[0];
rz(2.0509913) q[0];
rz(-1.0724267) q[1];
sx q[1];
rz(-1.6630395) q[1];
sx q[1];
rz(1.6987775) q[1];
rz(-0.49184797) q[2];
sx q[2];
rz(-1.576423) q[2];
sx q[2];
rz(1.5058422) q[2];
rz(-3.0092893) q[3];
sx q[3];
rz(-2.5322661) q[3];
sx q[3];
rz(0.54064396) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
