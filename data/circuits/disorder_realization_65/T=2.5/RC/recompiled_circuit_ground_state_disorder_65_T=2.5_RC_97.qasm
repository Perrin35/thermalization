OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.559691) q[0];
sx q[0];
rz(-0.087495916) q[0];
sx q[0];
rz(0.09905941) q[0];
rz(2.8506408) q[1];
sx q[1];
rz(-1.7506316) q[1];
sx q[1];
rz(-1.5150392) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45315021) q[0];
sx q[0];
rz(-2.6209798) q[0];
sx q[0];
rz(1.0941891) q[0];
rz(-2.6886772) q[2];
sx q[2];
rz(-2.0358815) q[2];
sx q[2];
rz(-2.1208491) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3063872) q[1];
sx q[1];
rz(-2.3162419) q[1];
sx q[1];
rz(-2.9236776) q[1];
rz(-pi) q[2];
rz(-1.8606676) q[3];
sx q[3];
rz(-1.8767281) q[3];
sx q[3];
rz(0.044222983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.060595162) q[2];
sx q[2];
rz(-0.011757714) q[2];
sx q[2];
rz(-3.0309929) q[2];
rz(-0.009875385) q[3];
sx q[3];
rz(-0.014009352) q[3];
sx q[3];
rz(2.2601155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9843543) q[0];
sx q[0];
rz(-3.053559) q[0];
sx q[0];
rz(-1.3798168) q[0];
rz(-3.1294322) q[1];
sx q[1];
rz(-0.98406839) q[1];
sx q[1];
rz(-1.5252569) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79560382) q[0];
sx q[0];
rz(-1.5143745) q[0];
sx q[0];
rz(-0.31375558) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3697769) q[2];
sx q[2];
rz(-1.9819385) q[2];
sx q[2];
rz(2.0697921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7213707) q[1];
sx q[1];
rz(-0.18850148) q[1];
sx q[1];
rz(-2.7367715) q[1];
rz(-2.8875927) q[3];
sx q[3];
rz(-0.45886654) q[3];
sx q[3];
rz(-0.34628689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.81185594) q[2];
sx q[2];
rz(-0.024710329) q[2];
sx q[2];
rz(-3.109566) q[2];
rz(-2.5305667) q[3];
sx q[3];
rz(-1.9910087) q[3];
sx q[3];
rz(-1.7648296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5830314) q[0];
sx q[0];
rz(-0.91201454) q[0];
sx q[0];
rz(-1.4645905) q[0];
rz(1.5088082) q[1];
sx q[1];
rz(-1.8784411) q[1];
sx q[1];
rz(0.060922932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446242) q[0];
sx q[0];
rz(-1.754027) q[0];
sx q[0];
rz(-0.38559648) q[0];
x q[1];
rz(0.10468049) q[2];
sx q[2];
rz(-1.4587683) q[2];
sx q[2];
rz(-2.688863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12184252) q[1];
sx q[1];
rz(-1.7438683) q[1];
sx q[1];
rz(-2.9815841) q[1];
rz(2.4931566) q[3];
sx q[3];
rz(-2.6649974) q[3];
sx q[3];
rz(2.4364382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1724725) q[2];
sx q[2];
rz(-3.1388333) q[2];
sx q[2];
rz(2.3976682) q[2];
rz(-0.38532358) q[3];
sx q[3];
rz(-0.0016366882) q[3];
sx q[3];
rz(0.90414944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1291523) q[0];
sx q[0];
rz(-3.0991982) q[0];
sx q[0];
rz(0.84268919) q[0];
rz(-0.95376247) q[1];
sx q[1];
rz(-1.5019491) q[1];
sx q[1];
rz(-1.5488254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212397) q[0];
sx q[0];
rz(-1.5505393) q[0];
sx q[0];
rz(-3.0244104) q[0];
x q[1];
rz(1.2090525) q[2];
sx q[2];
rz(-2.7186223) q[2];
sx q[2];
rz(-1.376938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4348244) q[1];
sx q[1];
rz(-0.93509669) q[1];
sx q[1];
rz(0.058556538) q[1];
rz(-pi) q[2];
rz(-2.3609045) q[3];
sx q[3];
rz(-1.4997291) q[3];
sx q[3];
rz(-0.64551693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.039975) q[2];
sx q[2];
rz(-0.016923252) q[2];
sx q[2];
rz(0.42538154) q[2];
rz(1.7915122) q[3];
sx q[3];
rz(-1.5451508) q[3];
sx q[3];
rz(0.30459705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4199583) q[0];
sx q[0];
rz(-3.1379134) q[0];
sx q[0];
rz(2.4256308) q[0];
rz(1.6250027) q[1];
sx q[1];
rz(-0.45418987) q[1];
sx q[1];
rz(-2.9862278) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6139406) q[0];
sx q[0];
rz(-0.10738534) q[0];
sx q[0];
rz(1.7201109) q[0];
x q[1];
rz(2.9174439) q[2];
sx q[2];
rz(-2.756065) q[2];
sx q[2];
rz(1.4003426) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1137312) q[1];
sx q[1];
rz(-1.3044954) q[1];
sx q[1];
rz(-1.4425502) q[1];
x q[2];
rz(2.7097059) q[3];
sx q[3];
rz(-2.1968522) q[3];
sx q[3];
rz(1.1736479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4510437) q[2];
sx q[2];
rz(-1.1210219) q[2];
sx q[2];
rz(-1.6070018) q[2];
rz(-2.0524041) q[3];
sx q[3];
rz(-3.1179699) q[3];
sx q[3];
rz(1.3634118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71517336) q[0];
sx q[0];
rz(-0.024113163) q[0];
sx q[0];
rz(2.5809848) q[0];
rz(0.86588612) q[1];
sx q[1];
rz(-3.1325603) q[1];
sx q[1];
rz(-0.83585709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27640057) q[0];
sx q[0];
rz(-1.5527343) q[0];
sx q[0];
rz(-2.8457321) q[0];
rz(1.4908815) q[2];
sx q[2];
rz(-2.7330988) q[2];
sx q[2];
rz(-1.6871883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33619704) q[1];
sx q[1];
rz(-1.8289036) q[1];
sx q[1];
rz(-1.4610402) q[1];
x q[2];
rz(1.3243598) q[3];
sx q[3];
rz(-1.8772238) q[3];
sx q[3];
rz(0.6014733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5777638) q[2];
sx q[2];
rz(-1.8356297) q[2];
sx q[2];
rz(1.5888265) q[2];
rz(-2.0896437) q[3];
sx q[3];
rz(-2.6237539) q[3];
sx q[3];
rz(-0.53798211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2688667) q[0];
sx q[0];
rz(-2.5476542) q[0];
sx q[0];
rz(-1.3237704) q[0];
rz(-0.84953228) q[1];
sx q[1];
rz(-3.140026) q[1];
sx q[1];
rz(0.82226396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13110833) q[0];
sx q[0];
rz(-2.9379651) q[0];
sx q[0];
rz(-1.078461) q[0];
x q[1];
rz(-1.6258114) q[2];
sx q[2];
rz(-1.233068) q[2];
sx q[2];
rz(-1.6134451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7436372) q[1];
sx q[1];
rz(-1.5642867) q[1];
sx q[1];
rz(-1.5012451) q[1];
rz(-pi) q[2];
rz(-1.7674462) q[3];
sx q[3];
rz(-1.3101642) q[3];
sx q[3];
rz(2.9109241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5620586) q[2];
sx q[2];
rz(-1.2349393) q[2];
sx q[2];
rz(2.298992) q[2];
rz(-2.5990727) q[3];
sx q[3];
rz(-0.018535651) q[3];
sx q[3];
rz(2.5491469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.056996718) q[0];
sx q[0];
rz(-1.5942986) q[0];
sx q[0];
rz(-1.0513167) q[0];
rz(1.7683138) q[1];
sx q[1];
rz(-3.1128502) q[1];
sx q[1];
rz(-1.9017259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9997885) q[0];
sx q[0];
rz(-1.3281137) q[0];
sx q[0];
rz(-1.61116) q[0];
rz(1.6592024) q[2];
sx q[2];
rz(-2.7475221) q[2];
sx q[2];
rz(-1.566377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6426601) q[1];
sx q[1];
rz(-0.98915425) q[1];
sx q[1];
rz(-10/(13*pi)) q[1];
rz(1.8375592) q[3];
sx q[3];
rz(-1.752231) q[3];
sx q[3];
rz(2.415588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54590589) q[2];
sx q[2];
rz(-3.1133856) q[2];
sx q[2];
rz(-0.14740454) q[2];
rz(1.6029415) q[3];
sx q[3];
rz(-1.4378865) q[3];
sx q[3];
rz(2.9586207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1786757) q[0];
sx q[0];
rz(-0.33483949) q[0];
sx q[0];
rz(1.8033002) q[0];
rz(1.4778522) q[1];
sx q[1];
rz(-1.1344974) q[1];
sx q[1];
rz(2.9683364) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4497183) q[0];
sx q[0];
rz(-2.1341748) q[0];
sx q[0];
rz(2.0823576) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7376027) q[2];
sx q[2];
rz(-2.1232018) q[2];
sx q[2];
rz(-0.45726038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45341982) q[1];
sx q[1];
rz(-1.416561) q[1];
sx q[1];
rz(-2.8919528) q[1];
rz(-1.5573172) q[3];
sx q[3];
rz(-1.5766249) q[3];
sx q[3];
rz(-1.7268501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.133539) q[2];
sx q[2];
rz(-0.011913813) q[2];
sx q[2];
rz(-1.0201721) q[2];
rz(-0.080848761) q[3];
sx q[3];
rz(-3.1404218) q[3];
sx q[3];
rz(-1.1297869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5197649) q[0];
sx q[0];
rz(-0.28525678) q[0];
sx q[0];
rz(-1.4813625) q[0];
rz(0.18352428) q[1];
sx q[1];
rz(-0.32569519) q[1];
sx q[1];
rz(0.087800177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38890281) q[0];
sx q[0];
rz(-1.6713033) q[0];
sx q[0];
rz(-2.697562) q[0];
rz(-pi) q[1];
rz(-1.4884639) q[2];
sx q[2];
rz(-1.8762686) q[2];
sx q[2];
rz(0.20037547) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1147881) q[1];
sx q[1];
rz(-0.27301306) q[1];
sx q[1];
rz(-2.6778912) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3751288) q[3];
sx q[3];
rz(-2.4371456) q[3];
sx q[3];
rz(2.7337027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8537019) q[2];
sx q[2];
rz(-0.012306865) q[2];
sx q[2];
rz(-2.2575209) q[2];
rz(-0.70841241) q[3];
sx q[3];
rz(-3.1413779) q[3];
sx q[3];
rz(0.29402548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52083279) q[0];
sx q[0];
rz(-1.5694869) q[0];
sx q[0];
rz(-1.6318305) q[0];
rz(-0.028342551) q[1];
sx q[1];
rz(-2.5181073) q[1];
sx q[1];
rz(0.070652031) q[1];
rz(0.44755837) q[2];
sx q[2];
rz(-1.9050514) q[2];
sx q[2];
rz(-1.5290268) q[2];
rz(0.72851758) q[3];
sx q[3];
rz(-1.3195724) q[3];
sx q[3];
rz(2.6825678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
