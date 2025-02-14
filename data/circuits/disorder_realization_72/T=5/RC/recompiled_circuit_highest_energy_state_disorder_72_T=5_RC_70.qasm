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
rz(0.25207818) q[0];
sx q[0];
rz(-2.5660388) q[0];
sx q[0];
rz(-0.12230305) q[0];
rz(1.1072493) q[1];
sx q[1];
rz(-0.9613494) q[1];
sx q[1];
rz(-0.038318757) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3275546) q[0];
sx q[0];
rz(-0.71056847) q[0];
sx q[0];
rz(-0.98301218) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11876583) q[2];
sx q[2];
rz(-1.6890172) q[2];
sx q[2];
rz(1.2615276) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9895175) q[1];
sx q[1];
rz(-0.65836009) q[1];
sx q[1];
rz(2.1170298) q[1];
x q[2];
rz(0.24738048) q[3];
sx q[3];
rz(-1.5508964) q[3];
sx q[3];
rz(0.54075356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4965839) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(-1.4487779) q[2];
rz(0.19041666) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(-2.9922488) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038302) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(-0.25800905) q[0];
rz(-1.5500655) q[1];
sx q[1];
rz(-1.240088) q[1];
sx q[1];
rz(0.16664997) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.907703) q[0];
sx q[0];
rz(-0.59032102) q[0];
sx q[0];
rz(0.17206828) q[0];
rz(-pi) q[1];
rz(1.2857976) q[2];
sx q[2];
rz(-1.9420331) q[2];
sx q[2];
rz(1.2165537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0936792) q[1];
sx q[1];
rz(-2.4166345) q[1];
sx q[1];
rz(-0.18904257) q[1];
rz(0.80662585) q[3];
sx q[3];
rz(-2.4203033) q[3];
sx q[3];
rz(-1.2794577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0698645) q[2];
sx q[2];
rz(-2.0191777) q[2];
sx q[2];
rz(-1.2321164) q[2];
rz(-0.92464906) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(-2.0360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6919747) q[0];
sx q[0];
rz(-1.469935) q[0];
sx q[0];
rz(3.1396507) q[0];
rz(-3.107403) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(-1.5431822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0372678) q[0];
sx q[0];
rz(-2.0195978) q[0];
sx q[0];
rz(-1.857313) q[0];
rz(-pi) q[1];
rz(0.73590455) q[2];
sx q[2];
rz(-1.8712448) q[2];
sx q[2];
rz(1.4020593) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.98108058) q[1];
sx q[1];
rz(-2.7857669) q[1];
sx q[1];
rz(-3.0291104) q[1];
x q[2];
rz(-0.24418025) q[3];
sx q[3];
rz(-0.93726087) q[3];
sx q[3];
rz(1.4118495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89692846) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-2.1042692) q[2];
rz(-1.1050998) q[3];
sx q[3];
rz(-1.6048071) q[3];
sx q[3];
rz(-1.1387811) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8741375) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(1.4111891) q[0];
rz(0.63938582) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(-3.0944518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1829202) q[0];
sx q[0];
rz(-0.65802461) q[0];
sx q[0];
rz(-0.47112314) q[0];
rz(-pi) q[1];
rz(0.36836715) q[2];
sx q[2];
rz(-2.6410111) q[2];
sx q[2];
rz(-0.49886242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5410536) q[1];
sx q[1];
rz(-2.4324634) q[1];
sx q[1];
rz(-2.1313271) q[1];
rz(-0.29508884) q[3];
sx q[3];
rz(-0.86936823) q[3];
sx q[3];
rz(1.6158582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.822927) q[2];
sx q[2];
rz(-1.9959799) q[2];
sx q[2];
rz(-1.9226496) q[2];
rz(-0.31560358) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(2.992673) q[3];
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
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8197935) q[0];
sx q[0];
rz(-0.55235523) q[0];
sx q[0];
rz(2.7929982) q[0];
rz(-1.8888585) q[1];
sx q[1];
rz(-1.9786973) q[1];
sx q[1];
rz(-1.3016275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4322223) q[0];
sx q[0];
rz(-0.8101058) q[0];
sx q[0];
rz(-0.10929426) q[0];
x q[1];
rz(-2.581004) q[2];
sx q[2];
rz(-1.9288256) q[2];
sx q[2];
rz(-2.5303417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.0685454) q[1];
sx q[1];
rz(-1.7918192) q[1];
sx q[1];
rz(-2.9566492) q[1];
rz(-1.0200649) q[3];
sx q[3];
rz(-1.8740843) q[3];
sx q[3];
rz(-1.0098977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18315135) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(1.4754254) q[2];
rz(-0.685855) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1118065) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(-1.4923805) q[0];
rz(-0.97914186) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(-1.917256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5021073) q[0];
sx q[0];
rz(-1.7687135) q[0];
sx q[0];
rz(-1.5396126) q[0];
x q[1];
rz(-0.54927214) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(0.21085462) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.853272) q[1];
sx q[1];
rz(-1.5832807) q[1];
sx q[1];
rz(-0.95697831) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94131366) q[3];
sx q[3];
rz(-1.286003) q[3];
sx q[3];
rz(2.3465057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41902038) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(2.2255619) q[2];
rz(-1.3845059) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(-0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0093805669) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(0.95426553) q[0];
rz(-2.7047899) q[1];
sx q[1];
rz(-1.5780459) q[1];
sx q[1];
rz(-3.0314235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1230135) q[0];
sx q[0];
rz(-1.8948104) q[0];
sx q[0];
rz(-1.3186245) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1071579) q[2];
sx q[2];
rz(-2.1948994) q[2];
sx q[2];
rz(2.1339061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4136849) q[1];
sx q[1];
rz(-1.2251108) q[1];
sx q[1];
rz(-3.137008) q[1];
x q[2];
rz(-2.0131936) q[3];
sx q[3];
rz(-2.5964649) q[3];
sx q[3];
rz(-2.9233208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7097912) q[2];
sx q[2];
rz(-0.0063889901) q[2];
sx q[2];
rz(0.94376454) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540045) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(0.34550825) q[0];
rz(2.0461138) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(1.5798205) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6248572) q[0];
sx q[0];
rz(-0.78291946) q[0];
sx q[0];
rz(-1.6153512) q[0];
rz(-pi) q[1];
rz(-1.4905632) q[2];
sx q[2];
rz(-2.4828665) q[2];
sx q[2];
rz(1.6541437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4862389) q[1];
sx q[1];
rz(-1.9353239) q[1];
sx q[1];
rz(0.93185394) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6027661) q[3];
sx q[3];
rz(-1.073146) q[3];
sx q[3];
rz(-0.15807334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3788508) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(1.4307107) q[2];
rz(2.6070969) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(-0.9526332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6119824) q[0];
sx q[0];
rz(-0.069267608) q[0];
sx q[0];
rz(1.3105422) q[0];
rz(0.88690859) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(1.2695674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7879155) q[0];
sx q[0];
rz(-1.6991109) q[0];
sx q[0];
rz(-1.2169891) q[0];
x q[1];
rz(-1.1940895) q[2];
sx q[2];
rz(-2.7393638) q[2];
sx q[2];
rz(2.5144983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0506057) q[1];
sx q[1];
rz(-2.7805353) q[1];
sx q[1];
rz(-2.2772842) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9501602) q[3];
sx q[3];
rz(-2.9342954) q[3];
sx q[3];
rz(1.8199004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5294007) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(-0.23294918) q[2];
rz(-1.285078) q[3];
sx q[3];
rz(-0.67684567) q[3];
sx q[3];
rz(2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1583629) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(-3.0443211) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(-0.21673094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9602338) q[0];
sx q[0];
rz(-0.17941071) q[0];
sx q[0];
rz(1.2563733) q[0];
rz(-pi) q[1];
rz(-1.6172269) q[2];
sx q[2];
rz(-1.7123607) q[2];
sx q[2];
rz(-0.19664779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72783732) q[1];
sx q[1];
rz(-1.0168795) q[1];
sx q[1];
rz(-2.3604849) q[1];
rz(-pi) q[2];
rz(-1.9166462) q[3];
sx q[3];
rz(-2.2035363) q[3];
sx q[3];
rz(0.94055292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8486166) q[2];
sx q[2];
rz(-0.28102195) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(-1.5589335) q[3];
sx q[3];
rz(-1.0844743) q[3];
sx q[3];
rz(-0.16286477) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611298) q[0];
sx q[0];
rz(-1.5123788) q[0];
sx q[0];
rz(2.6176591) q[0];
rz(1.0864661) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(3.1343717) q[2];
sx q[2];
rz(-1.5389969) q[2];
sx q[2];
rz(2.5015884) q[2];
rz(-0.85862715) q[3];
sx q[3];
rz(-1.5024019) q[3];
sx q[3];
rz(2.0988322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
