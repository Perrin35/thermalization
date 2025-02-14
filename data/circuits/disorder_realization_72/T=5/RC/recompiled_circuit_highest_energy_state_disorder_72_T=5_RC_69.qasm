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
rz(-2.8895145) q[0];
sx q[0];
rz(-0.57555389) q[0];
sx q[0];
rz(-3.0192896) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(-3.1032739) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2244683) q[0];
sx q[0];
rz(-1.2007133) q[0];
sx q[0];
rz(-0.94934772) q[0];
x q[1];
rz(2.3549809) q[2];
sx q[2];
rz(-0.16737882) q[2];
sx q[2];
rz(-0.47030631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8073544) q[1];
sx q[1];
rz(-2.1209201) q[1];
sx q[1];
rz(0.38205876) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0604912) q[3];
sx q[3];
rz(-0.24816324) q[3];
sx q[3];
rz(0.95141548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64500874) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(1.6928147) q[2];
rz(-0.19041666) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(-0.14934389) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038302) q[0];
sx q[0];
rz(-1.2685403) q[0];
sx q[0];
rz(0.25800905) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.240088) q[1];
sx q[1];
rz(0.16664997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9480707) q[0];
sx q[0];
rz(-1.4753454) q[0];
sx q[0];
rz(2.5581317) q[0];
rz(-0.62549297) q[2];
sx q[2];
rz(-2.6776367) q[2];
sx q[2];
rz(2.6044012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2981603) q[1];
sx q[1];
rz(-2.2800804) q[1];
sx q[1];
rz(-1.7357566) q[1];
x q[2];
rz(0.80662585) q[3];
sx q[3];
rz(-0.72128937) q[3];
sx q[3];
rz(-1.862135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0717281) q[2];
sx q[2];
rz(-2.0191777) q[2];
sx q[2];
rz(1.2321164) q[2];
rz(0.92464906) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(-1.1054976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.449618) q[0];
sx q[0];
rz(-1.469935) q[0];
sx q[0];
rz(0.0019419226) q[0];
rz(-0.034189668) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(1.5431822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7352075) q[0];
sx q[0];
rz(-1.8282561) q[0];
sx q[0];
rz(0.46528146) q[0];
x q[1];
rz(0.43242792) q[2];
sx q[2];
rz(-0.78410599) q[2];
sx q[2];
rz(-2.6570005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48422563) q[1];
sx q[1];
rz(-1.6099085) q[1];
sx q[1];
rz(2.7878321) q[1];
rz(-pi) q[2];
rz(0.24418025) q[3];
sx q[3];
rz(-2.2043318) q[3];
sx q[3];
rz(1.4118495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2446642) q[2];
sx q[2];
rz(-0.43411532) q[2];
sx q[2];
rz(2.1042692) q[2];
rz(-1.1050998) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(-2.0028116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.26745519) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(-1.4111891) q[0];
rz(2.5022068) q[1];
sx q[1];
rz(-1.4566028) q[1];
sx q[1];
rz(0.04714084) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706076) q[0];
sx q[0];
rz(-1.8520675) q[0];
sx q[0];
rz(2.5384643) q[0];
rz(-pi) q[1];
rz(0.36836715) q[2];
sx q[2];
rz(-0.50058156) q[2];
sx q[2];
rz(-2.6427302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47430965) q[1];
sx q[1];
rz(-1.9243) q[1];
sx q[1];
rz(-0.94236417) q[1];
x q[2];
rz(2.8465038) q[3];
sx q[3];
rz(-2.2722244) q[3];
sx q[3];
rz(-1.6158582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.822927) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(-1.9226496) q[2];
rz(0.31560358) q[3];
sx q[3];
rz(-0.054840755) q[3];
sx q[3];
rz(-0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3217992) q[0];
sx q[0];
rz(-0.55235523) q[0];
sx q[0];
rz(2.7929982) q[0];
rz(-1.8888585) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(-1.8399651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5515297) q[0];
sx q[0];
rz(-2.3746535) q[0];
sx q[0];
rz(-1.4566896) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61321494) q[2];
sx q[2];
rz(-2.4869031) q[2];
sx q[2];
rz(1.6729205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5036936) q[1];
sx q[1];
rz(-0.28721962) q[1];
sx q[1];
rz(0.88493185) q[1];
rz(-2.1215277) q[3];
sx q[3];
rz(-1.2675084) q[3];
sx q[3];
rz(2.1316949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18315135) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(-1.6661673) q[2];
rz(2.4557377) q[3];
sx q[3];
rz(-0.97969046) q[3];
sx q[3];
rz(2.6077152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1118065) q[0];
sx q[0];
rz(-0.15203467) q[0];
sx q[0];
rz(-1.6492122) q[0];
rz(2.1624508) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(-1.917256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5021073) q[0];
sx q[0];
rz(-1.7687135) q[0];
sx q[0];
rz(-1.60198) q[0];
rz(-pi) q[1];
rz(0.54927214) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(-0.21085462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2883207) q[1];
sx q[1];
rz(-1.558312) q[1];
sx q[1];
rz(0.95697831) q[1];
rz(-1.1093418) q[3];
sx q[3];
rz(-0.68285817) q[3];
sx q[3];
rz(2.7340555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41902038) q[2];
sx q[2];
rz(-1.4973065) q[2];
sx q[2];
rz(2.2255619) q[2];
rz(-1.7570868) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0093805669) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(0.95426553) q[0];
rz(-0.43680278) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(-3.0314235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6340652) q[0];
sx q[0];
rz(-1.332009) q[0];
sx q[0];
rz(0.33383835) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4442441) q[2];
sx q[2];
rz(-1.9983872) q[2];
sx q[2];
rz(-0.22874895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72790775) q[1];
sx q[1];
rz(-1.9164819) q[1];
sx q[1];
rz(0.0045846049) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8875868) q[3];
sx q[3];
rz(-2.0585103) q[3];
sx q[3];
rz(0.28764492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43180141) q[2];
sx q[2];
rz(-0.0063889901) q[2];
sx q[2];
rz(0.94376454) q[2];
rz(0.56378311) q[3];
sx q[3];
rz(-2.0457334) q[3];
sx q[3];
rz(2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540045) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(0.34550825) q[0];
rz(1.0954789) q[1];
sx q[1];
rz(-1.6637207) q[1];
sx q[1];
rz(-1.5617721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673543) q[0];
sx q[0];
rz(-2.3586732) q[0];
sx q[0];
rz(-1.5262414) q[0];
x q[1];
rz(-3.079633) q[2];
sx q[2];
rz(-2.2270348) q[2];
sx q[2];
rz(-1.7554754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17434828) q[1];
sx q[1];
rz(-2.1617608) q[1];
sx q[1];
rz(-2.6978542) q[1];
rz(-pi) q[2];
rz(0.058770725) q[3];
sx q[3];
rz(-2.6430025) q[3];
sx q[3];
rz(2.9166247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76274189) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(-1.7108819) q[2];
rz(2.6070969) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(2.1889595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6119824) q[0];
sx q[0];
rz(-0.069267608) q[0];
sx q[0];
rz(-1.8310504) q[0];
rz(0.88690859) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(1.2695674) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16989141) q[0];
sx q[0];
rz(-1.2200238) q[0];
sx q[0];
rz(-3.0049075) q[0];
x q[1];
rz(-1.1940895) q[2];
sx q[2];
rz(-0.40222886) q[2];
sx q[2];
rz(-2.5144983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0506057) q[1];
sx q[1];
rz(-2.7805353) q[1];
sx q[1];
rz(-0.86430849) q[1];
x q[2];
rz(1.6107913) q[3];
sx q[3];
rz(-1.3673395) q[3];
sx q[3];
rz(2.0154161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.612192) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(-2.9086435) q[2];
rz(1.8565146) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(0.3860093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9832298) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(-0.097271517) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-0.709788) q[1];
sx q[1];
rz(0.21673094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50057208) q[0];
sx q[0];
rz(-1.4002698) q[0];
sx q[0];
rz(0.056030355) q[0];
x q[1];
rz(0.14171504) q[2];
sx q[2];
rz(-1.6167621) q[2];
sx q[2];
rz(-1.7739997) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4137553) q[1];
sx q[1];
rz(-1.0168795) q[1];
sx q[1];
rz(-2.3604849) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4794934) q[3];
sx q[3];
rz(-1.2939014) q[3];
sx q[3];
rz(0.84018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29297605) q[2];
sx q[2];
rz(-0.28102195) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(1.5589335) q[3];
sx q[3];
rz(-1.0844743) q[3];
sx q[3];
rz(-2.9787279) q[3];
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
rz(pi/2) q[0];
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
rz(-2.0551266) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(0.0072210014) q[2];
sx q[2];
rz(-1.6025958) q[2];
sx q[2];
rz(-0.64000426) q[2];
rz(-0.090251017) q[3];
sx q[3];
rz(-0.8606438) q[3];
sx q[3];
rz(0.46910486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
