OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(-1.985328) q[1];
sx q[1];
rz(-1.1896689) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1253163) q[0];
sx q[0];
rz(-1.2287041) q[0];
sx q[0];
rz(-2.5972511) q[0];
x q[1];
rz(2.9615133) q[2];
sx q[2];
rz(-1.4305978) q[2];
sx q[2];
rz(-2.7498498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8813821) q[1];
sx q[1];
rz(-1.26169) q[1];
sx q[1];
rz(2.9865772) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0359997) q[3];
sx q[3];
rz(-0.70123226) q[3];
sx q[3];
rz(-2.0555156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4253915) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(-0.5973967) q[2];
rz(-1.3655837) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(-1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(3.0283668) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0229867) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(-0.042516275) q[0];
rz(-3.0683238) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(-2.427223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7624224) q[1];
sx q[1];
rz(-1.1602853) q[1];
sx q[1];
rz(-1.1344086) q[1];
rz(-pi) q[2];
rz(-1.9085957) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(2.256957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.018628) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(0.078331746) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(-1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-3.0139794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1193082) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(0.22247252) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67655501) q[2];
sx q[2];
rz(-0.99202079) q[2];
sx q[2];
rz(0.81625953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.322791) q[1];
sx q[1];
rz(-1.6530419) q[1];
sx q[1];
rz(1.3311885) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25281275) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(-1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(2.7691832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5126257) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(-0.15928282) q[0];
rz(-2.0580975) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(1.0207748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47632521) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(-1.1900351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1129205) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(1.8478912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(2.2327936) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.882778) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(-1.1313261) q[0];
rz(-0.5979901) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871183) q[0];
sx q[0];
rz(-1.1312383) q[0];
sx q[0];
rz(-1.6580824) q[0];
rz(-1.3555688) q[2];
sx q[2];
rz(-2.2903246) q[2];
sx q[2];
rz(2.1446251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.685806) q[1];
sx q[1];
rz(-1.5450486) q[1];
sx q[1];
rz(-0.2548494) q[1];
rz(-0.39354126) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(-0.81641203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.1452902) q[2];
rz(-0.14906135) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(-1.0173652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155708) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(-2.3964264) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3286896) q[0];
sx q[0];
rz(-2.0314386) q[0];
sx q[0];
rz(-0.33005379) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2007347) q[2];
sx q[2];
rz(-1.0002631) q[2];
sx q[2];
rz(0.42966118) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16490368) q[1];
sx q[1];
rz(-0.36786825) q[1];
sx q[1];
rz(1.5771754) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3429787) q[3];
sx q[3];
rz(-0.67043793) q[3];
sx q[3];
rz(-0.05712856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3559945) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(-0.28665001) q[2];
rz(-2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(-0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(1.0010304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46465835) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(1.2901165) q[0];
rz(-pi) q[1];
rz(0.34044388) q[2];
sx q[2];
rz(-1.921244) q[2];
sx q[2];
rz(3.0710789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66623464) q[1];
sx q[1];
rz(-1.959414) q[1];
sx q[1];
rz(0.16814853) q[1];
rz(2.4909199) q[3];
sx q[3];
rz(-1.7400029) q[3];
sx q[3];
rz(2.8229439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.004185685) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(3.0380847) q[2];
rz(0.59182709) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(-1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1471257) q[0];
sx q[0];
rz(-1.2110706) q[0];
sx q[0];
rz(-0.61867449) q[0];
rz(1.1578324) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(-2.4177891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0720313) q[1];
sx q[1];
rz(-1.8998635) q[1];
sx q[1];
rz(-2.0272994) q[1];
rz(-pi) q[2];
rz(2.249986) q[3];
sx q[3];
rz(-1.7779751) q[3];
sx q[3];
rz(-0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8395681) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(-0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(2.8651967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6733288) q[0];
sx q[0];
rz(-1.7305166) q[0];
sx q[0];
rz(2.7295652) q[0];
rz(1.7644291) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(-2.5335238) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8564057) q[1];
sx q[1];
rz(-1.4126084) q[1];
sx q[1];
rz(-1.3473947) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60864457) q[3];
sx q[3];
rz(-1.7749655) q[3];
sx q[3];
rz(0.58231402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(-3.0140871) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(1.030285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477662) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(-1.312027) q[0];
rz(-pi) q[1];
rz(1.2376386) q[2];
sx q[2];
rz(-2.2292456) q[2];
sx q[2];
rz(1.2308987) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.13223091) q[1];
sx q[1];
rz(-1.2538326) q[1];
sx q[1];
rz(0.46225458) q[1];
x q[2];
rz(2.685931) q[3];
sx q[3];
rz(-2.4767498) q[3];
sx q[3];
rz(-1.1841707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(2.4009005) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13070233) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(1.7651991) q[2];
sx q[2];
rz(-1.9971078) q[2];
sx q[2];
rz(2.8304921) q[2];
rz(-2.1445198) q[3];
sx q[3];
rz(-2.1763747) q[3];
sx q[3];
rz(2.5381266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
