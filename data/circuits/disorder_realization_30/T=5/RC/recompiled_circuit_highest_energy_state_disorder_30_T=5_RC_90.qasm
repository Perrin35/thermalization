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
rz(0.30000609) q[0];
sx q[0];
rz(-1.0862779) q[0];
sx q[0];
rz(2.3334184) q[0];
rz(-1.1440682) q[1];
sx q[1];
rz(-0.77163982) q[1];
sx q[1];
rz(-1.5789403) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15641744) q[0];
sx q[0];
rz(-2.1083207) q[0];
sx q[0];
rz(2.318105) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8380661) q[2];
sx q[2];
rz(-1.4216982) q[2];
sx q[2];
rz(1.2318512) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7464802) q[1];
sx q[1];
rz(-1.0653138) q[1];
sx q[1];
rz(-3.014747) q[1];
rz(-1.5366301) q[3];
sx q[3];
rz(-1.5251659) q[3];
sx q[3];
rz(2.6229317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19848862) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(-1.0342106) q[2];
rz(0.69283038) q[3];
sx q[3];
rz(-1.0537909) q[3];
sx q[3];
rz(1.0606162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3028054) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(-2.5651108) q[0];
rz(3.1209962) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(1.0284665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45555025) q[0];
sx q[0];
rz(-0.33838135) q[0];
sx q[0];
rz(-0.69940059) q[0];
rz(-0.89056252) q[2];
sx q[2];
rz(-0.46698585) q[2];
sx q[2];
rz(-0.61636965) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3334782) q[1];
sx q[1];
rz(-1.9030182) q[1];
sx q[1];
rz(1.1444886) q[1];
x q[2];
rz(-1.1787492) q[3];
sx q[3];
rz(-0.96091753) q[3];
sx q[3];
rz(1.4489685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91745201) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(1.8769598) q[2];
rz(0.97359109) q[3];
sx q[3];
rz(-1.4507989) q[3];
sx q[3];
rz(0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50315404) q[0];
sx q[0];
rz(-1.8373024) q[0];
sx q[0];
rz(-0.85357443) q[0];
rz(-2.0897934) q[1];
sx q[1];
rz(-0.97426668) q[1];
sx q[1];
rz(3.1265756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5334458) q[0];
sx q[0];
rz(-2.2813873) q[0];
sx q[0];
rz(3.1300302) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21785801) q[2];
sx q[2];
rz(-0.098420489) q[2];
sx q[2];
rz(2.276536) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89755631) q[1];
sx q[1];
rz(-1.6445531) q[1];
sx q[1];
rz(-2.4270127) q[1];
x q[2];
rz(-1.5068568) q[3];
sx q[3];
rz(-1.2815164) q[3];
sx q[3];
rz(2.6562986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1198279) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(-2.8601698) q[2];
rz(-1.4540539) q[3];
sx q[3];
rz(-1.2097996) q[3];
sx q[3];
rz(0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5948831) q[0];
sx q[0];
rz(-1.2458206) q[0];
sx q[0];
rz(-2.967714) q[0];
rz(1.3315382) q[1];
sx q[1];
rz(-1.4215819) q[1];
sx q[1];
rz(1.4283659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4177307) q[0];
sx q[0];
rz(-1.3784467) q[0];
sx q[0];
rz(0.014561166) q[0];
rz(-pi) q[1];
rz(-0.32555737) q[2];
sx q[2];
rz(-0.37806219) q[2];
sx q[2];
rz(-0.69272536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54879872) q[1];
sx q[1];
rz(-0.6585291) q[1];
sx q[1];
rz(-3.0786773) q[1];
rz(-pi) q[2];
rz(-1.9589728) q[3];
sx q[3];
rz(-1.2283162) q[3];
sx q[3];
rz(-1.6714718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4398769) q[2];
sx q[2];
rz(-2.3150847) q[2];
sx q[2];
rz(-0.41713777) q[2];
rz(0.16962984) q[3];
sx q[3];
rz(-3.0483584) q[3];
sx q[3];
rz(0.72181845) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193264) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(-2.8322423) q[0];
rz(0.7695235) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(1.5187029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5734539) q[0];
sx q[0];
rz(-0.72069695) q[0];
sx q[0];
rz(1.1842404) q[0];
rz(-pi) q[1];
rz(2.3974472) q[2];
sx q[2];
rz(-2.2501906) q[2];
sx q[2];
rz(1.8216004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9936693) q[1];
sx q[1];
rz(-1.6003612) q[1];
sx q[1];
rz(1.3815341) q[1];
rz(-pi) q[2];
rz(-0.098747323) q[3];
sx q[3];
rz(-1.4374466) q[3];
sx q[3];
rz(2.6790706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7421444) q[2];
sx q[2];
rz(-2.1592906) q[2];
sx q[2];
rz(2.1461416) q[2];
rz(-2.4230912) q[3];
sx q[3];
rz(-1.3281497) q[3];
sx q[3];
rz(-2.3464581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8566078) q[0];
sx q[0];
rz(-1.6866848) q[0];
sx q[0];
rz(-1.5307776) q[0];
rz(1.4467622) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(-1.4379427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79057825) q[0];
sx q[0];
rz(-1.5503128) q[0];
sx q[0];
rz(-1.162942) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1968568) q[2];
sx q[2];
rz(-2.1293921) q[2];
sx q[2];
rz(0.93511673) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29026689) q[1];
sx q[1];
rz(-1.0630634) q[1];
sx q[1];
rz(1.4749871) q[1];
x q[2];
rz(-1.5737278) q[3];
sx q[3];
rz(-2.3303707) q[3];
sx q[3];
rz(-2.7377759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1630254) q[2];
sx q[2];
rz(-1.359553) q[2];
sx q[2];
rz(1.2320409) q[2];
rz(1.9949404) q[3];
sx q[3];
rz(-1.3403284) q[3];
sx q[3];
rz(2.8310217) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6335886) q[0];
sx q[0];
rz(-2.3242943) q[0];
sx q[0];
rz(2.001413) q[0];
rz(-2.2445402) q[1];
sx q[1];
rz(-1.8042118) q[1];
sx q[1];
rz(0.030489771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94808364) q[0];
sx q[0];
rz(-1.7917669) q[0];
sx q[0];
rz(-3.0224209) q[0];
rz(1.8980128) q[2];
sx q[2];
rz(-0.95409617) q[2];
sx q[2];
rz(1.9685352) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5329689) q[1];
sx q[1];
rz(-0.86809671) q[1];
sx q[1];
rz(2.7547902) q[1];
rz(-pi) q[2];
rz(-2.0694437) q[3];
sx q[3];
rz(-0.79550084) q[3];
sx q[3];
rz(1.3971869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9474779) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(-1.6537559) q[3];
sx q[3];
rz(-2.1066809) q[3];
sx q[3];
rz(-1.3207159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3461935) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(-1.3854223) q[0];
rz(1.2162544) q[1];
sx q[1];
rz(-1.4477891) q[1];
sx q[1];
rz(0.48929712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3922244) q[0];
sx q[0];
rz(-2.4697692) q[0];
sx q[0];
rz(-3.1098614) q[0];
rz(-0.036542372) q[2];
sx q[2];
rz(-1.9068235) q[2];
sx q[2];
rz(-0.69118369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.341574) q[1];
sx q[1];
rz(-1.4765146) q[1];
sx q[1];
rz(-1.0977919) q[1];
x q[2];
rz(3.0396949) q[3];
sx q[3];
rz(-1.8188261) q[3];
sx q[3];
rz(2.0534648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3936002) q[2];
sx q[2];
rz(-1.9990498) q[2];
sx q[2];
rz(0.46621123) q[2];
rz(1.5318058) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(0.32061779) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521249) q[0];
sx q[0];
rz(-2.2418699) q[0];
sx q[0];
rz(3.0925282) q[0];
rz(-2.9009254) q[1];
sx q[1];
rz(-2.1235695) q[1];
sx q[1];
rz(3.1191471) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2769299) q[0];
sx q[0];
rz(-2.1500181) q[0];
sx q[0];
rz(3.1129709) q[0];
x q[1];
rz(0.38752611) q[2];
sx q[2];
rz(-1.8948613) q[2];
sx q[2];
rz(1.9112196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7514403) q[1];
sx q[1];
rz(-2.1933687) q[1];
sx q[1];
rz(-2.1639813) q[1];
x q[2];
rz(-2.0231699) q[3];
sx q[3];
rz(-2.2045772) q[3];
sx q[3];
rz(-2.3364802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(-0.9683041) q[2];
rz(0.83640313) q[3];
sx q[3];
rz(-0.29763779) q[3];
sx q[3];
rz(-1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7264929) q[0];
sx q[0];
rz(-1.8510171) q[0];
sx q[0];
rz(2.7628164) q[0];
rz(-3.1355766) q[1];
sx q[1];
rz(-2.6117987) q[1];
sx q[1];
rz(-2.5078497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260035) q[0];
sx q[0];
rz(-2.3657614) q[0];
sx q[0];
rz(2.6494527) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1442398) q[2];
sx q[2];
rz(-1.7055943) q[2];
sx q[2];
rz(-2.1950795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4937456) q[1];
sx q[1];
rz(-2.4253107) q[1];
sx q[1];
rz(-2.358308) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2692436) q[3];
sx q[3];
rz(-0.81065882) q[3];
sx q[3];
rz(-1.8332421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4358431) q[2];
sx q[2];
rz(-1.5219995) q[2];
sx q[2];
rz(-0.61326927) q[2];
rz(0.67522007) q[3];
sx q[3];
rz(-0.37874159) q[3];
sx q[3];
rz(0.81775445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262065) q[0];
sx q[0];
rz(-1.362726) q[0];
sx q[0];
rz(1.2363731) q[0];
rz(0.18648237) q[1];
sx q[1];
rz(-1.3874556) q[1];
sx q[1];
rz(1.5127771) q[1];
rz(-0.040493852) q[2];
sx q[2];
rz(-0.45999415) q[2];
sx q[2];
rz(-2.5754089) q[2];
rz(3.0499115) q[3];
sx q[3];
rz(-0.95060075) q[3];
sx q[3];
rz(-1.4178249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
