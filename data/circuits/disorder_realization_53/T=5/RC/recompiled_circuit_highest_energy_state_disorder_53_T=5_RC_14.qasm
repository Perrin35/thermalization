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
rz(-0.1102912) q[0];
sx q[0];
rz(-1.2854486) q[0];
sx q[0];
rz(-0.31874803) q[0];
rz(1.1822074) q[1];
sx q[1];
rz(-1.6637586) q[1];
sx q[1];
rz(-1.1629265) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6642374) q[0];
sx q[0];
rz(-1.2384982) q[0];
sx q[0];
rz(2.8927781) q[0];
rz(-3.0672795) q[2];
sx q[2];
rz(-0.19467672) q[2];
sx q[2];
rz(3.1271324) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.27027495) q[1];
sx q[1];
rz(-1.8757038) q[1];
sx q[1];
rz(2.8256846) q[1];
rz(0.031008677) q[3];
sx q[3];
rz(-2.5664799) q[3];
sx q[3];
rz(1.2367347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9691539) q[2];
sx q[2];
rz(-2.6814851) q[2];
sx q[2];
rz(-3.0058506) q[2];
rz(-2.8067449) q[3];
sx q[3];
rz(-1.2232774) q[3];
sx q[3];
rz(-0.39723435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.6927476) q[0];
sx q[0];
rz(-2.1362342) q[0];
sx q[0];
rz(-0.94648615) q[0];
rz(-1.954156) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(0.28396398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341742) q[0];
sx q[0];
rz(-1.7611146) q[0];
sx q[0];
rz(-1.756041) q[0];
x q[1];
rz(-2.2694964) q[2];
sx q[2];
rz(-1.7533894) q[2];
sx q[2];
rz(0.77699772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2614132) q[1];
sx q[1];
rz(-1.993426) q[1];
sx q[1];
rz(1.4201866) q[1];
x q[2];
rz(-1.2228187) q[3];
sx q[3];
rz(-1.2038411) q[3];
sx q[3];
rz(1.3104132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9339319) q[2];
sx q[2];
rz(-1.2958823) q[2];
sx q[2];
rz(-0.80061039) q[2];
rz(-1.3019531) q[3];
sx q[3];
rz(-2.0079565) q[3];
sx q[3];
rz(3.083526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.12450739) q[0];
sx q[0];
rz(-2.7094816) q[0];
sx q[0];
rz(0.79297638) q[0];
rz(-2.4555581) q[1];
sx q[1];
rz(-1.425309) q[1];
sx q[1];
rz(-2.3763903) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8898123) q[0];
sx q[0];
rz(-1.7818101) q[0];
sx q[0];
rz(3.0821256) q[0];
rz(-pi) q[1];
rz(3.0900218) q[2];
sx q[2];
rz(-1.5386536) q[2];
sx q[2];
rz(2.3598537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6531892) q[1];
sx q[1];
rz(-2.326366) q[1];
sx q[1];
rz(1.1865739) q[1];
rz(-pi) q[2];
rz(-0.92242494) q[3];
sx q[3];
rz(-2.2545345) q[3];
sx q[3];
rz(-0.91372638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56696314) q[2];
sx q[2];
rz(-2.4506863) q[2];
sx q[2];
rz(-1.8951269) q[2];
rz(1.1860819) q[3];
sx q[3];
rz(-1.3776774) q[3];
sx q[3];
rz(0.20868364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935811) q[0];
sx q[0];
rz(-1.3205386) q[0];
sx q[0];
rz(-3.0082974) q[0];
rz(-0.26946274) q[1];
sx q[1];
rz(-0.91454426) q[1];
sx q[1];
rz(-2.6752313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15233697) q[0];
sx q[0];
rz(-1.4513124) q[0];
sx q[0];
rz(-2.5463922) q[0];
rz(2.9507841) q[2];
sx q[2];
rz(-0.95170232) q[2];
sx q[2];
rz(1.9931249) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.129648) q[1];
sx q[1];
rz(-0.62800558) q[1];
sx q[1];
rz(2.8824214) q[1];
rz(-pi) q[2];
rz(1.6367988) q[3];
sx q[3];
rz(-2.0836692) q[3];
sx q[3];
rz(-2.6309225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93126297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(0.24410625) q[2];
rz(-1.7804451) q[3];
sx q[3];
rz(-0.94435349) q[3];
sx q[3];
rz(1.2066427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.72531438) q[0];
sx q[0];
rz(-0.84753528) q[0];
sx q[0];
rz(-0.086294802) q[0];
rz(1.3823973) q[1];
sx q[1];
rz(-2.081213) q[1];
sx q[1];
rz(-1.8550526) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2479682) q[0];
sx q[0];
rz(-2.7148406) q[0];
sx q[0];
rz(-1.990699) q[0];
rz(-2.3628144) q[2];
sx q[2];
rz(-1.4113103) q[2];
sx q[2];
rz(0.82568141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72107347) q[1];
sx q[1];
rz(-2.2700078) q[1];
sx q[1];
rz(-3.0341604) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1819918) q[3];
sx q[3];
rz(-0.7905851) q[3];
sx q[3];
rz(0.039856002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7354108) q[2];
sx q[2];
rz(-1.7291131) q[2];
sx q[2];
rz(1.1253051) q[2];
rz(-2.3451037) q[3];
sx q[3];
rz(-0.73553604) q[3];
sx q[3];
rz(-2.9430732) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35561246) q[0];
sx q[0];
rz(-0.2934083) q[0];
sx q[0];
rz(-0.32817131) q[0];
rz(2.7525821) q[1];
sx q[1];
rz(-1.5711454) q[1];
sx q[1];
rz(1.4555812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16598116) q[0];
sx q[0];
rz(-1.360438) q[0];
sx q[0];
rz(1.4882404) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7057538) q[2];
sx q[2];
rz(-1.4452057) q[2];
sx q[2];
rz(1.3940982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38992369) q[1];
sx q[1];
rz(-1.8018541) q[1];
sx q[1];
rz(0.91737813) q[1];
rz(1.9295272) q[3];
sx q[3];
rz(-1.609612) q[3];
sx q[3];
rz(-1.2286548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97633156) q[2];
sx q[2];
rz(-2.2165522) q[2];
sx q[2];
rz(0.46249214) q[2];
rz(2.1180604) q[3];
sx q[3];
rz(-1.0547124) q[3];
sx q[3];
rz(-2.3033477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3384712) q[0];
sx q[0];
rz(-1.7110889) q[0];
sx q[0];
rz(2.9679003) q[0];
rz(-2.9202785) q[1];
sx q[1];
rz(-0.68562713) q[1];
sx q[1];
rz(1.3632704) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1376138) q[0];
sx q[0];
rz(-0.44671808) q[0];
sx q[0];
rz(-1.1398619) q[0];
rz(-2.3214746) q[2];
sx q[2];
rz(-3.0228399) q[2];
sx q[2];
rz(1.135716) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89243556) q[1];
sx q[1];
rz(-1.1022864) q[1];
sx q[1];
rz(0.69312842) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2962988) q[3];
sx q[3];
rz(-2.7850683) q[3];
sx q[3];
rz(-0.58024065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26019105) q[2];
sx q[2];
rz(-1.9804201) q[2];
sx q[2];
rz(-1.4455522) q[2];
rz(-2.3675303) q[3];
sx q[3];
rz(-0.71662199) q[3];
sx q[3];
rz(-1.6707481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.562302) q[0];
sx q[0];
rz(-0.93920541) q[0];
sx q[0];
rz(2.8939409) q[0];
rz(0.66894764) q[1];
sx q[1];
rz(-1.9525783) q[1];
sx q[1];
rz(0.98562366) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4581873) q[0];
sx q[0];
rz(-0.56945812) q[0];
sx q[0];
rz(0.21413212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47565461) q[2];
sx q[2];
rz(-0.61810571) q[2];
sx q[2];
rz(0.8408159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7880747) q[1];
sx q[1];
rz(-1.9606143) q[1];
sx q[1];
rz(-2.8105763) q[1];
rz(-pi) q[2];
x q[2];
rz(3.139601) q[3];
sx q[3];
rz(-0.59496204) q[3];
sx q[3];
rz(1.686694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0011657) q[2];
sx q[2];
rz(-2.3900034) q[2];
sx q[2];
rz(1.8482194) q[2];
rz(-1.2398531) q[3];
sx q[3];
rz(-2.445502) q[3];
sx q[3];
rz(0.70824879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86466113) q[0];
sx q[0];
rz(-1.7129352) q[0];
sx q[0];
rz(-0.96555936) q[0];
rz(-0.37636617) q[1];
sx q[1];
rz(-1.8195567) q[1];
sx q[1];
rz(1.9115062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282283) q[0];
sx q[0];
rz(-1.8364062) q[0];
sx q[0];
rz(-1.6227386) q[0];
rz(-pi) q[1];
rz(-1.240991) q[2];
sx q[2];
rz(-1.1993186) q[2];
sx q[2];
rz(0.35201752) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.230229) q[1];
sx q[1];
rz(-1.5004284) q[1];
sx q[1];
rz(0.60649782) q[1];
rz(-pi) q[2];
rz(-3.0885124) q[3];
sx q[3];
rz(-0.81392589) q[3];
sx q[3];
rz(-2.2102714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2824715) q[2];
sx q[2];
rz(-0.66994795) q[2];
sx q[2];
rz(-3.0100789) q[2];
rz(-2.6604743) q[3];
sx q[3];
rz(-0.93004623) q[3];
sx q[3];
rz(2.6432945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0216825) q[0];
sx q[0];
rz(-0.57761884) q[0];
sx q[0];
rz(2.6525894) q[0];
rz(1.8148212) q[1];
sx q[1];
rz(-1.2811456) q[1];
sx q[1];
rz(2.9723523) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6833333) q[0];
sx q[0];
rz(-1.1348327) q[0];
sx q[0];
rz(3.0100214) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45648123) q[2];
sx q[2];
rz(-0.20083961) q[2];
sx q[2];
rz(1.072509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0161017) q[1];
sx q[1];
rz(-0.91052848) q[1];
sx q[1];
rz(-1.0339526) q[1];
x q[2];
rz(-1.2792688) q[3];
sx q[3];
rz(-1.8390553) q[3];
sx q[3];
rz(2.832178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5532316) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(-2.9300743) q[2];
rz(-1.616098) q[3];
sx q[3];
rz(-2.1414521) q[3];
sx q[3];
rz(2.8741527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.78032988) q[0];
sx q[0];
rz(-1.4365256) q[0];
sx q[0];
rz(-2.8509675) q[0];
rz(2.3296539) q[1];
sx q[1];
rz(-2.5135136) q[1];
sx q[1];
rz(1.5807349) q[1];
rz(0.042150368) q[2];
sx q[2];
rz(-2.7075645) q[2];
sx q[2];
rz(-0.27070207) q[2];
rz(2.2253401) q[3];
sx q[3];
rz(-0.8662681) q[3];
sx q[3];
rz(-2.7570799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
