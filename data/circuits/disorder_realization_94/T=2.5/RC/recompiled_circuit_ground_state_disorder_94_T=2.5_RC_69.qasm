OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2296978) q[0];
sx q[0];
rz(4.3877572) q[0];
sx q[0];
rz(10.945209) q[0];
rz(-2.7819832) q[1];
sx q[1];
rz(-0.26878992) q[1];
sx q[1];
rz(-2.6263045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775151) q[0];
sx q[0];
rz(-2.7381445) q[0];
sx q[0];
rz(0.53865428) q[0];
x q[1];
rz(0.47364606) q[2];
sx q[2];
rz(-0.59913991) q[2];
sx q[2];
rz(1.1205709) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1396619) q[1];
sx q[1];
rz(-1.3593622) q[1];
sx q[1];
rz(2.7359208) q[1];
rz(-0.73909594) q[3];
sx q[3];
rz(-2.1198273) q[3];
sx q[3];
rz(2.9573553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23060736) q[2];
sx q[2];
rz(-1.5585941) q[2];
sx q[2];
rz(0.67414635) q[2];
rz(-0.19787431) q[3];
sx q[3];
rz(-1.2400235) q[3];
sx q[3];
rz(1.6363293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3182217) q[0];
sx q[0];
rz(-0.30105337) q[0];
sx q[0];
rz(0.2163042) q[0];
rz(-0.11953106) q[1];
sx q[1];
rz(-1.1263589) q[1];
sx q[1];
rz(-1.4215887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5622373) q[0];
sx q[0];
rz(-0.52412141) q[0];
sx q[0];
rz(1.7299132) q[0];
rz(-pi) q[1];
rz(1.7825215) q[2];
sx q[2];
rz(-1.9634517) q[2];
sx q[2];
rz(0.77048877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9227289) q[1];
sx q[1];
rz(-2.0170171) q[1];
sx q[1];
rz(0.34901825) q[1];
rz(-pi) q[2];
rz(0.64038527) q[3];
sx q[3];
rz(-2.3043568) q[3];
sx q[3];
rz(-3.0471826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8372832) q[2];
sx q[2];
rz(-1.9459566) q[2];
sx q[2];
rz(-0.99011123) q[2];
rz(-0.75634161) q[3];
sx q[3];
rz(-2.4939311) q[3];
sx q[3];
rz(0.17624779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15683098) q[0];
sx q[0];
rz(-0.75060833) q[0];
sx q[0];
rz(-1.1358776) q[0];
rz(2.1319977) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(0.44581595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899781) q[0];
sx q[0];
rz(-1.5015242) q[0];
sx q[0];
rz(-1.4659856) q[0];
x q[1];
rz(3.1059044) q[2];
sx q[2];
rz(-0.70377398) q[2];
sx q[2];
rz(1.698871) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.517934) q[1];
sx q[1];
rz(-1.1594738) q[1];
sx q[1];
rz(-0.92075303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6246844) q[3];
sx q[3];
rz(-2.0175551) q[3];
sx q[3];
rz(1.7487546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3710215) q[2];
sx q[2];
rz(-1.0022481) q[2];
sx q[2];
rz(-0.64024964) q[2];
rz(-2.443743) q[3];
sx q[3];
rz(-1.6777439) q[3];
sx q[3];
rz(2.8016134) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3014389) q[0];
sx q[0];
rz(-2.2541663) q[0];
sx q[0];
rz(-1.812717) q[0];
rz(-1.2170732) q[1];
sx q[1];
rz(-2.4220059) q[1];
sx q[1];
rz(2.365239) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76809873) q[0];
sx q[0];
rz(-1.2006309) q[0];
sx q[0];
rz(-2.4959356) q[0];
x q[1];
rz(-2.0587875) q[2];
sx q[2];
rz(-1.0363724) q[2];
sx q[2];
rz(-0.39235175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66527077) q[1];
sx q[1];
rz(-1.5821536) q[1];
sx q[1];
rz(1.2922056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68459529) q[3];
sx q[3];
rz(-1.5719613) q[3];
sx q[3];
rz(0.62546414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2886469) q[2];
sx q[2];
rz(-2.8198346) q[2];
sx q[2];
rz(1.9256437) q[2];
rz(-1.1207885) q[3];
sx q[3];
rz(-1.2633163) q[3];
sx q[3];
rz(2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057366) q[0];
sx q[0];
rz(-0.97705066) q[0];
sx q[0];
rz(2.6781154) q[0];
rz(-1.0889168) q[1];
sx q[1];
rz(-1.8810279) q[1];
sx q[1];
rz(1.3190528) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.656502) q[0];
sx q[0];
rz(-2.7816371) q[0];
sx q[0];
rz(0.75004049) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9761808) q[2];
sx q[2];
rz(-2.3369419) q[2];
sx q[2];
rz(1.9555447) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12015776) q[1];
sx q[1];
rz(-1.48367) q[1];
sx q[1];
rz(-2.1539262) q[1];
rz(-pi) q[2];
rz(1.1811046) q[3];
sx q[3];
rz(-0.44868413) q[3];
sx q[3];
rz(-1.1758903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17135458) q[2];
sx q[2];
rz(-2.3562608) q[2];
sx q[2];
rz(0.63924092) q[2];
rz(0.81131896) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(1.607224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2163579) q[0];
sx q[0];
rz(-0.41330591) q[0];
sx q[0];
rz(-0.21892029) q[0];
rz(1.2159329) q[1];
sx q[1];
rz(-0.464012) q[1];
sx q[1];
rz(1.4930412) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0279044) q[0];
sx q[0];
rz(-2.0311715) q[0];
sx q[0];
rz(1.4186064) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0068917787) q[2];
sx q[2];
rz(-0.41756072) q[2];
sx q[2];
rz(1.2456976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56144372) q[1];
sx q[1];
rz(-0.83858788) q[1];
sx q[1];
rz(-1.5931908) q[1];
rz(-pi) q[2];
rz(-1.5739176) q[3];
sx q[3];
rz(-1.6699446) q[3];
sx q[3];
rz(-1.8120476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1047989) q[2];
sx q[2];
rz(-1.7918189) q[2];
sx q[2];
rz(1.3118504) q[2];
rz(-1.6364243) q[3];
sx q[3];
rz(-1.2994095) q[3];
sx q[3];
rz(-0.45984355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63603193) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(-1.7946515) q[0];
rz(-1.8147644) q[1];
sx q[1];
rz(-0.83419269) q[1];
sx q[1];
rz(2.7562275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1043804) q[0];
sx q[0];
rz(-1.9297486) q[0];
sx q[0];
rz(1.5370374) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1009675) q[2];
sx q[2];
rz(-2.0010002) q[2];
sx q[2];
rz(0.29597798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2750157) q[1];
sx q[1];
rz(-2.6948435) q[1];
sx q[1];
rz(-0.061468852) q[1];
rz(-0.47355424) q[3];
sx q[3];
rz(-1.9723168) q[3];
sx q[3];
rz(2.7724289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68606004) q[2];
sx q[2];
rz(-2.3251688) q[2];
sx q[2];
rz(-1.653999) q[2];
rz(-0.16573302) q[3];
sx q[3];
rz(-0.82560742) q[3];
sx q[3];
rz(-0.17416557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1626749) q[0];
sx q[0];
rz(-2.5161777) q[0];
sx q[0];
rz(0.22853525) q[0];
rz(2.8288815) q[1];
sx q[1];
rz(-2.2622908) q[1];
sx q[1];
rz(-1.3794587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9420825) q[0];
sx q[0];
rz(-1.4934015) q[0];
sx q[0];
rz(-0.40298526) q[0];
rz(0.32416043) q[2];
sx q[2];
rz(-1.7334786) q[2];
sx q[2];
rz(-1.206347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91614281) q[1];
sx q[1];
rz(-1.6054549) q[1];
sx q[1];
rz(-1.9198656) q[1];
rz(-pi) q[2];
rz(0.70693858) q[3];
sx q[3];
rz(-1.1683373) q[3];
sx q[3];
rz(-1.5618531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.032893) q[2];
sx q[2];
rz(-1.0445107) q[2];
sx q[2];
rz(-1.0408939) q[2];
rz(0.4852455) q[3];
sx q[3];
rz(-1.8294561) q[3];
sx q[3];
rz(2.7237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8682206) q[0];
sx q[0];
rz(-0.98216787) q[0];
sx q[0];
rz(0.81106538) q[0];
rz(1.2366933) q[1];
sx q[1];
rz(-2.2752454) q[1];
sx q[1];
rz(-0.19485697) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92835125) q[0];
sx q[0];
rz(-1.8272864) q[0];
sx q[0];
rz(2.8942906) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56320517) q[2];
sx q[2];
rz(-2.1426472) q[2];
sx q[2];
rz(0.88746136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5611546) q[1];
sx q[1];
rz(-1.913979) q[1];
sx q[1];
rz(-2.047206) q[1];
rz(-2.2236597) q[3];
sx q[3];
rz(-1.1102765) q[3];
sx q[3];
rz(-1.9916818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98086944) q[2];
sx q[2];
rz(-1.2474493) q[2];
sx q[2];
rz(2.5998739) q[2];
rz(1.3228275) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(2.8790348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59170359) q[0];
sx q[0];
rz(-1.5174958) q[0];
sx q[0];
rz(0.24895915) q[0];
rz(-1.9242363) q[1];
sx q[1];
rz(-1.267642) q[1];
sx q[1];
rz(-2.731146) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7792203) q[0];
sx q[0];
rz(-1.1817314) q[0];
sx q[0];
rz(2.7011407) q[0];
x q[1];
rz(-1.6617695) q[2];
sx q[2];
rz(-1.5814648) q[2];
sx q[2];
rz(-1.4271229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9767993) q[1];
sx q[1];
rz(-1.8135957) q[1];
sx q[1];
rz(2.040928) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4725222) q[3];
sx q[3];
rz(-0.72815547) q[3];
sx q[3];
rz(1.6403584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2534788) q[2];
sx q[2];
rz(-1.1682744) q[2];
sx q[2];
rz(1.1178364) q[2];
rz(-3.0228293) q[3];
sx q[3];
rz(-0.6898841) q[3];
sx q[3];
rz(-1.2004872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73040199) q[0];
sx q[0];
rz(-1.7350736) q[0];
sx q[0];
rz(-0.34179678) q[0];
rz(3.0939915) q[1];
sx q[1];
rz(-1.8879415) q[1];
sx q[1];
rz(-1.3157848) q[1];
rz(-2.32687) q[2];
sx q[2];
rz(-1.9152894) q[2];
sx q[2];
rz(0.49283129) q[2];
rz(2.1988395) q[3];
sx q[3];
rz(-2.2136064) q[3];
sx q[3];
rz(2.8436713) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
