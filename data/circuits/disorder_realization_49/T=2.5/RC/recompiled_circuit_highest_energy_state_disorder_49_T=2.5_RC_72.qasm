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
rz(0.024700392) q[0];
sx q[0];
rz(-1.2588809) q[0];
sx q[0];
rz(1.2727241) q[0];
rz(-1.7495573) q[1];
sx q[1];
rz(-1.7311544) q[1];
sx q[1];
rz(-2.2478204) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0451242) q[0];
sx q[0];
rz(-0.91334263) q[0];
sx q[0];
rz(-2.3923001) q[0];
rz(-0.29869334) q[2];
sx q[2];
rz(-1.325404) q[2];
sx q[2];
rz(-1.1358062) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3013348) q[1];
sx q[1];
rz(-0.71941676) q[1];
sx q[1];
rz(0.69566984) q[1];
rz(-pi) q[2];
rz(0.36811604) q[3];
sx q[3];
rz(-1.5874581) q[3];
sx q[3];
rz(-2.1429495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7749403) q[2];
sx q[2];
rz(-1.4619991) q[2];
sx q[2];
rz(2.6095663) q[2];
rz(0.73412791) q[3];
sx q[3];
rz(-0.2747772) q[3];
sx q[3];
rz(1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183913) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(-0.080408737) q[0];
rz(-0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(-0.72371662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88392576) q[0];
sx q[0];
rz(-2.2307255) q[0];
sx q[0];
rz(-1.0038478) q[0];
x q[1];
rz(-2.2369464) q[2];
sx q[2];
rz(-2.0331185) q[2];
sx q[2];
rz(-1.8682299) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5608175) q[1];
sx q[1];
rz(-2.5530949) q[1];
sx q[1];
rz(1.568041) q[1];
rz(-pi) q[2];
rz(-2.1172013) q[3];
sx q[3];
rz(-0.83163578) q[3];
sx q[3];
rz(1.7061159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(-2.2354324) q[2];
rz(0.3624889) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(-3.0947065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
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
rz(1.2694117) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(-1.0915225) q[0];
rz(3.0190234) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(1.3465808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4849834) q[0];
sx q[0];
rz(-0.35728982) q[0];
sx q[0];
rz(-1.4862509) q[0];
rz(-pi) q[1];
rz(2.7287219) q[2];
sx q[2];
rz(-0.70253583) q[2];
sx q[2];
rz(1.9091878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8217762) q[1];
sx q[1];
rz(-0.21142658) q[1];
sx q[1];
rz(1.3216622) q[1];
rz(2.6123206) q[3];
sx q[3];
rz(-1.4733757) q[3];
sx q[3];
rz(2.4306963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4284105) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(-2.909519) q[2];
rz(-1.170916) q[3];
sx q[3];
rz(-0.62058312) q[3];
sx q[3];
rz(-2.8569729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93037) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(-1.5027745) q[0];
rz(0.59208313) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(0.88964644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1574123) q[0];
sx q[0];
rz(-1.1129459) q[0];
sx q[0];
rz(-2.1628987) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2438898) q[2];
sx q[2];
rz(-2.0497344) q[2];
sx q[2];
rz(-0.51264352) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6742176) q[1];
sx q[1];
rz(-1.2548037) q[1];
sx q[1];
rz(0.18958474) q[1];
rz(-pi) q[2];
rz(2.3208404) q[3];
sx q[3];
rz(-2.6057557) q[3];
sx q[3];
rz(-1.7465296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64638102) q[2];
sx q[2];
rz(-1.6782328) q[2];
sx q[2];
rz(0.67592534) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(2.9743312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306358) q[0];
sx q[0];
rz(-0.3322596) q[0];
sx q[0];
rz(2.9462872) q[0];
rz(-3.0209814) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-2.9511071) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99283907) q[0];
sx q[0];
rz(-1.5061814) q[0];
sx q[0];
rz(3.0535327) q[0];
rz(-3.0342874) q[2];
sx q[2];
rz(-2.4621747) q[2];
sx q[2];
rz(-0.36919644) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0234649) q[1];
sx q[1];
rz(-1.6227229) q[1];
sx q[1];
rz(-3.0408188) q[1];
rz(-pi) q[2];
rz(-1.9917914) q[3];
sx q[3];
rz(-0.99830571) q[3];
sx q[3];
rz(1.5508625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6327989) q[2];
sx q[2];
rz(-1.4845279) q[2];
sx q[2];
rz(1.9602027) q[2];
rz(-1.7975636) q[3];
sx q[3];
rz(-2.400178) q[3];
sx q[3];
rz(0.77176362) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8089499) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(-2.0020265) q[0];
rz(2.471916) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(1.12961) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87605873) q[0];
sx q[0];
rz(-1.8647412) q[0];
sx q[0];
rz(-0.29220279) q[0];
rz(-pi) q[1];
rz(2.8390719) q[2];
sx q[2];
rz(-0.74069689) q[2];
sx q[2];
rz(-2.4940235) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.672513) q[1];
sx q[1];
rz(-1.7653078) q[1];
sx q[1];
rz(-2.7698293) q[1];
rz(-pi) q[2];
rz(-2.9608742) q[3];
sx q[3];
rz(-1.4152495) q[3];
sx q[3];
rz(0.29744086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42664042) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(0.40531522) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-2.5155641) q[3];
sx q[3];
rz(0.85957447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7760794) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(0.93589163) q[0];
rz(2.2731958) q[1];
sx q[1];
rz(-2.6135542) q[1];
sx q[1];
rz(-1.4194999) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4286242) q[0];
sx q[0];
rz(-1.0626864) q[0];
sx q[0];
rz(1.1626121) q[0];
rz(0.23250696) q[2];
sx q[2];
rz(-2.0662796) q[2];
sx q[2];
rz(-1.2459618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.81379997) q[1];
sx q[1];
rz(-2.2918502) q[1];
sx q[1];
rz(-2.2071597) q[1];
rz(-1.5704182) q[3];
sx q[3];
rz(-1.5725563) q[3];
sx q[3];
rz(-1.5906217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3203656) q[2];
sx q[2];
rz(-1.7419523) q[2];
sx q[2];
rz(-1.0199176) q[2];
rz(0.50926456) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34198636) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(1.1588143) q[0];
rz(-2.671057) q[1];
sx q[1];
rz(-1.4152941) q[1];
sx q[1];
rz(-1.6758957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7635135) q[0];
sx q[0];
rz(-1.5944885) q[0];
sx q[0];
rz(-3.0906584) q[0];
rz(0.69078867) q[2];
sx q[2];
rz(-1.8422519) q[2];
sx q[2];
rz(-1.1784306) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1294687) q[1];
sx q[1];
rz(-0.91255847) q[1];
sx q[1];
rz(-1.4736498) q[1];
rz(1.2837778) q[3];
sx q[3];
rz(-2.1987408) q[3];
sx q[3];
rz(2.7046741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1450242) q[2];
sx q[2];
rz(-1.8137167) q[2];
sx q[2];
rz(-0.81929755) q[2];
rz(-2.3451292) q[3];
sx q[3];
rz(-2.7345246) q[3];
sx q[3];
rz(-1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040937) q[0];
sx q[0];
rz(-0.096385328) q[0];
sx q[0];
rz(0.82164422) q[0];
rz(0.67283982) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(2.0988665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48548904) q[0];
sx q[0];
rz(-2.9239131) q[0];
sx q[0];
rz(-0.34164048) q[0];
rz(-pi) q[1];
rz(2.369528) q[2];
sx q[2];
rz(-2.0520718) q[2];
sx q[2];
rz(1.5636064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.753841) q[1];
sx q[1];
rz(-2.4341704) q[1];
sx q[1];
rz(-2.161987) q[1];
rz(-pi) q[2];
rz(-0.25298499) q[3];
sx q[3];
rz(-1.7182941) q[3];
sx q[3];
rz(1.1563039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3599856) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(2.2992415) q[2];
rz(0.53884566) q[3];
sx q[3];
rz(-1.6539961) q[3];
sx q[3];
rz(-0.52411383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.6738324) q[0];
sx q[0];
rz(-0.14586511) q[0];
sx q[0];
rz(-1.9061331) q[0];
rz(2.1927059) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(1.4804776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7296071) q[0];
sx q[0];
rz(-2.3959037) q[0];
sx q[0];
rz(1.1545769) q[0];
rz(-pi) q[1];
rz(-2.1820004) q[2];
sx q[2];
rz(-1.094154) q[2];
sx q[2];
rz(-0.075919064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93714156) q[1];
sx q[1];
rz(-1.4960663) q[1];
sx q[1];
rz(2.9765169) q[1];
rz(-pi) q[2];
rz(-1.8126902) q[3];
sx q[3];
rz(-1.8902361) q[3];
sx q[3];
rz(1.1238568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4296055) q[2];
sx q[2];
rz(-1.7179855) q[2];
sx q[2];
rz(-0.31022662) q[2];
rz(-0.45983908) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(-1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.9643758) q[0];
sx q[0];
rz(-2.0073267) q[0];
sx q[0];
rz(1.0649756) q[0];
rz(0.89827697) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(2.73478) q[2];
sx q[2];
rz(-2.1944703) q[2];
sx q[2];
rz(-2.2866972) q[2];
rz(-2.2326917) q[3];
sx q[3];
rz(-1.881897) q[3];
sx q[3];
rz(1.9492634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
