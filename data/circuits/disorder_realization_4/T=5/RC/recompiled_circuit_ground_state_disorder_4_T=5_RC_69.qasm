OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.088012785) q[0];
sx q[0];
rz(-2.8289284) q[0];
sx q[0];
rz(0.13392681) q[0];
rz(-0.0066095134) q[1];
sx q[1];
rz(4.7314965) q[1];
sx q[1];
rz(9.849698) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28025173) q[0];
sx q[0];
rz(-1.6854648) q[0];
sx q[0];
rz(3.0404509) q[0];
rz(2.0449465) q[2];
sx q[2];
rz(-0.72410781) q[2];
sx q[2];
rz(-2.6158138) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1462789) q[1];
sx q[1];
rz(-1.9094719) q[1];
sx q[1];
rz(-2.2168753) q[1];
x q[2];
rz(0.71436259) q[3];
sx q[3];
rz(-2.2788725) q[3];
sx q[3];
rz(-1.2618511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46140823) q[2];
sx q[2];
rz(-0.37698656) q[2];
sx q[2];
rz(-3.027463) q[2];
rz(0.38873172) q[3];
sx q[3];
rz(-2.8753493) q[3];
sx q[3];
rz(1.3346599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330133) q[0];
sx q[0];
rz(-2.7393434) q[0];
sx q[0];
rz(-0.026799686) q[0];
rz(1.1773479) q[1];
sx q[1];
rz(-0.64759308) q[1];
sx q[1];
rz(-3.1039544) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55521727) q[0];
sx q[0];
rz(-0.77459413) q[0];
sx q[0];
rz(0.019217773) q[0];
rz(-pi) q[1];
rz(2.9178647) q[2];
sx q[2];
rz(-1.7729974) q[2];
sx q[2];
rz(-3.1333609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85559713) q[1];
sx q[1];
rz(-1.0904795) q[1];
sx q[1];
rz(-0.91676401) q[1];
rz(2.3214171) q[3];
sx q[3];
rz(-1.7595152) q[3];
sx q[3];
rz(0.80126002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2521952) q[2];
sx q[2];
rz(-2.1571721) q[2];
sx q[2];
rz(2.9241015) q[2];
rz(2.7018231) q[3];
sx q[3];
rz(-1.35448) q[3];
sx q[3];
rz(0.10194889) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52461034) q[0];
sx q[0];
rz(-0.007402448) q[0];
sx q[0];
rz(-2.5819085) q[0];
rz(2.2451378) q[1];
sx q[1];
rz(-0.47210109) q[1];
sx q[1];
rz(0.40759531) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64227885) q[0];
sx q[0];
rz(-2.1526205) q[0];
sx q[0];
rz(-3.0034742) q[0];
rz(3.0733068) q[2];
sx q[2];
rz(-0.90203055) q[2];
sx q[2];
rz(-2.1055773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7771908) q[1];
sx q[1];
rz(-2.6620416) q[1];
sx q[1];
rz(2.9001791) q[1];
rz(-pi) q[2];
rz(1.0901426) q[3];
sx q[3];
rz(-0.60308057) q[3];
sx q[3];
rz(-1.9657976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1233623) q[2];
sx q[2];
rz(-1.1574278) q[2];
sx q[2];
rz(-2.1165712) q[2];
rz(1.4347264) q[3];
sx q[3];
rz(-0.44308174) q[3];
sx q[3];
rz(0.6492492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822815) q[0];
sx q[0];
rz(-2.870443) q[0];
sx q[0];
rz(-0.49736381) q[0];
rz(-1.2936032) q[1];
sx q[1];
rz(-2.135364) q[1];
sx q[1];
rz(-1.2327548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4068488) q[0];
sx q[0];
rz(-2.5355784) q[0];
sx q[0];
rz(2.7789609) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7035162) q[2];
sx q[2];
rz(-1.2428962) q[2];
sx q[2];
rz(-0.10645535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2805505) q[1];
sx q[1];
rz(-2.8192564) q[1];
sx q[1];
rz(-0.68962421) q[1];
rz(0.044305459) q[3];
sx q[3];
rz(-2.7409275) q[3];
sx q[3];
rz(1.6687499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5645912) q[2];
sx q[2];
rz(-2.5265054) q[2];
sx q[2];
rz(2.9992529) q[2];
rz(1.9650991) q[3];
sx q[3];
rz(-1.0005955) q[3];
sx q[3];
rz(-0.4308027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7560526) q[0];
sx q[0];
rz(-2.1591594) q[0];
sx q[0];
rz(-2.8681712) q[0];
rz(-2.7418819) q[1];
sx q[1];
rz(-2.3276261) q[1];
sx q[1];
rz(-2.8846557) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533511) q[0];
sx q[0];
rz(-1.0437766) q[0];
sx q[0];
rz(2.7668883) q[0];
rz(-1.1838205) q[2];
sx q[2];
rz(-1.3391025) q[2];
sx q[2];
rz(2.608172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2013984) q[1];
sx q[1];
rz(-1.640854) q[1];
sx q[1];
rz(-3.0959341) q[1];
x q[2];
rz(-1.2236274) q[3];
sx q[3];
rz(-0.70968443) q[3];
sx q[3];
rz(1.1598809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5895245) q[2];
sx q[2];
rz(-0.35577154) q[2];
sx q[2];
rz(0.75999981) q[2];
rz(0.99203569) q[3];
sx q[3];
rz(-1.9325247) q[3];
sx q[3];
rz(-1.9973283) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8511667) q[0];
sx q[0];
rz(-0.63153428) q[0];
sx q[0];
rz(0.60612154) q[0];
rz(2.0947314) q[1];
sx q[1];
rz(-1.650834) q[1];
sx q[1];
rz(-0.077233888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7116878) q[0];
sx q[0];
rz(-3.0992134) q[0];
sx q[0];
rz(-1.872657) q[0];
rz(0.31862835) q[2];
sx q[2];
rz(-2.0214404) q[2];
sx q[2];
rz(2.3261827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.162335) q[1];
sx q[1];
rz(-1.1282776) q[1];
sx q[1];
rz(-2.6975432) q[1];
x q[2];
rz(0.24155946) q[3];
sx q[3];
rz(-2.4643341) q[3];
sx q[3];
rz(-0.23158555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2250526) q[2];
sx q[2];
rz(-0.70987916) q[2];
sx q[2];
rz(0.27099657) q[2];
rz(2.5535876) q[3];
sx q[3];
rz(-0.80395144) q[3];
sx q[3];
rz(2.7325381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005017) q[0];
sx q[0];
rz(-1.611447) q[0];
sx q[0];
rz(2.8588168) q[0];
rz(-2.458789) q[1];
sx q[1];
rz(-2.2802201) q[1];
sx q[1];
rz(0.85321325) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7616939) q[0];
sx q[0];
rz(-0.81710862) q[0];
sx q[0];
rz(-2.0920112) q[0];
x q[1];
rz(0.81099895) q[2];
sx q[2];
rz(-1.2110707) q[2];
sx q[2];
rz(-0.0051509858) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.68487271) q[1];
sx q[1];
rz(-1.8032852) q[1];
sx q[1];
rz(0.16831919) q[1];
x q[2];
rz(2.1370579) q[3];
sx q[3];
rz(-2.9101203) q[3];
sx q[3];
rz(-1.8799409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0679438) q[2];
sx q[2];
rz(-2.6239519) q[2];
sx q[2];
rz(-0.29120564) q[2];
rz(1.2047042) q[3];
sx q[3];
rz(-2.5681345) q[3];
sx q[3];
rz(-1.5706435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7912927) q[0];
sx q[0];
rz(-1.5234103) q[0];
sx q[0];
rz(-2.5147901) q[0];
rz(-2.0794012) q[1];
sx q[1];
rz(-0.79963446) q[1];
sx q[1];
rz(-2.3041384) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7933374) q[0];
sx q[0];
rz(-1.1438864) q[0];
sx q[0];
rz(-0.12664533) q[0];
rz(-pi) q[1];
rz(1.8071354) q[2];
sx q[2];
rz(-1.9009942) q[2];
sx q[2];
rz(-0.2969674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22540846) q[1];
sx q[1];
rz(-1.7042158) q[1];
sx q[1];
rz(2.0853985) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4222048) q[3];
sx q[3];
rz(-2.0839543) q[3];
sx q[3];
rz(0.10529127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8219139) q[2];
sx q[2];
rz(-1.9210812) q[2];
sx q[2];
rz(-2.1004045) q[2];
rz(-2.239481) q[3];
sx q[3];
rz(-1.1646885) q[3];
sx q[3];
rz(-2.5743217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5015471) q[0];
sx q[0];
rz(-0.29842672) q[0];
sx q[0];
rz(3.0409467) q[0];
rz(-1.8177265) q[1];
sx q[1];
rz(-2.3034425) q[1];
sx q[1];
rz(-2.5460338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190191) q[0];
sx q[0];
rz(-0.77054502) q[0];
sx q[0];
rz(-2.3014803) q[0];
rz(-2.2329464) q[2];
sx q[2];
rz(-1.7344237) q[2];
sx q[2];
rz(-1.2589135) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42587316) q[1];
sx q[1];
rz(-1.8173479) q[1];
sx q[1];
rz(2.4049525) q[1];
x q[2];
rz(2.9902114) q[3];
sx q[3];
rz(-0.80946846) q[3];
sx q[3];
rz(-2.4224506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14463921) q[2];
sx q[2];
rz(-2.3914631) q[2];
sx q[2];
rz(-1.3422802) q[2];
rz(2.4053549) q[3];
sx q[3];
rz(-2.8056371) q[3];
sx q[3];
rz(0.098522447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017224273) q[0];
sx q[0];
rz(-2.4898744) q[0];
sx q[0];
rz(0.29618725) q[0];
rz(1.4172957) q[1];
sx q[1];
rz(-0.92767757) q[1];
sx q[1];
rz(-0.16253026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95537335) q[0];
sx q[0];
rz(-0.027039921) q[0];
sx q[0];
rz(1.3229738) q[0];
rz(0.98661454) q[2];
sx q[2];
rz(-0.59609813) q[2];
sx q[2];
rz(-1.7510027) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0386943) q[1];
sx q[1];
rz(-1.7084645) q[1];
sx q[1];
rz(0.18935151) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82270427) q[3];
sx q[3];
rz(-0.71559042) q[3];
sx q[3];
rz(-2.0377318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6235003) q[2];
sx q[2];
rz(-1.2178428) q[2];
sx q[2];
rz(2.6859542) q[2];
rz(-2.0677004) q[3];
sx q[3];
rz(-2.5116601) q[3];
sx q[3];
rz(2.5560801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940654) q[0];
sx q[0];
rz(-1.3205262) q[0];
sx q[0];
rz(2.5393215) q[0];
rz(0.31996721) q[1];
sx q[1];
rz(-1.1155557) q[1];
sx q[1];
rz(-1.3394578) q[1];
rz(1.6306277) q[2];
sx q[2];
rz(-1.2345805) q[2];
sx q[2];
rz(2.6474093) q[2];
rz(-1.4634168) q[3];
sx q[3];
rz(-0.79859514) q[3];
sx q[3];
rz(2.7664281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
