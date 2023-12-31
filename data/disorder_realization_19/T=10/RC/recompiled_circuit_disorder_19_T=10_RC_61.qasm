OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(4.5933525) q[0];
sx q[0];
rz(10.070355) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(-1.3728377) q[1];
sx q[1];
rz(1.6436613) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0092702) q[0];
sx q[0];
rz(-1.9267123) q[0];
sx q[0];
rz(-0.15412553) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2013118) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(-1.8688569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.363417) q[1];
sx q[1];
rz(-2.2925804) q[1];
sx q[1];
rz(1.282882) q[1];
x q[2];
rz(-0.50819355) q[3];
sx q[3];
rz(-2.4262706) q[3];
sx q[3];
rz(-0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-0.54498589) q[0];
rz(0.90822059) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-2.3166336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4175407) q[0];
sx q[0];
rz(-1.9736971) q[0];
sx q[0];
rz(2.9379662) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024436342) q[2];
sx q[2];
rz(-2.7385798) q[2];
sx q[2];
rz(1.8624061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3927041) q[1];
sx q[1];
rz(-0.57447937) q[1];
sx q[1];
rz(-0.21484612) q[1];
x q[2];
rz(-1.0042079) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(-2.2752938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-0.17641243) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(0.032657284) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(-1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(-3.1030531) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8124213) q[0];
sx q[0];
rz(-3.1338131) q[0];
sx q[0];
rz(1.4882985) q[0];
x q[1];
rz(0.58263393) q[2];
sx q[2];
rz(-2.8672672) q[2];
sx q[2];
rz(-2.3290616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.027027834) q[1];
sx q[1];
rz(-1.5092883) q[1];
sx q[1];
rz(-2.0608749) q[1];
x q[2];
rz(0.97088082) q[3];
sx q[3];
rz(-2.0754793) q[3];
sx q[3];
rz(1.8463299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0168266) q[0];
sx q[0];
rz(-0.15206465) q[0];
sx q[0];
rz(2.2269339) q[0];
rz(-pi) q[1];
rz(-2.6606576) q[2];
sx q[2];
rz(-2.0630884) q[2];
sx q[2];
rz(2.7155657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8198204) q[1];
sx q[1];
rz(-2.0211377) q[1];
sx q[1];
rz(-0.26178534) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5024662) q[3];
sx q[3];
rz(-1.5577321) q[3];
sx q[3];
rz(-0.55707219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9106456) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(-1.8738497) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(-0.1396133) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9169086) q[0];
sx q[0];
rz(-1.2901257) q[0];
sx q[0];
rz(-1.0018437) q[0];
rz(-pi) q[1];
rz(-0.80412229) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(0.59876195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2851583) q[1];
sx q[1];
rz(-2.563623) q[1];
sx q[1];
rz(0.85698378) q[1];
x q[2];
rz(-2.674621) q[3];
sx q[3];
rz(-1.3936371) q[3];
sx q[3];
rz(0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(2.1652083) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75509214) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5765502) q[0];
sx q[0];
rz(-2.8083907) q[0];
sx q[0];
rz(0.25175005) q[0];
x q[1];
rz(2.1377863) q[2];
sx q[2];
rz(-2.0096471) q[2];
sx q[2];
rz(1.7457419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9050161) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(-0.31355365) q[1];
x q[2];
rz(2.4470607) q[3];
sx q[3];
rz(-1.4834705) q[3];
sx q[3];
rz(0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(0.79052314) q[2];
rz(2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910949) q[0];
sx q[0];
rz(-2.0031345) q[0];
sx q[0];
rz(1.7763441) q[0];
rz(-pi) q[1];
rz(3.0076214) q[2];
sx q[2];
rz(-0.44951648) q[2];
sx q[2];
rz(1.1531032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.047445) q[1];
sx q[1];
rz(-0.74456753) q[1];
sx q[1];
rz(-0.57569699) q[1];
x q[2];
rz(3.1281934) q[3];
sx q[3];
rz(-1.0671167) q[3];
sx q[3];
rz(-0.79813938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(-3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84412557) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(-0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90826666) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(0.53091913) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2194527) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(2.5285072) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.016537746) q[1];
sx q[1];
rz(-2.9978831) q[1];
sx q[1];
rz(-2.2347666) q[1];
rz(-1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(-2.0549783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(1.0060271) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(-1.8267652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49004236) q[0];
sx q[0];
rz(-0.74768066) q[0];
sx q[0];
rz(0.3726532) q[0];
rz(-2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(3.0968015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.41373738) q[1];
sx q[1];
rz(-1.9836042) q[1];
sx q[1];
rz(-2.2933943) q[1];
rz(-pi) q[2];
rz(-2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(2.1237802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7018147) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(2.2626256) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4927917) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(0.0097983629) q[0];
x q[1];
rz(0.91498615) q[2];
sx q[2];
rz(-1.4631541) q[2];
sx q[2];
rz(0.3273302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7614903) q[1];
sx q[1];
rz(-2.3099766) q[1];
sx q[1];
rz(-0.2525316) q[1];
x q[2];
rz(1.8896905) q[3];
sx q[3];
rz(-1.0341757) q[3];
sx q[3];
rz(1.5434138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75523238) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(-1.8023087) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-2.3618868) q[2];
sx q[2];
rz(-1.4479326) q[2];
sx q[2];
rz(-0.89276199) q[2];
rz(-1.0380575) q[3];
sx q[3];
rz(-2.5235015) q[3];
sx q[3];
rz(-1.5550767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
