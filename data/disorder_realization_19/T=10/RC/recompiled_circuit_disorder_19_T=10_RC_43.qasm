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
rz(-1.6898328) q[0];
sx q[0];
rz(0.64557689) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(-1.3728377) q[1];
sx q[1];
rz(1.6436613) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71292114) q[0];
sx q[0];
rz(-0.38654583) q[0];
sx q[0];
rz(-1.9624233) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.475392) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(2.9205517) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9420535) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(0.3120504) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1708158) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(-0.67584544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(2.3085964) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67265636) q[0];
sx q[0];
rz(-2.7514973) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(2.3166336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23404113) q[0];
sx q[0];
rz(-1.3836765) q[0];
sx q[0];
rz(1.9812816) q[0];
rz(1.5603793) q[2];
sx q[2];
rz(-1.1679107) q[2];
sx q[2];
rz(-1.8889697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1383914) q[1];
sx q[1];
rz(-2.1304641) q[1];
sx q[1];
rz(-1.4336587) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1373848) q[3];
sx q[3];
rz(-1.405218) q[3];
sx q[3];
rz(2.2752938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-0.17641243) q[2];
rz(3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(-0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(3.1030531) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3241203) q[0];
sx q[0];
rz(-1.5701553) q[0];
sx q[0];
rz(-1.5785494) q[0];
rz(-pi) q[1];
rz(-1.4171717) q[2];
sx q[2];
rz(-1.3426174) q[2];
sx q[2];
rz(-1.7287901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4290532) q[1];
sx q[1];
rz(-0.4936115) q[1];
sx q[1];
rz(-1.7008971) q[1];
x q[2];
rz(-0.79629691) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(0.89087668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(1.3085261) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(2.9464338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.6911117) q[0];
sx q[0];
rz(3.0483732) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(1.3865711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26804231) q[1];
sx q[1];
rz(-2.6252384) q[1];
sx q[1];
rz(-1.079308) q[1];
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
rz(0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.8738497) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(-0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75487126) q[0];
sx q[0];
rz(-0.62749642) q[0];
sx q[0];
rz(-2.0621215) q[0];
x q[1];
rz(-2.3374704) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(2.5428307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7996019) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(2.0288543) q[1];
x q[2];
rz(-2.7630745) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-2.7639672) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-0.91517085) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.299303) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(-1.6567985) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85257951) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(0.41350565) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9050161) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(2.828039) q[1];
rz(1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(-2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(-0.36703584) q[0];
rz(-1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.4253915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050497748) q[0];
sx q[0];
rz(-1.1384581) q[0];
sx q[0];
rz(1.7763441) q[0];
rz(3.0076214) q[2];
sx q[2];
rz(-0.44951648) q[2];
sx q[2];
rz(-1.9884895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.047445) q[1];
sx q[1];
rz(-2.3970251) q[1];
sx q[1];
rz(2.5658957) q[1];
rz(-pi) q[2];
rz(1.5951049) q[3];
sx q[3];
rz(-2.6377502) q[3];
sx q[3];
rz(-2.3712096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9795064) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(-1.9308176) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(-2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(3.0650744) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-2.3892367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90826666) q[0];
sx q[0];
rz(-1.255863) q[0];
sx q[0];
rz(-0.53091913) q[0];
rz(-pi) q[1];
rz(-2.2194527) q[2];
sx q[2];
rz(-1.6862009) q[2];
sx q[2];
rz(-0.61308544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4560495) q[1];
sx q[1];
rz(-1.6838264) q[1];
sx q[1];
rz(3.0526524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30267834) q[3];
sx q[3];
rz(-2.5303826) q[3];
sx q[3];
rz(1.6906893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-1.0060271) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(-1.8267652) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7816759) q[0];
sx q[0];
rz(-1.3206375) q[0];
sx q[0];
rz(2.4292388) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36395145) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(-3.0968015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6445551) q[1];
sx q[1];
rz(-2.2215448) q[1];
sx q[1];
rz(-2.6130996) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1341348) q[3];
sx q[3];
rz(-0.82743001) q[3];
sx q[3];
rz(-2.5936562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7018147) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(-0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(-1.5022087) q[1];
sx q[1];
rz(0.39189664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4927917) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(0.0097983629) q[0];
rz(-pi) q[1];
rz(-1.7461807) q[2];
sx q[2];
rz(-0.66329623) q[2];
sx q[2];
rz(1.7593918) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3956986) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.8383154) q[1];
rz(2.6565353) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(-2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.77970589) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(2.7950381) q[3];
sx q[3];
rz(-2.0934436) q[3];
sx q[3];
rz(2.2128076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
