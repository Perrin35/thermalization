OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26389709) q[0];
sx q[0];
rz(2.1030302) q[0];
sx q[0];
rz(9.8197688) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(-1.7698987) q[1];
sx q[1];
rz(2.3514907) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34694865) q[0];
sx q[0];
rz(-0.99604368) q[0];
sx q[0];
rz(-1.6084421) q[0];
x q[1];
rz(2.3656125) q[2];
sx q[2];
rz(-1.8939514) q[2];
sx q[2];
rz(-2.6903642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8589646) q[1];
sx q[1];
rz(-0.55759768) q[1];
sx q[1];
rz(1.2931812) q[1];
rz(-pi) q[2];
rz(-2.4710841) q[3];
sx q[3];
rz(-0.6165579) q[3];
sx q[3];
rz(-2.3634865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52749741) q[2];
sx q[2];
rz(-0.082753269) q[2];
sx q[2];
rz(0.56438524) q[2];
rz(0.91974059) q[3];
sx q[3];
rz(-0.65110937) q[3];
sx q[3];
rz(1.1434327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.1564002) q[0];
sx q[0];
rz(-0.99034482) q[0];
sx q[0];
rz(-1.9957969) q[0];
rz(2.1339553) q[1];
sx q[1];
rz(-0.8877019) q[1];
sx q[1];
rz(0.064858286) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6332209) q[0];
sx q[0];
rz(-2.457058) q[0];
sx q[0];
rz(-0.9100911) q[0];
rz(-pi) q[1];
rz(-0.95352651) q[2];
sx q[2];
rz(-0.89260403) q[2];
sx q[2];
rz(-1.2645406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25003281) q[1];
sx q[1];
rz(-1.3866826) q[1];
sx q[1];
rz(1.5450983) q[1];
rz(-pi) q[2];
rz(2.3228538) q[3];
sx q[3];
rz(-0.68065879) q[3];
sx q[3];
rz(-0.37938947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2603944) q[2];
sx q[2];
rz(-0.10513267) q[2];
sx q[2];
rz(0.51327389) q[2];
rz(1.5848292) q[3];
sx q[3];
rz(-1.1791752) q[3];
sx q[3];
rz(-1.56196) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8455115) q[0];
sx q[0];
rz(-3.0873612) q[0];
sx q[0];
rz(-2.2965812) q[0];
rz(2.4301279) q[1];
sx q[1];
rz(-2.3399946) q[1];
sx q[1];
rz(2.5149288) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7539106) q[0];
sx q[0];
rz(-1.57608) q[0];
sx q[0];
rz(1.7668916) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0490295) q[2];
sx q[2];
rz(-2.704552) q[2];
sx q[2];
rz(2.4573987) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6271403) q[1];
sx q[1];
rz(-2.4042292) q[1];
sx q[1];
rz(1.3080578) q[1];
rz(-pi) q[2];
rz(-3.1342642) q[3];
sx q[3];
rz(-1.0798608) q[3];
sx q[3];
rz(1.0388663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5034379) q[2];
sx q[2];
rz(-1.2388836) q[2];
sx q[2];
rz(-0.16866355) q[2];
rz(2.9851798) q[3];
sx q[3];
rz(-2.5037239) q[3];
sx q[3];
rz(-0.28353459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6969358) q[0];
sx q[0];
rz(-2.038027) q[0];
sx q[0];
rz(-0.37100938) q[0];
rz(-1.8095398) q[1];
sx q[1];
rz(-1.5468372) q[1];
sx q[1];
rz(2.7807049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7545032) q[0];
sx q[0];
rz(-1.6564122) q[0];
sx q[0];
rz(1.4748917) q[0];
rz(-pi) q[1];
rz(-0.12207403) q[2];
sx q[2];
rz(-1.586682) q[2];
sx q[2];
rz(-2.0636914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.12965716) q[1];
sx q[1];
rz(-1.3016372) q[1];
sx q[1];
rz(0.98711826) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0293352) q[3];
sx q[3];
rz(-2.1730221) q[3];
sx q[3];
rz(-0.15439872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0648301) q[2];
sx q[2];
rz(-0.37134376) q[2];
sx q[2];
rz(1.4999464) q[2];
rz(2.0867945) q[3];
sx q[3];
rz(-1.0902371) q[3];
sx q[3];
rz(-1.1928026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.763279) q[0];
sx q[0];
rz(-1.167647) q[0];
sx q[0];
rz(1.8462697) q[0];
rz(0.12652215) q[1];
sx q[1];
rz(-0.69252068) q[1];
sx q[1];
rz(-1.2476547) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1039377) q[0];
sx q[0];
rz(-1.7154365) q[0];
sx q[0];
rz(-2.8580998) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40389668) q[2];
sx q[2];
rz(-2.3314709) q[2];
sx q[2];
rz(-0.4357117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8101471) q[1];
sx q[1];
rz(-1.1964015) q[1];
sx q[1];
rz(-1.7935171) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99301569) q[3];
sx q[3];
rz(-0.3330001) q[3];
sx q[3];
rz(1.1525994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7724472) q[2];
sx q[2];
rz(-0.36448604) q[2];
sx q[2];
rz(0.44949731) q[2];
rz(-0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(-2.0800169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19675572) q[0];
sx q[0];
rz(-2.0773092) q[0];
sx q[0];
rz(0.34360316) q[0];
rz(0.65394872) q[1];
sx q[1];
rz(-2.3908354) q[1];
sx q[1];
rz(0.59069815) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4784929) q[0];
sx q[0];
rz(-1.0859153) q[0];
sx q[0];
rz(2.8829178) q[0];
x q[1];
rz(2.3542064) q[2];
sx q[2];
rz(-0.91286589) q[2];
sx q[2];
rz(-1.7869365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0959582) q[1];
sx q[1];
rz(-1.5603298) q[1];
sx q[1];
rz(-1.4399308) q[1];
rz(-0.88149397) q[3];
sx q[3];
rz(-2.1432264) q[3];
sx q[3];
rz(1.6400169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5609453) q[2];
sx q[2];
rz(-0.86981589) q[2];
sx q[2];
rz(-0.84442863) q[2];
rz(2.1010418) q[3];
sx q[3];
rz(-1.8620164) q[3];
sx q[3];
rz(-1.1738663) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7548673) q[0];
sx q[0];
rz(-2.6193021) q[0];
sx q[0];
rz(-2.0544384) q[0];
rz(0.48671752) q[1];
sx q[1];
rz(-1.4954647) q[1];
sx q[1];
rz(-1.6513991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9947858) q[0];
sx q[0];
rz(-1.4658966) q[0];
sx q[0];
rz(-2.278743) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1370239) q[2];
sx q[2];
rz(-2.2360768) q[2];
sx q[2];
rz(2.4481004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5512276) q[1];
sx q[1];
rz(-1.9842973) q[1];
sx q[1];
rz(-1.9848003) q[1];
x q[2];
rz(0.20743117) q[3];
sx q[3];
rz(-1.279502) q[3];
sx q[3];
rz(-2.8136611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3460441) q[2];
sx q[2];
rz(-2.2717924) q[2];
sx q[2];
rz(0.33934936) q[2];
rz(-2.7782373) q[3];
sx q[3];
rz(-2.8715869) q[3];
sx q[3];
rz(1.5615162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5245847) q[0];
sx q[0];
rz(-2.0211077) q[0];
sx q[0];
rz(0.28550112) q[0];
rz(-0.78954804) q[1];
sx q[1];
rz(-1.5317081) q[1];
sx q[1];
rz(1.550536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.199478) q[0];
sx q[0];
rz(-1.3923595) q[0];
sx q[0];
rz(2.8392544) q[0];
x q[1];
rz(-1.8218846) q[2];
sx q[2];
rz(-2.6007915) q[2];
sx q[2];
rz(3.1201517) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4553677) q[1];
sx q[1];
rz(-1.9880725) q[1];
sx q[1];
rz(0.57699202) q[1];
x q[2];
rz(-1.490216) q[3];
sx q[3];
rz(-0.90513193) q[3];
sx q[3];
rz(-2.5192266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32621128) q[2];
sx q[2];
rz(-0.15924328) q[2];
sx q[2];
rz(0.86686575) q[2];
rz(-0.73831144) q[3];
sx q[3];
rz(-1.5840931) q[3];
sx q[3];
rz(-0.90832925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.065780491) q[0];
sx q[0];
rz(-2.3920238) q[0];
sx q[0];
rz(2.442389) q[0];
rz(2.5993644) q[1];
sx q[1];
rz(-1.3424073) q[1];
sx q[1];
rz(-0.29621616) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86481536) q[0];
sx q[0];
rz(-0.84057284) q[0];
sx q[0];
rz(3.0481799) q[0];
rz(-pi) q[1];
rz(-2.2454295) q[2];
sx q[2];
rz(-2.2632209) q[2];
sx q[2];
rz(1.9085753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6449127) q[1];
sx q[1];
rz(-2.1261922) q[1];
sx q[1];
rz(-0.72873656) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0551137) q[3];
sx q[3];
rz(-1.9473377) q[3];
sx q[3];
rz(-1.915177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88973033) q[2];
sx q[2];
rz(-0.44318649) q[2];
sx q[2];
rz(-2.6207793) q[2];
rz(3.1098747) q[3];
sx q[3];
rz(-2.7364276) q[3];
sx q[3];
rz(0.051137663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48145914) q[0];
sx q[0];
rz(-2.0429459) q[0];
sx q[0];
rz(-2.3302186) q[0];
rz(2.7157281) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(1.834555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.725106) q[0];
sx q[0];
rz(-1.567868) q[0];
sx q[0];
rz(3.1386969) q[0];
rz(2.0352896) q[2];
sx q[2];
rz(-1.5896322) q[2];
sx q[2];
rz(-2.1003758) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2073686) q[1];
sx q[1];
rz(-2.0133874) q[1];
sx q[1];
rz(1.9277444) q[1];
rz(-0.47932239) q[3];
sx q[3];
rz(-0.36277825) q[3];
sx q[3];
rz(2.5517705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9797152) q[2];
sx q[2];
rz(-0.76354176) q[2];
sx q[2];
rz(2.7757576) q[2];
rz(-0.93519917) q[3];
sx q[3];
rz(-1.5949944) q[3];
sx q[3];
rz(-0.38177761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314421) q[0];
sx q[0];
rz(-1.7773542) q[0];
sx q[0];
rz(-1.6553028) q[0];
rz(-1.4797795) q[1];
sx q[1];
rz(-1.9098837) q[1];
sx q[1];
rz(2.2178537) q[1];
rz(-2.1187028) q[2];
sx q[2];
rz(-1.6216424) q[2];
sx q[2];
rz(-2.4301823) q[2];
rz(0.26461919) q[3];
sx q[3];
rz(-1.7195279) q[3];
sx q[3];
rz(-1.4715094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
