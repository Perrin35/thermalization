OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9517684) q[0];
sx q[0];
rz(-2.4387852) q[0];
sx q[0];
rz(-1.2220609) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2015704) q[0];
sx q[0];
rz(-0.68212494) q[0];
sx q[0];
rz(-1.8322102) q[0];
x q[1];
rz(1.7403549) q[2];
sx q[2];
rz(-1.6785067) q[2];
sx q[2];
rz(-0.28125276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9492717) q[1];
sx q[1];
rz(-1.1484129) q[1];
sx q[1];
rz(1.3793089) q[1];
x q[2];
rz(2.7056115) q[3];
sx q[3];
rz(-2.0153305) q[3];
sx q[3];
rz(1.7874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9606955) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(-2.1107471) q[2];
rz(-0.24813063) q[3];
sx q[3];
rz(-1.3458359) q[3];
sx q[3];
rz(-0.26721755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(-1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(-2.7148278) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65157344) q[0];
sx q[0];
rz(-0.82758862) q[0];
sx q[0];
rz(3.1378967) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5019481) q[2];
sx q[2];
rz(-0.42176127) q[2];
sx q[2];
rz(-2.1157921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2093231) q[1];
sx q[1];
rz(-0.70581064) q[1];
sx q[1];
rz(0.16684823) q[1];
x q[2];
rz(1.9598448) q[3];
sx q[3];
rz(-1.7579683) q[3];
sx q[3];
rz(2.1120918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65853226) q[2];
sx q[2];
rz(-0.51653647) q[2];
sx q[2];
rz(0.10022441) q[2];
rz(2.0641067) q[3];
sx q[3];
rz(-1.5764377) q[3];
sx q[3];
rz(1.235435) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.539262) q[0];
sx q[0];
rz(-0.97719181) q[0];
sx q[0];
rz(1.4343028) q[0];
rz(0.82646787) q[1];
sx q[1];
rz(-1.6441328) q[1];
sx q[1];
rz(-0.43063393) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.193482) q[0];
sx q[0];
rz(-1.3010345) q[0];
sx q[0];
rz(-1.4007691) q[0];
x q[1];
rz(0.13027066) q[2];
sx q[2];
rz(-1.8415057) q[2];
sx q[2];
rz(-0.66616466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0525749) q[1];
sx q[1];
rz(-2.3529426) q[1];
sx q[1];
rz(2.1441205) q[1];
rz(1.0266193) q[3];
sx q[3];
rz(-1.2261571) q[3];
sx q[3];
rz(-1.0349864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87381252) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(0.38392797) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(-2.4887776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4189932) q[0];
sx q[0];
rz(-0.79177952) q[0];
sx q[0];
rz(1.0571887) q[0];
rz(-2.3253697) q[1];
sx q[1];
rz(-0.99333251) q[1];
sx q[1];
rz(0.23424558) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62158305) q[0];
sx q[0];
rz(-2.710973) q[0];
sx q[0];
rz(-1.8355811) q[0];
rz(-pi) q[1];
rz(-2.0969056) q[2];
sx q[2];
rz(-1.9434233) q[2];
sx q[2];
rz(-1.6372731) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87061674) q[1];
sx q[1];
rz(-0.91667316) q[1];
sx q[1];
rz(2.4821217) q[1];
x q[2];
rz(-0.47611634) q[3];
sx q[3];
rz(-2.0272539) q[3];
sx q[3];
rz(1.8511839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.28698841) q[2];
sx q[2];
rz(-1.333326) q[2];
sx q[2];
rz(-2.2387779) q[2];
rz(-1.758894) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(2.2949016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.46120241) q[0];
sx q[0];
rz(-1.2795376) q[0];
sx q[0];
rz(-2.4317106) q[0];
rz(-2.6704125) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(-0.064362854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854537) q[0];
sx q[0];
rz(-0.50489932) q[0];
sx q[0];
rz(-1.288815) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9079403) q[2];
sx q[2];
rz(-1.086768) q[2];
sx q[2];
rz(-2.7357227) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2379265) q[1];
sx q[1];
rz(-1.1519379) q[1];
sx q[1];
rz(2.5861857) q[1];
rz(-pi) q[2];
rz(-0.10897763) q[3];
sx q[3];
rz(-0.087317467) q[3];
sx q[3];
rz(0.84978896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.430696) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(-0.95332471) q[2];
rz(-0.55771762) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(0.52887708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3144658) q[0];
sx q[0];
rz(-1.0449266) q[0];
sx q[0];
rz(0.99660981) q[0];
rz(0.059545513) q[1];
sx q[1];
rz(-2.153502) q[1];
sx q[1];
rz(1.7893808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9080862) q[0];
sx q[0];
rz(-2.3829989) q[0];
sx q[0];
rz(2.3367466) q[0];
rz(0.77162051) q[2];
sx q[2];
rz(-1.3305656) q[2];
sx q[2];
rz(2.8738603) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.775573) q[1];
sx q[1];
rz(-1.8692353) q[1];
sx q[1];
rz(1.2828903) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91259802) q[3];
sx q[3];
rz(-1.293598) q[3];
sx q[3];
rz(0.67597055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.754564) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-2.582029) q[2];
rz(0.83545056) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(-1.3028418) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68592042) q[0];
sx q[0];
rz(-0.73021013) q[0];
sx q[0];
rz(-2.7749104) q[0];
rz(-2.2600251) q[1];
sx q[1];
rz(-2.5036948) q[1];
sx q[1];
rz(-1.7104023) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.139643) q[0];
sx q[0];
rz(-0.77660034) q[0];
sx q[0];
rz(0.5712255) q[0];
rz(-pi) q[1];
rz(-0.41421271) q[2];
sx q[2];
rz(-2.4117232) q[2];
sx q[2];
rz(3.102906) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.28364946) q[1];
sx q[1];
rz(-1.4542631) q[1];
sx q[1];
rz(1.5999419) q[1];
rz(-pi) q[2];
rz(1.335698) q[3];
sx q[3];
rz(-2.2155016) q[3];
sx q[3];
rz(0.39261445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7876579) q[2];
sx q[2];
rz(-2.6255609) q[2];
sx q[2];
rz(2.3616974) q[2];
rz(-2.6025313) q[3];
sx q[3];
rz(-1.5122248) q[3];
sx q[3];
rz(-0.61162925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270545) q[0];
sx q[0];
rz(-0.67033613) q[0];
sx q[0];
rz(-0.77914733) q[0];
rz(0.51891333) q[1];
sx q[1];
rz(-1.2823558) q[1];
sx q[1];
rz(2.7386477) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8810324) q[0];
sx q[0];
rz(-2.0391818) q[0];
sx q[0];
rz(0.82527501) q[0];
rz(-pi) q[1];
rz(0.31488089) q[2];
sx q[2];
rz(-1.1078896) q[2];
sx q[2];
rz(-2.7969282) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25499763) q[1];
sx q[1];
rz(-0.95980058) q[1];
sx q[1];
rz(1.081089) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99767606) q[3];
sx q[3];
rz(-1.3460288) q[3];
sx q[3];
rz(-0.59319699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0275823) q[2];
sx q[2];
rz(-1.9738013) q[2];
sx q[2];
rz(1.1172509) q[2];
rz(0.69636017) q[3];
sx q[3];
rz(-2.0317234) q[3];
sx q[3];
rz(1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6663412) q[0];
sx q[0];
rz(-2.9069558) q[0];
sx q[0];
rz(-2.7269205) q[0];
rz(-2.0798202) q[1];
sx q[1];
rz(-0.88645187) q[1];
sx q[1];
rz(-2.3019703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1621203) q[0];
sx q[0];
rz(-1.1651255) q[0];
sx q[0];
rz(-2.0590873) q[0];
rz(-2.1346852) q[2];
sx q[2];
rz(-2.208935) q[2];
sx q[2];
rz(1.2096805) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95123467) q[1];
sx q[1];
rz(-0.91427416) q[1];
sx q[1];
rz(1.7503529) q[1];
rz(0.59323816) q[3];
sx q[3];
rz(-2.7194033) q[3];
sx q[3];
rz(-1.4271023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7825369) q[2];
sx q[2];
rz(-1.3691207) q[2];
sx q[2];
rz(-3.0699406) q[2];
rz(-0.045529384) q[3];
sx q[3];
rz(-1.1979016) q[3];
sx q[3];
rz(0.45904747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.59938207) q[0];
sx q[0];
rz(-0.57294232) q[0];
sx q[0];
rz(-2.086916) q[0];
rz(-0.032546267) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(2.6194825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63037496) q[0];
sx q[0];
rz(-1.8763834) q[0];
sx q[0];
rz(0.60141984) q[0];
rz(-pi) q[1];
rz(2.1565151) q[2];
sx q[2];
rz(-0.92364468) q[2];
sx q[2];
rz(-0.070731846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13455277) q[1];
sx q[1];
rz(-1.947548) q[1];
sx q[1];
rz(-1.1621848) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9154112) q[3];
sx q[3];
rz(-1.0345248) q[3];
sx q[3];
rz(2.8482343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9944428) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(-1.7333376) q[2];
rz(-1.5857006) q[3];
sx q[3];
rz(-2.9097911) q[3];
sx q[3];
rz(2.2304992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1237352) q[0];
sx q[0];
rz(-1.8386848) q[0];
sx q[0];
rz(0.89717502) q[0];
rz(2.1625715) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(-2.945902) q[2];
sx q[2];
rz(-1.2441845) q[2];
sx q[2];
rz(1.2938538) q[2];
rz(-1.7579982) q[3];
sx q[3];
rz(-2.2254007) q[3];
sx q[3];
rz(2.7986516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
