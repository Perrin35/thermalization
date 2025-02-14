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
rz(0.99474466) q[0];
sx q[0];
rz(-2.5543307) q[0];
sx q[0];
rz(0.62556148) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(-2.151139) q[1];
sx q[1];
rz(0.78392309) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44255689) q[0];
sx q[0];
rz(-2.7074506) q[0];
sx q[0];
rz(-0.0084369466) q[0];
rz(-1.4874081) q[2];
sx q[2];
rz(-1.7094128) q[2];
sx q[2];
rz(1.0178138) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16456024) q[1];
sx q[1];
rz(-0.98233084) q[1];
sx q[1];
rz(-0.87314897) q[1];
rz(2.3638681) q[3];
sx q[3];
rz(-1.5081662) q[3];
sx q[3];
rz(2.3304813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94773942) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(0.6673153) q[2];
rz(-0.31467485) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(0.92710322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2941403) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(-0.2997998) q[0];
rz(2.1369797) q[1];
sx q[1];
rz(-0.71669465) q[1];
sx q[1];
rz(0.85003781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93545234) q[0];
sx q[0];
rz(-0.61203921) q[0];
sx q[0];
rz(-1.7382337) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17053594) q[2];
sx q[2];
rz(-1.6493622) q[2];
sx q[2];
rz(-2.6002392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5483173) q[1];
sx q[1];
rz(-1.9276345) q[1];
sx q[1];
rz(0.22290454) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1095457) q[3];
sx q[3];
rz(-0.40478313) q[3];
sx q[3];
rz(1.3758152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(0.66298318) q[2];
rz(-0.7524544) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.74333423) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(0.67984003) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-2.5847021) q[1];
sx q[1];
rz(-2.360875) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74405438) q[0];
sx q[0];
rz(-0.40073943) q[0];
sx q[0];
rz(-1.5903872) q[0];
x q[1];
rz(-0.0038006101) q[2];
sx q[2];
rz(-1.6569306) q[2];
sx q[2];
rz(-1.116556) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6284077) q[1];
sx q[1];
rz(-2.3806913) q[1];
sx q[1];
rz(0.55879467) q[1];
x q[2];
rz(-2.3060477) q[3];
sx q[3];
rz(-1.4073331) q[3];
sx q[3];
rz(-0.36629656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22689936) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(-0.66960382) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-2.1293631) q[3];
sx q[3];
rz(2.3646234) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25903073) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(2.2341109) q[0];
rz(-2.5471845) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(2.0118735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8639899) q[0];
sx q[0];
rz(-0.66501319) q[0];
sx q[0];
rz(-1.8825085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3200892) q[2];
sx q[2];
rz(-0.54566979) q[2];
sx q[2];
rz(0.16599338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4817244) q[1];
sx q[1];
rz(-2.4818083) q[1];
sx q[1];
rz(-0.45141141) q[1];
rz(1.8640737) q[3];
sx q[3];
rz(-1.8403111) q[3];
sx q[3];
rz(0.5130918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(-3.1224871) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(-3.0066971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.8371209) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(0.16803148) q[0];
rz(-1.9145603) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(-2.8438445) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79556161) q[0];
sx q[0];
rz(-0.77036422) q[0];
sx q[0];
rz(0.7839696) q[0];
rz(-2.1476168) q[2];
sx q[2];
rz(-1.0712396) q[2];
sx q[2];
rz(-3.1021038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49714607) q[1];
sx q[1];
rz(-2.9342817) q[1];
sx q[1];
rz(0.12122341) q[1];
x q[2];
rz(-0.00068126531) q[3];
sx q[3];
rz(-1.7295803) q[3];
sx q[3];
rz(0.53631594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(-2.3330073) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(1.3396858) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488778) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(0.46522796) q[0];
rz(2.6564927) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(-0.088265158) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53118491) q[0];
sx q[0];
rz(-2.0636798) q[0];
sx q[0];
rz(0.51778527) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2997238) q[2];
sx q[2];
rz(-1.5612023) q[2];
sx q[2];
rz(-0.95649715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0572964) q[1];
sx q[1];
rz(-0.7727333) q[1];
sx q[1];
rz(0.23717662) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8507455) q[3];
sx q[3];
rz(-0.29971443) q[3];
sx q[3];
rz(-2.7951473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(-0.49760231) q[2];
rz(2.6615214) q[3];
sx q[3];
rz(-0.92104715) q[3];
sx q[3];
rz(-0.068643071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29480544) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(1.4303327) q[0];
rz(1.3149186) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(-1.5865145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1436227) q[0];
sx q[0];
rz(-0.41150948) q[0];
sx q[0];
rz(1.2689778) q[0];
rz(0.58760037) q[2];
sx q[2];
rz(-0.94600224) q[2];
sx q[2];
rz(1.8163565) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58007121) q[1];
sx q[1];
rz(-2.5513756) q[1];
sx q[1];
rz(-1.3886222) q[1];
x q[2];
rz(2.2022543) q[3];
sx q[3];
rz(-1.5087481) q[3];
sx q[3];
rz(0.59988672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2985208) q[2];
sx q[2];
rz(-0.71121794) q[2];
sx q[2];
rz(2.2141875) q[2];
rz(-1.983042) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(0.096175171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(2.6655777) q[0];
rz(-0.99126518) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(0.083866619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80591969) q[0];
sx q[0];
rz(-2.0562547) q[0];
sx q[0];
rz(-0.40747633) q[0];
rz(-pi) q[1];
rz(0.4617347) q[2];
sx q[2];
rz(-1.8991915) q[2];
sx q[2];
rz(-2.6685614) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8955651) q[1];
sx q[1];
rz(-1.152134) q[1];
sx q[1];
rz(-3.0484867) q[1];
x q[2];
rz(1.3478892) q[3];
sx q[3];
rz(-0.37316868) q[3];
sx q[3];
rz(2.580433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91741651) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(0.64845294) q[2];
rz(2.2367541) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(-0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7998578) q[0];
sx q[0];
rz(-2.4315727) q[0];
sx q[0];
rz(-2.572686) q[0];
rz(-0.83373249) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(2.1368829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2434655) q[0];
sx q[0];
rz(-1.7633738) q[0];
sx q[0];
rz(0.63611998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0332074) q[2];
sx q[2];
rz(-1.7694607) q[2];
sx q[2];
rz(0.63420701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9819599) q[1];
sx q[1];
rz(-1.6430322) q[1];
sx q[1];
rz(-0.87149974) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4887722) q[3];
sx q[3];
rz(-0.98482212) q[3];
sx q[3];
rz(-0.88912933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4363165) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(-2.7936068) q[2];
rz(-2.3157388) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(-2.5364449) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572877) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(-2.9344946) q[0];
rz(0.53889489) q[1];
sx q[1];
rz(-0.62108827) q[1];
sx q[1];
rz(-0.60212392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3054863) q[0];
sx q[0];
rz(-1.0998816) q[0];
sx q[0];
rz(0.54022809) q[0];
rz(-pi) q[1];
rz(-0.58790284) q[2];
sx q[2];
rz(-1.9155115) q[2];
sx q[2];
rz(1.4522788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6844897) q[1];
sx q[1];
rz(-0.41375638) q[1];
sx q[1];
rz(1.3828418) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35134985) q[3];
sx q[3];
rz(-0.77882871) q[3];
sx q[3];
rz(0.035127775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78486097) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(1.3447364) q[2];
rz(-0.52062672) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(0.80210137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(2.3583892) q[1];
sx q[1];
rz(-1.7192817) q[1];
sx q[1];
rz(-1.5789938) q[1];
rz(-2.4827847) q[2];
sx q[2];
rz(-0.42429069) q[2];
sx q[2];
rz(-1.9784603) q[2];
rz(0.24012031) q[3];
sx q[3];
rz(-1.5821531) q[3];
sx q[3];
rz(1.5021363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
