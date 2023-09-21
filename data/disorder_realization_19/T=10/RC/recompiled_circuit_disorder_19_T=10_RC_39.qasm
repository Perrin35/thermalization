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
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0092702) q[0];
sx q[0];
rz(-1.2148804) q[0];
sx q[0];
rz(-2.9874671) q[0];
x q[1];
rz(-0.94028084) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(-1.2727357) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77817569) q[1];
sx q[1];
rz(-0.84901224) q[1];
sx q[1];
rz(1.282882) q[1];
x q[2];
rz(-1.1708158) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(-0.67584544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.67265636) q[0];
sx q[0];
rz(-2.7514973) q[0];
sx q[0];
rz(-0.54498589) q[0];
rz(-0.90822059) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(-2.3166336) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.724052) q[0];
sx q[0];
rz(-1.9736971) q[0];
sx q[0];
rz(2.9379662) q[0];
rz(-0.40290515) q[2];
sx q[2];
rz(-1.5612134) q[2];
sx q[2];
rz(0.314089) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7488885) q[1];
sx q[1];
rz(-2.5671133) q[1];
sx q[1];
rz(2.9267465) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9460658) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(-0.60002458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(0.17641243) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(3.1030531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32917133) q[0];
sx q[0];
rz(-3.1338131) q[0];
sx q[0];
rz(1.6532941) q[0];
rz(-pi) q[1];
rz(1.4171717) q[2];
sx q[2];
rz(-1.7989752) q[2];
sx q[2];
rz(-1.7287901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4290532) q[1];
sx q[1];
rz(-0.4936115) q[1];
sx q[1];
rz(-1.7008971) q[1];
rz(-pi) q[2];
rz(2.3452957) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(-2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5588351) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-0.19515881) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.450481) q[0];
sx q[0];
rz(-3.0483732) q[0];
rz(2.2825225) q[2];
sx q[2];
rz(-2.4675183) q[2];
sx q[2];
rz(2.7328343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8735503) q[1];
sx q[1];
rz(-0.51635427) q[1];
sx q[1];
rz(1.079308) q[1];
x q[2];
rz(-1.5545198) q[3];
sx q[3];
rz(-0.93173325) q[3];
sx q[3];
rz(-2.1181599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.8738497) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(2.0647678) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(0.37809125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9169086) q[0];
sx q[0];
rz(-1.2901257) q[0];
sx q[0];
rz(1.0018437) q[0];
x q[1];
rz(0.73363186) q[2];
sx q[2];
rz(-0.94978226) q[2];
sx q[2];
rz(-0.48895198) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0542478) q[1];
sx q[1];
rz(-1.9964881) q[1];
sx q[1];
rz(2.737983) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.674621) q[3];
sx q[3];
rz(-1.7479556) q[3];
sx q[3];
rz(2.5339479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
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
rz(-2.2264218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422896) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(-1.4847941) q[0];
rz(-1.0038063) q[2];
sx q[2];
rz(-1.1319455) q[2];
sx q[2];
rz(1.3958508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23657654) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(-0.31355365) q[1];
rz(-pi) q[2];
rz(2.4470607) q[3];
sx q[3];
rz(-1.6581222) q[3];
sx q[3];
rz(-0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.4253915) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050497748) q[0];
sx q[0];
rz(-2.0031345) q[0];
sx q[0];
rz(1.3652486) q[0];
rz(-pi) q[1];
rz(-2.695589) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(-2.8446978) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2196301) q[1];
sx q[1];
rz(-1.1929409) q[1];
sx q[1];
rz(2.4835543) q[1];
rz(1.5951049) q[3];
sx q[3];
rz(-2.6377502) q[3];
sx q[3];
rz(-2.3712096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84412557) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9936179) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(2.5698635) q[0];
rz(2.997142) q[2];
sx q[2];
rz(-2.2144197) q[2];
sx q[2];
rz(-0.87063906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4560495) q[1];
sx q[1];
rz(-1.4577663) q[1];
sx q[1];
rz(-3.0526524) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(2.0549783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2032808) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(1.1095307) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(-2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(1.0060271) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(-1.8267652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3599167) q[0];
sx q[0];
rz(-1.8209551) q[0];
sx q[0];
rz(2.4292388) q[0];
x q[1];
rz(1.7039653) q[2];
sx q[2];
rz(-1.2097934) q[2];
sx q[2];
rz(1.5732869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4970376) q[1];
sx q[1];
rz(-2.2215448) q[1];
sx q[1];
rz(2.6130996) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0074578961) q[3];
sx q[3];
rz(-0.82743001) q[3];
sx q[3];
rz(2.5936562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(2.8961199) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(-0.47732863) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7219287) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(0.28668177) q[0];
rz(2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(0.39189664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.648801) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(-0.0097983629) q[0];
rz(-pi) q[1];
rz(-3.0060843) q[2];
sx q[2];
rz(-2.2221609) q[2];
sx q[2];
rz(-1.9806005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3956986) q[1];
sx q[1];
rz(-0.77334009) q[1];
sx q[1];
rz(-1.3032773) q[1];
x q[2];
rz(-0.48505731) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(0.96998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.9840476) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(-1.339284) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(0.77970589) q[2];
sx q[2];
rz(-1.4479326) q[2];
sx q[2];
rz(-0.89276199) q[2];
rz(2.1035351) q[3];
sx q[3];
rz(-2.5235015) q[3];
sx q[3];
rz(-1.5550767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
