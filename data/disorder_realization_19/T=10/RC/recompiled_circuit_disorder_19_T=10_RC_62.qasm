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
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(1.4979314) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6490373) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(1.9306246) q[0];
rz(-2.2013118) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(1.2727357) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.363417) q[1];
sx q[1];
rz(-2.2925804) q[1];
sx q[1];
rz(-1.8587106) q[1];
rz(-2.492339) q[3];
sx q[3];
rz(-1.2459727) q[3];
sx q[3];
rz(1.1326102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-0.34040889) q[2];
rz(2.3085964) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.814099) q[1];
sx q[1];
rz(-2.3166336) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23404113) q[0];
sx q[0];
rz(-1.7579161) q[0];
sx q[0];
rz(-1.1603111) q[0];
rz(-pi) q[1];
rz(-2.7386875) q[2];
sx q[2];
rz(-1.5612134) q[2];
sx q[2];
rz(2.8275037) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1383914) q[1];
sx q[1];
rz(-2.1304641) q[1];
sx q[1];
rz(-1.707934) q[1];
rz(-pi) q[2];
rz(1.8726146) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-0.17641243) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(-1.6167971) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-0.038539561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3241203) q[0];
sx q[0];
rz(-1.5701553) q[0];
sx q[0];
rz(1.5630432) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5589587) q[2];
sx q[2];
rz(-2.8672672) q[2];
sx q[2];
rz(2.3290616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5765499) q[1];
sx q[1];
rz(-2.0598663) q[1];
sx q[1];
rz(-0.069688571) q[1];
rz(-pi) q[2];
rz(2.3452957) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(-2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38347605) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-0.34969774) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(0.19515881) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12476607) q[0];
sx q[0];
rz(-2.989528) q[0];
sx q[0];
rz(2.2269339) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1148557) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(-1.3865711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2745167) q[1];
sx q[1];
rz(-1.80596) q[1];
sx q[1];
rz(2.0348674) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1196932) q[3];
sx q[3];
rz(-2.5023513) q[3];
sx q[3];
rz(-0.99614776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.267743) q[2];
rz(2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-2.7635014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.224684) q[0];
sx q[0];
rz(-1.8514669) q[0];
sx q[0];
rz(-1.0018437) q[0];
rz(-pi) q[1];
rz(2.4079608) q[2];
sx q[2];
rz(-2.1918104) q[2];
sx q[2];
rz(2.6526407) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8564344) q[1];
sx q[1];
rz(-0.57796961) q[1];
sx q[1];
rz(2.2846089) q[1];
rz(2.7630745) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(-2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-2.7639672) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897402) q[0];
sx q[0];
rz(-1.6523598) q[0];
sx q[0];
rz(0.32342644) q[0];
rz(1.0038063) q[2];
sx q[2];
rz(-1.1319455) q[2];
sx q[2];
rz(1.7457419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8296216) q[1];
sx q[1];
rz(-1.258007) q[1];
sx q[1];
rz(-1.4986067) q[1];
rz(-3.0056474) q[3];
sx q[3];
rz(-2.4424986) q[3];
sx q[3];
rz(1.618315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6043828) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-2.3510695) q[2];
rz(2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(-2.7745568) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.7162011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910949) q[0];
sx q[0];
rz(-1.1384581) q[0];
sx q[0];
rz(1.3652486) q[0];
rz(-pi) q[1];
rz(-0.44600365) q[2];
sx q[2];
rz(-1.6288695) q[2];
sx q[2];
rz(0.29689483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.047445) q[1];
sx q[1];
rz(-2.3970251) q[1];
sx q[1];
rz(-2.5658957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5951049) q[3];
sx q[3];
rz(-0.50384249) q[3];
sx q[3];
rz(0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84412557) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233326) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(0.53091913) q[0];
rz(-1.3812177) q[2];
sx q[2];
rz(-0.65738064) q[2];
sx q[2];
rz(-2.0331403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2462818) q[1];
sx q[1];
rz(-1.4824251) q[1];
sx q[1];
rz(1.6842711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30267834) q[3];
sx q[3];
rz(-0.61121002) q[3];
sx q[3];
rz(-1.4509033) q[3];
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
rz(-2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(-2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(1.0060271) q[0];
rz(-0.26793119) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.3148274) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3599167) q[0];
sx q[0];
rz(-1.3206375) q[0];
sx q[0];
rz(-2.4292388) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7039653) q[2];
sx q[2];
rz(-1.2097934) q[2];
sx q[2];
rz(1.5683057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6445551) q[1];
sx q[1];
rz(-2.2215448) q[1];
sx q[1];
rz(0.52849309) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1341348) q[3];
sx q[3];
rz(-0.82743001) q[3];
sx q[3];
rz(0.5479365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(-2.8961199) q[2];
rz(2.7108575) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(-2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-2.749696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083188699) q[0];
sx q[0];
rz(-1.5624816) q[0];
sx q[0];
rz(-2.1283098) q[0];
rz(-pi) q[1];
rz(-2.2266065) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(2.8142625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74589409) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.3032773) q[1];
rz(-pi) q[2];
rz(0.55962555) q[3];
sx q[3];
rz(-1.2979753) q[3];
sx q[3];
rz(-2.9469957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(-1.8023087) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(1.3988613) q[2];
sx q[2];
rz(-0.79851616) q[2];
sx q[2];
rz(-2.584138) q[2];
rz(-2.1035351) q[3];
sx q[3];
rz(-0.61809117) q[3];
sx q[3];
rz(1.5865159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
