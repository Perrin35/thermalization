OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1146381) q[0];
sx q[0];
rz(-1.4517598) q[0];
sx q[0];
rz(-0.64557689) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(-1.3728377) q[1];
sx q[1];
rz(1.6436613) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71292114) q[0];
sx q[0];
rz(-2.7550468) q[0];
sx q[0];
rz(1.1791694) q[0];
rz(0.94028084) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(-1.8688569) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9420535) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(2.8295423) q[1];
rz(-pi) q[2];
rz(-1.1708158) q[3];
sx q[3];
rz(-0.96066374) q[3];
sx q[3];
rz(0.67584544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(-2.3085964) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(0.82495904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4175407) q[0];
sx q[0];
rz(-1.1678956) q[0];
sx q[0];
rz(0.20362644) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1171563) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(1.2791866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0032012) q[1];
sx q[1];
rz(-2.1304641) q[1];
sx q[1];
rz(-1.4336587) q[1];
rz(-pi) q[2];
rz(-1.2689781) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(0.038539561) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32917133) q[0];
sx q[0];
rz(-3.1338131) q[0];
sx q[0];
rz(-1.6532941) q[0];
rz(-pi) q[1];
rz(2.9107895) q[2];
sx q[2];
rz(-1.421184) q[2];
sx q[2];
rz(2.9485867) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5765499) q[1];
sx q[1];
rz(-2.0598663) q[1];
sx q[1];
rz(-3.0719041) q[1];
rz(-pi) q[2];
rz(0.58979804) q[3];
sx q[3];
rz(-2.0876948) q[3];
sx q[3];
rz(-0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(1.3085261) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-0.19515881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965658) q[0];
sx q[0];
rz(-1.4782527) q[0];
sx q[0];
rz(1.4499614) q[0];
rz(-pi) q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(-1.7550215) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8198204) q[1];
sx q[1];
rz(-2.0211377) q[1];
sx q[1];
rz(-2.8798073) q[1];
x q[2];
rz(-1.5870729) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(-2.1181599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1068263) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(2.7635014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.970801) q[0];
sx q[0];
rz(-1.0266725) q[0];
sx q[0];
rz(2.8118954) q[0];
rz(-pi) q[1];
rz(2.4079608) q[2];
sx q[2];
rz(-2.1918104) q[2];
sx q[2];
rz(-0.48895198) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0873449) q[1];
sx q[1];
rz(-1.9964881) q[1];
sx q[1];
rz(2.737983) q[1];
rz(-pi) q[2];
rz(2.674621) q[3];
sx q[3];
rz(-1.7479556) q[3];
sx q[3];
rz(0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-0.37762541) q[0];
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
rz(0.24419063) q[0];
sx q[0];
rz(-1.4892329) q[0];
sx q[0];
rz(-0.32342644) q[0];
x q[1];
rz(0.50778163) q[2];
sx q[2];
rz(-2.0785329) q[2];
sx q[2];
rz(-0.43916647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8296216) q[1];
sx q[1];
rz(-1.8835856) q[1];
sx q[1];
rz(-1.4986067) q[1];
rz(-pi) q[2];
rz(2.4470607) q[3];
sx q[3];
rz(-1.4834705) q[3];
sx q[3];
rz(0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(-0.36703584) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4331626) q[0];
sx q[0];
rz(-1.7571974) q[0];
sx q[0];
rz(2.7011046) q[0];
x q[1];
rz(0.13397127) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(1.1531032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2196301) q[1];
sx q[1];
rz(-1.9486517) q[1];
sx q[1];
rz(-0.65803836) q[1];
rz(-pi) q[2];
rz(1.0670788) q[3];
sx q[3];
rz(-1.5825315) q[3];
sx q[3];
rz(-2.3624682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(-1.9308176) q[2];
rz(0.0018421729) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(-0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(0.076518245) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589688) q[0];
sx q[0];
rz(-2.0730744) q[0];
sx q[0];
rz(1.2095905) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1444507) q[2];
sx q[2];
rz(-0.92717294) q[2];
sx q[2];
rz(-2.2709536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4560495) q[1];
sx q[1];
rz(-1.4577663) q[1];
sx q[1];
rz(0.088940253) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(-1.0866144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.8267652) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1413404) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(1.2452026) q[0];
rz(-0.33816955) q[2];
sx q[2];
rz(-0.38376946) q[2];
sx q[2];
rz(-1.9308117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.41373738) q[1];
sx q[1];
rz(-1.1579885) q[1];
sx q[1];
rz(-0.84819838) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7018147) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058404) q[0];
sx q[0];
rz(-1.579111) q[0];
sx q[0];
rz(-1.0132829) q[0];
rz(-2.2266065) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(2.8142625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3956986) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(-1.3032773) q[1];
rz(-pi) q[2];
rz(-2.6565353) q[3];
sx q[3];
rz(-2.5254446) q[3];
sx q[3];
rz(0.96998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(1.7427313) q[2];
sx q[2];
rz(-2.3430765) q[2];
sx q[2];
rz(0.55745468) q[2];
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