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
rz(-0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(1.9306246) q[0];
rz(-0.94028084) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(-1.2727357) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5421996) q[1];
sx q[1];
rz(-1.7855872) q[1];
sx q[1];
rz(0.7426803) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.492339) q[3];
sx q[3];
rz(-1.2459727) q[3];
sx q[3];
rz(1.1326102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(2.3085964) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.32749367) q[1];
sx q[1];
rz(-0.82495904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.724052) q[0];
sx q[0];
rz(-1.9736971) q[0];
sx q[0];
rz(-2.9379662) q[0];
rz(2.7386875) q[2];
sx q[2];
rz(-1.5803792) q[2];
sx q[2];
rz(2.8275037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7488885) q[1];
sx q[1];
rz(-0.57447937) q[1];
sx q[1];
rz(-2.9267465) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2689781) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(-2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(-0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(0.038539561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8949216) q[0];
sx q[0];
rz(-1.5785494) q[0];
sx q[0];
rz(3.1409516) q[0];
rz(0.58263393) q[2];
sx q[2];
rz(-0.27432549) q[2];
sx q[2];
rz(2.3290616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5650428) q[1];
sx q[1];
rz(-2.0598663) q[1];
sx q[1];
rz(-0.069688571) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79629691) q[3];
sx q[3];
rz(-2.3781804) q[3];
sx q[3];
rz(2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58275756) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(0.91252404) q[2];
rz(1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(0.19515881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12476607) q[0];
sx q[0];
rz(-2.989528) q[0];
sx q[0];
rz(2.2269339) q[0];
rz(-pi) q[1];
rz(-0.85907016) q[2];
sx q[2];
rz(-0.67407437) q[2];
sx q[2];
rz(-2.7328343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.86707592) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(1.1067252) q[1];
rz(-pi) q[2];
rz(-3.1196932) q[3];
sx q[3];
rz(-2.5023513) q[3];
sx q[3];
rz(0.99614776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.8738497) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-0.1396133) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-2.7635014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.970801) q[0];
sx q[0];
rz(-1.0266725) q[0];
sx q[0];
rz(-0.32969726) q[0];
x q[1];
rz(-2.3374704) q[2];
sx q[2];
rz(-0.99493775) q[2];
sx q[2];
rz(-2.5428307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34199076) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(1.1127383) q[1];
rz(-0.37851815) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(-2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(2.2536229) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-0.37762541) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-0.91517085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24419063) q[0];
sx q[0];
rz(-1.4892329) q[0];
sx q[0];
rz(2.8181662) q[0];
rz(-pi) q[1];
rz(-0.50778163) q[2];
sx q[2];
rz(-1.0630597) q[2];
sx q[2];
rz(-0.43916647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23657654) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(-2.828039) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0056474) q[3];
sx q[3];
rz(-0.69909401) q[3];
sx q[3];
rz(1.5232777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(-2.524232) q[3];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(-0.36703584) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.7162011) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7084301) q[0];
sx q[0];
rz(-1.7571974) q[0];
sx q[0];
rz(0.44048803) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13397127) q[2];
sx q[2];
rz(-0.44951648) q[2];
sx q[2];
rz(-1.1531032) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5126309) q[1];
sx q[1];
rz(-0.96620622) q[1];
sx q[1];
rz(1.1058034) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5951049) q[3];
sx q[3];
rz(-0.50384249) q[3];
sx q[3];
rz(-0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1620862) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(-3.1397505) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(2.4842998) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
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
rz(1.9936179) q[0];
sx q[0];
rz(-2.5320842) q[0];
sx q[0];
rz(0.57172914) q[0];
rz(-1.760375) q[2];
sx q[2];
rz(-2.484212) q[2];
sx q[2];
rz(1.1084523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2462818) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(-1.6842711) q[1];
rz(-pi) q[2];
rz(-2.5520677) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(0.37026438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9922441) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(2.1355656) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(1.3148274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00025230322) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(-1.2452026) q[0];
x q[1];
rz(1.4376274) q[2];
sx q[2];
rz(-1.2097934) q[2];
sx q[2];
rz(-1.5732869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73003522) q[1];
sx q[1];
rz(-2.3282603) q[1];
sx q[1];
rz(0.98585339) q[1];
rz(-pi) q[2];
rz(-1.5626838) q[3];
sx q[3];
rz(-0.74339657) q[3];
sx q[3];
rz(0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(2.7108575) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(-2.664264) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4196639) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(2.749696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083188699) q[0];
sx q[0];
rz(-1.5624816) q[0];
sx q[0];
rz(-1.0132829) q[0];
rz(-pi) q[1];
rz(-0.91498615) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(-2.8142625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1230159) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(-0.8155483) q[1];
x q[2];
rz(1.8896905) q[3];
sx q[3];
rz(-2.1074169) q[3];
sx q[3];
rz(1.5981789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(2.3618868) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(0.34655456) q[3];
sx q[3];
rz(-1.0481491) q[3];
sx q[3];
rz(-0.92878503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];