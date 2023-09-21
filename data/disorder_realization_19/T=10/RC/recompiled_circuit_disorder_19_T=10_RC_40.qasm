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
rz(-1.6436613) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(1.210968) q[0];
rz(-1.6662007) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(0.22104095) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5421996) q[1];
sx q[1];
rz(-1.3560055) q[1];
sx q[1];
rz(0.7426803) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9707768) q[3];
sx q[3];
rz(-0.96066374) q[3];
sx q[3];
rz(2.4657472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92007414) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.39009538) q[0];
sx q[0];
rz(0.54498589) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(-2.3166336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2087814) q[0];
sx q[0];
rz(-2.6926846) q[0];
sx q[0];
rz(-2.0138028) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024436342) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(1.2791866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3927041) q[1];
sx q[1];
rz(-2.5671133) q[1];
sx q[1];
rz(2.9267465) q[1];
x q[2];
rz(-2.1373848) q[3];
sx q[3];
rz(-1.405218) q[3];
sx q[3];
rz(0.86629888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(-0.038539561) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8949216) q[0];
sx q[0];
rz(-1.5785494) q[0];
sx q[0];
rz(-0.00064108032) q[0];
rz(-0.23080319) q[2];
sx q[2];
rz(-1.421184) q[2];
sx q[2];
rz(2.9485867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5650428) q[1];
sx q[1];
rz(-2.0598663) q[1];
sx q[1];
rz(3.0719041) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58979804) q[3];
sx q[3];
rz(-2.0876948) q[3];
sx q[3];
rz(-0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(0.19515881) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0168266) q[0];
sx q[0];
rz(-2.989528) q[0];
sx q[0];
rz(0.91465871) q[0];
rz(2.6606576) q[2];
sx q[2];
rz(-2.0630884) q[2];
sx q[2];
rz(-2.7155657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8735503) q[1];
sx q[1];
rz(-0.51635427) q[1];
sx q[1];
rz(2.0622846) q[1];
rz(-2.5024662) q[3];
sx q[3];
rz(-1.5577321) q[3];
sx q[3];
rz(-2.5845205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.8738497) q[2];
rz(-1.0686482) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(2.7635014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75487126) q[0];
sx q[0];
rz(-2.5140962) q[0];
sx q[0];
rz(2.0621215) q[0];
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
rz(-2.0873449) q[1];
sx q[1];
rz(-1.9964881) q[1];
sx q[1];
rz(-2.737983) q[1];
x q[2];
rz(-0.37851815) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(-2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-2.2264218) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-2.2890131) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(2.728087) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9050161) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(2.828039) q[1];
rz(-pi) q[2];
rz(-0.69453199) q[3];
sx q[3];
rz(-1.4834705) q[3];
sx q[3];
rz(0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(1.7162011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6293684) q[0];
sx q[0];
rz(-0.47591305) q[0];
sx q[0];
rz(0.41643629) q[0];
rz(2.695589) q[2];
sx q[2];
rz(-1.6288695) q[2];
sx q[2];
rz(0.29689483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.094147625) q[1];
sx q[1];
rz(-2.3970251) q[1];
sx q[1];
rz(0.57569699) q[1];
x q[2];
rz(-1.5951049) q[3];
sx q[3];
rz(-0.50384249) q[3];
sx q[3];
rz(0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(-1.9308176) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(3.0650744) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589688) q[0];
sx q[0];
rz(-2.0730744) q[0];
sx q[0];
rz(1.2095905) q[0];
rz(-pi) q[1];
rz(0.92213995) q[2];
sx q[2];
rz(-1.6862009) q[2];
sx q[2];
rz(-0.61308544) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2462818) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(-1.4573216) q[1];
x q[2];
rz(-0.58952491) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(2.7713283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2032808) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(2.032062) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-2.1355656) q[0];
rz(2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.8267652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3599167) q[0];
sx q[0];
rz(-1.8209551) q[0];
sx q[0];
rz(-0.71235384) q[0];
x q[1];
rz(-0.33816955) q[2];
sx q[2];
rz(-0.38376946) q[2];
sx q[2];
rz(1.2107809) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4970376) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(-2.6130996) q[1];
x q[2];
rz(-1.5626838) q[3];
sx q[3];
rz(-2.3981961) q[3];
sx q[3];
rz(-0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(-1.5022087) q[1];
sx q[1];
rz(-0.39189664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.648801) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(-0.0097983629) q[0];
x q[1];
rz(0.91498615) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(-0.3273302) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3956986) q[1];
sx q[1];
rz(-0.77334009) q[1];
sx q[1];
rz(-1.3032773) q[1];
x q[2];
rz(-1.2519022) q[3];
sx q[3];
rz(-1.0341757) q[3];
sx q[3];
rz(-1.5981789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3863603) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-2.9677283) q[2];
sx q[2];
rz(-0.78730135) q[2];
sx q[2];
rz(0.80136328) q[2];
rz(-0.34655456) q[3];
sx q[3];
rz(-2.0934436) q[3];
sx q[3];
rz(2.2128076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];