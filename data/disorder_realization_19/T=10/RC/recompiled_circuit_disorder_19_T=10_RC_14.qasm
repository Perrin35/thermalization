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
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71292114) q[0];
sx q[0];
rz(-2.7550468) q[0];
sx q[0];
rz(-1.9624233) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0718578) q[2];
sx q[2];
rz(-2.2001534) q[2];
sx q[2];
rz(-0.33915181) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.363417) q[1];
sx q[1];
rz(-2.2925804) q[1];
sx q[1];
rz(-1.8587106) q[1];
rz(-0.64925361) q[3];
sx q[3];
rz(-1.89562) q[3];
sx q[3];
rz(-2.0089825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(2.3085964) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(0.54498589) q[0];
rz(0.90822059) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-2.3166336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.724052) q[0];
sx q[0];
rz(-1.1678956) q[0];
sx q[0];
rz(2.9379662) q[0];
rz(3.1171563) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(1.8624061) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3927041) q[1];
sx q[1];
rz(-2.5671133) q[1];
sx q[1];
rz(2.9267465) q[1];
rz(2.9460658) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(-2.5415681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(3.0107064) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(3.1030531) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8949216) q[0];
sx q[0];
rz(-1.5785494) q[0];
sx q[0];
rz(-0.00064108032) q[0];
rz(2.9107895) q[2];
sx q[2];
rz(-1.7204086) q[2];
sx q[2];
rz(0.19300592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1145648) q[1];
sx q[1];
rz(-1.5092883) q[1];
sx q[1];
rz(2.0608749) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97088082) q[3];
sx q[3];
rz(-1.0661134) q[3];
sx q[3];
rz(-1.8463299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(0.91252404) q[2];
rz(1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.450481) q[0];
sx q[0];
rz(0.093219482) q[0];
rz(-pi) q[1];
rz(-2.6606576) q[2];
sx q[2];
rz(-2.0630884) q[2];
sx q[2];
rz(-0.42602691) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3217722) q[1];
sx q[1];
rz(-1.1204549) q[1];
sx q[1];
rz(0.26178534) q[1];
rz(-1.5870729) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.267743) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-0.37809125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3867214) q[0];
sx q[0];
rz(-0.62749642) q[0];
sx q[0];
rz(-1.0794712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80412229) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(-2.5428307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7996019) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(-1.1127383) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37851815) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(0.62687031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(0.97638431) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75509214) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(-0.37762541) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5650425) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(0.25175005) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1377863) q[2];
sx q[2];
rz(-1.1319455) q[2];
sx q[2];
rz(1.7457419) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9050161) q[1];
sx q[1];
rz(-1.502115) q[1];
sx q[1];
rz(-2.828039) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(0.79052314) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73615605) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(-1.4253915) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51222425) q[0];
sx q[0];
rz(-2.6656796) q[0];
sx q[0];
rz(0.41643629) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6351498) q[2];
sx q[2];
rz(-2.0159941) q[2];
sx q[2];
rz(1.839947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62896171) q[1];
sx q[1];
rz(-2.1753864) q[1];
sx q[1];
rz(-1.1058034) q[1];
rz(-0.013399259) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(0.79813938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90826666) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(-2.6106735) q[0];
rz(2.997142) q[2];
sx q[2];
rz(-2.2144197) q[2];
sx q[2];
rz(2.2709536) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(1.6842711) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.3148274) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7816759) q[0];
sx q[0];
rz(-1.8209551) q[0];
sx q[0];
rz(-0.71235384) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7776412) q[2];
sx q[2];
rz(-1.4462573) q[2];
sx q[2];
rz(3.0968015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6445551) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(-2.6130996) q[1];
x q[2];
rz(1.5789088) q[3];
sx q[3];
rz(-2.3981961) q[3];
sx q[3];
rz(-0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-2.7108575) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(-0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(0.39189664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.648801) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(-3.1317943) q[0];
rz(3.0060843) q[2];
sx q[2];
rz(-2.2221609) q[2];
sx q[2];
rz(-1.1609921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7614903) q[1];
sx q[1];
rz(-0.83161608) q[1];
sx q[1];
rz(0.2525316) q[1];
rz(-2.6565353) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(2.5907497) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(-1.339284) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(1.3988613) q[2];
sx q[2];
rz(-0.79851616) q[2];
sx q[2];
rz(-2.584138) q[2];
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