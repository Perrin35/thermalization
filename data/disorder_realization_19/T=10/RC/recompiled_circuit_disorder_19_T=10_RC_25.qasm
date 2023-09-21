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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71292114) q[0];
sx q[0];
rz(-0.38654583) q[0];
sx q[0];
rz(-1.9624233) q[0];
x q[1];
rz(2.2013118) q[2];
sx q[2];
rz(-1.5144381) q[2];
sx q[2];
rz(-1.8688569) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1995391) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(0.3120504) q[1];
x q[2];
rz(-1.1708158) q[3];
sx q[3];
rz(-0.96066374) q[3];
sx q[3];
rz(-2.4657472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(0.82495904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2087814) q[0];
sx q[0];
rz(-0.44890807) q[0];
sx q[0];
rz(-2.0138028) q[0];
rz(-pi) q[1];
rz(0.024436342) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(1.2791866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0032012) q[1];
sx q[1];
rz(-2.1304641) q[1];
sx q[1];
rz(1.4336587) q[1];
x q[2];
rz(-1.2689781) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(0.17641243) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(0.038539561) q[1];
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
rz(-2.5589587) q[2];
sx q[2];
rz(-2.8672672) q[2];
sx q[2];
rz(-2.3290616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4290532) q[1];
sx q[1];
rz(-2.6479811) q[1];
sx q[1];
rz(1.4406956) q[1];
rz(0.58979804) q[3];
sx q[3];
rz(-1.0538978) q[3];
sx q[3];
rz(0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-0.91252404) q[2];
rz(-1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168266) q[0];
sx q[0];
rz(-0.15206465) q[0];
sx q[0];
rz(-2.2269339) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(1.3865711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8735503) q[1];
sx q[1];
rz(-0.51635427) q[1];
sx q[1];
rz(-1.079308) q[1];
rz(-0.021899453) q[3];
sx q[3];
rz(-2.5023513) q[3];
sx q[3];
rz(2.1454449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9106456) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1068263) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(2.7635014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.224684) q[0];
sx q[0];
rz(-1.8514669) q[0];
sx q[0];
rz(-2.1397489) q[0];
x q[1];
rz(-0.81850448) q[2];
sx q[2];
rz(-2.2193925) q[2];
sx q[2];
rz(-1.4865781) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34199076) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(2.0288543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3729172) q[3];
sx q[3];
rz(-2.0298925) q[3];
sx q[3];
rz(-1.0517694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2698764) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(2.1652083) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-0.91517085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8422896) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(1.4847941) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2890131) q[2];
sx q[2];
rz(-2.4396509) q[2];
sx q[2];
rz(-2.728087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.542791) q[1];
sx q[1];
rz(-2.8208477) q[1];
sx q[1];
rz(2.9221605) q[1];
x q[2];
rz(-2.4470607) q[3];
sx q[3];
rz(-1.6581222) q[3];
sx q[3];
rz(-2.9897523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-2.3510695) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4331626) q[0];
sx q[0];
rz(-1.3843952) q[0];
sx q[0];
rz(-2.7011046) q[0];
x q[1];
rz(1.6351498) q[2];
sx q[2];
rz(-1.1255985) q[2];
sx q[2];
rz(-1.839947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.094147625) q[1];
sx q[1];
rz(-0.74456753) q[1];
sx q[1];
rz(-2.5658957) q[1];
rz(-pi) q[2];
rz(2.0745139) q[3];
sx q[3];
rz(-1.5825315) q[3];
sx q[3];
rz(2.3624682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(-1.210775) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(-2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(-2.3892367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233326) q[0];
sx q[0];
rz(-1.255863) q[0];
sx q[0];
rz(-2.6106735) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.760375) q[2];
sx q[2];
rz(-2.484212) q[2];
sx q[2];
rz(1.1084523) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1250549) q[1];
sx q[1];
rz(-0.14370951) q[1];
sx q[1];
rz(-0.90682604) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30267834) q[3];
sx q[3];
rz(-0.61121002) q[3];
sx q[3];
rz(-1.6906893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2032808) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(-1.8267652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00025230322) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(-1.89639) q[0];
x q[1];
rz(-1.7039653) q[2];
sx q[2];
rz(-1.2097934) q[2];
sx q[2];
rz(-1.5732869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73003522) q[1];
sx q[1];
rz(-2.3282603) q[1];
sx q[1];
rz(0.98585339) q[1];
rz(-pi) q[2];
rz(-2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(-1.0178125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-2.749696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4742728) q[0];
sx q[0];
rz(-0.55756888) q[0];
sx q[0];
rz(-1.5550818) q[0];
rz(-3.0060843) q[2];
sx q[2];
rz(-0.91943179) q[2];
sx q[2];
rz(1.9806005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1230159) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(-2.3260444) q[1];
rz(-2.5819671) q[3];
sx q[3];
rz(-1.8436173) q[3];
sx q[3];
rz(2.9469957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.9840476) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(0.77970589) q[2];
sx q[2];
rz(-1.4479326) q[2];
sx q[2];
rz(-0.89276199) q[2];
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
