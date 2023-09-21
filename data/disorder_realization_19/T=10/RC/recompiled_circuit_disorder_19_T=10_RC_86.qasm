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
rz(0.37880701) q[1];
sx q[1];
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71292114) q[0];
sx q[0];
rz(-2.7550468) q[0];
sx q[0];
rz(1.9624233) q[0];
x q[1];
rz(-1.6662007) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(-2.9205517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.363417) q[1];
sx q[1];
rz(-0.84901224) q[1];
sx q[1];
rz(1.282882) q[1];
rz(2.6333991) q[3];
sx q[3];
rz(-2.4262706) q[3];
sx q[3];
rz(-0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(0.54498589) q[0];
rz(-0.90822059) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(2.3166336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2087814) q[0];
sx q[0];
rz(-0.44890807) q[0];
sx q[0];
rz(-2.0138028) q[0];
rz(3.1171563) q[2];
sx q[2];
rz(-2.7385798) q[2];
sx q[2];
rz(-1.8624061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7488885) q[1];
sx q[1];
rz(-2.5671133) q[1];
sx q[1];
rz(-0.21484612) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0042079) q[3];
sx q[3];
rz(-1.405218) q[3];
sx q[3];
rz(-2.2752938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(0.17641243) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(0.46974716) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(-0.038539561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24667106) q[0];
sx q[0];
rz(-1.5630432) q[0];
sx q[0];
rz(-0.00064108032) q[0];
x q[1];
rz(-2.9107895) q[2];
sx q[2];
rz(-1.421184) q[2];
sx q[2];
rz(0.19300592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5650428) q[1];
sx q[1];
rz(-1.0817263) q[1];
sx q[1];
rz(-3.0719041) q[1];
rz(-pi) q[2];
rz(-0.58979804) q[3];
sx q[3];
rz(-2.0876948) q[3];
sx q[3];
rz(-3.0977109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58275756) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(-0.91252404) q[2];
rz(-1.3085261) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(1.4413888) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168266) q[0];
sx q[0];
rz(-2.989528) q[0];
sx q[0];
rz(-0.91465871) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.9907021) q[2];
sx q[2];
rz(1.3865711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2745167) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(-1.1067252) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5545198) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(-1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(2.0647678) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(0.37809125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.970801) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(-2.8118954) q[0];
rz(2.4079608) q[2];
sx q[2];
rz(-0.94978226) q[2];
sx q[2];
rz(-2.6526407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34199076) q[1];
sx q[1];
rz(-1.9365891) q[1];
sx q[1];
rz(-1.1127383) q[1];
rz(2.7630745) q[3];
sx q[3];
rz(-2.6444728) q[3];
sx q[3];
rz(2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-2.7639672) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(2.2264218) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24419063) q[0];
sx q[0];
rz(-1.6523598) q[0];
sx q[0];
rz(-0.32342644) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0038063) q[2];
sx q[2];
rz(-2.0096471) q[2];
sx q[2];
rz(1.3958508) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5988016) q[1];
sx q[1];
rz(-0.32074499) q[1];
sx q[1];
rz(0.21943211) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1359453) q[3];
sx q[3];
rz(-2.4424986) q[3];
sx q[3];
rz(1.5232777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(-2.8516155) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050497748) q[0];
sx q[0];
rz(-1.1384581) q[0];
sx q[0];
rz(-1.3652486) q[0];
rz(-2.695589) q[2];
sx q[2];
rz(-1.6288695) q[2];
sx q[2];
rz(-0.29689483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9219626) q[1];
sx q[1];
rz(-1.9486517) q[1];
sx q[1];
rz(0.65803836) q[1];
rz(-pi) q[2];
rz(3.1281934) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(-2.3434533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(3.0650744) q[0];
rz(-0.67529768) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(-2.3892367) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589688) q[0];
sx q[0];
rz(-1.0685182) q[0];
sx q[0];
rz(1.9320022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.760375) q[2];
sx q[2];
rz(-2.484212) q[2];
sx q[2];
rz(2.0331403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(-1.4573216) q[1];
rz(-0.58952491) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(2.7713283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2032808) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(0.00014649815) q[2];
rz(2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(2.1355656) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.3148274) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49004236) q[0];
sx q[0];
rz(-2.393912) q[0];
sx q[0];
rz(-0.3726532) q[0];
x q[1];
rz(-2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(-0.044791128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6445551) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(2.6130996) q[1];
rz(-pi) q[2];
rz(1.5789088) q[3];
sx q[3];
rz(-0.74339657) q[3];
sx q[3];
rz(-2.5826366) q[3];
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
rz(-2.8961199) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
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
rz(0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-2.749696) q[1];
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
rz(0.91498615) q[2];
sx q[2];
rz(-1.4631541) q[2];
sx q[2];
rz(0.3273302) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3956986) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(-1.3032773) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8896905) q[3];
sx q[3];
rz(-2.1074169) q[3];
sx q[3];
rz(1.5434138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(-1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.17386439) q[2];
sx q[2];
rz(-2.3542913) q[2];
sx q[2];
rz(-2.3402294) q[2];
rz(-1.0212392) q[3];
sx q[3];
rz(-1.869535) q[3];
sx q[3];
rz(0.46365999) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
