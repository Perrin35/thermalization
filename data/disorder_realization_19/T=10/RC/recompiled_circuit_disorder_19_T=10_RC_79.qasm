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
rz(2.4960158) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(-1.6436613) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49255532) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(1.210968) q[0];
rz(0.94028084) q[2];
sx q[2];
rz(-1.5144381) q[2];
sx q[2];
rz(-1.2727357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77817569) q[1];
sx q[1];
rz(-0.84901224) q[1];
sx q[1];
rz(1.282882) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6333991) q[3];
sx q[3];
rz(-2.4262706) q[3];
sx q[3];
rz(-3.1014266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(-0.34040889) q[2];
rz(-2.3085964) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.67265636) q[0];
sx q[0];
rz(-2.7514973) q[0];
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
rz(-1.724052) q[0];
sx q[0];
rz(-1.1678956) q[0];
sx q[0];
rz(0.20362644) q[0];
rz(3.1171563) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(-1.2791866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7823239) q[1];
sx q[1];
rz(-1.6869079) q[1];
sx q[1];
rz(0.5639204) q[1];
rz(-pi) q[2];
rz(1.2689781) q[3];
sx q[3];
rz(-0.58773731) q[3];
sx q[3];
rz(-0.9580108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(0.46974716) q[0];
rz(-1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(0.038539561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8174724) q[0];
sx q[0];
rz(-1.5714374) q[0];
sx q[0];
rz(1.5785494) q[0];
rz(2.9107895) q[2];
sx q[2];
rz(-1.421184) q[2];
sx q[2];
rz(-0.19300592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7125394) q[1];
sx q[1];
rz(-2.6479811) q[1];
sx q[1];
rz(1.4406956) q[1];
x q[2];
rz(-0.79629691) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(-2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38347605) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-0.19515881) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965658) q[0];
sx q[0];
rz(-1.66334) q[0];
sx q[0];
rz(-1.4499614) q[0];
x q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.3865711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8198204) q[1];
sx q[1];
rz(-1.1204549) q[1];
sx q[1];
rz(2.8798073) q[1];
rz(-pi) q[2];
rz(-1.5545198) q[3];
sx q[3];
rz(-0.93173325) q[3];
sx q[3];
rz(-2.1181599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.970801) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(0.32969726) q[0];
rz(-pi) q[1];
rz(2.3230882) q[2];
sx q[2];
rz(-2.2193925) q[2];
sx q[2];
rz(1.6550145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0542478) q[1];
sx q[1];
rz(-1.9964881) q[1];
sx q[1];
rz(-2.737983) q[1];
rz(-pi) q[2];
rz(-2.674621) q[3];
sx q[3];
rz(-1.7479556) q[3];
sx q[3];
rz(-0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(0.97638431) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(-0.91517085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8422896) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(-1.6567985) q[0];
rz(-pi) q[1];
rz(-1.0038063) q[2];
sx q[2];
rz(-2.0096471) q[2];
sx q[2];
rz(-1.3958508) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8296216) q[1];
sx q[1];
rz(-1.258007) q[1];
sx q[1];
rz(-1.642986) q[1];
x q[2];
rz(1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(-2.8516155) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.4253915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4331626) q[0];
sx q[0];
rz(-1.3843952) q[0];
sx q[0];
rz(-2.7011046) q[0];
x q[1];
rz(0.13397127) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(-1.9884895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2196301) q[1];
sx q[1];
rz(-1.1929409) q[1];
sx q[1];
rz(-0.65803836) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1281934) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(0.79813938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.210775) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84412557) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(0.076518245) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90826666) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(2.6106735) q[0];
x q[1];
rz(2.2194527) q[2];
sx q[2];
rz(-1.6862009) q[2];
sx q[2];
rz(0.61308544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.4824251) q[1];
sx q[1];
rz(-1.6842711) q[1];
rz(-1.3648871) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(2.0549783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2032808) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(0.00014649815) q[2];
rz(1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9922441) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-1.0060271) q[0];
rz(-0.26793119) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.8267652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3599167) q[0];
sx q[0];
rz(-1.3206375) q[0];
sx q[0];
rz(-0.71235384) q[0];
rz(-pi) q[1];
rz(2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(0.044791128) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6445551) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(-2.6130996) q[1];
rz(1.5626838) q[3];
sx q[3];
rz(-2.3981961) q[3];
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
rz(2.8961199) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(0.47732863) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(0.39189664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.648801) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(-0.0097983629) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7461807) q[2];
sx q[2];
rz(-2.4782964) q[2];
sx q[2];
rz(-1.3822008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0185768) q[1];
sx q[1];
rz(-1.7565109) q[1];
sx q[1];
rz(0.8155483) q[1];
rz(-2.6565353) q[3];
sx q[3];
rz(-2.5254446) q[3];
sx q[3];
rz(-2.171606) q[3];
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
rz(0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.75523238) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(0.17386439) q[2];
sx q[2];
rz(-0.78730135) q[2];
sx q[2];
rz(0.80136328) q[2];
rz(2.1203534) q[3];
sx q[3];
rz(-1.869535) q[3];
sx q[3];
rz(0.46365999) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];