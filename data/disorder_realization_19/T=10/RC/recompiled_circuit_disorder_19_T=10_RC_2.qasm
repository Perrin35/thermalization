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
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(1.4979314) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6490373) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(-1.9306246) q[0];
rz(-pi) q[1];
rz(-1.6662007) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(0.22104095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5421996) q[1];
sx q[1];
rz(-1.7855872) q[1];
sx q[1];
rz(-0.7426803) q[1];
x q[2];
rz(0.50819355) q[3];
sx q[3];
rz(-2.4262706) q[3];
sx q[3];
rz(0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(0.90822059) q[1];
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
rz(-2.2087814) q[0];
sx q[0];
rz(-0.44890807) q[0];
sx q[0];
rz(1.1277898) q[0];
rz(-pi) q[1];
rz(1.5812133) q[2];
sx q[2];
rz(-1.1679107) q[2];
sx q[2];
rz(1.8889697) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35926871) q[1];
sx q[1];
rz(-1.4546848) q[1];
sx q[1];
rz(0.5639204) q[1];
rz(1.8726146) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(-0.9580108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-0.17641243) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(3.1030531) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32917133) q[0];
sx q[0];
rz(-0.0077795452) q[0];
sx q[0];
rz(-1.4882985) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23080319) q[2];
sx q[2];
rz(-1.7204086) q[2];
sx q[2];
rz(-0.19300592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4290532) q[1];
sx q[1];
rz(-2.6479811) q[1];
sx q[1];
rz(1.4406956) q[1];
x q[2];
rz(-0.58979804) q[3];
sx q[3];
rz(-1.0538978) q[3];
sx q[3];
rz(-0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(0.19515881) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0965658) q[0];
sx q[0];
rz(-1.4782527) q[0];
sx q[0];
rz(-1.4499614) q[0];
rz(-0.85907016) q[2];
sx q[2];
rz(-0.67407437) q[2];
sx q[2];
rz(0.40875834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
x q[2];
rz(3.1196932) q[3];
sx q[3];
rz(-2.5023513) q[3];
sx q[3];
rz(-0.99614776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.0347663) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(-0.37809125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17079167) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(-2.8118954) q[0];
rz(-pi) q[1];
rz(0.73363186) q[2];
sx q[2];
rz(-2.1918104) q[2];
sx q[2];
rz(-2.6526407) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0542478) q[1];
sx q[1];
rz(-1.1451045) q[1];
sx q[1];
rz(-0.40360968) q[1];
rz(-pi) q[2];
rz(0.46697163) q[3];
sx q[3];
rz(-1.7479556) q[3];
sx q[3];
rz(-0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(-1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75509214) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-2.7639672) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(2.2264218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5650425) q[0];
sx q[0];
rz(-2.8083907) q[0];
sx q[0];
rz(2.8898426) q[0];
x q[1];
rz(2.2890131) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(-2.728087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.542791) q[1];
sx q[1];
rz(-0.32074499) q[1];
sx q[1];
rz(0.21943211) q[1];
rz(-pi) q[2];
rz(0.1359453) q[3];
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
rz(0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(2.3510695) q[2];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(-2.7745568) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(1.4253915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51222425) q[0];
sx q[0];
rz(-2.6656796) q[0];
sx q[0];
rz(-2.7251564) q[0];
rz(-0.44600365) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(-0.29689483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.047445) q[1];
sx q[1];
rz(-2.3970251) q[1];
sx q[1];
rz(0.57569699) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0745139) q[3];
sx q[3];
rz(-1.5825315) q[3];
sx q[3];
rz(2.3624682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.210775) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84412557) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(3.0650744) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-2.3892367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1479748) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(-0.57172914) q[0];
x q[1];
rz(0.92213995) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(-2.5285072) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.3648871) q[3];
sx q[3];
rz(-2.1504953) q[3];
sx q[3];
rz(2.0549783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(0.00014649815) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00025230322) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(-1.89639) q[0];
x q[1];
rz(2.7776412) q[2];
sx q[2];
rz(-1.4462573) q[2];
sx q[2];
rz(-0.044791128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41373738) q[1];
sx q[1];
rz(-1.9836042) q[1];
sx q[1];
rz(0.84819838) q[1];
rz(-pi) q[2];
rz(1.5626838) q[3];
sx q[3];
rz(-0.74339657) q[3];
sx q[3];
rz(2.5826366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(-0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7219287) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.648801) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(0.0097983629) q[0];
rz(-0.91498615) q[2];
sx q[2];
rz(-1.4631541) q[2];
sx q[2];
rz(-0.3273302) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74589409) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.8383154) q[1];
rz(-pi) q[2];
rz(-1.8896905) q[3];
sx q[3];
rz(-2.1074169) q[3];
sx q[3];
rz(1.5434138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-0.77970589) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(1.0212392) q[3];
sx q[3];
rz(-1.2720576) q[3];
sx q[3];
rz(-2.6779327) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
