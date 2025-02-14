OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9775554) q[0];
sx q[0];
rz(-1.2241751) q[0];
sx q[0];
rz(1.2807711) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(0.26689902) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42519841) q[0];
sx q[0];
rz(-3.0073637) q[0];
sx q[0];
rz(-2.3421046) q[0];
rz(-pi) q[1];
x q[1];
rz(2.884567) q[2];
sx q[2];
rz(-2.2429209) q[2];
sx q[2];
rz(2.4804403) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1303838) q[1];
sx q[1];
rz(-0.88700548) q[1];
sx q[1];
rz(2.5793772) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21453114) q[3];
sx q[3];
rz(-2.6900468) q[3];
sx q[3];
rz(0.2487693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12307564) q[2];
sx q[2];
rz(-0.96394959) q[2];
sx q[2];
rz(-0.71259552) q[2];
rz(0.16768843) q[3];
sx q[3];
rz(-2.2919787) q[3];
sx q[3];
rz(-2.6970421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.0406822) q[0];
sx q[0];
rz(-1.7829144) q[0];
sx q[0];
rz(1.4605301) q[0];
rz(-1.679861) q[1];
sx q[1];
rz(-2.0748383) q[1];
sx q[1];
rz(1.1163968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8921601) q[0];
sx q[0];
rz(-1.5276225) q[0];
sx q[0];
rz(-0.024726111) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1266534) q[2];
sx q[2];
rz(-1.6995322) q[2];
sx q[2];
rz(-1.5644846) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4907995) q[1];
sx q[1];
rz(-1.3099346) q[1];
sx q[1];
rz(-0.38968226) q[1];
x q[2];
rz(1.7061911) q[3];
sx q[3];
rz(-1.97768) q[3];
sx q[3];
rz(-3.1171006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98550335) q[2];
sx q[2];
rz(-2.6458793) q[2];
sx q[2];
rz(2.8893341) q[2];
rz(-2.0856527) q[3];
sx q[3];
rz(-2.2838433) q[3];
sx q[3];
rz(-1.0004388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5591739) q[0];
sx q[0];
rz(-0.90455872) q[0];
sx q[0];
rz(0.82758033) q[0];
rz(2.6170392) q[1];
sx q[1];
rz(-1.2453715) q[1];
sx q[1];
rz(0.96447271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23339749) q[0];
sx q[0];
rz(-2.1475773) q[0];
sx q[0];
rz(-0.080349313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9461412) q[2];
sx q[2];
rz(-1.4895194) q[2];
sx q[2];
rz(0.92513212) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8207404) q[1];
sx q[1];
rz(-2.2005282) q[1];
sx q[1];
rz(-1.6321502) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12196864) q[3];
sx q[3];
rz(-2.6971574) q[3];
sx q[3];
rz(-1.1201064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39782897) q[2];
sx q[2];
rz(-0.22481329) q[2];
sx q[2];
rz(1.942983) q[2];
rz(-1.648929) q[3];
sx q[3];
rz(-0.97380939) q[3];
sx q[3];
rz(-0.57884136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156133) q[0];
sx q[0];
rz(-0.15376832) q[0];
sx q[0];
rz(1.2157259) q[0];
rz(-2.1513596) q[1];
sx q[1];
rz(-2.4001887) q[1];
sx q[1];
rz(0.56891099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7051681) q[0];
sx q[0];
rz(-1.2878083) q[0];
sx q[0];
rz(0.97375906) q[0];
x q[1];
rz(-1.2014451) q[2];
sx q[2];
rz(-1.8180038) q[2];
sx q[2];
rz(-1.9245242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7987631) q[1];
sx q[1];
rz(-1.1278648) q[1];
sx q[1];
rz(2.8196067) q[1];
x q[2];
rz(-2.4510108) q[3];
sx q[3];
rz(-0.8775228) q[3];
sx q[3];
rz(-2.1335672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8625921) q[2];
sx q[2];
rz(-0.96070868) q[2];
sx q[2];
rz(-0.29435364) q[2];
rz(-2.477008) q[3];
sx q[3];
rz(-2.3817101) q[3];
sx q[3];
rz(1.6871281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.45254529) q[0];
sx q[0];
rz(-1.4280467) q[0];
sx q[0];
rz(0.31387615) q[0];
rz(-2.2398056) q[1];
sx q[1];
rz(-1.3548464) q[1];
sx q[1];
rz(-1.4656167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70756992) q[0];
sx q[0];
rz(-1.5450053) q[0];
sx q[0];
rz(-1.4517734) q[0];
rz(0.11326934) q[2];
sx q[2];
rz(-1.3217032) q[2];
sx q[2];
rz(-1.3048362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61205381) q[1];
sx q[1];
rz(-1.1837237) q[1];
sx q[1];
rz(-1.4899292) q[1];
x q[2];
rz(-0.63057804) q[3];
sx q[3];
rz(-1.6007716) q[3];
sx q[3];
rz(1.1512299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19447154) q[2];
sx q[2];
rz(-0.7496382) q[2];
sx q[2];
rz(-1.0844024) q[2];
rz(0.32367745) q[3];
sx q[3];
rz(-2.4249228) q[3];
sx q[3];
rz(-0.22423854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913469) q[0];
sx q[0];
rz(-2.8960189) q[0];
sx q[0];
rz(-2.7299951) q[0];
rz(-2.4161074) q[1];
sx q[1];
rz(-1.2440224) q[1];
sx q[1];
rz(1.4010319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0593911) q[0];
sx q[0];
rz(-0.93125805) q[0];
sx q[0];
rz(-1.3277505) q[0];
rz(1.696642) q[2];
sx q[2];
rz(-2.3072444) q[2];
sx q[2];
rz(3.1395314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7159375) q[1];
sx q[1];
rz(-1.9129941) q[1];
sx q[1];
rz(1.9075431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4947004) q[3];
sx q[3];
rz(-1.1206566) q[3];
sx q[3];
rz(1.3967747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1650042) q[2];
sx q[2];
rz(-0.62825957) q[2];
sx q[2];
rz(-1.4264301) q[2];
rz(-0.12380883) q[3];
sx q[3];
rz(-1.0694458) q[3];
sx q[3];
rz(0.4933221) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2774778) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(1.8051099) q[0];
rz(2.2026964) q[1];
sx q[1];
rz(-2.0195596) q[1];
sx q[1];
rz(2.8575361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8495318) q[0];
sx q[0];
rz(-2.1725328) q[0];
sx q[0];
rz(-1.9868149) q[0];
x q[1];
rz(-0.16347537) q[2];
sx q[2];
rz(-1.3473841) q[2];
sx q[2];
rz(2.7976409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3160243) q[1];
sx q[1];
rz(-0.033016769) q[1];
sx q[1];
rz(1.2570639) q[1];
x q[2];
rz(0.18064336) q[3];
sx q[3];
rz(-1.4339281) q[3];
sx q[3];
rz(1.6302072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9629918) q[2];
sx q[2];
rz(-0.68135771) q[2];
sx q[2];
rz(-0.94092384) q[2];
rz(-0.00087794463) q[3];
sx q[3];
rz(-1.2255171) q[3];
sx q[3];
rz(-3.0555449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451999) q[0];
sx q[0];
rz(-0.89458507) q[0];
sx q[0];
rz(-2.9659502) q[0];
rz(-1.9215709) q[1];
sx q[1];
rz(-2.2717387) q[1];
sx q[1];
rz(-0.87179914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052446121) q[0];
sx q[0];
rz(-1.8348357) q[0];
sx q[0];
rz(-2.9857078) q[0];
rz(-pi) q[1];
rz(2.7705338) q[2];
sx q[2];
rz(-0.81766978) q[2];
sx q[2];
rz(-2.121986) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6260274) q[1];
sx q[1];
rz(-2.5142418) q[1];
sx q[1];
rz(-0.27317843) q[1];
rz(-3.0230247) q[3];
sx q[3];
rz(-2.5817079) q[3];
sx q[3];
rz(-0.15794755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4153727) q[2];
sx q[2];
rz(-1.2803593) q[2];
sx q[2];
rz(1.0225164) q[2];
rz(-0.38698777) q[3];
sx q[3];
rz(-1.3849247) q[3];
sx q[3];
rz(2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33685327) q[0];
sx q[0];
rz(-1.9751208) q[0];
sx q[0];
rz(0.26300305) q[0];
rz(2.8291342) q[1];
sx q[1];
rz(-2.0461021) q[1];
sx q[1];
rz(1.9591263) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0175041) q[0];
sx q[0];
rz(-2.5944983) q[0];
sx q[0];
rz(-1.9518243) q[0];
rz(0.55652852) q[2];
sx q[2];
rz(-0.91224837) q[2];
sx q[2];
rz(1.5040894) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39315614) q[1];
sx q[1];
rz(-1.9416205) q[1];
sx q[1];
rz(-1.7840339) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.617917) q[3];
sx q[3];
rz(-1.0645691) q[3];
sx q[3];
rz(-0.62577137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.075228127) q[2];
sx q[2];
rz(-1.3508537) q[2];
sx q[2];
rz(2.8933375) q[2];
rz(-2.9634641) q[3];
sx q[3];
rz(-1.0823366) q[3];
sx q[3];
rz(0.78290141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42661509) q[0];
sx q[0];
rz(-2.7474032) q[0];
sx q[0];
rz(2.07975) q[0];
rz(1.3007523) q[1];
sx q[1];
rz(-1.5251093) q[1];
sx q[1];
rz(1.1311857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12515629) q[0];
sx q[0];
rz(-1.1068692) q[0];
sx q[0];
rz(2.6317276) q[0];
rz(-pi) q[1];
rz(0.71604095) q[2];
sx q[2];
rz(-2.5114905) q[2];
sx q[2];
rz(1.014647) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50000514) q[1];
sx q[1];
rz(-1.4785071) q[1];
sx q[1];
rz(-2.7422043) q[1];
x q[2];
rz(0.33343306) q[3];
sx q[3];
rz(-2.5134183) q[3];
sx q[3];
rz(1.6510568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0807557) q[2];
sx q[2];
rz(-0.26641014) q[2];
sx q[2];
rz(2.5401435) q[2];
rz(-2.1001749) q[3];
sx q[3];
rz(-1.4789378) q[3];
sx q[3];
rz(-1.7634332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1489442) q[0];
sx q[0];
rz(-2.5288378) q[0];
sx q[0];
rz(0.80768325) q[0];
rz(0.07769892) q[1];
sx q[1];
rz(-0.54010375) q[1];
sx q[1];
rz(1.7355951) q[1];
rz(-1.0098828) q[2];
sx q[2];
rz(-2.461888) q[2];
sx q[2];
rz(1.6585435) q[2];
rz(-2.5004456) q[3];
sx q[3];
rz(-2.2401886) q[3];
sx q[3];
rz(-0.2284579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
