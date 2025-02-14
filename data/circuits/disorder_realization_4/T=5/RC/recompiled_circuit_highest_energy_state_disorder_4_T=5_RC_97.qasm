OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8430804) q[0];
sx q[0];
rz(-1.2895583) q[0];
sx q[0];
rz(-0.02331743) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(7.0206849) q[1];
sx q[1];
rz(7.7590396) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97853547) q[0];
sx q[0];
rz(-0.082162372) q[0];
sx q[0];
rz(-2.1294247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.453773) q[2];
sx q[2];
rz(-0.3323148) q[2];
sx q[2];
rz(1.3261283) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3838073) q[1];
sx q[1];
rz(-1.9646891) q[1];
sx q[1];
rz(0.76683137) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9528927) q[3];
sx q[3];
rz(-0.77761071) q[3];
sx q[3];
rz(2.4617406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4618571) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(-0.30065817) q[2];
rz(-0.65827185) q[3];
sx q[3];
rz(-0.39977795) q[3];
sx q[3];
rz(2.5532653) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268426) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(0.17587371) q[0];
rz(0.56100065) q[1];
sx q[1];
rz(-2.439552) q[1];
sx q[1];
rz(2.9452513) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6628274) q[0];
sx q[0];
rz(-1.2365554) q[0];
sx q[0];
rz(-1.1007376) q[0];
rz(-pi) q[1];
rz(-0.11949338) q[2];
sx q[2];
rz(-1.4548848) q[2];
sx q[2];
rz(1.8279861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2113116) q[1];
sx q[1];
rz(-2.086211) q[1];
sx q[1];
rz(0.0023018259) q[1];
rz(-pi) q[2];
rz(1.4370055) q[3];
sx q[3];
rz(-1.2302629) q[3];
sx q[3];
rz(-1.9386558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.089513) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(1.6717795) q[2];
rz(-2.6164264) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(-0.66450459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1216275) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(-0.63051939) q[0];
rz(1.9969253) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-3.0847881) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6358179) q[0];
sx q[0];
rz(-2.5774101) q[0];
sx q[0];
rz(2.3057372) q[0];
rz(0.47614758) q[2];
sx q[2];
rz(-2.2000929) q[2];
sx q[2];
rz(2.4193633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2395798) q[1];
sx q[1];
rz(-0.48643349) q[1];
sx q[1];
rz(0.3063267) q[1];
rz(-pi) q[2];
rz(-1.1278387) q[3];
sx q[3];
rz(-2.6747799) q[3];
sx q[3];
rz(-1.6396322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28222617) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(-0.43914208) q[2];
rz(-1.7450843) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.0948931) q[0];
sx q[0];
rz(-0.72143227) q[0];
sx q[0];
rz(2.1266345) q[0];
rz(-0.062648423) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(0.28422022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0089909) q[0];
sx q[0];
rz(-1.9997445) q[0];
sx q[0];
rz(-1.1852077) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8880803) q[2];
sx q[2];
rz(-0.51157839) q[2];
sx q[2];
rz(1.6588841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9584459) q[1];
sx q[1];
rz(-1.5588405) q[1];
sx q[1];
rz(0.35838106) q[1];
x q[2];
rz(-1.6246068) q[3];
sx q[3];
rz(-2.887886) q[3];
sx q[3];
rz(-1.0347317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3011498) q[2];
sx q[2];
rz(-1.0089077) q[2];
sx q[2];
rz(0.76048771) q[2];
rz(0.31989756) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(-2.1130051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.29192057) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(-0.40529761) q[0];
rz(-0.48794508) q[1];
sx q[1];
rz(-2.0351724) q[1];
sx q[1];
rz(0.36283666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10754171) q[0];
sx q[0];
rz(-1.1404317) q[0];
sx q[0];
rz(1.910188) q[0];
rz(-1.9486729) q[2];
sx q[2];
rz(-1.0185756) q[2];
sx q[2];
rz(-3.0966126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.88089534) q[1];
sx q[1];
rz(-2.915307) q[1];
sx q[1];
rz(-0.21842167) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7801301) q[3];
sx q[3];
rz(-2.4130954) q[3];
sx q[3];
rz(-1.1379153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1399347) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(0.72251594) q[2];
rz(-0.51269382) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2672511) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-2.8668218) q[0];
rz(2.784506) q[1];
sx q[1];
rz(-1.5029224) q[1];
sx q[1];
rz(-1.9836609) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0561075) q[0];
sx q[0];
rz(-0.72626136) q[0];
sx q[0];
rz(-1.1994491) q[0];
x q[1];
rz(2.2311796) q[2];
sx q[2];
rz(-1.0779625) q[2];
sx q[2];
rz(-0.016066859) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14640644) q[1];
sx q[1];
rz(-2.1651092) q[1];
sx q[1];
rz(1.3239149) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66863184) q[3];
sx q[3];
rz(-2.2205846) q[3];
sx q[3];
rz(1.970552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.57924119) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-2.5518899) q[2];
rz(1.3127182) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(0.5635128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24484672) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(-0.27031159) q[0];
rz(2.7145794) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(-0.40831533) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6231096) q[0];
sx q[0];
rz(-0.38869959) q[0];
sx q[0];
rz(2.8220909) q[0];
x q[1];
rz(1.6913658) q[2];
sx q[2];
rz(-0.70115108) q[2];
sx q[2];
rz(-0.98159678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22173126) q[1];
sx q[1];
rz(-1.1581005) q[1];
sx q[1];
rz(-2.2151715) q[1];
rz(-2.4189878) q[3];
sx q[3];
rz(-1.6480903) q[3];
sx q[3];
rz(2.3376436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2030877) q[2];
sx q[2];
rz(-2.0577343) q[2];
sx q[2];
rz(2.2173524) q[2];
rz(-2.3067394) q[3];
sx q[3];
rz(-2.3122841) q[3];
sx q[3];
rz(-0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860109) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(-2.7741449) q[0];
rz(-2.5732749) q[1];
sx q[1];
rz(-1.3325997) q[1];
sx q[1];
rz(-0.41499358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2872754) q[0];
sx q[0];
rz(-1.1181896) q[0];
sx q[0];
rz(-0.29083473) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6122323) q[2];
sx q[2];
rz(-1.0870013) q[2];
sx q[2];
rz(-0.31365487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7447978) q[1];
sx q[1];
rz(-0.84229453) q[1];
sx q[1];
rz(0.70835967) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0051963768) q[3];
sx q[3];
rz(-1.2250161) q[3];
sx q[3];
rz(0.41692641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35855287) q[2];
sx q[2];
rz(-2.2542605) q[2];
sx q[2];
rz(2.3084194) q[2];
rz(-0.35946515) q[3];
sx q[3];
rz(-1.0901674) q[3];
sx q[3];
rz(1.6370714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9678765) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-3.1267082) q[0];
rz(1.124565) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(3.0344149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8421321) q[0];
sx q[0];
rz(-2.3981574) q[0];
sx q[0];
rz(1.3495665) q[0];
rz(2.6769019) q[2];
sx q[2];
rz(-0.87848262) q[2];
sx q[2];
rz(-1.2274881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1746782) q[1];
sx q[1];
rz(-2.0474829) q[1];
sx q[1];
rz(1.4033474) q[1];
rz(-2.8802683) q[3];
sx q[3];
rz(-1.3991941) q[3];
sx q[3];
rz(1.2444519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(0.85961071) q[2];
rz(0.25935069) q[3];
sx q[3];
rz(-1.7813464) q[3];
sx q[3];
rz(2.7249961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49941007) q[0];
sx q[0];
rz(-0.12633093) q[0];
sx q[0];
rz(-1.1580178) q[0];
rz(-0.49531373) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(2.8445629) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659781) q[0];
sx q[0];
rz(-2.8225464) q[0];
sx q[0];
rz(-1.5837529) q[0];
rz(-0.13003243) q[2];
sx q[2];
rz(-2.8388925) q[2];
sx q[2];
rz(1.5720194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9969284) q[1];
sx q[1];
rz(-1.4132309) q[1];
sx q[1];
rz(2.8466606) q[1];
x q[2];
rz(0.90822753) q[3];
sx q[3];
rz(-1.5836827) q[3];
sx q[3];
rz(-2.9290309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26934066) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(2.2807109) q[2];
rz(2.6779029) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(-1.1369107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9853482) q[0];
sx q[0];
rz(-1.5736268) q[0];
sx q[0];
rz(-1.9851984) q[0];
rz(1.8078177) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(1.0063408) q[2];
sx q[2];
rz(-2.153572) q[2];
sx q[2];
rz(-1.5700271) q[2];
rz(0.64846497) q[3];
sx q[3];
rz(-0.59352661) q[3];
sx q[3];
rz(1.2980579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
