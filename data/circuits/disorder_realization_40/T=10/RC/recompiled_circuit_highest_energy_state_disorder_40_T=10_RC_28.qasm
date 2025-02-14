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
rz(2.1514819) q[0];
sx q[0];
rz(2.4408542) q[0];
sx q[0];
rz(10.182575) q[0];
rz(-2.4611729) q[1];
sx q[1];
rz(-1.0935723) q[1];
sx q[1];
rz(3.0419066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62698088) q[0];
sx q[0];
rz(-2.520769) q[0];
sx q[0];
rz(2.6326535) q[0];
rz(-0.52844471) q[2];
sx q[2];
rz(-2.4917654) q[2];
sx q[2];
rz(-0.53534269) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0792011) q[1];
sx q[1];
rz(-1.3729368) q[1];
sx q[1];
rz(1.3049091) q[1];
rz(-pi) q[2];
rz(-1.4138347) q[3];
sx q[3];
rz(-2.779134) q[3];
sx q[3];
rz(-0.80998224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1216571) q[2];
sx q[2];
rz(-1.1172373) q[2];
sx q[2];
rz(3.115444) q[2];
rz(-1.3610241) q[3];
sx q[3];
rz(-1.0739505) q[3];
sx q[3];
rz(1.1092383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015198135) q[0];
sx q[0];
rz(-0.69271883) q[0];
sx q[0];
rz(-1.3808845) q[0];
rz(0.2805925) q[1];
sx q[1];
rz(-2.3255489) q[1];
sx q[1];
rz(2.7516344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0637012) q[0];
sx q[0];
rz(-2.3305463) q[0];
sx q[0];
rz(1.6249932) q[0];
rz(3.105279) q[2];
sx q[2];
rz(-0.89082375) q[2];
sx q[2];
rz(2.3009686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.094495) q[1];
sx q[1];
rz(-0.19097155) q[1];
sx q[1];
rz(0.32778962) q[1];
rz(-pi) q[2];
rz(1.5481351) q[3];
sx q[3];
rz(-2.0840328) q[3];
sx q[3];
rz(-3.0340241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0661162) q[2];
sx q[2];
rz(-0.37049299) q[2];
sx q[2];
rz(2.9541435) q[2];
rz(0.97505331) q[3];
sx q[3];
rz(-1.2504028) q[3];
sx q[3];
rz(-1.3343078) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27075574) q[0];
sx q[0];
rz(-2.2360531) q[0];
sx q[0];
rz(1.723473) q[0];
rz(-2.1899147) q[1];
sx q[1];
rz(-2.8680809) q[1];
sx q[1];
rz(0.66114122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9959629) q[0];
sx q[0];
rz(-1.055777) q[0];
sx q[0];
rz(1.4670232) q[0];
rz(-0.69532891) q[2];
sx q[2];
rz(-1.2711356) q[2];
sx q[2];
rz(-1.9315182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3959655) q[1];
sx q[1];
rz(-2.9955609) q[1];
sx q[1];
rz(-1.795799) q[1];
rz(0.16942008) q[3];
sx q[3];
rz(-0.5656618) q[3];
sx q[3];
rz(-3.0974714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4989138) q[2];
sx q[2];
rz(-1.4475409) q[2];
sx q[2];
rz(-1.0399237) q[2];
rz(-2.6652279) q[3];
sx q[3];
rz(-0.56291181) q[3];
sx q[3];
rz(-2.3001455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845402) q[0];
sx q[0];
rz(-0.81890506) q[0];
sx q[0];
rz(-0.99405974) q[0];
rz(2.0301863) q[1];
sx q[1];
rz(-1.2329654) q[1];
sx q[1];
rz(-2.9434189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4497714) q[0];
sx q[0];
rz(-1.9110962) q[0];
sx q[0];
rz(-1.6671851) q[0];
rz(-pi) q[1];
rz(-1.092574) q[2];
sx q[2];
rz(-1.0727777) q[2];
sx q[2];
rz(-0.39229506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6804643) q[1];
sx q[1];
rz(-1.2733508) q[1];
sx q[1];
rz(-0.72338413) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1311225) q[3];
sx q[3];
rz(-2.098985) q[3];
sx q[3];
rz(-1.8757602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9626832) q[2];
sx q[2];
rz(-1.6390025) q[2];
sx q[2];
rz(-2.7965014) q[2];
rz(-2.1033449) q[3];
sx q[3];
rz(-2.0227261) q[3];
sx q[3];
rz(-3.069675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2800696) q[0];
sx q[0];
rz(-0.11174209) q[0];
sx q[0];
rz(0.34568632) q[0];
rz(0.12738906) q[1];
sx q[1];
rz(-0.40475875) q[1];
sx q[1];
rz(1.60188) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87563959) q[0];
sx q[0];
rz(-1.709076) q[0];
sx q[0];
rz(-1.7931275) q[0];
rz(-2.0284838) q[2];
sx q[2];
rz(-2.1808592) q[2];
sx q[2];
rz(-1.8230566) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8495765) q[1];
sx q[1];
rz(-0.56835163) q[1];
sx q[1];
rz(-0.59225099) q[1];
x q[2];
rz(0.57347371) q[3];
sx q[3];
rz(-1.6454564) q[3];
sx q[3];
rz(0.71686137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6973528) q[2];
sx q[2];
rz(-1.3109861) q[2];
sx q[2];
rz(0.41651192) q[2];
rz(0.087827772) q[3];
sx q[3];
rz(-0.4952687) q[3];
sx q[3];
rz(-2.7145568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86513615) q[0];
sx q[0];
rz(-2.3910531) q[0];
sx q[0];
rz(0.99391341) q[0];
rz(-1.7991426) q[1];
sx q[1];
rz(-1.3748704) q[1];
sx q[1];
rz(1.3329175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48644984) q[0];
sx q[0];
rz(-2.764439) q[0];
sx q[0];
rz(2.7906453) q[0];
rz(-pi) q[1];
rz(-0.41282515) q[2];
sx q[2];
rz(-1.7805837) q[2];
sx q[2];
rz(-2.0693463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.7018334) q[1];
sx q[1];
rz(-1.5343795) q[1];
sx q[1];
rz(-1.4324051) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52458472) q[3];
sx q[3];
rz(-0.3915873) q[3];
sx q[3];
rz(-0.064849555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6028613) q[2];
sx q[2];
rz(-1.8659325) q[2];
sx q[2];
rz(-0.71962792) q[2];
rz(2.1384278) q[3];
sx q[3];
rz(-1.7374618) q[3];
sx q[3];
rz(-2.7145332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.3908983) q[0];
sx q[0];
rz(-2.7321132) q[0];
sx q[0];
rz(0.77792186) q[0];
rz(-1.3748417) q[1];
sx q[1];
rz(-2.1781616) q[1];
sx q[1];
rz(0.28819293) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88633986) q[0];
sx q[0];
rz(-1.1561516) q[0];
sx q[0];
rz(-2.9575081) q[0];
rz(-pi) q[1];
rz(-2.8164045) q[2];
sx q[2];
rz(-0.033043229) q[2];
sx q[2];
rz(0.059800241) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.035305939) q[1];
sx q[1];
rz(-1.0540153) q[1];
sx q[1];
rz(1.27716) q[1];
x q[2];
rz(1.8882264) q[3];
sx q[3];
rz(-0.39690382) q[3];
sx q[3];
rz(-1.1008747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9162743) q[2];
sx q[2];
rz(-2.306873) q[2];
sx q[2];
rz(1.6342596) q[2];
rz(-2.7507014) q[3];
sx q[3];
rz(-1.8796128) q[3];
sx q[3];
rz(1.5220557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759719) q[0];
sx q[0];
rz(-0.93637744) q[0];
sx q[0];
rz(0.17882624) q[0];
rz(1.7504182) q[1];
sx q[1];
rz(-1.6988674) q[1];
sx q[1];
rz(0.19405445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3349105) q[0];
sx q[0];
rz(-1.9437978) q[0];
sx q[0];
rz(-1.0399517) q[0];
rz(-0.75066815) q[2];
sx q[2];
rz(-1.0714415) q[2];
sx q[2];
rz(-2.4908092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0240101) q[1];
sx q[1];
rz(-0.76723209) q[1];
sx q[1];
rz(-1.3831286) q[1];
x q[2];
rz(0.91839478) q[3];
sx q[3];
rz(-1.4271591) q[3];
sx q[3];
rz(-1.4220765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27851823) q[2];
sx q[2];
rz(-2.9926509) q[2];
sx q[2];
rz(-1.2433368) q[2];
rz(2.8362823) q[3];
sx q[3];
rz(-1.52933) q[3];
sx q[3];
rz(-2.7436411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0447926) q[0];
sx q[0];
rz(-0.069792062) q[0];
sx q[0];
rz(-3.0332562) q[0];
rz(2.4414252) q[1];
sx q[1];
rz(-1.6077176) q[1];
sx q[1];
rz(-2.7549423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413992) q[0];
sx q[0];
rz(-0.94111573) q[0];
sx q[0];
rz(-2.5085418) q[0];
rz(-pi) q[1];
rz(0.79603802) q[2];
sx q[2];
rz(-2.926119) q[2];
sx q[2];
rz(2.3448159) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0476847) q[1];
sx q[1];
rz(-2.3833774) q[1];
sx q[1];
rz(0.42102937) q[1];
rz(-2.2240102) q[3];
sx q[3];
rz(-2.4303103) q[3];
sx q[3];
rz(-0.68071738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.391523) q[2];
sx q[2];
rz(-2.0400901) q[2];
sx q[2];
rz(-2.2186685) q[2];
rz(0.37724885) q[3];
sx q[3];
rz(-2.173285) q[3];
sx q[3];
rz(1.7678461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7467576) q[0];
sx q[0];
rz(-2.6834798) q[0];
sx q[0];
rz(0.75570345) q[0];
rz(-2.34756) q[1];
sx q[1];
rz(-2.463412) q[1];
sx q[1];
rz(-1.979535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1080429) q[0];
sx q[0];
rz(-0.4539412) q[0];
sx q[0];
rz(2.0583711) q[0];
x q[1];
rz(-0.65067567) q[2];
sx q[2];
rz(-1.2066368) q[2];
sx q[2];
rz(-2.3641567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5291497) q[1];
sx q[1];
rz(-0.72436391) q[1];
sx q[1];
rz(2.0606478) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28323876) q[3];
sx q[3];
rz(-1.9593718) q[3];
sx q[3];
rz(2.376414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0572759) q[2];
sx q[2];
rz(-0.82547775) q[2];
sx q[2];
rz(-3.0153073) q[2];
rz(-2.5430191) q[3];
sx q[3];
rz(-1.5263661) q[3];
sx q[3];
rz(1.9954782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3446658) q[0];
sx q[0];
rz(-1.666259) q[0];
sx q[0];
rz(1.689612) q[0];
rz(-1.7563734) q[1];
sx q[1];
rz(-1.9671556) q[1];
sx q[1];
rz(3.0671493) q[1];
rz(-2.47592) q[2];
sx q[2];
rz(-2.4187928) q[2];
sx q[2];
rz(-1.7841675) q[2];
rz(-2.606288) q[3];
sx q[3];
rz(-0.44620958) q[3];
sx q[3];
rz(-0.37887497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
