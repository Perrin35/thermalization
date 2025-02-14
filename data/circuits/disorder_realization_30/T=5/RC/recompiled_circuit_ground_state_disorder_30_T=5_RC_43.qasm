OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(-0.51751408) q[0];
rz(-1.8139047) q[1];
sx q[1];
rz(-0.89193901) q[1];
sx q[1];
rz(2.5995624) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5919246) q[0];
sx q[0];
rz(-0.79093638) q[0];
sx q[0];
rz(0.24178453) q[0];
x q[1];
rz(-2.3752604) q[2];
sx q[2];
rz(-1.3856941) q[2];
sx q[2];
rz(1.9232149) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0412647) q[1];
sx q[1];
rz(-0.81721837) q[1];
sx q[1];
rz(-2.8689117) q[1];
x q[2];
rz(0.23302257) q[3];
sx q[3];
rz(-1.4015504) q[3];
sx q[3];
rz(-0.47692933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0077670495) q[2];
sx q[2];
rz(-2.0962891) q[2];
sx q[2];
rz(2.2110151) q[2];
rz(-3.0104356) q[3];
sx q[3];
rz(-0.37808642) q[3];
sx q[3];
rz(1.650943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5432878) q[0];
sx q[0];
rz(-2.2853993) q[0];
sx q[0];
rz(-1.244586) q[0];
rz(-0.92567956) q[1];
sx q[1];
rz(-1.2558179) q[1];
sx q[1];
rz(-1.6474887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0605436) q[0];
sx q[0];
rz(-0.88969066) q[0];
sx q[0];
rz(2.7174755) q[0];
rz(-pi) q[1];
rz(1.235512) q[2];
sx q[2];
rz(-2.2163894) q[2];
sx q[2];
rz(2.6123227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32116649) q[1];
sx q[1];
rz(-1.120472) q[1];
sx q[1];
rz(1.9558332) q[1];
rz(-2.8374313) q[3];
sx q[3];
rz(-1.5790325) q[3];
sx q[3];
rz(-1.3502163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99574789) q[2];
sx q[2];
rz(-0.44439849) q[2];
sx q[2];
rz(2.2399529) q[2];
rz(-1.3580648) q[3];
sx q[3];
rz(-1.3172904) q[3];
sx q[3];
rz(-0.54734126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667592) q[0];
sx q[0];
rz(-1.1792553) q[0];
sx q[0];
rz(-1.1460079) q[0];
rz(0.92578069) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(-2.944223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58008798) q[0];
sx q[0];
rz(-2.5338123) q[0];
sx q[0];
rz(-1.32919) q[0];
x q[1];
rz(-1.8756798) q[2];
sx q[2];
rz(-1.291687) q[2];
sx q[2];
rz(0.32372083) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92835411) q[1];
sx q[1];
rz(-0.23121195) q[1];
sx q[1];
rz(-1.7417045) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8160964) q[3];
sx q[3];
rz(-2.5739658) q[3];
sx q[3];
rz(2.5095255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4027412) q[2];
sx q[2];
rz(-0.98029843) q[2];
sx q[2];
rz(-2.951238) q[2];
rz(-3.0800152) q[3];
sx q[3];
rz(-0.96934167) q[3];
sx q[3];
rz(-2.0891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15678081) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(-2.5132827) q[0];
rz(-1.5208987) q[1];
sx q[1];
rz(-1.2077121) q[1];
sx q[1];
rz(-0.77888387) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91858868) q[0];
sx q[0];
rz(-2.1445334) q[0];
sx q[0];
rz(2.6911435) q[0];
rz(-2.7812296) q[2];
sx q[2];
rz(-1.833263) q[2];
sx q[2];
rz(-1.9356464) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35793909) q[1];
sx q[1];
rz(-3.0920368) q[1];
sx q[1];
rz(0.2752343) q[1];
rz(-pi) q[2];
rz(-0.67653894) q[3];
sx q[3];
rz(-2.153844) q[3];
sx q[3];
rz(0.18932115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4817619) q[2];
sx q[2];
rz(-2.063282) q[2];
sx q[2];
rz(3.1415494) q[2];
rz(1.0846042) q[3];
sx q[3];
rz(-1.1559887) q[3];
sx q[3];
rz(-2.675975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(2.4701009) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(-1.2985562) q[0];
rz(2.5647054) q[1];
sx q[1];
rz(-1.8678317) q[1];
sx q[1];
rz(2.6947122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80432207) q[0];
sx q[0];
rz(-1.7511586) q[0];
sx q[0];
rz(2.4867663) q[0];
rz(1.8700897) q[2];
sx q[2];
rz(-1.4993877) q[2];
sx q[2];
rz(1.7091319) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.943534) q[1];
sx q[1];
rz(-1.3683649) q[1];
sx q[1];
rz(-1.1620896) q[1];
rz(1.5112259) q[3];
sx q[3];
rz(-1.7559663) q[3];
sx q[3];
rz(2.6310754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4345066) q[2];
sx q[2];
rz(-0.27927566) q[2];
sx q[2];
rz(-2.9217829) q[2];
rz(1.6541727) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(2.4421104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952482) q[0];
sx q[0];
rz(-0.43287745) q[0];
sx q[0];
rz(-2.2440946) q[0];
rz(1.7706361) q[1];
sx q[1];
rz(-1.0400583) q[1];
sx q[1];
rz(0.8955566) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61293759) q[0];
sx q[0];
rz(-2.2870758) q[0];
sx q[0];
rz(-1.3837293) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3562548) q[2];
sx q[2];
rz(-0.78162748) q[2];
sx q[2];
rz(1.8434032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59437737) q[1];
sx q[1];
rz(-1.9037694) q[1];
sx q[1];
rz(2.7463288) q[1];
rz(-pi) q[2];
rz(-2.343781) q[3];
sx q[3];
rz(-0.63029248) q[3];
sx q[3];
rz(-0.52784656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7728277) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(-2.9844798) q[2];
rz(-0.67240063) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(3.0290643) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2812578) q[0];
sx q[0];
rz(-0.66547886) q[0];
sx q[0];
rz(-1.6819287) q[0];
rz(0.36059391) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(1.8416539) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4111709) q[0];
sx q[0];
rz(-1.9864559) q[0];
sx q[0];
rz(1.4372197) q[0];
rz(-pi) q[1];
rz(-0.62076135) q[2];
sx q[2];
rz(-1.1859535) q[2];
sx q[2];
rz(2.8892725) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5438788) q[1];
sx q[1];
rz(-1.2744441) q[1];
sx q[1];
rz(1.4696025) q[1];
rz(2.5266493) q[3];
sx q[3];
rz(-1.9087026) q[3];
sx q[3];
rz(0.29454052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.931687) q[2];
sx q[2];
rz(-2.1114712) q[2];
sx q[2];
rz(-0.21200655) q[2];
rz(2.0906406) q[3];
sx q[3];
rz(-1.3764328) q[3];
sx q[3];
rz(-0.088137805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5365005) q[0];
sx q[0];
rz(-0.78992805) q[0];
sx q[0];
rz(-2.9686046) q[0];
rz(1.1146924) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(-3.0317422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45013406) q[0];
sx q[0];
rz(-0.5930674) q[0];
sx q[0];
rz(2.791138) q[0];
rz(0.82288701) q[2];
sx q[2];
rz(-1.7639065) q[2];
sx q[2];
rz(0.88179526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50625769) q[1];
sx q[1];
rz(-1.3999363) q[1];
sx q[1];
rz(1.3116763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6906711) q[3];
sx q[3];
rz(-1.0195512) q[3];
sx q[3];
rz(0.84037432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7383808) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(2.476725) q[2];
rz(-2.6730149) q[3];
sx q[3];
rz(-1.6178308) q[3];
sx q[3];
rz(-2.328228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68798962) q[0];
sx q[0];
rz(-0.60992321) q[0];
sx q[0];
rz(2.4000121) q[0];
rz(1.954151) q[1];
sx q[1];
rz(-1.5267742) q[1];
sx q[1];
rz(0.56914079) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65780902) q[0];
sx q[0];
rz(-1.8632194) q[0];
sx q[0];
rz(0.34088366) q[0];
rz(-2.6879758) q[2];
sx q[2];
rz(-1.5065644) q[2];
sx q[2];
rz(-1.0421696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7267159) q[1];
sx q[1];
rz(-1.1099713) q[1];
sx q[1];
rz(1.0970835) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69066234) q[3];
sx q[3];
rz(-0.96513018) q[3];
sx q[3];
rz(0.31177917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6471214) q[2];
sx q[2];
rz(-2.1529866) q[2];
sx q[2];
rz(0.76623255) q[2];
rz(-1.9689485) q[3];
sx q[3];
rz(-1.5394883) q[3];
sx q[3];
rz(-0.22181454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93429339) q[0];
sx q[0];
rz(-2.793029) q[0];
sx q[0];
rz(2.9497414) q[0];
rz(2.3244997) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(0.70768913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1199214) q[0];
sx q[0];
rz(-1.7448398) q[0];
sx q[0];
rz(-0.24876033) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0470263) q[2];
sx q[2];
rz(-2.4942538) q[2];
sx q[2];
rz(-0.29905427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0470815) q[1];
sx q[1];
rz(-0.99267497) q[1];
sx q[1];
rz(2.6972527) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5557184) q[3];
sx q[3];
rz(-1.0535686) q[3];
sx q[3];
rz(-0.45211238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2231458) q[2];
sx q[2];
rz(-2.196849) q[2];
sx q[2];
rz(2.4165912) q[2];
rz(2.9265192) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(0.71896368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216777) q[0];
sx q[0];
rz(-2.324993) q[0];
sx q[0];
rz(-2.8391229) q[0];
rz(2.0401781) q[1];
sx q[1];
rz(-1.2322203) q[1];
sx q[1];
rz(0.89422918) q[1];
rz(2.106059) q[2];
sx q[2];
rz(-1.4695833) q[2];
sx q[2];
rz(1.5493456) q[2];
rz(-1.5433031) q[3];
sx q[3];
rz(-2.0381654) q[3];
sx q[3];
rz(0.5034133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
