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
rz(-1.7582769) q[0];
sx q[0];
rz(4.5780616) q[0];
sx q[0];
rz(8.4599001) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(-0.42958346) q[1];
sx q[1];
rz(2.092195) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8937738) q[0];
sx q[0];
rz(-0.6914833) q[0];
sx q[0];
rz(-1.0214424) q[0];
rz(-1.2464748) q[2];
sx q[2];
rz(-1.4359546) q[2];
sx q[2];
rz(-2.2244649) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2856372) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(0.4680856) q[1];
rz(-pi) q[2];
rz(-1.4780294) q[3];
sx q[3];
rz(-2.1481107) q[3];
sx q[3];
rz(2.4924459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6875978) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(0.018608658) q[2];
rz(0.54443693) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740771) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(-2.6065705) q[0];
rz(0.53994838) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(-1.7417057) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8606621) q[0];
sx q[0];
rz(-1.9561523) q[0];
sx q[0];
rz(-1.0866665) q[0];
x q[1];
rz(-0.56324701) q[2];
sx q[2];
rz(-1.0457102) q[2];
sx q[2];
rz(1.4525177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8230086) q[1];
sx q[1];
rz(-0.46097791) q[1];
sx q[1];
rz(2.5557842) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0008282) q[3];
sx q[3];
rz(-1.0335575) q[3];
sx q[3];
rz(-1.1972103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0743559) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(-2.2155217) q[2];
rz(-2.1615084) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(-1.5128304) q[0];
rz(0.28383645) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(-2.0294752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29942) q[0];
sx q[0];
rz(-1.5636235) q[0];
sx q[0];
rz(0.0070863574) q[0];
x q[1];
rz(-2.1826964) q[2];
sx q[2];
rz(-0.13223967) q[2];
sx q[2];
rz(-2.7041777) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2123472) q[1];
sx q[1];
rz(-0.51873365) q[1];
sx q[1];
rz(-2.7580845) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9341957) q[3];
sx q[3];
rz(-1.8135241) q[3];
sx q[3];
rz(0.12871615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6185559) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(-2.0396566) q[2];
rz(-2.2526422) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76982826) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(0.68921047) q[0];
rz(0.31461942) q[1];
sx q[1];
rz(-2.2210821) q[1];
sx q[1];
rz(-1.4124195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310881) q[0];
sx q[0];
rz(-2.304739) q[0];
sx q[0];
rz(-0.42829163) q[0];
rz(-pi) q[1];
rz(0.41823776) q[2];
sx q[2];
rz(-2.9091638) q[2];
sx q[2];
rz(-1.4168036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3312348) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(-2.9879301) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.054677) q[3];
sx q[3];
rz(-1.7431419) q[3];
sx q[3];
rz(-2.2215171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41669258) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(-0.55595428) q[2];
rz(1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(-0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81412643) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(-2.5471174) q[0];
rz(0.56198436) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(2.1655703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9851889) q[0];
sx q[0];
rz(-1.9984584) q[0];
sx q[0];
rz(-0.39910631) q[0];
x q[1];
rz(3.1107509) q[2];
sx q[2];
rz(-0.75428666) q[2];
sx q[2];
rz(-0.13067836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4415293) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(-1.1928012) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30140437) q[3];
sx q[3];
rz(-1.7782974) q[3];
sx q[3];
rz(-2.1077771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(0.61069926) q[2];
rz(-2.8357909) q[3];
sx q[3];
rz(-2.1761201) q[3];
sx q[3];
rz(1.8122199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139451) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(-1.6754643) q[0];
rz(-0.77955359) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-2.4028042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1818016) q[0];
sx q[0];
rz(-2.7914146) q[0];
sx q[0];
rz(2.2018593) q[0];
rz(-pi) q[1];
rz(-0.30679484) q[2];
sx q[2];
rz(-1.5606797) q[2];
sx q[2];
rz(-0.9900569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2704084) q[1];
sx q[1];
rz(-1.021402) q[1];
sx q[1];
rz(2.710538) q[1];
rz(0.36130623) q[3];
sx q[3];
rz(-1.2194467) q[3];
sx q[3];
rz(1.460399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61591992) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(-2.2789148) q[2];
rz(-2.681813) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(1.9988683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2343242) q[0];
sx q[0];
rz(-2.7643804) q[0];
sx q[0];
rz(-2.1211076) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(-0.78757706) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4754935) q[0];
sx q[0];
rz(-1.8084053) q[0];
sx q[0];
rz(0.65639021) q[0];
rz(2.8767005) q[2];
sx q[2];
rz(-0.26703003) q[2];
sx q[2];
rz(-1.9445813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0246525) q[1];
sx q[1];
rz(-1.7491915) q[1];
sx q[1];
rz(-2.0703038) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5191139) q[3];
sx q[3];
rz(-2.5530836) q[3];
sx q[3];
rz(-1.9543598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82137498) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(0.58464948) q[2];
rz(0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24340165) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(-2.8133494) q[0];
rz(-2.7499061) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(1.6627056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97349629) q[0];
sx q[0];
rz(-0.35051051) q[0];
sx q[0];
rz(-2.4309733) q[0];
x q[1];
rz(-2.0441891) q[2];
sx q[2];
rz(-2.0461296) q[2];
sx q[2];
rz(1.2021827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7110243) q[1];
sx q[1];
rz(-1.4757753) q[1];
sx q[1];
rz(1.8501758) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9401266) q[3];
sx q[3];
rz(-1.374326) q[3];
sx q[3];
rz(-0.35614355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(-2.3528986) q[2];
rz(-0.24104077) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(-2.327976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64860827) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(0.96187821) q[0];
rz(-2.9073763) q[1];
sx q[1];
rz(-1.1385671) q[1];
sx q[1];
rz(-0.73807565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3320112) q[0];
sx q[0];
rz(-1.9652848) q[0];
sx q[0];
rz(1.8277192) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5066824) q[2];
sx q[2];
rz(-1.2037841) q[2];
sx q[2];
rz(2.5653332) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9752561) q[1];
sx q[1];
rz(-1.280002) q[1];
sx q[1];
rz(-1.1329805) q[1];
rz(-1.8515737) q[3];
sx q[3];
rz(-2.2178136) q[3];
sx q[3];
rz(-0.10296497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48478475) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(1.8812995) q[2];
rz(-1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.242908) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(2.2743478) q[0];
rz(1.0137089) q[1];
sx q[1];
rz(-1.3288682) q[1];
sx q[1];
rz(-2.0713461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263171) q[0];
sx q[0];
rz(-1.7219825) q[0];
sx q[0];
rz(0.12355208) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4325244) q[2];
sx q[2];
rz(-0.96154172) q[2];
sx q[2];
rz(-0.23808646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6079191) q[1];
sx q[1];
rz(-1.9551139) q[1];
sx q[1];
rz(-2.584105) q[1];
rz(-2.4730014) q[3];
sx q[3];
rz(-0.1771268) q[3];
sx q[3];
rz(1.9884381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(-1.2316068) q[2];
rz(-1.6407137) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83723849) q[0];
sx q[0];
rz(-1.1501034) q[0];
sx q[0];
rz(1.3788086) q[0];
rz(2.7267743) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(-1.5157221) q[2];
sx q[2];
rz(-1.1255506) q[2];
sx q[2];
rz(-0.19565565) q[2];
rz(0.98607705) q[3];
sx q[3];
rz(-0.90860962) q[3];
sx q[3];
rz(-1.2398401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
