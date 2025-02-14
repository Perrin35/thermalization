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
rz(2.6296339) q[0];
sx q[0];
rz(-0.40357885) q[0];
sx q[0];
rz(-3.1374748) q[0];
rz(-2.4601958) q[1];
sx q[1];
rz(-1.8713142) q[1];
sx q[1];
rz(-1.4407925) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4096298) q[0];
sx q[0];
rz(-2.8131631) q[0];
sx q[0];
rz(-1.0148125) q[0];
x q[1];
rz(-2.7782562) q[2];
sx q[2];
rz(-0.78509313) q[2];
sx q[2];
rz(1.3183644) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3057249) q[1];
sx q[1];
rz(-0.9192217) q[1];
sx q[1];
rz(2.0209794) q[1];
x q[2];
rz(-3.1295092) q[3];
sx q[3];
rz(-1.7444495) q[3];
sx q[3];
rz(2.1326667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9422841) q[2];
sx q[2];
rz(-1.0465304) q[2];
sx q[2];
rz(0.47041565) q[2];
rz(2.2453902) q[3];
sx q[3];
rz(-1.6737409) q[3];
sx q[3];
rz(-1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6386221) q[0];
sx q[0];
rz(-0.48415411) q[0];
sx q[0];
rz(0.36889398) q[0];
rz(-2.6781354) q[1];
sx q[1];
rz(-1.2838485) q[1];
sx q[1];
rz(-0.2643815) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437816) q[0];
sx q[0];
rz(-0.66566313) q[0];
sx q[0];
rz(1.3206318) q[0];
rz(-pi) q[1];
rz(-0.55230852) q[2];
sx q[2];
rz(-2.8495516) q[2];
sx q[2];
rz(2.8949182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40603134) q[1];
sx q[1];
rz(-1.2792865) q[1];
sx q[1];
rz(-2.8245165) q[1];
x q[2];
rz(0.12118487) q[3];
sx q[3];
rz(-0.78599343) q[3];
sx q[3];
rz(-1.7517029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5276864) q[2];
sx q[2];
rz(-2.9759585) q[2];
sx q[2];
rz(-1.8370834) q[2];
rz(-2.8218609) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(-0.50362292) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9673135) q[0];
sx q[0];
rz(-0.19198424) q[0];
sx q[0];
rz(-2.018003) q[0];
rz(2.7536821) q[1];
sx q[1];
rz(-1.8759517) q[1];
sx q[1];
rz(-2.8111615) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1242535) q[0];
sx q[0];
rz(-1.3934008) q[0];
sx q[0];
rz(-1.3040845) q[0];
rz(-pi) q[1];
rz(1.4777884) q[2];
sx q[2];
rz(-2.1762117) q[2];
sx q[2];
rz(0.85697848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6092005) q[1];
sx q[1];
rz(-0.48621854) q[1];
sx q[1];
rz(2.4052252) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.167424) q[3];
sx q[3];
rz(-1.0181352) q[3];
sx q[3];
rz(1.0018347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72184163) q[2];
sx q[2];
rz(-1.9998877) q[2];
sx q[2];
rz(-1.9640131) q[2];
rz(-1.8925331) q[3];
sx q[3];
rz(-1.306059) q[3];
sx q[3];
rz(0.038399847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053442001) q[0];
sx q[0];
rz(-1.4295239) q[0];
sx q[0];
rz(-3.0528659) q[0];
rz(0.42452043) q[1];
sx q[1];
rz(-2.4219234) q[1];
sx q[1];
rz(1.4094062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85195711) q[0];
sx q[0];
rz(-1.9047613) q[0];
sx q[0];
rz(-0.1969239) q[0];
rz(-2.652651) q[2];
sx q[2];
rz(-1.0283111) q[2];
sx q[2];
rz(0.20636339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1115422) q[1];
sx q[1];
rz(-2.4218028) q[1];
sx q[1];
rz(-0.58776249) q[1];
rz(-pi) q[2];
rz(2.0921246) q[3];
sx q[3];
rz(-1.2884669) q[3];
sx q[3];
rz(1.2897593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7523664) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(2.7092095) q[2];
rz(-2.6835175) q[3];
sx q[3];
rz(-1.658344) q[3];
sx q[3];
rz(1.0079916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78419375) q[0];
sx q[0];
rz(-0.14506871) q[0];
sx q[0];
rz(0.31132895) q[0];
rz(-1.9642824) q[1];
sx q[1];
rz(-1.9639587) q[1];
sx q[1];
rz(0.48431531) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949849) q[0];
sx q[0];
rz(-1.7292882) q[0];
sx q[0];
rz(0.096317795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3179146) q[2];
sx q[2];
rz(-0.30159471) q[2];
sx q[2];
rz(-2.4689134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2023485) q[1];
sx q[1];
rz(-0.40196291) q[1];
sx q[1];
rz(-1.1999997) q[1];
x q[2];
rz(2.43552) q[3];
sx q[3];
rz(-2.9459402) q[3];
sx q[3];
rz(-2.5062989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3207265) q[2];
sx q[2];
rz(-1.4297239) q[2];
sx q[2];
rz(-2.7166264) q[2];
rz(3.037437) q[3];
sx q[3];
rz(-0.36915532) q[3];
sx q[3];
rz(2.4764376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686309) q[0];
sx q[0];
rz(-2.5017128) q[0];
sx q[0];
rz(0.1567008) q[0];
rz(2.8796097) q[1];
sx q[1];
rz(-1.1527088) q[1];
sx q[1];
rz(-1.6798457) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7279434) q[0];
sx q[0];
rz(-2.5143412) q[0];
sx q[0];
rz(1.2722737) q[0];
x q[1];
rz(-1.842672) q[2];
sx q[2];
rz(-0.71496925) q[2];
sx q[2];
rz(-0.78503099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20355454) q[1];
sx q[1];
rz(-2.9239836) q[1];
sx q[1];
rz(-0.59534351) q[1];
x q[2];
rz(-0.42243345) q[3];
sx q[3];
rz(-2.3452873) q[3];
sx q[3];
rz(2.5811206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67019749) q[2];
sx q[2];
rz(-1.4963701) q[2];
sx q[2];
rz(0.71693286) q[2];
rz(-0.21643058) q[3];
sx q[3];
rz(-1.8559034) q[3];
sx q[3];
rz(-1.1855679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4099429) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(1.2662079) q[0];
rz(2.3833497) q[1];
sx q[1];
rz(-1.9634602) q[1];
sx q[1];
rz(1.0479124) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3558725) q[0];
sx q[0];
rz(-1.8685172) q[0];
sx q[0];
rz(0.96426086) q[0];
x q[1];
rz(2.6897383) q[2];
sx q[2];
rz(-2.0819217) q[2];
sx q[2];
rz(-3.1060807) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0430773) q[1];
sx q[1];
rz(-1.2081573) q[1];
sx q[1];
rz(-3.0671544) q[1];
rz(-pi) q[2];
rz(2.6909184) q[3];
sx q[3];
rz(-1.3503805) q[3];
sx q[3];
rz(-2.8251951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3332112) q[2];
sx q[2];
rz(-1.0071249) q[2];
sx q[2];
rz(3.132931) q[2];
rz(-2.1374785) q[3];
sx q[3];
rz(-2.6486371) q[3];
sx q[3];
rz(2.139411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36326161) q[0];
sx q[0];
rz(-0.14334981) q[0];
sx q[0];
rz(0.95247254) q[0];
rz(1.5337503) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(1.0362157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2619963) q[0];
sx q[0];
rz(-2.41925) q[0];
sx q[0];
rz(-0.68601261) q[0];
rz(-pi) q[1];
rz(-0.63663738) q[2];
sx q[2];
rz(-1.271406) q[2];
sx q[2];
rz(2.4518397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0360603) q[1];
sx q[1];
rz(-2.1478473) q[1];
sx q[1];
rz(2.5358776) q[1];
x q[2];
rz(2.4149706) q[3];
sx q[3];
rz(-1.22261) q[3];
sx q[3];
rz(-0.49648703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8448392) q[2];
sx q[2];
rz(-0.71594683) q[2];
sx q[2];
rz(2.1412264) q[2];
rz(-1.865546) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(1.067151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.100383) q[0];
sx q[0];
rz(-0.90827933) q[0];
sx q[0];
rz(-2.8874183) q[0];
rz(-1.3379478) q[1];
sx q[1];
rz(-2.281052) q[1];
sx q[1];
rz(-2.3635704) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8786273) q[0];
sx q[0];
rz(-0.74889442) q[0];
sx q[0];
rz(0.48567943) q[0];
x q[1];
rz(0.67963375) q[2];
sx q[2];
rz(-2.8850515) q[2];
sx q[2];
rz(-0.40479615) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58147704) q[1];
sx q[1];
rz(-0.72598493) q[1];
sx q[1];
rz(-2.5037161) q[1];
x q[2];
rz(-1.144514) q[3];
sx q[3];
rz(-0.23565764) q[3];
sx q[3];
rz(0.2453177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0622056) q[2];
sx q[2];
rz(-1.5176682) q[2];
sx q[2];
rz(-0.13266955) q[2];
rz(-2.0343871) q[3];
sx q[3];
rz(-0.58834326) q[3];
sx q[3];
rz(1.4440822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74141136) q[0];
sx q[0];
rz(-1.9453229) q[0];
sx q[0];
rz(2.3924526) q[0];
rz(-0.49508849) q[1];
sx q[1];
rz(-1.2352713) q[1];
sx q[1];
rz(-1.7503768) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.063569) q[0];
sx q[0];
rz(-1.1207523) q[0];
sx q[0];
rz(0.70722945) q[0];
rz(2.3074564) q[2];
sx q[2];
rz(-2.5747402) q[2];
sx q[2];
rz(-1.2421654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0616845) q[1];
sx q[1];
rz(-0.40929738) q[1];
sx q[1];
rz(3.0319935) q[1];
rz(-2.0999317) q[3];
sx q[3];
rz(-1.1253998) q[3];
sx q[3];
rz(-0.22280773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64080015) q[2];
sx q[2];
rz(-1.4747138) q[2];
sx q[2];
rz(-0.87727171) q[2];
rz(-0.52214617) q[3];
sx q[3];
rz(-2.3418661) q[3];
sx q[3];
rz(-1.4570025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40182879) q[0];
sx q[0];
rz(-2.5943828) q[0];
sx q[0];
rz(1.7964969) q[0];
rz(-1.8632035) q[1];
sx q[1];
rz(-2.8291193) q[1];
sx q[1];
rz(3.1183174) q[1];
rz(-0.042365778) q[2];
sx q[2];
rz(-2.5103217) q[2];
sx q[2];
rz(-1.1787189) q[2];
rz(2.4458281) q[3];
sx q[3];
rz(-0.56752612) q[3];
sx q[3];
rz(-0.023434536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
