OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1385652) q[0];
sx q[0];
rz(-0.87831098) q[0];
sx q[0];
rz(2.3105829) q[0];
rz(2.5660958) q[1];
sx q[1];
rz(-0.73605186) q[1];
sx q[1];
rz(-0.73986685) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66798009) q[0];
sx q[0];
rz(-1.4987317) q[0];
sx q[0];
rz(2.9524809) q[0];
rz(2.1388571) q[2];
sx q[2];
rz(-1.4502) q[2];
sx q[2];
rz(1.2029778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5062472) q[1];
sx q[1];
rz(-0.50134778) q[1];
sx q[1];
rz(1.8789199) q[1];
x q[2];
rz(1.7629729) q[3];
sx q[3];
rz(-1.9488244) q[3];
sx q[3];
rz(0.7188294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2089219) q[2];
sx q[2];
rz(-0.17503665) q[2];
sx q[2];
rz(-0.33898655) q[2];
rz(-2.637376) q[3];
sx q[3];
rz(-1.0176858) q[3];
sx q[3];
rz(-2.0174111) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15892383) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(0.50952953) q[0];
rz(1.6391899) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(-0.94430077) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4128542) q[0];
sx q[0];
rz(-1.7704238) q[0];
sx q[0];
rz(2.9908604) q[0];
x q[1];
rz(0.048205094) q[2];
sx q[2];
rz(-1.4314326) q[2];
sx q[2];
rz(-2.4037619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2118476) q[1];
sx q[1];
rz(-2.3073688) q[1];
sx q[1];
rz(0.75726189) q[1];
x q[2];
rz(-0.89280309) q[3];
sx q[3];
rz(-1.2965186) q[3];
sx q[3];
rz(2.1712042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79537359) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-0.53595558) q[2];
rz(2.0761944) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(0.65142256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(-2.405622) q[0];
rz(-0.53572267) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(0.89964286) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9743703) q[0];
sx q[0];
rz(-1.6851387) q[0];
sx q[0];
rz(0.13916441) q[0];
rz(-1.2399142) q[2];
sx q[2];
rz(-2.5106259) q[2];
sx q[2];
rz(-3.0598989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8232628) q[1];
sx q[1];
rz(-2.0644651) q[1];
sx q[1];
rz(3.0521293) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1944481) q[3];
sx q[3];
rz(-1.8032142) q[3];
sx q[3];
rz(-2.2853394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54124093) q[2];
sx q[2];
rz(-1.9466126) q[2];
sx q[2];
rz(-0.80292732) q[2];
rz(2.3927355) q[3];
sx q[3];
rz(-1.7436946) q[3];
sx q[3];
rz(-2.9822541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.4117389) q[0];
sx q[0];
rz(-0.13036048) q[0];
rz(-1.2654001) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(-3.0317543) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8396436) q[0];
sx q[0];
rz(-0.57900864) q[0];
sx q[0];
rz(0.92424519) q[0];
rz(-pi) q[1];
x q[1];
rz(0.420094) q[2];
sx q[2];
rz(-0.60852988) q[2];
sx q[2];
rz(-2.9836754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59200724) q[1];
sx q[1];
rz(-1.9535258) q[1];
sx q[1];
rz(0.81603433) q[1];
x q[2];
rz(0.63427744) q[3];
sx q[3];
rz(-1.1697672) q[3];
sx q[3];
rz(2.4672535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.639223) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(2.7184674) q[2];
rz(2.1028178) q[3];
sx q[3];
rz(-0.3813425) q[3];
sx q[3];
rz(2.1151306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410685) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(-2.8084843) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(2.8531029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4216327) q[0];
sx q[0];
rz(-0.47556092) q[0];
sx q[0];
rz(-1.4325761) q[0];
rz(-pi) q[1];
rz(-1.1568858) q[2];
sx q[2];
rz(-0.26662808) q[2];
sx q[2];
rz(-1.3150584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0049952) q[1];
sx q[1];
rz(-1.4815287) q[1];
sx q[1];
rz(-3.0725293) q[1];
x q[2];
rz(1.1909423) q[3];
sx q[3];
rz(-2.4831746) q[3];
sx q[3];
rz(-1.0059716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5491817) q[2];
sx q[2];
rz(-1.9224527) q[2];
sx q[2];
rz(0.16331095) q[2];
rz(-1.9773989) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(-1.7827079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.0613681) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(-1.226271) q[0];
rz(-1.3310883) q[1];
sx q[1];
rz(-0.93003479) q[1];
sx q[1];
rz(0.61940449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6339267) q[0];
sx q[0];
rz(-1.0867501) q[0];
sx q[0];
rz(2.593301) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8885968) q[2];
sx q[2];
rz(-1.3770475) q[2];
sx q[2];
rz(-0.27744833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2575469) q[1];
sx q[1];
rz(-2.9830898) q[1];
sx q[1];
rz(-1.201215) q[1];
rz(-pi) q[2];
rz(-0.19120817) q[3];
sx q[3];
rz(-1.9421028) q[3];
sx q[3];
rz(-3.1131668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2254534) q[2];
sx q[2];
rz(-1.658354) q[2];
sx q[2];
rz(-2.7062866) q[2];
rz(-3.0127323) q[3];
sx q[3];
rz(-1.3110327) q[3];
sx q[3];
rz(-2.784957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72047609) q[0];
sx q[0];
rz(-2.5832376) q[0];
sx q[0];
rz(2.5893353) q[0];
rz(0.61739677) q[1];
sx q[1];
rz(-1.9300902) q[1];
sx q[1];
rz(1.5026106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033806697) q[0];
sx q[0];
rz(-1.5272041) q[0];
sx q[0];
rz(-1.1826452) q[0];
rz(-pi) q[1];
rz(0.1467948) q[2];
sx q[2];
rz(-1.0705612) q[2];
sx q[2];
rz(3.0584832) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4070081) q[1];
sx q[1];
rz(-2.609715) q[1];
sx q[1];
rz(1.8854669) q[1];
rz(0.076528744) q[3];
sx q[3];
rz(-0.26004836) q[3];
sx q[3];
rz(-1.5299152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5611394) q[2];
sx q[2];
rz(-0.20623198) q[2];
sx q[2];
rz(1.0882161) q[2];
rz(-3.1214118) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(1.5599686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7839171) q[0];
sx q[0];
rz(-0.38991424) q[0];
sx q[0];
rz(-2.1141323) q[0];
rz(-1.1032392) q[1];
sx q[1];
rz(-1.5733066) q[1];
sx q[1];
rz(-1.9267513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2347082) q[0];
sx q[0];
rz(-0.32433332) q[0];
sx q[0];
rz(1.2417481) q[0];
rz(-0.45502624) q[2];
sx q[2];
rz(-1.1385673) q[2];
sx q[2];
rz(2.2042008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5826792) q[1];
sx q[1];
rz(-1.2238992) q[1];
sx q[1];
rz(2.2668462) q[1];
rz(-2.4818871) q[3];
sx q[3];
rz(-0.72835975) q[3];
sx q[3];
rz(-1.212478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8814298) q[2];
sx q[2];
rz(-1.9378928) q[2];
sx q[2];
rz(-1.0677968) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(1.1394181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.41336173) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(0.2555787) q[0];
rz(2.3981587) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(1.5240634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7937318) q[0];
sx q[0];
rz(-1.1208236) q[0];
sx q[0];
rz(-1.0211358) q[0];
rz(-pi) q[1];
rz(2.5211224) q[2];
sx q[2];
rz(-1.3201642) q[2];
sx q[2];
rz(2.5567101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7453654) q[1];
sx q[1];
rz(-2.2364625) q[1];
sx q[1];
rz(-2.971632) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9202523) q[3];
sx q[3];
rz(-0.67103681) q[3];
sx q[3];
rz(-2.2411335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1114379) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(-1.5343522) q[2];
rz(-2.0700908) q[3];
sx q[3];
rz(-1.6000308) q[3];
sx q[3];
rz(-1.967954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8906422) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(0.51837921) q[0];
rz(-0.80825835) q[1];
sx q[1];
rz(-1.6812485) q[1];
sx q[1];
rz(-1.4498651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5980412) q[0];
sx q[0];
rz(-1.8596453) q[0];
sx q[0];
rz(0.01057118) q[0];
rz(2.4956365) q[2];
sx q[2];
rz(-2.5017512) q[2];
sx q[2];
rz(1.3825939) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3634285) q[1];
sx q[1];
rz(-1.320667) q[1];
sx q[1];
rz(2.5787337) q[1];
x q[2];
rz(-2.4042481) q[3];
sx q[3];
rz(-2.4196845) q[3];
sx q[3];
rz(1.9487716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0290252) q[2];
sx q[2];
rz(-1.2351278) q[2];
sx q[2];
rz(-0.9355363) q[2];
rz(-1.4714636) q[3];
sx q[3];
rz(-1.6664489) q[3];
sx q[3];
rz(1.0940301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4569296) q[0];
sx q[0];
rz(-0.59964947) q[0];
sx q[0];
rz(-0.82053091) q[0];
rz(-2.6928071) q[1];
sx q[1];
rz(-1.9955336) q[1];
sx q[1];
rz(1.21036) q[1];
rz(1.6937428) q[2];
sx q[2];
rz(-2.7715383) q[2];
sx q[2];
rz(-0.72889974) q[2];
rz(0.56002496) q[3];
sx q[3];
rz(-1.8631794) q[3];
sx q[3];
rz(1.2003492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
