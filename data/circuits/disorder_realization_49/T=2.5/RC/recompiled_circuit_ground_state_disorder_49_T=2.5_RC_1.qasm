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
rz(-0.83100975) q[0];
rz(-0.57549685) q[1];
sx q[1];
rz(3.8776445) q[1];
sx q[1];
rz(10.164645) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4736126) q[0];
sx q[0];
rz(-1.6428609) q[0];
sx q[0];
rz(2.9524809) q[0];
rz(-1.0027356) q[2];
sx q[2];
rz(-1.4502) q[2];
sx q[2];
rz(1.2029778) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.98348599) q[1];
sx q[1];
rz(-2.0465104) q[1];
sx q[1];
rz(2.9768894) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3786198) q[3];
sx q[3];
rz(-1.1927682) q[3];
sx q[3];
rz(-2.4227633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9326707) q[2];
sx q[2];
rz(-0.17503665) q[2];
sx q[2];
rz(-0.33898655) q[2];
rz(-2.637376) q[3];
sx q[3];
rz(-2.1239069) q[3];
sx q[3];
rz(2.0174111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9826688) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(-0.50952953) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-2.8645611) q[1];
sx q[1];
rz(-0.94430077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0137607) q[0];
sx q[0];
rz(-1.7185129) q[0];
sx q[0];
rz(-1.7726519) q[0];
x q[1];
rz(-1.71032) q[2];
sx q[2];
rz(-1.6185337) q[2];
sx q[2];
rz(0.83966694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1183853) q[1];
sx q[1];
rz(-2.1392576) q[1];
sx q[1];
rz(2.2190907) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7949132) q[3];
sx q[3];
rz(-2.2190385) q[3];
sx q[3];
rz(0.38564206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79537359) q[2];
sx q[2];
rz(-1.5920762) q[2];
sx q[2];
rz(-0.53595558) q[2];
rz(2.0761944) q[3];
sx q[3];
rz(-2.8841116) q[3];
sx q[3];
rz(-0.65142256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(-2.405622) q[0];
rz(-2.60587) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(-0.89964286) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7539983) q[0];
sx q[0];
rz(-1.7090461) q[0];
sx q[0];
rz(1.686245) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96615919) q[2];
sx q[2];
rz(-1.7636429) q[2];
sx q[2];
rz(1.9230587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8232628) q[1];
sx q[1];
rz(-2.0644651) q[1];
sx q[1];
rz(-0.089463316) q[1];
rz(0.24921649) q[3];
sx q[3];
rz(-1.936541) q[3];
sx q[3];
rz(-2.5178227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6003517) q[2];
sx q[2];
rz(-1.1949801) q[2];
sx q[2];
rz(-2.3386653) q[2];
rz(0.7488572) q[3];
sx q[3];
rz(-1.7436946) q[3];
sx q[3];
rz(-0.15933855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51711851) q[0];
sx q[0];
rz(-1.7298537) q[0];
sx q[0];
rz(-0.13036048) q[0];
rz(-1.2654001) q[1];
sx q[1];
rz(-1.9644901) q[1];
sx q[1];
rz(-0.1098384) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097561) q[0];
sx q[0];
rz(-2.0227814) q[0];
sx q[0];
rz(0.37518895) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.420094) q[2];
sx q[2];
rz(-0.60852988) q[2];
sx q[2];
rz(2.9836754) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.500587) q[1];
sx q[1];
rz(-2.2595123) q[1];
sx q[1];
rz(0.50488774) q[1];
x q[2];
rz(-2.5205344) q[3];
sx q[3];
rz(-0.73535669) q[3];
sx q[3];
rz(-1.7572973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.639223) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(0.42312527) q[2];
rz(-1.0387748) q[3];
sx q[3];
rz(-0.3813425) q[3];
sx q[3];
rz(2.1151306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410685) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(2.8084843) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(-2.8531029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56474876) q[0];
sx q[0];
rz(-2.0414519) q[0];
sx q[0];
rz(-0.0708357) q[0];
rz(-pi) q[1];
rz(1.1568858) q[2];
sx q[2];
rz(-0.26662808) q[2];
sx q[2];
rz(1.3150584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5696249) q[1];
sx q[1];
rz(-1.5020084) q[1];
sx q[1];
rz(-1.6602762) q[1];
rz(-1.1909423) q[3];
sx q[3];
rz(-2.4831746) q[3];
sx q[3];
rz(1.0059716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5491817) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(0.16331095) q[2];
rz(1.9773989) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(1.7827079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50766599) q[0];
sx q[0];
rz(-2.0548425) q[0];
sx q[0];
rz(-2.593301) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2529958) q[2];
sx q[2];
rz(-1.7645451) q[2];
sx q[2];
rz(0.27744833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0521026) q[1];
sx q[1];
rz(-1.627843) q[1];
sx q[1];
rz(-1.7187579) q[1];
rz(-0.19120817) q[3];
sx q[3];
rz(-1.9421028) q[3];
sx q[3];
rz(-3.1131668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2254534) q[2];
sx q[2];
rz(-1.658354) q[2];
sx q[2];
rz(2.7062866) q[2];
rz(0.12886038) q[3];
sx q[3];
rz(-1.3110327) q[3];
sx q[3];
rz(0.35663566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211166) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(-2.5893353) q[0];
rz(2.5241959) q[1];
sx q[1];
rz(-1.2115024) q[1];
sx q[1];
rz(-1.6389821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7107781) q[0];
sx q[0];
rz(-2.7511247) q[0];
sx q[0];
rz(1.6855408) q[0];
rz(-pi) q[1];
rz(-2.9947979) q[2];
sx q[2];
rz(-1.0705612) q[2];
sx q[2];
rz(-0.083109476) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43728033) q[1];
sx q[1];
rz(-1.7284135) q[1];
sx q[1];
rz(-1.0606517) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8822717) q[3];
sx q[3];
rz(-1.590456) q[3];
sx q[3];
rz(0.033084083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5611394) q[2];
sx q[2];
rz(-0.20623198) q[2];
sx q[2];
rz(1.0882161) q[2];
rz(-0.020180833) q[3];
sx q[3];
rz(-1.8686649) q[3];
sx q[3];
rz(-1.5816241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35767558) q[0];
sx q[0];
rz(-0.38991424) q[0];
sx q[0];
rz(-1.0274603) q[0];
rz(-2.0383535) q[1];
sx q[1];
rz(-1.5682861) q[1];
sx q[1];
rz(1.2148414) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2347082) q[0];
sx q[0];
rz(-2.8172593) q[0];
sx q[0];
rz(-1.8998446) q[0];
rz(-pi) q[1];
rz(0.45502624) q[2];
sx q[2];
rz(-1.1385673) q[2];
sx q[2];
rz(-2.2042008) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.852927) q[1];
sx q[1];
rz(-0.92354362) q[1];
sx q[1];
rz(-2.701328) q[1];
rz(-pi) q[2];
rz(-2.4818871) q[3];
sx q[3];
rz(-0.72835975) q[3];
sx q[3];
rz(1.9291147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2601629) q[2];
sx q[2];
rz(-1.2036999) q[2];
sx q[2];
rz(2.0737958) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(1.1394181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282309) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(2.886014) q[0];
rz(-0.74343395) q[1];
sx q[1];
rz(-1.5579222) q[1];
sx q[1];
rz(-1.5240634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96252464) q[0];
sx q[0];
rz(-1.0810548) q[0];
sx q[0];
rz(2.6262002) q[0];
x q[1];
rz(-0.62047024) q[2];
sx q[2];
rz(-1.3201642) q[2];
sx q[2];
rz(-0.58488256) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39622722) q[1];
sx q[1];
rz(-2.2364625) q[1];
sx q[1];
rz(-2.971632) q[1];
x q[2];
rz(2.2117046) q[3];
sx q[3];
rz(-1.7853338) q[3];
sx q[3];
rz(-0.39232871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1114379) q[2];
sx q[2];
rz(-1.7058426) q[2];
sx q[2];
rz(-1.6072404) q[2];
rz(1.0715019) q[3];
sx q[3];
rz(-1.5415618) q[3];
sx q[3];
rz(-1.1736386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8906422) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(-0.80825835) q[1];
sx q[1];
rz(-1.6812485) q[1];
sx q[1];
rz(1.6917276) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302561) q[0];
sx q[0];
rz(-1.5809296) q[0];
sx q[0];
rz(1.8596605) q[0];
rz(-2.4956365) q[2];
sx q[2];
rz(-0.63984144) q[2];
sx q[2];
rz(1.3825939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7781642) q[1];
sx q[1];
rz(-1.320667) q[1];
sx q[1];
rz(2.5787337) q[1];
rz(-pi) q[2];
rz(-0.57761044) q[3];
sx q[3];
rz(-1.1104212) q[3];
sx q[3];
rz(-2.1652997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1125674) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(-0.9355363) q[2];
rz(1.6701291) q[3];
sx q[3];
rz(-1.4751438) q[3];
sx q[3];
rz(2.0475625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6846631) q[0];
sx q[0];
rz(-0.59964947) q[0];
sx q[0];
rz(-0.82053091) q[0];
rz(-0.44878557) q[1];
sx q[1];
rz(-1.1460591) q[1];
sx q[1];
rz(-1.9312327) q[1];
rz(1.4478499) q[2];
sx q[2];
rz(-0.37005432) q[2];
sx q[2];
rz(2.4126929) q[2];
rz(-2.5815677) q[3];
sx q[3];
rz(-1.8631794) q[3];
sx q[3];
rz(1.2003492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
