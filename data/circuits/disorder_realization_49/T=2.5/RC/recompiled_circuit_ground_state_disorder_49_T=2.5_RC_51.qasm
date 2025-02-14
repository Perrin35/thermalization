OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0030274) q[0];
sx q[0];
rz(4.0199036) q[0];
sx q[0];
rz(10.255788) q[0];
rz(2.5660958) q[1];
sx q[1];
rz(-0.73605186) q[1];
sx q[1];
rz(2.4017258) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4736126) q[0];
sx q[0];
rz(-1.4987317) q[0];
sx q[0];
rz(0.18911171) q[0];
rz(-1.3492435) q[2];
sx q[2];
rz(-2.5622517) q[2];
sx q[2];
rz(-2.960083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4783096) q[1];
sx q[1];
rz(-1.7170719) q[1];
sx q[1];
rz(-1.0895132) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3786198) q[3];
sx q[3];
rz(-1.9488244) q[3];
sx q[3];
rz(0.7188294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2089219) q[2];
sx q[2];
rz(-2.966556) q[2];
sx q[2];
rz(-0.33898655) q[2];
rz(0.50421667) q[3];
sx q[3];
rz(-2.1239069) q[3];
sx q[3];
rz(2.0174111) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9826688) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(-2.6320631) q[0];
rz(-1.5024028) q[1];
sx q[1];
rz(-2.8645611) q[1];
sx q[1];
rz(-2.1972919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751289) q[0];
sx q[0];
rz(-2.8920566) q[0];
sx q[0];
rz(0.93231045) q[0];
rz(-pi) q[1];
rz(-1.239907) q[2];
sx q[2];
rz(-0.14741405) q[2];
sx q[2];
rz(2.7380163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2118476) q[1];
sx q[1];
rz(-2.3073688) q[1];
sx q[1];
rz(-0.75726189) q[1];
rz(-pi) q[2];
rz(-0.89280309) q[3];
sx q[3];
rz(-1.2965186) q[3];
sx q[3];
rz(2.1712042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79537359) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-0.53595558) q[2];
rz(2.0761944) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(-2.4901701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029243328) q[0];
sx q[0];
rz(-1.1493244) q[0];
sx q[0];
rz(-2.405622) q[0];
rz(0.53572267) q[1];
sx q[1];
rz(-1.842272) q[1];
sx q[1];
rz(-2.2419498) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0544708) q[0];
sx q[0];
rz(-0.17987862) q[0];
sx q[0];
rz(2.4500671) q[0];
x q[1];
rz(1.9016784) q[2];
sx q[2];
rz(-0.63096672) q[2];
sx q[2];
rz(3.0598989) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6361743) q[1];
sx q[1];
rz(-2.6405425) q[1];
sx q[1];
rz(-1.4062642) q[1];
x q[2];
rz(0.24921649) q[3];
sx q[3];
rz(-1.2050516) q[3];
sx q[3];
rz(2.5178227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6003517) q[2];
sx q[2];
rz(-1.9466126) q[2];
sx q[2];
rz(0.80292732) q[2];
rz(-0.7488572) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(2.9822541) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.7298537) q[0];
sx q[0];
rz(-3.0112322) q[0];
rz(1.8761926) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(-3.0317543) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8322873) q[0];
sx q[0];
rz(-1.9067295) q[0];
sx q[0];
rz(1.0898587) q[0];
x q[1];
rz(-0.420094) q[2];
sx q[2];
rz(-0.60852988) q[2];
sx q[2];
rz(-0.15791721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7848282) q[1];
sx q[1];
rz(-2.3127529) q[1];
sx q[1];
rz(2.1020562) q[1];
rz(-2.5205344) q[3];
sx q[3];
rz(-2.406236) q[3];
sx q[3];
rz(-1.3842954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5023697) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(-2.7184674) q[2];
rz(1.0387748) q[3];
sx q[3];
rz(-2.7602502) q[3];
sx q[3];
rz(2.1151306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005242) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(0.27960676) q[0];
rz(-0.33310834) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(0.28848973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97388291) q[0];
sx q[0];
rz(-1.6339193) q[0];
sx q[0];
rz(-2.0424675) q[0];
rz(-pi) q[1];
rz(1.8158378) q[2];
sx q[2];
rz(-1.4646272) q[2];
sx q[2];
rz(2.9965056) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0049952) q[1];
sx q[1];
rz(-1.660064) q[1];
sx q[1];
rz(0.069063314) q[1];
rz(-2.1937859) q[3];
sx q[3];
rz(-1.3419328) q[3];
sx q[3];
rz(2.8826437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5491817) q[2];
sx q[2];
rz(-1.9224527) q[2];
sx q[2];
rz(0.16331095) q[2];
rz(-1.1641938) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(-1.3588847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0802245) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(1.9153216) q[0];
rz(-1.3310883) q[1];
sx q[1];
rz(-0.93003479) q[1];
sx q[1];
rz(0.61940449) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41202711) q[0];
sx q[0];
rz(-2.4270227) q[0];
sx q[0];
rz(0.78972915) q[0];
x q[1];
rz(-0.66483562) q[2];
sx q[2];
rz(-2.8241984) q[2];
sx q[2];
rz(-1.2081255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51019564) q[1];
sx q[1];
rz(-1.7185154) q[1];
sx q[1];
rz(3.0839171) q[1];
x q[2];
rz(-2.0248687) q[3];
sx q[3];
rz(-2.7259856) q[3];
sx q[3];
rz(0.46166438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2254534) q[2];
sx q[2];
rz(-1.4832387) q[2];
sx q[2];
rz(-0.43530604) q[2];
rz(3.0127323) q[3];
sx q[3];
rz(-1.83056) q[3];
sx q[3];
rz(-2.784957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211166) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(-2.5893353) q[0];
rz(0.61739677) q[1];
sx q[1];
rz(-1.2115024) q[1];
sx q[1];
rz(-1.5026106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7107781) q[0];
sx q[0];
rz(-0.39046791) q[0];
sx q[0];
rz(1.6855408) q[0];
rz(1.065997) q[2];
sx q[2];
rz(-1.4420954) q[2];
sx q[2];
rz(-1.7247049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7345846) q[1];
sx q[1];
rz(-0.53187766) q[1];
sx q[1];
rz(1.8854669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.076528744) q[3];
sx q[3];
rz(-0.26004836) q[3];
sx q[3];
rz(1.5299152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5804533) q[2];
sx q[2];
rz(-0.20623198) q[2];
sx q[2];
rz(-2.0533766) q[2];
rz(3.1214118) q[3];
sx q[3];
rz(-1.8686649) q[3];
sx q[3];
rz(-1.5816241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35767558) q[0];
sx q[0];
rz(-0.38991424) q[0];
sx q[0];
rz(1.0274603) q[0];
rz(1.1032392) q[1];
sx q[1];
rz(-1.5733066) q[1];
sx q[1];
rz(-1.2148414) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4924859) q[0];
sx q[0];
rz(-1.6739573) q[0];
sx q[0];
rz(1.2627559) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0452477) q[2];
sx q[2];
rz(-1.1602959) q[2];
sx q[2];
rz(-0.83555789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62544981) q[1];
sx q[1];
rz(-2.3770077) q[1];
sx q[1];
rz(-1.0574052) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.070511) q[3];
sx q[3];
rz(-1.0169344) q[3];
sx q[3];
rz(2.0171693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2601629) q[2];
sx q[2];
rz(-1.2036999) q[2];
sx q[2];
rz(2.0737958) q[2];
rz(2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(-1.1394181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282309) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(-0.2555787) q[0];
rz(2.3981587) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(-1.6175293) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3012863) q[0];
sx q[0];
rz(-0.69536007) q[0];
sx q[0];
rz(-2.3170503) q[0];
rz(2.5211224) q[2];
sx q[2];
rz(-1.3201642) q[2];
sx q[2];
rz(2.5567101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12518203) q[1];
sx q[1];
rz(-0.68380721) q[1];
sx q[1];
rz(1.7829624) q[1];
rz(2.2117046) q[3];
sx q[3];
rz(-1.7853338) q[3];
sx q[3];
rz(2.7492639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1114379) q[2];
sx q[2];
rz(-1.7058426) q[2];
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
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-1.2509505) q[0];
sx q[0];
rz(-2.8548456) q[0];
sx q[0];
rz(0.51837921) q[0];
rz(-2.3333343) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(-1.4498651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5980412) q[0];
sx q[0];
rz(-1.2819474) q[0];
sx q[0];
rz(-0.01057118) q[0];
x q[1];
rz(0.6459562) q[2];
sx q[2];
rz(-2.5017512) q[2];
sx q[2];
rz(-1.3825939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5603578) q[1];
sx q[1];
rz(-2.5311845) q[1];
sx q[1];
rz(-2.6950652) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57761044) q[3];
sx q[3];
rz(-2.0311714) q[3];
sx q[3];
rz(-2.1652997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1125674) q[2];
sx q[2];
rz(-1.2351278) q[2];
sx q[2];
rz(0.9355363) q[2];
rz(1.6701291) q[3];
sx q[3];
rz(-1.4751438) q[3];
sx q[3];
rz(2.0475625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6846631) q[0];
sx q[0];
rz(-2.5419432) q[0];
sx q[0];
rz(2.3210617) q[0];
rz(2.6928071) q[1];
sx q[1];
rz(-1.1460591) q[1];
sx q[1];
rz(-1.9312327) q[1];
rz(-1.4478499) q[2];
sx q[2];
rz(-2.7715383) q[2];
sx q[2];
rz(-0.72889974) q[2];
rz(1.2294235) q[3];
sx q[3];
rz(-2.1044272) q[3];
sx q[3];
rz(2.9499346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
