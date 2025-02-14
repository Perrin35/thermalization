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
rz(-2.2632817) q[0];
sx q[0];
rz(0.83100975) q[0];
rz(-0.57549685) q[1];
sx q[1];
rz(3.8776445) q[1];
sx q[1];
rz(10.164645) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66798009) q[0];
sx q[0];
rz(-1.6428609) q[0];
sx q[0];
rz(0.18911171) q[0];
x q[1];
rz(1.0027356) q[2];
sx q[2];
rz(-1.6913927) q[2];
sx q[2];
rz(-1.9386148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5062472) q[1];
sx q[1];
rz(-0.50134778) q[1];
sx q[1];
rz(1.8789199) q[1];
rz(-pi) q[2];
rz(-1.7629729) q[3];
sx q[3];
rz(-1.1927682) q[3];
sx q[3];
rz(0.7188294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2089219) q[2];
sx q[2];
rz(-0.17503665) q[2];
sx q[2];
rz(0.33898655) q[2];
rz(0.50421667) q[3];
sx q[3];
rz(-2.1239069) q[3];
sx q[3];
rz(2.0174111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9826688) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(0.50952953) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(0.94430077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751289) q[0];
sx q[0];
rz(-0.24953609) q[0];
sx q[0];
rz(-0.93231045) q[0];
x q[1];
rz(-1.239907) q[2];
sx q[2];
rz(-2.9941786) q[2];
sx q[2];
rz(0.40357631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2118476) q[1];
sx q[1];
rz(-2.3073688) q[1];
sx q[1];
rz(-0.75726189) q[1];
rz(-1.1491165) q[3];
sx q[3];
rz(-2.4184368) q[3];
sx q[3];
rz(2.2167517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79537359) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-0.53595558) q[2];
rz(-2.0761944) q[3];
sx q[3];
rz(-2.8841116) q[3];
sx q[3];
rz(0.65142256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029243328) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(0.73597062) q[0];
rz(-2.60587) q[1];
sx q[1];
rz(-1.842272) q[1];
sx q[1];
rz(-2.2419498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0544708) q[0];
sx q[0];
rz(-0.17987862) q[0];
sx q[0];
rz(2.4500671) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96615919) q[2];
sx q[2];
rz(-1.3779497) q[2];
sx q[2];
rz(-1.218534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31832987) q[1];
sx q[1];
rz(-2.0644651) q[1];
sx q[1];
rz(-3.0521293) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1429575) q[3];
sx q[3];
rz(-2.7021926) q[3];
sx q[3];
rz(1.8993401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54124093) q[2];
sx q[2];
rz(-1.1949801) q[2];
sx q[2];
rz(0.80292732) q[2];
rz(2.3927355) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(2.9822541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
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
rz(-1.9644901) q[1];
sx q[1];
rz(-0.1098384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3093053) q[0];
sx q[0];
rz(-1.9067295) q[0];
sx q[0];
rz(-1.0898587) q[0];
rz(-pi) q[1];
rz(-0.56657882) q[2];
sx q[2];
rz(-1.3354805) q[2];
sx q[2];
rz(1.7641774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.500587) q[1];
sx q[1];
rz(-0.88208032) q[1];
sx q[1];
rz(-2.6367049) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5073152) q[3];
sx q[3];
rz(-1.1697672) q[3];
sx q[3];
rz(2.4672535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5023697) q[2];
sx q[2];
rz(-1.9618192) q[2];
sx q[2];
rz(-0.42312527) q[2];
rz(2.1028178) q[3];
sx q[3];
rz(-2.7602502) q[3];
sx q[3];
rz(-2.1151306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005242) q[0];
sx q[0];
rz(-1.2491913) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(-0.33310834) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(0.28848973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4216327) q[0];
sx q[0];
rz(-0.47556092) q[0];
sx q[0];
rz(1.4325761) q[0];
x q[1];
rz(-1.9847068) q[2];
sx q[2];
rz(-2.8749646) q[2];
sx q[2];
rz(1.8265343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1365974) q[1];
sx q[1];
rz(-1.660064) q[1];
sx q[1];
rz(-0.069063314) q[1];
rz(-pi) q[2];
rz(0.94780677) q[3];
sx q[3];
rz(-1.3419328) q[3];
sx q[3];
rz(2.8826437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5491817) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(2.9782817) q[2];
rz(1.9773989) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(-1.3588847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(1.3310883) q[1];
sx q[1];
rz(-0.93003479) q[1];
sx q[1];
rz(-0.61940449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41202711) q[0];
sx q[0];
rz(-2.4270227) q[0];
sx q[0];
rz(2.3518635) q[0];
rz(-pi) q[1];
rz(-2.476757) q[2];
sx q[2];
rz(-2.8241984) q[2];
sx q[2];
rz(-1.9334671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51019564) q[1];
sx q[1];
rz(-1.7185154) q[1];
sx q[1];
rz(-3.0839171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.193229) q[3];
sx q[3];
rz(-1.7488297) q[3];
sx q[3];
rz(-1.61249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9161393) q[2];
sx q[2];
rz(-1.658354) q[2];
sx q[2];
rz(0.43530604) q[2];
rz(-3.0127323) q[3];
sx q[3];
rz(-1.3110327) q[3];
sx q[3];
rz(0.35663566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211166) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(0.55225736) q[0];
rz(-0.61739677) q[1];
sx q[1];
rz(-1.2115024) q[1];
sx q[1];
rz(-1.6389821) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4308145) q[0];
sx q[0];
rz(-2.7511247) q[0];
sx q[0];
rz(1.6855408) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8322629) q[2];
sx q[2];
rz(-0.51957031) q[2];
sx q[2];
rz(-0.38214035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.43728033) q[1];
sx q[1];
rz(-1.7284135) q[1];
sx q[1];
rz(-2.0809409) q[1];
rz(-pi) q[2];
rz(1.5504567) q[3];
sx q[3];
rz(-1.830066) q[3];
sx q[3];
rz(-1.532497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5611394) q[2];
sx q[2];
rz(-2.9353607) q[2];
sx q[2];
rz(1.0882161) q[2];
rz(-0.020180833) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(1.5816241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35767558) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(-2.1141323) q[0];
rz(-1.1032392) q[1];
sx q[1];
rz(-1.5733066) q[1];
sx q[1];
rz(-1.9267513) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64910674) q[0];
sx q[0];
rz(-1.6739573) q[0];
sx q[0];
rz(1.8788368) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0452477) q[2];
sx q[2];
rz(-1.1602959) q[2];
sx q[2];
rz(-0.83555789) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5161428) q[1];
sx q[1];
rz(-2.3770077) q[1];
sx q[1];
rz(1.0574052) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4818871) q[3];
sx q[3];
rz(-0.72835975) q[3];
sx q[3];
rz(1.212478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2601629) q[2];
sx q[2];
rz(-1.2036999) q[2];
sx q[2];
rz(2.0737958) q[2];
rz(2.2504375) q[3];
sx q[3];
rz(-0.46025899) q[3];
sx q[3];
rz(-2.0021745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-2.886014) q[0];
rz(-2.3981587) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(1.6175293) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96252464) q[0];
sx q[0];
rz(-1.0810548) q[0];
sx q[0];
rz(0.51539246) q[0];
rz(-pi) q[1];
rz(-1.8756549) q[2];
sx q[2];
rz(-0.97248021) q[2];
sx q[2];
rz(1.9802633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39622722) q[1];
sx q[1];
rz(-2.2364625) q[1];
sx q[1];
rz(0.16996064) q[1];
rz(1.9202523) q[3];
sx q[3];
rz(-0.67103681) q[3];
sx q[3];
rz(-0.90045917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1114379) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(1.5343522) q[2];
rz(2.0700908) q[3];
sx q[3];
rz(-1.6000308) q[3];
sx q[3];
rz(1.967954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8906422) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(2.3333343) q[1];
sx q[1];
rz(-1.6812485) q[1];
sx q[1];
rz(-1.4498651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5064556) q[0];
sx q[0];
rz(-0.28903693) q[0];
sx q[0];
rz(-1.6063547) q[0];
rz(-pi) q[1];
rz(2.6053455) q[2];
sx q[2];
rz(-1.203158) q[2];
sx q[2];
rz(-2.4095031) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0891493) q[1];
sx q[1];
rz(-1.027453) q[1];
sx q[1];
rz(1.2774317) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0363186) q[3];
sx q[3];
rz(-2.0819398) q[3];
sx q[3];
rz(2.8289464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1125674) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(2.2060564) q[2];
rz(-1.4714636) q[3];
sx q[3];
rz(-1.6664489) q[3];
sx q[3];
rz(1.0940301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6846631) q[0];
sx q[0];
rz(-0.59964947) q[0];
sx q[0];
rz(-0.82053091) q[0];
rz(-2.6928071) q[1];
sx q[1];
rz(-1.9955336) q[1];
sx q[1];
rz(1.21036) q[1];
rz(-1.4478499) q[2];
sx q[2];
rz(-2.7715383) q[2];
sx q[2];
rz(-0.72889974) q[2];
rz(-0.56002496) q[3];
sx q[3];
rz(-1.2784132) q[3];
sx q[3];
rz(-1.9412435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
