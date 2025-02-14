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
rz(0.91659651) q[0];
sx q[0];
rz(-1.7594113) q[0];
sx q[0];
rz(1.4974282) q[0];
rz(-2.1388571) q[2];
sx q[2];
rz(-1.4502) q[2];
sx q[2];
rz(-1.2029778) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5062472) q[1];
sx q[1];
rz(-0.50134778) q[1];
sx q[1];
rz(-1.2626727) q[1];
rz(-2.7571477) q[3];
sx q[3];
rz(-1.3923402) q[3];
sx q[3];
rz(-0.78027356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2089219) q[2];
sx q[2];
rz(-2.966556) q[2];
sx q[2];
rz(0.33898655) q[2];
rz(0.50421667) q[3];
sx q[3];
rz(-1.0176858) q[3];
sx q[3];
rz(1.1241815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9826688) q[0];
sx q[0];
rz(-0.6944387) q[0];
sx q[0];
rz(-2.6320631) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-2.8645611) q[1];
sx q[1];
rz(-0.94430077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12783192) q[0];
sx q[0];
rz(-1.4230797) q[0];
sx q[0];
rz(-1.7726519) q[0];
x q[1];
rz(-3.0933876) q[2];
sx q[2];
rz(-1.71016) q[2];
sx q[2];
rz(-0.73783079) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9348) q[1];
sx q[1];
rz(-2.1046608) q[1];
sx q[1];
rz(2.4660048) q[1];
x q[2];
rz(0.89280309) q[3];
sx q[3];
rz(-1.8450741) q[3];
sx q[3];
rz(-0.97038847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3462191) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-2.6056371) q[2];
rz(-1.0653982) q[3];
sx q[3];
rz(-2.8841116) q[3];
sx q[3];
rz(-0.65142256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(0.73597062) q[0];
rz(-0.53572267) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(0.89964286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1672223) q[0];
sx q[0];
rz(-1.6851387) q[0];
sx q[0];
rz(-3.0024282) q[0];
rz(-0.23304184) q[2];
sx q[2];
rz(-2.1626806) q[2];
sx q[2];
rz(0.48392236) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9316072) q[1];
sx q[1];
rz(-1.6495541) q[1];
sx q[1];
rz(-1.0754536) q[1];
x q[2];
rz(0.24921649) q[3];
sx q[3];
rz(-1.2050516) q[3];
sx q[3];
rz(-0.62376991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54124093) q[2];
sx q[2];
rz(-1.9466126) q[2];
sx q[2];
rz(2.3386653) q[2];
rz(2.3927355) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(-0.15933855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.7298537) q[0];
sx q[0];
rz(-3.0112322) q[0];
rz(-1.2654001) q[1];
sx q[1];
rz(-1.9644901) q[1];
sx q[1];
rz(-0.1098384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30194908) q[0];
sx q[0];
rz(-0.57900864) q[0];
sx q[0];
rz(-2.2173475) q[0];
rz(-pi) q[1];
rz(0.56657882) q[2];
sx q[2];
rz(-1.8061122) q[2];
sx q[2];
rz(1.7641774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7848282) q[1];
sx q[1];
rz(-2.3127529) q[1];
sx q[1];
rz(2.1020562) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5205344) q[3];
sx q[3];
rz(-0.73535669) q[3];
sx q[3];
rz(1.7572973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5023697) q[2];
sx q[2];
rz(-1.9618192) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005242) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(-0.27960676) q[0];
rz(0.33310834) q[1];
sx q[1];
rz(-0.50222841) q[1];
sx q[1];
rz(0.28848973) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97388291) q[0];
sx q[0];
rz(-1.6339193) q[0];
sx q[0];
rz(-1.0991251) q[0];
rz(-pi) q[1];
rz(-3.0321799) q[2];
sx q[2];
rz(-1.3271626) q[2];
sx q[2];
rz(1.7423767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0049952) q[1];
sx q[1];
rz(-1.4815287) q[1];
sx q[1];
rz(-3.0725293) q[1];
rz(-1.9506504) q[3];
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
rz(-1.592411) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(2.9782817) q[2];
rz(-1.9773989) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(1.3588847) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613681) q[0];
sx q[0];
rz(-0.017266406) q[0];
sx q[0];
rz(1.9153216) q[0];
rz(1.3310883) q[1];
sx q[1];
rz(-2.2115579) q[1];
sx q[1];
rz(-2.5221882) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50766599) q[0];
sx q[0];
rz(-2.0548425) q[0];
sx q[0];
rz(-2.593301) q[0];
x q[1];
rz(-2.476757) q[2];
sx q[2];
rz(-0.31739429) q[2];
sx q[2];
rz(1.9334671) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.631397) q[1];
sx q[1];
rz(-1.7185154) q[1];
sx q[1];
rz(3.0839171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9483637) q[3];
sx q[3];
rz(-1.7488297) q[3];
sx q[3];
rz(-1.5291027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72047609) q[0];
sx q[0];
rz(-2.5832376) q[0];
sx q[0];
rz(-0.55225736) q[0];
rz(-0.61739677) q[1];
sx q[1];
rz(-1.9300902) q[1];
sx q[1];
rz(-1.5026106) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5548067) q[0];
sx q[0];
rz(-1.1830336) q[0];
sx q[0];
rz(-3.094502) q[0];
rz(1.8322629) q[2];
sx q[2];
rz(-0.51957031) q[2];
sx q[2];
rz(0.38214035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0459111) q[1];
sx q[1];
rz(-2.0740182) q[1];
sx q[1];
rz(2.9614425) q[1];
x q[2];
rz(-0.25932094) q[3];
sx q[3];
rz(-1.5511366) q[3];
sx q[3];
rz(-3.1085086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5611394) q[2];
sx q[2];
rz(-2.9353607) q[2];
sx q[2];
rz(-2.0533766) q[2];
rz(3.1214118) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(1.5816241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7839171) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(-1.0274603) q[0];
rz(-2.0383535) q[1];
sx q[1];
rz(-1.5682861) q[1];
sx q[1];
rz(-1.9267513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64910674) q[0];
sx q[0];
rz(-1.4676354) q[0];
sx q[0];
rz(1.8788368) q[0];
rz(-pi) q[1];
rz(1.096345) q[2];
sx q[2];
rz(-1.1602959) q[2];
sx q[2];
rz(-0.83555789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2886656) q[1];
sx q[1];
rz(-2.218049) q[1];
sx q[1];
rz(-0.44026466) q[1];
rz(-pi) q[2];
rz(2.4818871) q[3];
sx q[3];
rz(-2.4132329) q[3];
sx q[3];
rz(1.9291147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8814298) q[2];
sx q[2];
rz(-1.2036999) q[2];
sx q[2];
rz(-1.0677968) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(-2.0021745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282309) q[0];
sx q[0];
rz(-0.75465337) q[0];
sx q[0];
rz(0.2555787) q[0];
rz(-2.3981587) q[1];
sx q[1];
rz(-1.5579222) q[1];
sx q[1];
rz(-1.6175293) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34786087) q[0];
sx q[0];
rz(-1.1208236) q[0];
sx q[0];
rz(2.1204569) q[0];
rz(-pi) q[1];
rz(-2.7268098) q[2];
sx q[2];
rz(-2.4786502) q[2];
sx q[2];
rz(-2.489733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0164106) q[1];
sx q[1];
rz(-2.4577854) q[1];
sx q[1];
rz(-1.7829624) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92988807) q[3];
sx q[3];
rz(-1.7853338) q[3];
sx q[3];
rz(-0.39232871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0301547) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(1.6072404) q[2];
rz(1.0715019) q[3];
sx q[3];
rz(-1.5415618) q[3];
sx q[3];
rz(1.967954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2509505) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(2.3333343) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(1.4498651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54355147) q[0];
sx q[0];
rz(-1.2819474) q[0];
sx q[0];
rz(-0.01057118) q[0];
rz(-1.9920182) q[2];
sx q[2];
rz(-1.0738157) q[2];
sx q[2];
rz(0.62825655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7781642) q[1];
sx q[1];
rz(-1.320667) q[1];
sx q[1];
rz(-2.5787337) q[1];
rz(2.1052741) q[3];
sx q[3];
rz(-1.0596529) q[3];
sx q[3];
rz(-0.3126463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(2.2060564) q[2];
rz(1.6701291) q[3];
sx q[3];
rz(-1.6664489) q[3];
sx q[3];
rz(1.0940301) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.2032897) q[2];
sx q[2];
rz(-1.5264282) q[2];
sx q[2];
rz(-2.1849968) q[2];
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
