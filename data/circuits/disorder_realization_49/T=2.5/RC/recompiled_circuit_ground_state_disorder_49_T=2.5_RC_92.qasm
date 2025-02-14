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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5985896) q[0];
sx q[0];
rz(-2.9393688) q[0];
sx q[0];
rz(0.36665066) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0027356) q[2];
sx q[2];
rz(-1.4502) q[2];
sx q[2];
rz(1.2029778) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4783096) q[1];
sx q[1];
rz(-1.4245207) q[1];
sx q[1];
rz(-1.0895132) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7571477) q[3];
sx q[3];
rz(-1.3923402) q[3];
sx q[3];
rz(2.3613191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9326707) q[2];
sx q[2];
rz(-2.966556) q[2];
sx q[2];
rz(-2.8026061) q[2];
rz(-2.637376) q[3];
sx q[3];
rz(-2.1239069) q[3];
sx q[3];
rz(-1.1241815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9826688) q[0];
sx q[0];
rz(-0.6944387) q[0];
sx q[0];
rz(2.6320631) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(0.94430077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12783192) q[0];
sx q[0];
rz(-1.7185129) q[0];
sx q[0];
rz(-1.3689408) q[0];
x q[1];
rz(1.239907) q[2];
sx q[2];
rz(-0.14741405) q[2];
sx q[2];
rz(0.40357631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92974508) q[1];
sx q[1];
rz(-2.3073688) q[1];
sx q[1];
rz(-0.75726189) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1491165) q[3];
sx q[3];
rz(-0.72315589) q[3];
sx q[3];
rz(-2.2167517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3462191) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-0.53595558) q[2];
rz(-1.0653982) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(0.65142256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.029243328) q[0];
sx q[0];
rz(-1.1493244) q[0];
sx q[0];
rz(2.405622) q[0];
rz(0.53572267) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(-0.89964286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0544708) q[0];
sx q[0];
rz(-2.961714) q[0];
sx q[0];
rz(0.69152559) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9016784) q[2];
sx q[2];
rz(-0.63096672) q[2];
sx q[2];
rz(-0.081693782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9316072) q[1];
sx q[1];
rz(-1.6495541) q[1];
sx q[1];
rz(1.0754536) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(0.54124093) q[2];
sx q[2];
rz(-1.9466126) q[2];
sx q[2];
rz(-2.3386653) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.4117389) q[0];
sx q[0];
rz(3.0112322) q[0];
rz(1.2654001) q[1];
sx q[1];
rz(-1.9644901) q[1];
sx q[1];
rz(-3.0317543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30194908) q[0];
sx q[0];
rz(-0.57900864) q[0];
sx q[0];
rz(-2.2173475) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5750138) q[2];
sx q[2];
rz(-1.3354805) q[2];
sx q[2];
rz(1.7641774) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.500587) q[1];
sx q[1];
rz(-0.88208032) q[1];
sx q[1];
rz(0.50488774) q[1];
x q[2];
rz(2.5205344) q[3];
sx q[3];
rz(-0.73535669) q[3];
sx q[3];
rz(-1.3842954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.639223) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(-0.42312527) q[2];
rz(-2.1028178) q[3];
sx q[3];
rz(-2.7602502) q[3];
sx q[3];
rz(-1.0264621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005242) q[0];
sx q[0];
rz(-1.2491913) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(-2.8084843) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(-0.28848973) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71995994) q[0];
sx q[0];
rz(-0.47556092) q[0];
sx q[0];
rz(-1.4325761) q[0];
x q[1];
rz(0.10941271) q[2];
sx q[2];
rz(-1.81443) q[2];
sx q[2];
rz(1.399216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34510558) q[1];
sx q[1];
rz(-3.0287841) q[1];
sx q[1];
rz(-2.2276001) q[1];
rz(-pi) q[2];
rz(0.94780677) q[3];
sx q[3];
rz(-1.3419328) q[3];
sx q[3];
rz(-0.2589489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5491817) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(-0.16331095) q[2];
rz(1.9773989) q[3];
sx q[3];
rz(-2.4653698) q[3];
sx q[3];
rz(1.3588847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.0613681) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(1.226271) q[0];
rz(-1.3310883) q[1];
sx q[1];
rz(-0.93003479) q[1];
sx q[1];
rz(0.61940449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50766599) q[0];
sx q[0];
rz(-2.0548425) q[0];
sx q[0];
rz(-0.54829161) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66483562) q[2];
sx q[2];
rz(-0.31739429) q[2];
sx q[2];
rz(1.9334671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.08949) q[1];
sx q[1];
rz(-1.5137496) q[1];
sx q[1];
rz(1.4228348) q[1];
rz(-pi) q[2];
rz(1.1167239) q[3];
sx q[3];
rz(-0.41560706) q[3];
sx q[3];
rz(2.6799283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2254534) q[2];
sx q[2];
rz(-1.4832387) q[2];
sx q[2];
rz(2.7062866) q[2];
rz(-3.0127323) q[3];
sx q[3];
rz(-1.83056) q[3];
sx q[3];
rz(2.784957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
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
rz(-2.5241959) q[1];
sx q[1];
rz(-1.2115024) q[1];
sx q[1];
rz(1.6389821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7107781) q[0];
sx q[0];
rz(-2.7511247) q[0];
sx q[0];
rz(1.6855408) q[0];
x q[1];
rz(2.9947979) q[2];
sx q[2];
rz(-2.0710315) q[2];
sx q[2];
rz(3.0584832) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7345846) q[1];
sx q[1];
rz(-0.53187766) q[1];
sx q[1];
rz(-1.8854669) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8822717) q[3];
sx q[3];
rz(-1.590456) q[3];
sx q[3];
rz(-0.033084083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5804533) q[2];
sx q[2];
rz(-0.20623198) q[2];
sx q[2];
rz(-1.0882161) q[2];
rz(-3.1214118) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(1.5599686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35767558) q[0];
sx q[0];
rz(-0.38991424) q[0];
sx q[0];
rz(2.1141323) q[0];
rz(1.1032392) q[1];
sx q[1];
rz(-1.5733066) q[1];
sx q[1];
rz(1.9267513) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64910674) q[0];
sx q[0];
rz(-1.6739573) q[0];
sx q[0];
rz(-1.2627559) q[0];
rz(-pi) q[1];
rz(-0.45502624) q[2];
sx q[2];
rz(-1.1385673) q[2];
sx q[2];
rz(2.2042008) q[2];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.72835975) q[3];
sx q[3];
rz(-1.9291147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2601629) q[2];
sx q[2];
rz(-1.9378928) q[2];
sx q[2];
rz(-2.0737958) q[2];
rz(2.2504375) q[3];
sx q[3];
rz(-0.46025899) q[3];
sx q[3];
rz(1.1394181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282309) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(2.886014) q[0];
rz(2.3981587) q[1];
sx q[1];
rz(-1.5579222) q[1];
sx q[1];
rz(-1.5240634) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7937318) q[0];
sx q[0];
rz(-1.1208236) q[0];
sx q[0];
rz(-1.0211358) q[0];
rz(-0.62047024) q[2];
sx q[2];
rz(-1.3201642) q[2];
sx q[2];
rz(2.5567101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2801622) q[1];
sx q[1];
rz(-1.4373684) q[1];
sx q[1];
rz(0.89806865) q[1];
rz(-0.26542191) q[3];
sx q[3];
rz(-0.94688225) q[3];
sx q[3];
rz(1.8056295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0301547) q[2];
sx q[2];
rz(-1.7058426) q[2];
sx q[2];
rz(1.5343522) q[2];
rz(2.0700908) q[3];
sx q[3];
rz(-1.5415618) q[3];
sx q[3];
rz(1.1736386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2509505) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(0.51837921) q[0];
rz(-2.3333343) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(1.6917276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5980412) q[0];
sx q[0];
rz(-1.8596453) q[0];
sx q[0];
rz(-3.1310215) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53624714) q[2];
sx q[2];
rz(-1.203158) q[2];
sx q[2];
rz(-2.4095031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5603578) q[1];
sx q[1];
rz(-0.61040813) q[1];
sx q[1];
rz(-0.44652744) q[1];
rz(-pi) q[2];
rz(-0.57761044) q[3];
sx q[3];
rz(-2.0311714) q[3];
sx q[3];
rz(-0.97629298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.2351278) q[2];
sx q[2];
rz(-2.2060564) q[2];
rz(-1.6701291) q[3];
sx q[3];
rz(-1.4751438) q[3];
sx q[3];
rz(1.0940301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.6937428) q[2];
sx q[2];
rz(-0.37005432) q[2];
sx q[2];
rz(2.4126929) q[2];
rz(-1.9121691) q[3];
sx q[3];
rz(-2.1044272) q[3];
sx q[3];
rz(2.9499346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
