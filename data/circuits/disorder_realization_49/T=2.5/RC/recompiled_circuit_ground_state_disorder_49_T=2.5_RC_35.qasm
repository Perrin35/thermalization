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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66798009) q[0];
sx q[0];
rz(-1.6428609) q[0];
sx q[0];
rz(-2.9524809) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14278463) q[2];
sx q[2];
rz(-1.0073642) q[2];
sx q[2];
rz(0.44445064) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1581067) q[1];
sx q[1];
rz(-2.0465104) q[1];
sx q[1];
rz(2.9768894) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38444491) q[3];
sx q[3];
rz(-1.3923402) q[3];
sx q[3];
rz(-2.3613191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2089219) q[2];
sx q[2];
rz(-2.966556) q[2];
sx q[2];
rz(-2.8026061) q[2];
rz(2.637376) q[3];
sx q[3];
rz(-1.0176858) q[3];
sx q[3];
rz(2.0174111) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15892383) q[0];
sx q[0];
rz(-0.6944387) q[0];
sx q[0];
rz(0.50952953) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(-2.1972919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4128542) q[0];
sx q[0];
rz(-1.7704238) q[0];
sx q[0];
rz(0.15073225) q[0];
rz(-pi) q[1];
rz(-1.9016857) q[2];
sx q[2];
rz(-2.9941786) q[2];
sx q[2];
rz(-0.40357631) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92974508) q[1];
sx q[1];
rz(-2.3073688) q[1];
sx q[1];
rz(-0.75726189) q[1];
rz(-0.89280309) q[3];
sx q[3];
rz(-1.8450741) q[3];
sx q[3];
rz(-2.1712042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3462191) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-2.6056371) q[2];
rz(1.0653982) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(2.4901701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029243328) q[0];
sx q[0];
rz(-1.1493244) q[0];
sx q[0];
rz(0.73597062) q[0];
rz(2.60587) q[1];
sx q[1];
rz(-1.842272) q[1];
sx q[1];
rz(-0.89964286) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0871219) q[0];
sx q[0];
rz(-2.961714) q[0];
sx q[0];
rz(2.4500671) q[0];
rz(1.9016784) q[2];
sx q[2];
rz(-0.63096672) q[2];
sx q[2];
rz(3.0598989) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8232628) q[1];
sx q[1];
rz(-2.0644651) q[1];
sx q[1];
rz(-0.089463316) q[1];
x q[2];
rz(-0.99863515) q[3];
sx q[3];
rz(-2.7021926) q[3];
sx q[3];
rz(-1.2422526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6003517) q[2];
sx q[2];
rz(-1.1949801) q[2];
sx q[2];
rz(-0.80292732) q[2];
rz(2.3927355) q[3];
sx q[3];
rz(-1.7436946) q[3];
sx q[3];
rz(0.15933855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.4117389) q[0];
sx q[0];
rz(-0.13036048) q[0];
rz(-1.8761926) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(-0.1098384) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3093053) q[0];
sx q[0];
rz(-1.9067295) q[0];
sx q[0];
rz(-2.0517339) q[0];
x q[1];
rz(0.56657882) q[2];
sx q[2];
rz(-1.3354805) q[2];
sx q[2];
rz(1.3774152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7848282) q[1];
sx q[1];
rz(-0.82883976) q[1];
sx q[1];
rz(1.0395365) q[1];
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
rz(-pi/2) q[1];
rz(1.639223) q[2];
sx q[2];
rz(-1.9618192) q[2];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1410685) q[0];
sx q[0];
rz(-1.2491913) q[0];
sx q[0];
rz(-0.27960676) q[0];
rz(0.33310834) q[1];
sx q[1];
rz(-0.50222841) q[1];
sx q[1];
rz(0.28848973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56474876) q[0];
sx q[0];
rz(-2.0414519) q[0];
sx q[0];
rz(3.070757) q[0];
rz(-0.10941271) q[2];
sx q[2];
rz(-1.3271626) q[2];
sx q[2];
rz(1.399216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5696249) q[1];
sx q[1];
rz(-1.5020084) q[1];
sx q[1];
rz(-1.4813165) q[1];
rz(-pi) q[2];
rz(1.9506504) q[3];
sx q[3];
rz(-0.65841802) q[3];
sx q[3];
rz(2.135621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5491817) q[2];
sx q[2];
rz(-1.9224527) q[2];
sx q[2];
rz(-0.16331095) q[2];
rz(-1.9773989) q[3];
sx q[3];
rz(-2.4653698) q[3];
sx q[3];
rz(1.7827079) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613681) q[0];
sx q[0];
rz(-0.017266406) q[0];
sx q[0];
rz(-1.226271) q[0];
rz(-1.8105043) q[1];
sx q[1];
rz(-0.93003479) q[1];
sx q[1];
rz(2.5221882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41202711) q[0];
sx q[0];
rz(-0.71456996) q[0];
sx q[0];
rz(0.78972915) q[0];
rz(-pi) q[1];
rz(-1.7707497) q[2];
sx q[2];
rz(-1.8189578) q[2];
sx q[2];
rz(-1.2436109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0521026) q[1];
sx q[1];
rz(-1.5137496) q[1];
sx q[1];
rz(-1.4228348) q[1];
x q[2];
rz(-1.9483637) q[3];
sx q[3];
rz(-1.7488297) q[3];
sx q[3];
rz(-1.5291027) q[3];
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
rz(-pi/2) q[1];
x q[3];
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
rz(-0.72047609) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(2.5893353) q[0];
rz(-0.61739677) q[1];
sx q[1];
rz(-1.9300902) q[1];
sx q[1];
rz(-1.5026106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5548067) q[0];
sx q[0];
rz(-1.1830336) q[0];
sx q[0];
rz(-3.094502) q[0];
x q[1];
rz(1.065997) q[2];
sx q[2];
rz(-1.6994972) q[2];
sx q[2];
rz(-1.4168878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4070081) q[1];
sx q[1];
rz(-2.609715) q[1];
sx q[1];
rz(1.2561258) q[1];
rz(1.5504567) q[3];
sx q[3];
rz(-1.830066) q[3];
sx q[3];
rz(1.6090956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5804533) q[2];
sx q[2];
rz(-0.20623198) q[2];
sx q[2];
rz(-2.0533766) q[2];
rz(-0.020180833) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(1.5816241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7839171) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(2.1141323) q[0];
rz(-1.1032392) q[1];
sx q[1];
rz(-1.5733066) q[1];
sx q[1];
rz(-1.9267513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2526557) q[0];
sx q[0];
rz(-1.8771457) q[0];
sx q[0];
rz(-3.0333748) q[0];
rz(-2.0452477) q[2];
sx q[2];
rz(-1.9812968) q[2];
sx q[2];
rz(-2.3060348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62544981) q[1];
sx q[1];
rz(-0.76458496) q[1];
sx q[1];
rz(2.0841875) q[1];
rz(2.5276466) q[3];
sx q[3];
rz(-1.1505732) q[3];
sx q[3];
rz(0.16638923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2601629) q[2];
sx q[2];
rz(-1.2036999) q[2];
sx q[2];
rz(2.0737958) q[2];
rz(2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(2.0021745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.7282309) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(-0.2555787) q[0];
rz(-0.74343395) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(-1.6175293) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3012863) q[0];
sx q[0];
rz(-0.69536007) q[0];
sx q[0];
rz(0.82454234) q[0];
rz(-pi) q[1];
rz(2.7268098) q[2];
sx q[2];
rz(-0.66294248) q[2];
sx q[2];
rz(0.65185968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12518203) q[1];
sx q[1];
rz(-2.4577854) q[1];
sx q[1];
rz(-1.7829624) q[1];
rz(1.9202523) q[3];
sx q[3];
rz(-0.67103681) q[3];
sx q[3];
rz(-0.90045917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1114379) q[2];
sx q[2];
rz(-1.7058426) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2509505) q[0];
sx q[0];
rz(-2.8548456) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(-2.3333343) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(1.6917276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302561) q[0];
sx q[0];
rz(-1.5606631) q[0];
sx q[0];
rz(-1.8596605) q[0];
rz(-pi) q[1];
rz(1.9920182) q[2];
sx q[2];
rz(-1.0738157) q[2];
sx q[2];
rz(-0.62825655) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58123484) q[1];
sx q[1];
rz(-2.5311845) q[1];
sx q[1];
rz(2.6950652) q[1];
x q[2];
rz(2.5639822) q[3];
sx q[3];
rz(-1.1104212) q[3];
sx q[3];
rz(0.97629298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.2351278) q[2];
sx q[2];
rz(0.9355363) q[2];
rz(-1.4714636) q[3];
sx q[3];
rz(-1.4751438) q[3];
sx q[3];
rz(2.0475625) q[3];
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
x q[0];
x q[1];
sx q[2];
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
