OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0177512) q[0];
sx q[0];
rz(-3.0525644) q[0];
sx q[0];
rz(-2.428399) q[0];
x q[1];
rz(0.97857742) q[2];
sx q[2];
rz(-2.7263612) q[2];
sx q[2];
rz(2.4821971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9477387) q[1];
sx q[1];
rz(-0.96809909) q[1];
sx q[1];
rz(-1.0268474) q[1];
rz(-2.7636823) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(-0.4370673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1258939) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(-1.4651441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629048) q[0];
sx q[0];
rz(-1.573274) q[0];
sx q[0];
rz(0.72420995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0287839) q[2];
sx q[2];
rz(-0.53456842) q[2];
sx q[2];
rz(-0.88876681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84491731) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(-2.1293473) q[1];
rz(-pi) q[2];
rz(-1.2259237) q[3];
sx q[3];
rz(-2.2861835) q[3];
sx q[3];
rz(-1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(2.3410472) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8974509) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(-1.98154) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86513743) q[2];
sx q[2];
rz(-2.1102998) q[2];
sx q[2];
rz(1.7973961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86299455) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(1.5763271) q[1];
x q[2];
rz(2.839746) q[3];
sx q[3];
rz(-2.2491124) q[3];
sx q[3];
rz(-0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(1.7819972) q[2];
rz(-0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-1.0850614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0264498) q[0];
sx q[0];
rz(-2.0949674) q[0];
sx q[0];
rz(-1.3036222) q[0];
rz(1.9374574) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(-2.1249078) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3738721) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(2.9994681) q[1];
rz(0.55235858) q[3];
sx q[3];
rz(-2.4529152) q[3];
sx q[3];
rz(2.3104582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(0.98168215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9156993) q[0];
sx q[0];
rz(-0.75076538) q[0];
sx q[0];
rz(2.4322926) q[0];
rz(-pi) q[1];
rz(-1.328674) q[2];
sx q[2];
rz(-1.4171346) q[2];
sx q[2];
rz(-2.9850933) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2213858) q[1];
sx q[1];
rz(-2.033794) q[1];
sx q[1];
rz(1.7358001) q[1];
x q[2];
rz(-0.97335191) q[3];
sx q[3];
rz(-1.6320328) q[3];
sx q[3];
rz(-2.6254568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(2.7691675) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-2.9763124) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0342456) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(0.29883595) q[0];
rz(-pi) q[1];
rz(-0.68768244) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(0.24535594) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0300671) q[1];
sx q[1];
rz(-2.1038052) q[1];
sx q[1];
rz(-2.0433321) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1214439) q[3];
sx q[3];
rz(-1.3580139) q[3];
sx q[3];
rz(1.7920115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(0.2120367) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(-2.4051037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6155646) q[0];
sx q[0];
rz(-3.0897339) q[0];
sx q[0];
rz(-1.9066914) q[0];
x q[1];
rz(0.94524224) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(-1.9117102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6474825) q[1];
sx q[1];
rz(-2.3586015) q[1];
sx q[1];
rz(-0.25581911) q[1];
rz(-0.69865366) q[3];
sx q[3];
rz(-2.5210288) q[3];
sx q[3];
rz(-0.26649775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-0.53156701) q[2];
rz(-2.452204) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078995973) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-2.2198026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(0.69358967) q[0];
rz(-pi) q[1];
rz(-0.28618257) q[2];
sx q[2];
rz(-0.85368644) q[2];
sx q[2];
rz(-1.66695) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8340048) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(-0.25968857) q[1];
rz(-1.7796302) q[3];
sx q[3];
rz(-2.4224671) q[3];
sx q[3];
rz(-1.5987087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-0.67970651) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70521351) q[0];
sx q[0];
rz(-1.0562911) q[0];
sx q[0];
rz(1.3374469) q[0];
rz(-0.87551261) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(1.8347486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19787994) q[1];
sx q[1];
rz(-0.76849741) q[1];
sx q[1];
rz(-2.6430623) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71984843) q[3];
sx q[3];
rz(-2.3206629) q[3];
sx q[3];
rz(1.5035226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(-2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.1766599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3424073) q[0];
sx q[0];
rz(-0.17051324) q[0];
sx q[0];
rz(-1.2810345) q[0];
rz(-pi) q[1];
rz(-2.2628225) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(0.15904418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9709819) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(1.7112205) q[1];
x q[2];
rz(2.9858399) q[3];
sx q[3];
rz(-0.51625801) q[3];
sx q[3];
rz(-1.2092276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(1.127355) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-2.7453616) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(0.30504967) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(3.1191961) q[3];
sx q[3];
rz(-2.7907484) q[3];
sx q[3];
rz(-0.69805935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];