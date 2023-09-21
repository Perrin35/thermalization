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
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(2.4584682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8389991) q[0];
sx q[0];
rz(-1.6380881) q[0];
sx q[0];
rz(-1.5124613) q[0];
rz(-pi) q[1];
rz(-0.97857742) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(2.4821971) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5176757) q[1];
sx q[1];
rz(-2.3530934) q[1];
sx q[1];
rz(2.4967525) q[1];
rz(-pi) q[2];
rz(1.0055553) q[3];
sx q[3];
rz(-1.2473277) q[3];
sx q[3];
rz(-2.2068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.6764486) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4629048) q[0];
sx q[0];
rz(-1.573274) q[0];
sx q[0];
rz(2.4173827) q[0];
rz(1.1128088) q[2];
sx q[2];
rz(-2.6070242) q[2];
sx q[2];
rz(2.2528258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.022545594) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(-0.70336282) q[1];
rz(-pi) q[2];
rz(1.2259237) q[3];
sx q[3];
rz(-2.2861835) q[3];
sx q[3];
rz(-1.8638924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(2.3410472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3886867) q[0];
sx q[0];
rz(-1.1645002) q[0];
sx q[0];
rz(2.9857062) q[0];
rz(-pi) q[1];
rz(0.82528798) q[2];
sx q[2];
rz(-0.85916677) q[2];
sx q[2];
rz(2.3724144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87677466) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(-0.012621565) q[1];
rz(-pi) q[2];
rz(-1.9242026) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.7819972) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-0.46491369) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(1.0850614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1151428) q[0];
sx q[0];
rz(-1.0466252) q[0];
sx q[0];
rz(1.8379704) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1539139) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(1.093986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76772056) q[1];
sx q[1];
rz(-1.4704736) q[1];
sx q[1];
rz(-0.14212455) q[1];
x q[2];
rz(0.61120175) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(1.1838278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(0.51182169) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22589332) q[0];
sx q[0];
rz(-0.75076538) q[0];
sx q[0];
rz(-0.70930003) q[0];
x q[1];
rz(1.328674) q[2];
sx q[2];
rz(-1.4171346) q[2];
sx q[2];
rz(-0.15649934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2770734) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(2.8237052) q[1];
rz(-pi) q[2];
rz(-2.1682407) q[3];
sx q[3];
rz(-1.6320328) q[3];
sx q[3];
rz(-0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-0.16528027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0342456) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(2.8427567) q[0];
x q[1];
rz(-0.084947305) q[2];
sx q[2];
rz(-0.68945486) q[2];
sx q[2];
rz(1.2598318) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.936441) q[1];
sx q[1];
rz(-1.9736104) q[1];
sx q[1];
rz(-0.58516296) q[1];
rz(-1.1214439) q[3];
sx q[3];
rz(-1.7835788) q[3];
sx q[3];
rz(1.7920115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(0.2120367) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(2.4051037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6155646) q[0];
sx q[0];
rz(-3.0897339) q[0];
sx q[0];
rz(-1.9066914) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8924106) q[2];
sx q[2];
rz(-1.2417925) q[2];
sx q[2];
rz(0.56275425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6474825) q[1];
sx q[1];
rz(-0.78299114) q[1];
sx q[1];
rz(-0.25581911) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.442939) q[3];
sx q[3];
rz(-2.5210288) q[3];
sx q[3];
rz(-2.8750949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-0.53156701) q[2];
rz(2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(-1.2980365) q[0];
rz(2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-0.92179006) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4066276) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(0.69358967) q[0];
rz(-pi) q[1];
rz(2.3085262) q[2];
sx q[2];
rz(-1.7852011) q[2];
sx q[2];
rz(-0.094878541) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31729749) q[1];
sx q[1];
rz(-1.824914) q[1];
sx q[1];
rz(-1.3586678) q[1];
rz(-pi) q[2];
rz(-2.2790518) q[3];
sx q[3];
rz(-1.7077912) q[3];
sx q[3];
rz(-0.18602895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8885324) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-0.67970651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25526991) q[0];
sx q[0];
rz(-2.5810044) q[0];
sx q[0];
rz(0.38829304) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7858511) q[2];
sx q[2];
rz(-2.4348223) q[2];
sx q[2];
rz(-0.099345318) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9437127) q[1];
sx q[1];
rz(-2.3730952) q[1];
sx q[1];
rz(2.6430623) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4623975) q[3];
sx q[3];
rz(-2.0742356) q[3];
sx q[3];
rz(-2.5356843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-2.2909686) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.95497447) q[1];
sx q[1];
rz(1.9649327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.048621) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(-3.0924348) q[0];
x q[1];
rz(0.93675905) q[2];
sx q[2];
rz(-2.0271795) q[2];
sx q[2];
rz(0.87181834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51917045) q[1];
sx q[1];
rz(-1.6455489) q[1];
sx q[1];
rz(1.012411) q[1];
x q[2];
rz(-0.51104607) q[3];
sx q[3];
rz(-1.4941477) q[3];
sx q[3];
rz(-2.6443036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8979793) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(1.127355) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-1.0424785) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(1.8489807) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(1.5789923) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
