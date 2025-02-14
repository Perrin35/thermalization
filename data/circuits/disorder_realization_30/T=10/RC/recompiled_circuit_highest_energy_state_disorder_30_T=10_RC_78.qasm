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
rz(-2.2089145) q[0];
sx q[0];
rz(-0.39251602) q[0];
sx q[0];
rz(-2.0445332) q[0];
rz(-1.8707844) q[1];
sx q[1];
rz(4.3391736) q[1];
sx q[1];
rz(14.270562) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559833) q[0];
sx q[0];
rz(-2.3893111) q[0];
sx q[0];
rz(-0.52623574) q[0];
rz(-1.6439866) q[2];
sx q[2];
rz(-1.3393667) q[2];
sx q[2];
rz(0.65999962) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7858251) q[1];
sx q[1];
rz(-0.70097979) q[1];
sx q[1];
rz(2.088083) q[1];
x q[2];
rz(0.64793628) q[3];
sx q[3];
rz(-1.5625232) q[3];
sx q[3];
rz(-0.054516597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8561594) q[2];
sx q[2];
rz(-0.23468748) q[2];
sx q[2];
rz(-2.6402546) q[2];
rz(-1.3242138) q[3];
sx q[3];
rz(-2.5960077) q[3];
sx q[3];
rz(-1.6747624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45601869) q[0];
sx q[0];
rz(-2.7360003) q[0];
sx q[0];
rz(1.6803886) q[0];
rz(-0.17924084) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(-2.5040212) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9623827) q[0];
sx q[0];
rz(-1.80733) q[0];
sx q[0];
rz(-1.3358227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2540206) q[2];
sx q[2];
rz(-1.5272123) q[2];
sx q[2];
rz(2.2769059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3481787) q[1];
sx q[1];
rz(-1.2169588) q[1];
sx q[1];
rz(-0.89576141) q[1];
rz(-2.6735503) q[3];
sx q[3];
rz(-0.7193102) q[3];
sx q[3];
rz(-1.5908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3011938) q[2];
sx q[2];
rz(-1.0777148) q[2];
sx q[2];
rz(-3.0974498) q[2];
rz(-0.61331493) q[3];
sx q[3];
rz(-1.3215348) q[3];
sx q[3];
rz(2.3348552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089207) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(-1.3432304) q[0];
rz(-1.7228458) q[1];
sx q[1];
rz(-0.90444618) q[1];
sx q[1];
rz(2.338063) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4889824) q[0];
sx q[0];
rz(-1.6296415) q[0];
sx q[0];
rz(2.4613358) q[0];
x q[1];
rz(-2.5577689) q[2];
sx q[2];
rz(-1.1808504) q[2];
sx q[2];
rz(-2.4175298) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0002001) q[1];
sx q[1];
rz(-1.2931001) q[1];
sx q[1];
rz(0.11118366) q[1];
x q[2];
rz(-0.50602956) q[3];
sx q[3];
rz(-1.2057736) q[3];
sx q[3];
rz(2.6377691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.558305) q[2];
sx q[2];
rz(-2.1896157) q[2];
sx q[2];
rz(1.0179016) q[2];
rz(-1.2244276) q[3];
sx q[3];
rz(-1.8095576) q[3];
sx q[3];
rz(2.0126191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.1430436) q[0];
sx q[0];
rz(-2.3426549) q[0];
sx q[0];
rz(0.96027389) q[0];
rz(0.1943365) q[1];
sx q[1];
rz(-1.4413709) q[1];
sx q[1];
rz(3.1365373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086297) q[0];
sx q[0];
rz(-0.89187183) q[0];
sx q[0];
rz(2.4358948) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6089392) q[2];
sx q[2];
rz(-1.7542217) q[2];
sx q[2];
rz(-2.9596552) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27144602) q[1];
sx q[1];
rz(-0.24695858) q[1];
sx q[1];
rz(2.6629763) q[1];
rz(2.9937135) q[3];
sx q[3];
rz(-2.085146) q[3];
sx q[3];
rz(-1.0837931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2527577) q[2];
sx q[2];
rz(-1.5450666) q[2];
sx q[2];
rz(-1.6952391) q[2];
rz(-1.3085922) q[3];
sx q[3];
rz(-1.3819709) q[3];
sx q[3];
rz(-2.5368209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62110353) q[0];
sx q[0];
rz(-0.72431722) q[0];
sx q[0];
rz(-0.58897585) q[0];
rz(-0.0079872459) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(-1.084682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81122855) q[0];
sx q[0];
rz(-0.34872751) q[0];
sx q[0];
rz(0.74411157) q[0];
rz(-pi) q[1];
rz(0.029726278) q[2];
sx q[2];
rz(-0.27966248) q[2];
sx q[2];
rz(3.0936808) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8950348) q[1];
sx q[1];
rz(-2.2854574) q[1];
sx q[1];
rz(-1.0393875) q[1];
rz(-pi) q[2];
rz(1.6570377) q[3];
sx q[3];
rz(-2.4771593) q[3];
sx q[3];
rz(-2.2905615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7320554) q[2];
sx q[2];
rz(-2.1246702) q[2];
sx q[2];
rz(1.3746877) q[2];
rz(-3.0459259) q[3];
sx q[3];
rz(-1.5868264) q[3];
sx q[3];
rz(-0.63053757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365527) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(-1.8010944) q[0];
rz(-0.66572491) q[1];
sx q[1];
rz(-1.3434854) q[1];
sx q[1];
rz(2.417477) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63987565) q[0];
sx q[0];
rz(-2.6611009) q[0];
sx q[0];
rz(-2.2697422) q[0];
x q[1];
rz(0.99704717) q[2];
sx q[2];
rz(-2.0553751) q[2];
sx q[2];
rz(-0.73198971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.648702) q[1];
sx q[1];
rz(-0.48149432) q[1];
sx q[1];
rz(0.29701155) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0047896623) q[3];
sx q[3];
rz(-1.4744722) q[3];
sx q[3];
rz(2.9504058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2035227) q[2];
sx q[2];
rz(-1.6911493) q[2];
sx q[2];
rz(-0.098048992) q[2];
rz(-0.34052643) q[3];
sx q[3];
rz(-0.90172684) q[3];
sx q[3];
rz(-0.19671973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53511867) q[0];
sx q[0];
rz(-2.4381194) q[0];
sx q[0];
rz(3.086049) q[0];
rz(-2.7659888) q[1];
sx q[1];
rz(-0.77508488) q[1];
sx q[1];
rz(1.6966049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42471573) q[0];
sx q[0];
rz(-0.77936855) q[0];
sx q[0];
rz(-2.9911441) q[0];
rz(-1.3705424) q[2];
sx q[2];
rz(-2.8828388) q[2];
sx q[2];
rz(3.0508397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58126175) q[1];
sx q[1];
rz(-0.37788299) q[1];
sx q[1];
rz(3.1043917) q[1];
rz(-pi) q[2];
rz(1.1678004) q[3];
sx q[3];
rz(-1.1366362) q[3];
sx q[3];
rz(0.94769615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.043341788) q[2];
sx q[2];
rz(-1.3037553) q[2];
sx q[2];
rz(1.8275758) q[2];
rz(-3.1116389) q[3];
sx q[3];
rz(-1.5104537) q[3];
sx q[3];
rz(-2.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444262) q[0];
sx q[0];
rz(-0.63969669) q[0];
sx q[0];
rz(1.2514914) q[0];
rz(-0.80998069) q[1];
sx q[1];
rz(-2.4046343) q[1];
sx q[1];
rz(-1.4553778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43622323) q[0];
sx q[0];
rz(-1.3443265) q[0];
sx q[0];
rz(-2.3821077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6258442) q[2];
sx q[2];
rz(-1.0491129) q[2];
sx q[2];
rz(-1.59984) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6941713) q[1];
sx q[1];
rz(-2.2024367) q[1];
sx q[1];
rz(2.3775435) q[1];
rz(-pi) q[2];
rz(-2.496893) q[3];
sx q[3];
rz(-0.94333269) q[3];
sx q[3];
rz(2.6312466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84888419) q[2];
sx q[2];
rz(-1.7851189) q[2];
sx q[2];
rz(-1.9902309) q[2];
rz(3.1274318) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(1.8748698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8779811) q[0];
sx q[0];
rz(-0.71313715) q[0];
sx q[0];
rz(-2.1942595) q[0];
rz(-2.5454648) q[1];
sx q[1];
rz(-1.7201741) q[1];
sx q[1];
rz(1.5625585) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2939908) q[0];
sx q[0];
rz(-0.85968535) q[0];
sx q[0];
rz(-1.343665) q[0];
x q[1];
rz(-1.4100227) q[2];
sx q[2];
rz(-1.1005745) q[2];
sx q[2];
rz(-1.9660814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.642338) q[1];
sx q[1];
rz(-2.0216536) q[1];
sx q[1];
rz(0.87085215) q[1];
rz(-0.72266717) q[3];
sx q[3];
rz(-0.30084601) q[3];
sx q[3];
rz(-0.85153714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8922213) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(-0.48386827) q[2];
rz(-2.5441235) q[3];
sx q[3];
rz(-1.6941083) q[3];
sx q[3];
rz(-0.57815236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396962) q[0];
sx q[0];
rz(-1.5902436) q[0];
sx q[0];
rz(0.36648146) q[0];
rz(-1.81555) q[1];
sx q[1];
rz(-2.1791024) q[1];
sx q[1];
rz(-2.0749626) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410927) q[0];
sx q[0];
rz(-1.0783195) q[0];
sx q[0];
rz(-0.89333138) q[0];
rz(-pi) q[1];
rz(2.5353955) q[2];
sx q[2];
rz(-1.6087039) q[2];
sx q[2];
rz(0.67680046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8588094) q[1];
sx q[1];
rz(-1.8121392) q[1];
sx q[1];
rz(-0.96940746) q[1];
rz(-pi) q[2];
rz(-1.6347061) q[3];
sx q[3];
rz(-0.47807594) q[3];
sx q[3];
rz(1.4538387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2610953) q[2];
sx q[2];
rz(-2.3636621) q[2];
sx q[2];
rz(-1.0658537) q[2];
rz(2.8737658) q[3];
sx q[3];
rz(-1.4466176) q[3];
sx q[3];
rz(-0.48521391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601226) q[0];
sx q[0];
rz(-2.1385834) q[0];
sx q[0];
rz(-2.9271097) q[0];
rz(2.8196234) q[1];
sx q[1];
rz(-1.8245158) q[1];
sx q[1];
rz(1.239924) q[1];
rz(-2.8602045) q[2];
sx q[2];
rz(-1.5631957) q[2];
sx q[2];
rz(-1.4917628) q[2];
rz(-2.0044708) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(1.6655123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
