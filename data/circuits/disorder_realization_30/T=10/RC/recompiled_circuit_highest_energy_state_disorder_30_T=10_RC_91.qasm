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
rz(1.0970595) q[0];
rz(1.2708083) q[1];
sx q[1];
rz(-1.1975809) q[1];
sx q[1];
rz(1.4374011) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88560933) q[0];
sx q[0];
rz(-0.75228158) q[0];
sx q[0];
rz(2.6153569) q[0];
rz(-pi) q[1];
rz(0.3008879) q[2];
sx q[2];
rz(-0.2425293) q[2];
sx q[2];
rz(0.96939847) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4258051) q[1];
sx q[1];
rz(-2.1658848) q[1];
sx q[1];
rz(2.7462105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5811722) q[3];
sx q[3];
rz(-2.2187067) q[3];
sx q[3];
rz(-1.5100175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2854332) q[2];
sx q[2];
rz(-2.9069052) q[2];
sx q[2];
rz(-2.6402546) q[2];
rz(1.3242138) q[3];
sx q[3];
rz(-0.54558498) q[3];
sx q[3];
rz(-1.6747624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.685574) q[0];
sx q[0];
rz(-0.40559232) q[0];
sx q[0];
rz(-1.6803886) q[0];
rz(-0.17924084) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(-2.5040212) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1792099) q[0];
sx q[0];
rz(-1.80733) q[0];
sx q[0];
rz(-1.80577) q[0];
rz(-1.887572) q[2];
sx q[2];
rz(-1.6143804) q[2];
sx q[2];
rz(-2.2769059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.04794807) q[1];
sx q[1];
rz(-0.94442299) q[1];
sx q[1];
rz(-0.44194024) q[1];
x q[2];
rz(-2.4781391) q[3];
sx q[3];
rz(-1.8725978) q[3];
sx q[3];
rz(-2.7581907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8403988) q[2];
sx q[2];
rz(-2.0638778) q[2];
sx q[2];
rz(-3.0974498) q[2];
rz(-2.5282777) q[3];
sx q[3];
rz(-1.8200579) q[3];
sx q[3];
rz(-0.80673748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089207) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(1.3432304) q[0];
rz(-1.7228458) q[1];
sx q[1];
rz(-0.90444618) q[1];
sx q[1];
rz(-0.80352965) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99074686) q[0];
sx q[0];
rz(-0.68239337) q[0];
sx q[0];
rz(3.0482024) q[0];
rz(-pi) q[1];
rz(2.0284925) q[2];
sx q[2];
rz(-2.105793) q[2];
sx q[2];
rz(2.0488103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5403998) q[1];
sx q[1];
rz(-1.6777039) q[1];
sx q[1];
rz(-1.8501297) q[1];
rz(-pi) q[2];
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
rz(-0.58328763) q[2];
sx q[2];
rz(-2.1896157) q[2];
sx q[2];
rz(1.0179016) q[2];
rz(-1.2244276) q[3];
sx q[3];
rz(-1.3320351) q[3];
sx q[3];
rz(-2.0126191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.99854904) q[0];
sx q[0];
rz(-0.7989378) q[0];
sx q[0];
rz(-2.1813188) q[0];
rz(2.9472561) q[1];
sx q[1];
rz(-1.7002218) q[1];
sx q[1];
rz(-0.0050553102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1881205) q[0];
sx q[0];
rz(-2.0998828) q[0];
sx q[0];
rz(0.75624589) q[0];
rz(2.9580367) q[2];
sx q[2];
rz(-1.608299) q[2];
sx q[2];
rz(1.3958193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21987629) q[1];
sx q[1];
rz(-1.3520693) q[1];
sx q[1];
rz(1.686386) q[1];
x q[2];
rz(-1.3157337) q[3];
sx q[3];
rz(-2.6082468) q[3];
sx q[3];
rz(1.7637788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2527577) q[2];
sx q[2];
rz(-1.5965261) q[2];
sx q[2];
rz(1.6952391) q[2];
rz(-1.3085922) q[3];
sx q[3];
rz(-1.7596217) q[3];
sx q[3];
rz(-0.60477177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62110353) q[0];
sx q[0];
rz(-0.72431722) q[0];
sx q[0];
rz(2.5526168) q[0];
rz(3.1336054) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(-1.084682) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046291489) q[0];
sx q[0];
rz(-1.3372375) q[0];
sx q[0];
rz(-2.8802242) q[0];
rz(-3.1118664) q[2];
sx q[2];
rz(-0.27966248) q[2];
sx q[2];
rz(3.0936808) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95653043) q[1];
sx q[1];
rz(-1.9635728) q[1];
sx q[1];
rz(-0.78861936) q[1];
x q[2];
rz(-1.6570377) q[3];
sx q[3];
rz(-0.6644333) q[3];
sx q[3];
rz(0.85103112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40953723) q[2];
sx q[2];
rz(-1.0169225) q[2];
sx q[2];
rz(-1.766905) q[2];
rz(-3.0459259) q[3];
sx q[3];
rz(-1.5547662) q[3];
sx q[3];
rz(0.63053757) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.0365527) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(-1.3404982) q[0];
rz(-2.4758677) q[1];
sx q[1];
rz(-1.3434854) q[1];
sx q[1];
rz(-2.417477) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042302) q[0];
sx q[0];
rz(-1.2688338) q[0];
sx q[0];
rz(-1.1911376) q[0];
rz(-pi) q[1];
rz(-0.55990368) q[2];
sx q[2];
rz(-1.0698058) q[2];
sx q[2];
rz(-2.5952114) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3427975) q[1];
sx q[1];
rz(-1.7067486) q[1];
sx q[1];
rz(2.6782381) q[1];
x q[2];
rz(-1.4744711) q[3];
sx q[3];
rz(-1.5755638) q[3];
sx q[3];
rz(-1.7624439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2035227) q[2];
sx q[2];
rz(-1.6911493) q[2];
sx q[2];
rz(3.0435437) q[2];
rz(2.8010662) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(0.19671973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.606474) q[0];
sx q[0];
rz(-0.70347324) q[0];
sx q[0];
rz(0.055543609) q[0];
rz(2.7659888) q[1];
sx q[1];
rz(-0.77508488) q[1];
sx q[1];
rz(1.4449878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21475131) q[0];
sx q[0];
rz(-0.80252778) q[0];
sx q[0];
rz(-1.7178128) q[0];
rz(-1.7710502) q[2];
sx q[2];
rz(-0.25875388) q[2];
sx q[2];
rz(-0.090752964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58126175) q[1];
sx q[1];
rz(-0.37788299) q[1];
sx q[1];
rz(0.037200971) q[1];
rz(-pi) q[2];
rz(0.7020601) q[3];
sx q[3];
rz(-2.5581048) q[3];
sx q[3];
rz(1.4017915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.043341788) q[2];
sx q[2];
rz(-1.3037553) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(-0.029953778) q[3];
sx q[3];
rz(-1.6311389) q[3];
sx q[3];
rz(-2.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.39716649) q[0];
sx q[0];
rz(-2.501896) q[0];
sx q[0];
rz(-1.8901012) q[0];
rz(-0.80998069) q[1];
sx q[1];
rz(-0.73695838) q[1];
sx q[1];
rz(1.4553778) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7969709) q[0];
sx q[0];
rz(-0.83528548) q[0];
sx q[0];
rz(-1.263144) q[0];
x q[1];
rz(-3.0461629) q[2];
sx q[2];
rz(-0.52431267) q[2];
sx q[2];
rz(-1.4316259) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4659459) q[1];
sx q[1];
rz(-0.94871229) q[1];
sx q[1];
rz(-2.3282566) q[1];
rz(2.2627566) q[3];
sx q[3];
rz(-2.2744826) q[3];
sx q[3];
rz(1.7233985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84888419) q[2];
sx q[2];
rz(-1.7851189) q[2];
sx q[2];
rz(-1.1513618) q[2];
rz(-0.014160841) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(-1.2667228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8779811) q[0];
sx q[0];
rz(-0.71313715) q[0];
sx q[0];
rz(0.94733316) q[0];
rz(0.5961279) q[1];
sx q[1];
rz(-1.4214186) q[1];
sx q[1];
rz(1.5790342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5681155) q[0];
sx q[0];
rz(-1.7422424) q[0];
sx q[0];
rz(2.417592) q[0];
x q[1];
rz(0.47548465) q[2];
sx q[2];
rz(-1.4275996) q[2];
sx q[2];
rz(-2.6729613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5925929) q[1];
sx q[1];
rz(-2.3301417) q[1];
sx q[1];
rz(0.92632067) q[1];
x q[2];
rz(-0.72266717) q[3];
sx q[3];
rz(-2.8407466) q[3];
sx q[3];
rz(0.85153714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8922213) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(-2.6577244) q[2];
rz(2.5441235) q[3];
sx q[3];
rz(-1.6941083) q[3];
sx q[3];
rz(-2.5634403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4018965) q[0];
sx q[0];
rz(-1.551349) q[0];
sx q[0];
rz(0.36648146) q[0];
rz(1.81555) q[1];
sx q[1];
rz(-0.96249023) q[1];
sx q[1];
rz(-2.0749626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13895282) q[0];
sx q[0];
rz(-2.3275597) q[0];
sx q[0];
rz(0.86281438) q[0];
x q[1];
rz(-3.0751247) q[2];
sx q[2];
rz(-2.53436) q[2];
sx q[2];
rz(0.83938423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.047160867) q[1];
sx q[1];
rz(-0.64241132) q[1];
sx q[1];
rz(1.1604527) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0935547) q[3];
sx q[3];
rz(-1.6001836) q[3];
sx q[3];
rz(-0.060197006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2610953) q[2];
sx q[2];
rz(-0.77793056) q[2];
sx q[2];
rz(1.0658537) q[2];
rz(2.8737658) q[3];
sx q[3];
rz(-1.6949751) q[3];
sx q[3];
rz(-2.6563787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5814701) q[0];
sx q[0];
rz(-1.0030092) q[0];
sx q[0];
rz(0.21448294) q[0];
rz(-0.32196925) q[1];
sx q[1];
rz(-1.8245158) q[1];
sx q[1];
rz(1.239924) q[1];
rz(0.027364597) q[2];
sx q[2];
rz(-0.28148808) q[2];
sx q[2];
rz(0.1053216) q[2];
rz(2.7291344) q[3];
sx q[3];
rz(-1.9720244) q[3];
sx q[3];
rz(0.26396863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
