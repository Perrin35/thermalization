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
rz(0.95712823) q[0];
sx q[0];
rz(-1.4883833) q[0];
sx q[0];
rz(1.1487577) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(1.6869071) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4948954) q[0];
sx q[0];
rz(-2.2410164) q[0];
sx q[0];
rz(-0.47166069) q[0];
rz(-pi) q[1];
rz(-0.11496719) q[2];
sx q[2];
rz(-0.76412725) q[2];
sx q[2];
rz(0.72944356) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0599938) q[1];
sx q[1];
rz(-0.76548701) q[1];
sx q[1];
rz(-0.19123921) q[1];
rz(-0.1249073) q[3];
sx q[3];
rz(-2.4797399) q[3];
sx q[3];
rz(0.13368363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.181432) q[2];
sx q[2];
rz(-2.1655607) q[2];
sx q[2];
rz(1.206548) q[2];
rz(1.9801961) q[3];
sx q[3];
rz(-0.56777081) q[3];
sx q[3];
rz(1.5709706) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777765) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(0.97208446) q[0];
rz(0.26845911) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(2.6867902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1546611) q[0];
sx q[0];
rz(-1.6973746) q[0];
sx q[0];
rz(-1.2722227) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9399547) q[2];
sx q[2];
rz(-0.41339407) q[2];
sx q[2];
rz(1.8231376) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.090524448) q[1];
sx q[1];
rz(-1.6220434) q[1];
sx q[1];
rz(-0.4649051) q[1];
rz(-pi) q[2];
x q[2];
rz(2.11449) q[3];
sx q[3];
rz(-2.0683401) q[3];
sx q[3];
rz(3.080201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1099757) q[2];
sx q[2];
rz(-0.72828186) q[2];
sx q[2];
rz(-2.1419683) q[2];
rz(-3.055618) q[3];
sx q[3];
rz(-1.5558259) q[3];
sx q[3];
rz(-0.011367817) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40083945) q[0];
sx q[0];
rz(-1.6735621) q[0];
sx q[0];
rz(2.9174347) q[0];
rz(1.5466746) q[1];
sx q[1];
rz(-1.5432576) q[1];
sx q[1];
rz(1.8348144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9973361) q[0];
sx q[0];
rz(-1.4486509) q[0];
sx q[0];
rz(-2.5886378) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3575338) q[2];
sx q[2];
rz(-0.90558115) q[2];
sx q[2];
rz(-0.61675081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3114709) q[1];
sx q[1];
rz(-1.1563627) q[1];
sx q[1];
rz(0.45763514) q[1];
rz(0.31587028) q[3];
sx q[3];
rz(-2.4474553) q[3];
sx q[3];
rz(2.4296076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.811502) q[2];
sx q[2];
rz(-1.5218488) q[2];
sx q[2];
rz(2.1211993) q[2];
rz(1.329782) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(1.5248732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.735585) q[0];
sx q[0];
rz(-1.0574874) q[0];
sx q[0];
rz(-1.1757346) q[0];
rz(2.8658087) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(0.63794678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36786554) q[0];
sx q[0];
rz(-1.5380895) q[0];
sx q[0];
rz(-1.5601394) q[0];
rz(0.3023382) q[2];
sx q[2];
rz(-2.4983642) q[2];
sx q[2];
rz(-1.1272205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0093858) q[1];
sx q[1];
rz(-1.8957545) q[1];
sx q[1];
rz(-0.64909776) q[1];
rz(2.8218587) q[3];
sx q[3];
rz(-1.481062) q[3];
sx q[3];
rz(0.054320099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63567579) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(-1.4471311) q[2];
rz(0.21008374) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(2.213018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90784812) q[0];
sx q[0];
rz(-0.15234983) q[0];
sx q[0];
rz(2.0727169) q[0];
rz(1.2999889) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(1.2058421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17888363) q[0];
sx q[0];
rz(-0.41085748) q[0];
sx q[0];
rz(-1.370048) q[0];
rz(2.2236125) q[2];
sx q[2];
rz(-1.5783572) q[2];
sx q[2];
rz(-0.078469097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1002378) q[1];
sx q[1];
rz(-2.2861028) q[1];
sx q[1];
rz(-3.0042786) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30773325) q[3];
sx q[3];
rz(-2.640297) q[3];
sx q[3];
rz(-3.1395819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24870366) q[2];
sx q[2];
rz(-0.37962571) q[2];
sx q[2];
rz(-0.068566337) q[2];
rz(0.93962234) q[3];
sx q[3];
rz(-1.4087804) q[3];
sx q[3];
rz(1.2287963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3786316) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(2.6201541) q[0];
rz(1.5154845) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(0.62561402) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9245878) q[0];
sx q[0];
rz(-0.5408322) q[0];
sx q[0];
rz(-1.5103673) q[0];
rz(-pi) q[1];
rz(-0.31452532) q[2];
sx q[2];
rz(-2.7447766) q[2];
sx q[2];
rz(-2.0443969) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82078104) q[1];
sx q[1];
rz(-0.50034517) q[1];
sx q[1];
rz(-1.9920177) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42194052) q[3];
sx q[3];
rz(-1.5152389) q[3];
sx q[3];
rz(1.0797323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1336512) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(0.1055183) q[2];
rz(0.96406913) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(-2.156637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725726) q[0];
sx q[0];
rz(-2.9264086) q[0];
sx q[0];
rz(-1.7286812) q[0];
rz(-0.37881306) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(-2.5884195) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2755097) q[0];
sx q[0];
rz(-1.9635154) q[0];
sx q[0];
rz(-2.2662972) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.000053) q[2];
sx q[2];
rz(-1.4813378) q[2];
sx q[2];
rz(-2.3613561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87482086) q[1];
sx q[1];
rz(-1.8541341) q[1];
sx q[1];
rz(-1.5288407) q[1];
x q[2];
rz(-2.9384261) q[3];
sx q[3];
rz(-1.3152514) q[3];
sx q[3];
rz(1.1145962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(0.41213948) q[2];
rz(-0.66017094) q[3];
sx q[3];
rz(-2.6001402) q[3];
sx q[3];
rz(0.49303833) q[3];
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
rz(-2.2039345) q[0];
sx q[0];
rz(-2.0964607) q[0];
sx q[0];
rz(-0.20189051) q[0];
rz(-1.0578602) q[1];
sx q[1];
rz(-0.36483929) q[1];
sx q[1];
rz(0.26022628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8956586) q[0];
sx q[0];
rz(-1.9484451) q[0];
sx q[0];
rz(2.5833355) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43739983) q[2];
sx q[2];
rz(-1.0964616) q[2];
sx q[2];
rz(-2.2651644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3263798) q[1];
sx q[1];
rz(-1.2252439) q[1];
sx q[1];
rz(1.8680845) q[1];
x q[2];
rz(-1.4020355) q[3];
sx q[3];
rz(-0.91723727) q[3];
sx q[3];
rz(0.26168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95320931) q[2];
sx q[2];
rz(-1.4346432) q[2];
sx q[2];
rz(-0.58237135) q[2];
rz(-0.44238704) q[3];
sx q[3];
rz(-1.235639) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22170947) q[0];
sx q[0];
rz(-1.6136805) q[0];
sx q[0];
rz(1.7145994) q[0];
rz(1.7859979) q[1];
sx q[1];
rz(-1.4878863) q[1];
sx q[1];
rz(1.7048763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6528931) q[0];
sx q[0];
rz(-0.32036361) q[0];
sx q[0];
rz(-1.6344749) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6609909) q[2];
sx q[2];
rz(-1.3029104) q[2];
sx q[2];
rz(2.783423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7683803) q[1];
sx q[1];
rz(-1.5850968) q[1];
sx q[1];
rz(1.8119903) q[1];
x q[2];
rz(2.0974227) q[3];
sx q[3];
rz(-1.4344525) q[3];
sx q[3];
rz(2.6201893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9068678) q[2];
sx q[2];
rz(-1.2805254) q[2];
sx q[2];
rz(1.6024164) q[2];
rz(-2.5336044) q[3];
sx q[3];
rz(-1.7323078) q[3];
sx q[3];
rz(2.4517945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36668396) q[0];
sx q[0];
rz(-2.4803949) q[0];
sx q[0];
rz(-2.3181424) q[0];
rz(-0.89556328) q[1];
sx q[1];
rz(-1.3500373) q[1];
sx q[1];
rz(-0.38111883) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8117763) q[0];
sx q[0];
rz(-1.0346221) q[0];
sx q[0];
rz(-2.0945805) q[0];
rz(-pi) q[1];
rz(2.7398749) q[2];
sx q[2];
rz(-1.7251996) q[2];
sx q[2];
rz(0.32333514) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5504073) q[1];
sx q[1];
rz(-1.0246312) q[1];
sx q[1];
rz(0.3355432) q[1];
rz(-pi) q[2];
rz(-2.5140829) q[3];
sx q[3];
rz(-1.1063647) q[3];
sx q[3];
rz(3.0919558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33619189) q[2];
sx q[2];
rz(-2.7298584) q[2];
sx q[2];
rz(2.1545048) q[2];
rz(1.5385212) q[3];
sx q[3];
rz(-2.5542732) q[3];
sx q[3];
rz(-2.9015598) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22454746) q[0];
sx q[0];
rz(-1.9018835) q[0];
sx q[0];
rz(-1.6214669) q[0];
rz(0.70980258) q[1];
sx q[1];
rz(-2.5820422) q[1];
sx q[1];
rz(-1.6834264) q[1];
rz(1.456121) q[2];
sx q[2];
rz(-1.9330474) q[2];
sx q[2];
rz(-1.025552) q[2];
rz(2.1556606) q[3];
sx q[3];
rz(-2.2227047) q[3];
sx q[3];
rz(-0.97490132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
