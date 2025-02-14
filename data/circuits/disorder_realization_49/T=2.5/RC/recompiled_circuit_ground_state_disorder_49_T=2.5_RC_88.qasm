OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1385652) q[0];
sx q[0];
rz(-0.87831098) q[0];
sx q[0];
rz(-0.83100975) q[0];
rz(-0.57549685) q[1];
sx q[1];
rz(3.8776445) q[1];
sx q[1];
rz(10.164645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66798009) q[0];
sx q[0];
rz(-1.4987317) q[0];
sx q[0];
rz(0.18911171) q[0];
rz(-1.7923492) q[2];
sx q[2];
rz(-2.5622517) q[2];
sx q[2];
rz(2.960083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4783096) q[1];
sx q[1];
rz(-1.4245207) q[1];
sx q[1];
rz(-2.0520794) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44828592) q[3];
sx q[3];
rz(-0.42196754) q[3];
sx q[3];
rz(1.2039876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9326707) q[2];
sx q[2];
rz(-2.966556) q[2];
sx q[2];
rz(0.33898655) q[2];
rz(-2.637376) q[3];
sx q[3];
rz(-1.0176858) q[3];
sx q[3];
rz(-2.0174111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15892383) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(-2.6320631) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-2.8645611) q[1];
sx q[1];
rz(-0.94430077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4128542) q[0];
sx q[0];
rz(-1.3711689) q[0];
sx q[0];
rz(-2.9908604) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0933876) q[2];
sx q[2];
rz(-1.4314326) q[2];
sx q[2];
rz(-2.4037619) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9348) q[1];
sx q[1];
rz(-1.0369318) q[1];
sx q[1];
rz(2.4660048) q[1];
rz(-0.34667947) q[3];
sx q[3];
rz(-2.2190385) q[3];
sx q[3];
rz(2.7559506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3462191) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(2.6056371) q[2];
rz(1.0653982) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(2.4901701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(2.405622) q[0];
rz(-2.60587) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(2.2419498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38759432) q[0];
sx q[0];
rz(-1.4325465) q[0];
sx q[0];
rz(1.686245) q[0];
rz(-1.9016784) q[2];
sx q[2];
rz(-2.5106259) q[2];
sx q[2];
rz(-0.081693782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31832987) q[1];
sx q[1];
rz(-2.0644651) q[1];
sx q[1];
rz(3.0521293) q[1];
x q[2];
rz(2.1429575) q[3];
sx q[3];
rz(-0.43940001) q[3];
sx q[3];
rz(-1.8993401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6003517) q[2];
sx q[2];
rz(-1.1949801) q[2];
sx q[2];
rz(-2.3386653) q[2];
rz(0.7488572) q[3];
sx q[3];
rz(-1.7436946) q[3];
sx q[3];
rz(-0.15933855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.4117389) q[0];
sx q[0];
rz(-0.13036048) q[0];
rz(1.2654001) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(3.0317543) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8396436) q[0];
sx q[0];
rz(-2.562584) q[0];
sx q[0];
rz(0.92424519) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7214987) q[2];
sx q[2];
rz(-0.60852988) q[2];
sx q[2];
rz(0.15791721) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59200724) q[1];
sx q[1];
rz(-1.9535258) q[1];
sx q[1];
rz(2.3255583) q[1];
x q[2];
rz(2.5205344) q[3];
sx q[3];
rz(-2.406236) q[3];
sx q[3];
rz(-1.7572973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.639223) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(-0.42312527) q[2];
rz(-1.0387748) q[3];
sx q[3];
rz(-0.3813425) q[3];
sx q[3];
rz(-1.0264621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410685) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(-0.27960676) q[0];
rz(-0.33310834) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(-2.8531029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4216327) q[0];
sx q[0];
rz(-0.47556092) q[0];
sx q[0];
rz(-1.7090165) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9847068) q[2];
sx q[2];
rz(-2.8749646) q[2];
sx q[2];
rz(1.3150584) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34510558) q[1];
sx q[1];
rz(-0.11280858) q[1];
sx q[1];
rz(2.2276001) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1937859) q[3];
sx q[3];
rz(-1.7996598) q[3];
sx q[3];
rz(-0.2589489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.592411) q[2];
sx q[2];
rz(-1.9224527) q[2];
sx q[2];
rz(0.16331095) q[2];
rz(1.9773989) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(1.7827079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802245) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(-1.9153216) q[0];
rz(1.8105043) q[1];
sx q[1];
rz(-2.2115579) q[1];
sx q[1];
rz(-0.61940449) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7295655) q[0];
sx q[0];
rz(-2.4270227) q[0];
sx q[0];
rz(0.78972915) q[0];
rz(-pi) q[1];
rz(-2.8885968) q[2];
sx q[2];
rz(-1.7645451) q[2];
sx q[2];
rz(0.27744833) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2575469) q[1];
sx q[1];
rz(-0.15850286) q[1];
sx q[1];
rz(-1.201215) q[1];
rz(2.0248687) q[3];
sx q[3];
rz(-2.7259856) q[3];
sx q[3];
rz(-0.46166438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2254534) q[2];
sx q[2];
rz(-1.658354) q[2];
sx q[2];
rz(2.7062866) q[2];
rz(0.12886038) q[3];
sx q[3];
rz(-1.3110327) q[3];
sx q[3];
rz(-2.784957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72047609) q[0];
sx q[0];
rz(-2.5832376) q[0];
sx q[0];
rz(-0.55225736) q[0];
rz(2.5241959) q[1];
sx q[1];
rz(-1.9300902) q[1];
sx q[1];
rz(-1.5026106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5867859) q[0];
sx q[0];
rz(-1.9585591) q[0];
sx q[0];
rz(-0.047090637) q[0];
rz(-pi) q[1];
rz(-1.3093298) q[2];
sx q[2];
rz(-0.51957031) q[2];
sx q[2];
rz(-2.7594523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0459111) q[1];
sx q[1];
rz(-1.0675745) q[1];
sx q[1];
rz(-2.9614425) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25932094) q[3];
sx q[3];
rz(-1.5511366) q[3];
sx q[3];
rz(0.033084083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5804533) q[2];
sx q[2];
rz(-2.9353607) q[2];
sx q[2];
rz(-2.0533766) q[2];
rz(3.1214118) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(-1.5599686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7839171) q[0];
sx q[0];
rz(-0.38991424) q[0];
sx q[0];
rz(-2.1141323) q[0];
rz(1.1032392) q[1];
sx q[1];
rz(-1.5682861) q[1];
sx q[1];
rz(1.2148414) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9068844) q[0];
sx q[0];
rz(-2.8172593) q[0];
sx q[0];
rz(-1.2417481) q[0];
rz(2.0452477) q[2];
sx q[2];
rz(-1.1602959) q[2];
sx q[2];
rz(0.83555789) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5161428) q[1];
sx q[1];
rz(-2.3770077) q[1];
sx q[1];
rz(-2.0841875) q[1];
x q[2];
rz(-0.61394604) q[3];
sx q[3];
rz(-1.1505732) q[3];
sx q[3];
rz(0.16638923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2601629) q[2];
sx q[2];
rz(-1.9378928) q[2];
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
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282309) q[0];
sx q[0];
rz(-0.75465337) q[0];
sx q[0];
rz(-0.2555787) q[0];
rz(-0.74343395) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(1.5240634) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3012863) q[0];
sx q[0];
rz(-2.4462326) q[0];
sx q[0];
rz(-0.82454234) q[0];
rz(2.5211224) q[2];
sx q[2];
rz(-1.8214284) q[2];
sx q[2];
rz(0.58488256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12518203) q[1];
sx q[1];
rz(-2.4577854) q[1];
sx q[1];
rz(1.7829624) q[1];
rz(-1.2213403) q[3];
sx q[3];
rz(-2.4705558) q[3];
sx q[3];
rz(-2.2411335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0301547) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(-1.5343522) q[2];
rz(1.0715019) q[3];
sx q[3];
rz(-1.5415618) q[3];
sx q[3];
rz(1.967954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8906422) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(0.80825835) q[1];
sx q[1];
rz(-1.6812485) q[1];
sx q[1];
rz(1.4498651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0302561) q[0];
sx q[0];
rz(-1.5606631) q[0];
sx q[0];
rz(-1.2819321) q[0];
rz(-1.9920182) q[2];
sx q[2];
rz(-2.067777) q[2];
sx q[2];
rz(2.5133361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7781642) q[1];
sx q[1];
rz(-1.320667) q[1];
sx q[1];
rz(2.5787337) q[1];
x q[2];
rz(2.1052741) q[3];
sx q[3];
rz(-1.0596529) q[3];
sx q[3];
rz(-0.3126463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(2.2060564) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.6846631) q[0];
sx q[0];
rz(-2.5419432) q[0];
sx q[0];
rz(2.3210617) q[0];
rz(0.44878557) q[1];
sx q[1];
rz(-1.9955336) q[1];
sx q[1];
rz(1.21036) q[1];
rz(1.2032897) q[2];
sx q[2];
rz(-1.5264282) q[2];
sx q[2];
rz(-2.1849968) q[2];
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
