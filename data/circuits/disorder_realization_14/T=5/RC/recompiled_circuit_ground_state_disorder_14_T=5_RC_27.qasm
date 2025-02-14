OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8811964) q[0];
sx q[0];
rz(-1.908778) q[0];
sx q[0];
rz(2.8137299) q[0];
rz(-0.39363632) q[1];
sx q[1];
rz(4.7596158) q[1];
sx q[1];
rz(14.16852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1122788) q[0];
sx q[0];
rz(-1.3697222) q[0];
sx q[0];
rz(1.5758833) q[0];
rz(-pi) q[1];
rz(0.27199409) q[2];
sx q[2];
rz(-1.0390721) q[2];
sx q[2];
rz(-0.069892757) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2214927) q[1];
sx q[1];
rz(-1.3650238) q[1];
sx q[1];
rz(-2.747072) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39253916) q[3];
sx q[3];
rz(-1.2219567) q[3];
sx q[3];
rz(2.4667645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9615122) q[2];
sx q[2];
rz(-2.1891258) q[2];
sx q[2];
rz(-2.3251779) q[2];
rz(-0.69711971) q[3];
sx q[3];
rz(-2.7643047) q[3];
sx q[3];
rz(-2.2009946) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546251) q[0];
sx q[0];
rz(-1.0579728) q[0];
sx q[0];
rz(-0.70607287) q[0];
rz(0.1708897) q[1];
sx q[1];
rz(-2.6306174) q[1];
sx q[1];
rz(-0.99542803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083849523) q[0];
sx q[0];
rz(-2.665976) q[0];
sx q[0];
rz(0.75288478) q[0];
rz(0.49071838) q[2];
sx q[2];
rz(-2.2123537) q[2];
sx q[2];
rz(0.12841719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66227205) q[1];
sx q[1];
rz(-2.8779128) q[1];
sx q[1];
rz(0.94698833) q[1];
rz(-pi) q[2];
rz(-0.63748116) q[3];
sx q[3];
rz(-1.5439171) q[3];
sx q[3];
rz(1.0589375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55676111) q[2];
sx q[2];
rz(-2.2274667) q[2];
sx q[2];
rz(-2.8625028) q[2];
rz(-0.68785214) q[3];
sx q[3];
rz(-0.22356859) q[3];
sx q[3];
rz(2.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85182178) q[0];
sx q[0];
rz(-0.91854799) q[0];
sx q[0];
rz(2.7626792) q[0];
rz(1.1907578) q[1];
sx q[1];
rz(-2.7749116) q[1];
sx q[1];
rz(-2.3653638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.533894) q[0];
sx q[0];
rz(-1.4149117) q[0];
sx q[0];
rz(-2.0034231) q[0];
rz(-pi) q[1];
rz(-1.4188781) q[2];
sx q[2];
rz(-1.3891274) q[2];
sx q[2];
rz(0.40975299) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4632193) q[1];
sx q[1];
rz(-1.3161945) q[1];
sx q[1];
rz(0.71779743) q[1];
rz(2.0950731) q[3];
sx q[3];
rz(-1.5308497) q[3];
sx q[3];
rz(1.0429045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1227405) q[2];
sx q[2];
rz(-1.2383702) q[2];
sx q[2];
rz(0.59993258) q[2];
rz(-1.2395202) q[3];
sx q[3];
rz(-1.7299088) q[3];
sx q[3];
rz(2.5304573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.83830738) q[0];
sx q[0];
rz(-2.4628283) q[0];
sx q[0];
rz(1.432206) q[0];
rz(-2.2749061) q[1];
sx q[1];
rz(-1.7984093) q[1];
sx q[1];
rz(1.0923045) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286273) q[0];
sx q[0];
rz(-1.019729) q[0];
sx q[0];
rz(-2.4222213) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42756568) q[2];
sx q[2];
rz(-2.1852772) q[2];
sx q[2];
rz(-0.63860047) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6553147) q[1];
sx q[1];
rz(-2.0351145) q[1];
sx q[1];
rz(-1.5686591) q[1];
x q[2];
rz(-0.74104167) q[3];
sx q[3];
rz(-2.1730086) q[3];
sx q[3];
rz(-2.9054537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0563353) q[2];
sx q[2];
rz(-1.9838355) q[2];
sx q[2];
rz(-1.866327) q[2];
rz(-1.3750252) q[3];
sx q[3];
rz(-1.9019889) q[3];
sx q[3];
rz(1.0331155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.105724) q[0];
sx q[0];
rz(-2.2585456) q[0];
sx q[0];
rz(-1.7778273) q[0];
rz(0.56963244) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(0.5035351) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0346876) q[0];
sx q[0];
rz(-1.4464753) q[0];
sx q[0];
rz(-2.9479545) q[0];
x q[1];
rz(-2.4120008) q[2];
sx q[2];
rz(-2.1520832) q[2];
sx q[2];
rz(0.94129291) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9787748) q[1];
sx q[1];
rz(-1.4895338) q[1];
sx q[1];
rz(0.012465076) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8721759) q[3];
sx q[3];
rz(-0.82357025) q[3];
sx q[3];
rz(0.034589442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2786431) q[2];
sx q[2];
rz(-0.38662616) q[2];
sx q[2];
rz(-1.6031727) q[2];
rz(0.1362416) q[3];
sx q[3];
rz(-1.4038266) q[3];
sx q[3];
rz(1.2374102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513026) q[0];
sx q[0];
rz(-0.61838111) q[0];
sx q[0];
rz(-2.8001617) q[0];
rz(0.731172) q[1];
sx q[1];
rz(-2.5656504) q[1];
sx q[1];
rz(-2.0700571) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8971469) q[0];
sx q[0];
rz(-1.9715152) q[0];
sx q[0];
rz(0.44735502) q[0];
x q[1];
rz(-0.086929079) q[2];
sx q[2];
rz(-0.61682103) q[2];
sx q[2];
rz(1.4035743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.711763) q[1];
sx q[1];
rz(-2.2380658) q[1];
sx q[1];
rz(2.1416452) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5626181) q[3];
sx q[3];
rz(-2.6371752) q[3];
sx q[3];
rz(2.2498508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7858481) q[2];
sx q[2];
rz(-1.637849) q[2];
sx q[2];
rz(0.49622932) q[2];
rz(-1.548454) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(2.1350071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9736495) q[0];
sx q[0];
rz(-1.4730467) q[0];
sx q[0];
rz(0.014852511) q[0];
rz(1.4133833) q[1];
sx q[1];
rz(-2.035886) q[1];
sx q[1];
rz(2.3766439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0891554) q[0];
sx q[0];
rz(-0.90298072) q[0];
sx q[0];
rz(-2.3994156) q[0];
rz(-0.68335376) q[2];
sx q[2];
rz(-1.8477447) q[2];
sx q[2];
rz(-2.5954557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1582163) q[1];
sx q[1];
rz(-1.9675273) q[1];
sx q[1];
rz(-1.5755754) q[1];
x q[2];
rz(1.9946003) q[3];
sx q[3];
rz(-1.7693345) q[3];
sx q[3];
rz(-2.1295223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0726274) q[2];
sx q[2];
rz(-1.5921389) q[2];
sx q[2];
rz(0.7507945) q[2];
rz(2.1454861) q[3];
sx q[3];
rz(-0.80544296) q[3];
sx q[3];
rz(-2.6944845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873416) q[0];
sx q[0];
rz(-0.60289201) q[0];
sx q[0];
rz(2.531429) q[0];
rz(0.24972406) q[1];
sx q[1];
rz(-1.4796673) q[1];
sx q[1];
rz(-3.0009559) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8630757) q[0];
sx q[0];
rz(-0.91916537) q[0];
sx q[0];
rz(2.5550575) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9570661) q[2];
sx q[2];
rz(-0.51287309) q[2];
sx q[2];
rz(-1.41092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9156981) q[1];
sx q[1];
rz(-1.4984301) q[1];
sx q[1];
rz(-0.27727978) q[1];
rz(0.76042995) q[3];
sx q[3];
rz(-1.3730341) q[3];
sx q[3];
rz(-2.9341079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3934624) q[2];
sx q[2];
rz(-0.40337864) q[2];
sx q[2];
rz(-1.0900137) q[2];
rz(-1.2202834) q[3];
sx q[3];
rz(-1.0459432) q[3];
sx q[3];
rz(0.88855851) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1423993) q[0];
sx q[0];
rz(-1.1141454) q[0];
sx q[0];
rz(2.2166369) q[0];
rz(-1.6389182) q[1];
sx q[1];
rz(-2.131856) q[1];
sx q[1];
rz(0.22770539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0696237) q[0];
sx q[0];
rz(-1.6771131) q[0];
sx q[0];
rz(-2.1031455) q[0];
rz(2.7289189) q[2];
sx q[2];
rz(-2.5209171) q[2];
sx q[2];
rz(0.35457573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5081577) q[1];
sx q[1];
rz(-2.7351008) q[1];
sx q[1];
rz(0.50199957) q[1];
rz(-2.5277522) q[3];
sx q[3];
rz(-1.4110397) q[3];
sx q[3];
rz(1.0626807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0324675) q[2];
sx q[2];
rz(-0.6296857) q[2];
sx q[2];
rz(-2.6510748) q[2];
rz(1.0444752) q[3];
sx q[3];
rz(-0.46382469) q[3];
sx q[3];
rz(0.0096376816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.57546416) q[0];
sx q[0];
rz(-0.45970356) q[0];
sx q[0];
rz(-1.8756961) q[0];
rz(1.5220386) q[1];
sx q[1];
rz(-1.8892037) q[1];
sx q[1];
rz(0.44580805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3193008) q[0];
sx q[0];
rz(-2.5133488) q[0];
sx q[0];
rz(-0.99530812) q[0];
rz(-pi) q[1];
rz(3.0580373) q[2];
sx q[2];
rz(-0.5131459) q[2];
sx q[2];
rz(-1.6864849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3054786) q[1];
sx q[1];
rz(-2.1784601) q[1];
sx q[1];
rz(-1.7377678) q[1];
rz(-pi) q[2];
rz(1.494094) q[3];
sx q[3];
rz(-2.8674403) q[3];
sx q[3];
rz(2.2488684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.012099115) q[2];
sx q[2];
rz(-2.1603277) q[2];
sx q[2];
rz(1.0274461) q[2];
rz(-1.3198613) q[3];
sx q[3];
rz(-2.8362995) q[3];
sx q[3];
rz(-2.6125438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8158067) q[0];
sx q[0];
rz(-1.3675714) q[0];
sx q[0];
rz(2.0969781) q[0];
rz(-0.927399) q[1];
sx q[1];
rz(-1.8721885) q[1];
sx q[1];
rz(2.4354557) q[1];
rz(2.8140284) q[2];
sx q[2];
rz(-1.8675928) q[2];
sx q[2];
rz(1.7525385) q[2];
rz(0.61946034) q[3];
sx q[3];
rz(-1.3347698) q[3];
sx q[3];
rz(2.6556849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
