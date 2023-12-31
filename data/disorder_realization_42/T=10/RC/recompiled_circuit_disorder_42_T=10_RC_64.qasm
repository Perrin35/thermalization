OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6289829) q[0];
sx q[0];
rz(4.1689685) q[0];
sx q[0];
rz(12.203759) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61603123) q[0];
sx q[0];
rz(-1.6934868) q[0];
sx q[0];
rz(-0.37567715) q[0];
x q[1];
rz(-1.2149493) q[2];
sx q[2];
rz(-2.1321745) q[2];
sx q[2];
rz(-1.9422308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8056148) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(-0.34259818) q[1];
rz(-pi) q[2];
rz(-2.0658675) q[3];
sx q[3];
rz(-2.6945811) q[3];
sx q[3];
rz(1.8398374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(2.1851052) q[2];
rz(0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(-1.1454426) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024591) q[0];
sx q[0];
rz(-0.67507889) q[0];
sx q[0];
rz(-2.549987) q[0];
rz(-0.59402324) q[2];
sx q[2];
rz(-1.1813287) q[2];
sx q[2];
rz(2.4246852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5686381) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(3.1152578) q[1];
x q[2];
rz(-2.8511091) q[3];
sx q[3];
rz(-2.4986914) q[3];
sx q[3];
rz(3.052352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(2.144311) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230351) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(-0.9054786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3072976) q[0];
sx q[0];
rz(-0.30954888) q[0];
sx q[0];
rz(2.4718667) q[0];
x q[1];
rz(-2.3891874) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(-2.4706555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14721522) q[1];
sx q[1];
rz(-2.4265687) q[1];
sx q[1];
rz(1.4579525) q[1];
rz(-pi) q[2];
x q[2];
rz(1.521103) q[3];
sx q[3];
rz(-0.71411413) q[3];
sx q[3];
rz(0.36220887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7604312) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(2.1543489) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-2.9614017) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-3.0155244) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(-1.9741845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8465189) q[0];
sx q[0];
rz(-2.3210038) q[0];
sx q[0];
rz(-1.2942765) q[0];
rz(-pi) q[1];
rz(-0.31801571) q[2];
sx q[2];
rz(-0.52581767) q[2];
sx q[2];
rz(-1.1434681) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0427525) q[1];
sx q[1];
rz(-2.4272857) q[1];
sx q[1];
rz(0.40978281) q[1];
x q[2];
rz(-0.79951841) q[3];
sx q[3];
rz(-0.82516731) q[3];
sx q[3];
rz(-2.094401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1642078) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(2.6389129) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-1.0031507) q[1];
sx q[1];
rz(-0.46708333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.607476) q[0];
sx q[0];
rz(-2.7152938) q[0];
sx q[0];
rz(-2.0712453) q[0];
x q[1];
rz(1.3181114) q[2];
sx q[2];
rz(-2.1420896) q[2];
sx q[2];
rz(-0.90734446) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0568697) q[1];
sx q[1];
rz(-1.1493756) q[1];
sx q[1];
rz(0.23878581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15848666) q[3];
sx q[3];
rz(-2.2500258) q[3];
sx q[3];
rz(2.2331626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46665329) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-0.49218991) q[2];
rz(0.28218937) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(-1.496398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85747) q[0];
sx q[0];
rz(-0.5623445) q[0];
sx q[0];
rz(-1.9447295) q[0];
rz(-1.1804579) q[2];
sx q[2];
rz(-0.57363735) q[2];
sx q[2];
rz(-2.8611285) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.306957) q[1];
sx q[1];
rz(-0.35528696) q[1];
sx q[1];
rz(-1.5026232) q[1];
x q[2];
rz(-3.1268678) q[3];
sx q[3];
rz(-2.428599) q[3];
sx q[3];
rz(1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4795586) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.5007639) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783766) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(1.0923882) q[0];
rz(1.7842402) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(-3.1299652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8788293) q[0];
sx q[0];
rz(-1.8352404) q[0];
sx q[0];
rz(-2.7468365) q[0];
rz(-pi) q[1];
rz(-1.5426703) q[2];
sx q[2];
rz(-1.5962432) q[2];
sx q[2];
rz(-2.5487119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36794084) q[1];
sx q[1];
rz(-2.5950187) q[1];
sx q[1];
rz(1.5252588) q[1];
rz(-pi) q[2];
x q[2];
rz(2.743268) q[3];
sx q[3];
rz(-2.4279804) q[3];
sx q[3];
rz(0.45809612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(0.28798506) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31036723) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.134682) q[0];
sx q[0];
rz(-1.3575166) q[0];
sx q[0];
rz(-2.8365305) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0817501) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(-0.062054141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89477506) q[1];
sx q[1];
rz(-0.90118876) q[1];
sx q[1];
rz(-2.0991904) q[1];
x q[2];
rz(1.3573523) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(-1.3444927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68547738) q[0];
sx q[0];
rz(-2.5427186) q[0];
sx q[0];
rz(-0.53703888) q[0];
x q[1];
rz(0.91782848) q[2];
sx q[2];
rz(-2.4036916) q[2];
sx q[2];
rz(2.3248621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8930248) q[1];
sx q[1];
rz(-0.5859468) q[1];
sx q[1];
rz(-2.663661) q[1];
x q[2];
rz(2.0719028) q[3];
sx q[3];
rz(-1.9074829) q[3];
sx q[3];
rz(1.5130373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-2.9163196) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63554791) q[0];
sx q[0];
rz(-1.5061) q[0];
sx q[0];
rz(-0.062775469) q[0];
rz(2.225869) q[2];
sx q[2];
rz(-1.7172126) q[2];
sx q[2];
rz(2.5953948) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8577733) q[1];
sx q[1];
rz(-0.76446629) q[1];
sx q[1];
rz(-2.43999) q[1];
x q[2];
rz(2.9519134) q[3];
sx q[3];
rz(-1.9072201) q[3];
sx q[3];
rz(1.7419025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-0.069996746) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.6488279) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(0.4060612) q[2];
sx q[2];
rz(-2.1880697) q[2];
sx q[2];
rz(-0.38068117) q[2];
rz(1.6886961) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
