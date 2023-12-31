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
rz(-2.1142168) q[0];
sx q[0];
rz(2.7789814) q[0];
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1385961) q[0];
sx q[0];
rz(-1.9435104) q[0];
sx q[0];
rz(1.7025823) q[0];
x q[1];
rz(-0.50589675) q[2];
sx q[2];
rz(-2.4873173) q[2];
sx q[2];
rz(1.8088532) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29400723) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(2.1104382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9699691) q[3];
sx q[3];
rz(-1.3639543) q[3];
sx q[3];
rz(-2.9575461) q[3];
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
rz(-0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(1.6199002) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.99615) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94978588) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(1.1508862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1109353) q[2];
sx q[2];
rz(-2.1150555) q[2];
sx q[2];
rz(-1.1049251) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4821856) q[1];
sx q[1];
rz(-0.30989753) q[1];
sx q[1];
rz(-1.6531497) q[1];
x q[2];
rz(-2.8511091) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(0.089240616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48250616) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(2.144311) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230351) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-0.9054786) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3072976) q[0];
sx q[0];
rz(-2.8320438) q[0];
sx q[0];
rz(0.66972591) q[0];
rz(-2.3891874) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(0.67093713) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3382197) q[1];
sx q[1];
rz(-1.6446911) q[1];
sx q[1];
rz(2.282663) q[1];
rz(-pi) q[2];
rz(0.85729349) q[3];
sx q[3];
rz(-1.5382574) q[3];
sx q[3];
rz(1.9705704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3811615) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(3.0958214) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(0.12606829) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.9741845) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6895034) q[0];
sx q[0];
rz(-2.3514682) q[0];
sx q[0];
rz(-2.8566314) q[0];
rz(-1.3912958) q[2];
sx q[2];
rz(-2.0677535) q[2];
sx q[2];
rz(-2.3617982) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0988401) q[1];
sx q[1];
rz(-2.4272857) q[1];
sx q[1];
rz(-2.7318098) q[1];
rz(-pi) q[2];
rz(-0.91058369) q[3];
sx q[3];
rz(-2.108421) q[3];
sx q[3];
rz(1.1066574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(-1.1133105) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-2.6745093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665849) q[0];
sx q[0];
rz(-1.1995312) q[0];
sx q[0];
rz(0.21455547) q[0];
rz(0.37093016) q[2];
sx q[2];
rz(-2.5226454) q[2];
sx q[2];
rz(-1.3528454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52006066) q[1];
sx q[1];
rz(-2.660775) q[1];
sx q[1];
rz(2.0562999) q[1];
x q[2];
rz(0.88539447) q[3];
sx q[3];
rz(-1.4476895) q[3];
sx q[3];
rz(2.3791594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(2.6494027) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(-1.7135235) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(-1.496398) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8498358) q[0];
sx q[0];
rz(-2.09016) q[0];
sx q[0];
rz(0.22626466) q[0];
x q[1];
rz(1.9611347) q[2];
sx q[2];
rz(-2.5679553) q[2];
sx q[2];
rz(2.8611285) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8346356) q[1];
sx q[1];
rz(-2.7863057) q[1];
sx q[1];
rz(1.5026232) q[1];
rz(-pi) q[2];
rz(1.5835285) q[3];
sx q[3];
rz(-2.2836962) q[3];
sx q[3];
rz(-1.7604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(1.0923882) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(-0.011627442) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25225885) q[0];
sx q[0];
rz(-2.6703435) q[0];
sx q[0];
rz(2.528119) q[0];
rz(1.5989223) q[2];
sx q[2];
rz(-1.5453494) q[2];
sx q[2];
rz(2.5487119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9776476) q[1];
sx q[1];
rz(-1.5471336) q[1];
sx q[1];
rz(-1.0246828) q[1];
x q[2];
rz(-2.743268) q[3];
sx q[3];
rz(-0.71361226) q[3];
sx q[3];
rz(0.45809612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(-1.9912432) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(-2.7638226) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9858915) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(-0.62472384) q[0];
x q[1];
rz(-1.703007) q[2];
sx q[2];
rz(-1.6301179) q[2];
sx q[2];
rz(1.5166264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32840604) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(2.3996668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7842403) q[3];
sx q[3];
rz(-2.6297914) q[3];
sx q[3];
rz(1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(0.10929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(2.9633203) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423791) q[0];
sx q[0];
rz(-1.2782492) q[0];
sx q[0];
rz(-0.53036687) q[0];
rz(2.2237642) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(2.3248621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8930248) q[1];
sx q[1];
rz(-2.5556459) q[1];
sx q[1];
rz(0.47793169) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7618802) q[3];
sx q[3];
rz(-1.100193) q[3];
sx q[3];
rz(-0.1212561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69958413) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(1.3423086) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(-1.8348947) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(1.0888938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63554791) q[0];
sx q[0];
rz(-1.6354927) q[0];
sx q[0];
rz(-3.0788172) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8082952) q[2];
sx q[2];
rz(-2.4727159) q[2];
sx q[2];
rz(2.3046658) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1354462) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(-1.0165434) q[1];
x q[2];
rz(2.9519134) q[3];
sx q[3];
rz(-1.2343725) q[3];
sx q[3];
rz(1.3996901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(-0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(1.5169253) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(-0.4060612) q[2];
sx q[2];
rz(-0.95352298) q[2];
sx q[2];
rz(2.7609115) q[2];
rz(1.9061448) q[3];
sx q[3];
rz(-1.609758) q[3];
sx q[3];
rz(2.5555425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
