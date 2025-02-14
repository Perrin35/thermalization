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
rz(0.23586805) q[0];
sx q[0];
rz(-1.1455102) q[0];
sx q[0];
rz(0.048576485) q[0];
rz(1.5161169) q[1];
sx q[1];
rz(4.2344344) q[1];
sx q[1];
rz(8.7090127) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5465821) q[0];
sx q[0];
rz(-1.5413741) q[0];
sx q[0];
rz(-1.3999697) q[0];
rz(-pi) q[1];
rz(1.2907998) q[2];
sx q[2];
rz(-1.3716619) q[2];
sx q[2];
rz(3.0906349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6977344) q[1];
sx q[1];
rz(-0.83699534) q[1];
sx q[1];
rz(-2.5910282) q[1];
rz(-pi) q[2];
rz(-0.51718398) q[3];
sx q[3];
rz(-1.2150914) q[3];
sx q[3];
rz(-2.9773447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6526661) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(-2.0860591) q[2];
rz(-2.7355898) q[3];
sx q[3];
rz(-1.8279653) q[3];
sx q[3];
rz(-0.82160151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838298) q[0];
sx q[0];
rz(-2.1021748) q[0];
sx q[0];
rz(2.5854172) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(-0.083273085) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10262415) q[0];
sx q[0];
rz(-2.4100465) q[0];
sx q[0];
rz(2.6464173) q[0];
rz(-pi) q[1];
rz(-0.57751285) q[2];
sx q[2];
rz(-1.0339875) q[2];
sx q[2];
rz(-2.2978738) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1412388) q[1];
sx q[1];
rz(-1.3508995) q[1];
sx q[1];
rz(2.6882437) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2733685) q[3];
sx q[3];
rz(-0.69528841) q[3];
sx q[3];
rz(2.8657058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1978153) q[2];
sx q[2];
rz(-1.3894206) q[2];
sx q[2];
rz(-0.49276349) q[2];
rz(2.1150151) q[3];
sx q[3];
rz(-0.30288282) q[3];
sx q[3];
rz(-1.7290285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76348412) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(0.43877959) q[0];
rz(-1.4525061) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(-0.39295331) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2344491) q[0];
sx q[0];
rz(-0.4919855) q[0];
sx q[0];
rz(-1.5116631) q[0];
rz(-pi) q[1];
rz(2.0991481) q[2];
sx q[2];
rz(-1.5425041) q[2];
sx q[2];
rz(-1.229778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90653268) q[1];
sx q[1];
rz(-1.1533914) q[1];
sx q[1];
rz(-0.077666186) q[1];
rz(-1.4912823) q[3];
sx q[3];
rz(-0.92303716) q[3];
sx q[3];
rz(0.11707725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4517333) q[2];
sx q[2];
rz(-2.2971575) q[2];
sx q[2];
rz(-0.27352697) q[2];
rz(-2.5290329) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(-0.85649049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3278811) q[0];
sx q[0];
rz(-0.59309816) q[0];
sx q[0];
rz(-2.3408422) q[0];
rz(-2.4294991) q[1];
sx q[1];
rz(-1.1199896) q[1];
sx q[1];
rz(2.6584279) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2321076) q[0];
sx q[0];
rz(-1.5714025) q[0];
sx q[0];
rz(0.74269791) q[0];
rz(-1.0739378) q[2];
sx q[2];
rz(-1.0845079) q[2];
sx q[2];
rz(2.4118347) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3222683) q[1];
sx q[1];
rz(-1.5996278) q[1];
sx q[1];
rz(-2.4954456) q[1];
x q[2];
rz(1.9075059) q[3];
sx q[3];
rz(-0.90970618) q[3];
sx q[3];
rz(-0.31240901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98372769) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(-1.857081) q[2];
rz(0.48786783) q[3];
sx q[3];
rz(-1.4372466) q[3];
sx q[3];
rz(1.454486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294285) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(0.46464768) q[0];
rz(-2.1931785) q[1];
sx q[1];
rz(-2.3034838) q[1];
sx q[1];
rz(1.4533739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4522355) q[0];
sx q[0];
rz(-3.0984555) q[0];
sx q[0];
rz(-2.1998911) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74044944) q[2];
sx q[2];
rz(-1.3794823) q[2];
sx q[2];
rz(-0.96093169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.067277597) q[1];
sx q[1];
rz(-0.96100531) q[1];
sx q[1];
rz(1.6526755) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0539758) q[3];
sx q[3];
rz(-2.2272416) q[3];
sx q[3];
rz(2.2051728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1571618) q[2];
sx q[2];
rz(-2.0941907) q[2];
sx q[2];
rz(2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406469) q[0];
sx q[0];
rz(-0.41579682) q[0];
sx q[0];
rz(-1.0619324) q[0];
rz(0.95147079) q[1];
sx q[1];
rz(-0.87363344) q[1];
sx q[1];
rz(-1.6208614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5976819) q[0];
sx q[0];
rz(-1.2939699) q[0];
sx q[0];
rz(3.0863161) q[0];
rz(-pi) q[1];
rz(-0.64487793) q[2];
sx q[2];
rz(-1.0928003) q[2];
sx q[2];
rz(2.1410112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1870059) q[1];
sx q[1];
rz(-2.329064) q[1];
sx q[1];
rz(-1.7962667) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6470408) q[3];
sx q[3];
rz(-0.44122094) q[3];
sx q[3];
rz(2.0989204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.025853) q[2];
sx q[2];
rz(-1.226475) q[2];
sx q[2];
rz(-0.24570492) q[2];
rz(-1.4541516) q[3];
sx q[3];
rz(-1.5116296) q[3];
sx q[3];
rz(0.83034849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4858953) q[0];
sx q[0];
rz(-0.73807722) q[0];
sx q[0];
rz(-2.4468716) q[0];
rz(-0.10063902) q[1];
sx q[1];
rz(-2.1058319) q[1];
sx q[1];
rz(-0.86520854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9964906) q[0];
sx q[0];
rz(-1.9887513) q[0];
sx q[0];
rz(1.5254024) q[0];
rz(-pi) q[1];
rz(-0.018888028) q[2];
sx q[2];
rz(-2.2311418) q[2];
sx q[2];
rz(0.65062269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6859436) q[1];
sx q[1];
rz(-0.63243094) q[1];
sx q[1];
rz(-1.8258105) q[1];
x q[2];
rz(2.1221625) q[3];
sx q[3];
rz(-2.1688281) q[3];
sx q[3];
rz(0.74265146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9575017) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(2.2313879) q[2];
rz(1.4522067) q[3];
sx q[3];
rz(-1.9096749) q[3];
sx q[3];
rz(2.6773101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24377395) q[0];
sx q[0];
rz(-1.047736) q[0];
sx q[0];
rz(0.62328231) q[0];
rz(0.15577623) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(0.26330858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1415188) q[0];
sx q[0];
rz(-0.21528582) q[0];
sx q[0];
rz(1.6208642) q[0];
rz(-pi) q[1];
rz(-1.510235) q[2];
sx q[2];
rz(-1.4944296) q[2];
sx q[2];
rz(1.5983943) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7548482) q[1];
sx q[1];
rz(-1.4939018) q[1];
sx q[1];
rz(-2.5645178) q[1];
rz(-pi) q[2];
rz(-1.0271163) q[3];
sx q[3];
rz(-0.86302084) q[3];
sx q[3];
rz(2.2563427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89173633) q[2];
sx q[2];
rz(-0.24954924) q[2];
sx q[2];
rz(2.2535394) q[2];
rz(0.92787162) q[3];
sx q[3];
rz(-2.0136621) q[3];
sx q[3];
rz(-2.1685062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230098) q[0];
sx q[0];
rz(-0.50964481) q[0];
sx q[0];
rz(2.6408559) q[0];
rz(-1.1205193) q[1];
sx q[1];
rz(-2.3605533) q[1];
sx q[1];
rz(2.3401071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028541738) q[0];
sx q[0];
rz(-0.071001425) q[0];
sx q[0];
rz(-0.027790471) q[0];
rz(-pi) q[1];
rz(0.5653462) q[2];
sx q[2];
rz(-2.0259116) q[2];
sx q[2];
rz(2.3836977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3858955) q[1];
sx q[1];
rz(-0.79078005) q[1];
sx q[1];
rz(-2.5189931) q[1];
x q[2];
rz(2.9467907) q[3];
sx q[3];
rz(-1.8240989) q[3];
sx q[3];
rz(1.4089438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72006172) q[2];
sx q[2];
rz(-1.3214107) q[2];
sx q[2];
rz(-0.23019543) q[2];
rz(3.1318393) q[3];
sx q[3];
rz(-1.516927) q[3];
sx q[3];
rz(-2.9856288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.6941876) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(2.3719924) q[0];
rz(-2.1790478) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(0.98709551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0113163) q[0];
sx q[0];
rz(-0.53142953) q[0];
sx q[0];
rz(-1.7087144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0133621) q[2];
sx q[2];
rz(-0.28403541) q[2];
sx q[2];
rz(1.7839873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35008123) q[1];
sx q[1];
rz(-1.6120076) q[1];
sx q[1];
rz(0.033153127) q[1];
rz(-pi) q[2];
rz(-0.77642595) q[3];
sx q[3];
rz(-1.9382432) q[3];
sx q[3];
rz(2.020379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7139682) q[2];
sx q[2];
rz(-2.2043113) q[2];
sx q[2];
rz(-0.06812185) q[2];
rz(2.6918329) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(-1.7033887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9925256) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(0.29009157) q[1];
sx q[1];
rz(-0.42872226) q[1];
sx q[1];
rz(0.52679481) q[1];
rz(-0.49911015) q[2];
sx q[2];
rz(-1.2600095) q[2];
sx q[2];
rz(2.5799203) q[2];
rz(-1.5896593) q[3];
sx q[3];
rz(-1.9692148) q[3];
sx q[3];
rz(-2.0773902) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
