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
rz(-2.9057246) q[0];
sx q[0];
rz(-1.9960825) q[0];
sx q[0];
rz(3.0930162) q[0];
rz(1.5161169) q[1];
sx q[1];
rz(-2.0487509) q[1];
sx q[1];
rz(-0.71576524) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19312035) q[0];
sx q[0];
rz(-0.17331757) q[0];
sx q[0];
rz(1.3993708) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2907998) q[2];
sx q[2];
rz(-1.3716619) q[2];
sx q[2];
rz(-3.0906349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44385829) q[1];
sx q[1];
rz(-2.3045973) q[1];
sx q[1];
rz(-0.55056449) q[1];
x q[2];
rz(2.6244087) q[3];
sx q[3];
rz(-1.2150914) q[3];
sx q[3];
rz(0.16424792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48892659) q[2];
sx q[2];
rz(-2.581614) q[2];
sx q[2];
rz(-1.0555335) q[2];
rz(-2.7355898) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(0.82160151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15776289) q[0];
sx q[0];
rz(-2.1021748) q[0];
sx q[0];
rz(2.5854172) q[0];
rz(-1.8122541) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(3.0583196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52516925) q[0];
sx q[0];
rz(-0.94248191) q[0];
sx q[0];
rz(1.1675906) q[0];
rz(-pi) q[1];
rz(-0.82845344) q[2];
sx q[2];
rz(-0.76702416) q[2];
sx q[2];
rz(-1.3924862) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.992112) q[1];
sx q[1];
rz(-0.50051033) q[1];
sx q[1];
rz(2.6697199) q[1];
x q[2];
rz(-0.86822416) q[3];
sx q[3];
rz(-0.69528841) q[3];
sx q[3];
rz(-0.27588683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9437774) q[2];
sx q[2];
rz(-1.3894206) q[2];
sx q[2];
rz(2.6488292) q[2];
rz(-1.0265776) q[3];
sx q[3];
rz(-2.8387098) q[3];
sx q[3];
rz(-1.4125642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76348412) q[0];
sx q[0];
rz(-2.5788827) q[0];
sx q[0];
rz(-0.43877959) q[0];
rz(-1.6890866) q[1];
sx q[1];
rz(-0.32183281) q[1];
sx q[1];
rz(-0.39295331) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753112) q[0];
sx q[0];
rz(-1.5428758) q[0];
sx q[0];
rz(1.0795388) q[0];
rz(-pi) q[1];
rz(1.0424445) q[2];
sx q[2];
rz(-1.5990886) q[2];
sx q[2];
rz(1.9118146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6958017) q[1];
sx q[1];
rz(-1.6417827) q[1];
sx q[1];
rz(1.1522715) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6503104) q[3];
sx q[3];
rz(-2.2185555) q[3];
sx q[3];
rz(0.11707725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4517333) q[2];
sx q[2];
rz(-0.84443513) q[2];
sx q[2];
rz(2.8680657) q[2];
rz(0.61255974) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(2.2851022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3278811) q[0];
sx q[0];
rz(-2.5484945) q[0];
sx q[0];
rz(-0.80075049) q[0];
rz(-0.71209359) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(2.6584279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33813223) q[0];
sx q[0];
rz(-0.82809859) q[0];
sx q[0];
rz(-1.5716193) q[0];
x q[1];
rz(-2.6001873) q[2];
sx q[2];
rz(-2.0057937) q[2];
sx q[2];
rz(-2.5487399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.871328) q[1];
sx q[1];
rz(-2.21663) q[1];
sx q[1];
rz(-1.5346908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4523101) q[3];
sx q[3];
rz(-1.8346255) q[3];
sx q[3];
rz(-2.0949013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.157865) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(-1.857081) q[2];
rz(0.48786783) q[3];
sx q[3];
rz(-1.4372466) q[3];
sx q[3];
rz(-1.6871066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81216413) q[0];
sx q[0];
rz(-2.282113) q[0];
sx q[0];
rz(0.46464768) q[0];
rz(2.1931785) q[1];
sx q[1];
rz(-2.3034838) q[1];
sx q[1];
rz(1.6882187) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4522355) q[0];
sx q[0];
rz(-3.0984555) q[0];
sx q[0];
rz(2.1998911) q[0];
rz(-pi) q[1];
rz(-1.8273962) q[2];
sx q[2];
rz(-0.84689665) q[2];
sx q[2];
rz(0.43780299) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.075045295) q[1];
sx q[1];
rz(-0.61457026) q[1];
sx q[1];
rz(-0.11654186) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5713163) q[3];
sx q[3];
rz(-2.3304984) q[3];
sx q[3];
rz(-0.18660422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98443085) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(0.39443991) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30094576) q[0];
sx q[0];
rz(-0.41579682) q[0];
sx q[0];
rz(2.0796602) q[0];
rz(-2.1901219) q[1];
sx q[1];
rz(-2.2679592) q[1];
sx q[1];
rz(1.6208614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3979231) q[0];
sx q[0];
rz(-2.8594404) q[0];
sx q[0];
rz(1.762853) q[0];
x q[1];
rz(-2.4967147) q[2];
sx q[2];
rz(-2.0487924) q[2];
sx q[2];
rz(-1.0005815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54023149) q[1];
sx q[1];
rz(-1.733832) q[1];
sx q[1];
rz(-2.3705179) q[1];
x q[2];
rz(0.39395515) q[3];
sx q[3];
rz(-1.7749014) q[3];
sx q[3];
rz(-2.1597854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1157397) q[2];
sx q[2];
rz(-1.226475) q[2];
sx q[2];
rz(2.8958877) q[2];
rz(1.687441) q[3];
sx q[3];
rz(-1.629963) q[3];
sx q[3];
rz(2.3112442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65569735) q[0];
sx q[0];
rz(-2.4035154) q[0];
sx q[0];
rz(-0.6947211) q[0];
rz(-0.10063902) q[1];
sx q[1];
rz(-2.1058319) q[1];
sx q[1];
rz(2.2763841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14510205) q[0];
sx q[0];
rz(-1.1528413) q[0];
sx q[0];
rz(1.5254024) q[0];
rz(-0.018888028) q[2];
sx q[2];
rz(-2.2311418) q[2];
sx q[2];
rz(0.65062269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45564902) q[1];
sx q[1];
rz(-0.63243094) q[1];
sx q[1];
rz(-1.8258105) q[1];
rz(-pi) q[2];
rz(-0.65552222) q[3];
sx q[3];
rz(-2.351774) q[3];
sx q[3];
rz(1.5721377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18409099) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(0.91020477) q[2];
rz(1.4522067) q[3];
sx q[3];
rz(-1.2319177) q[3];
sx q[3];
rz(-2.6773101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24377395) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(-2.5183103) q[0];
rz(-0.15577623) q[1];
sx q[1];
rz(-0.77729762) q[1];
sx q[1];
rz(0.26330858) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52180852) q[0];
sx q[0];
rz(-1.5601047) q[0];
sx q[0];
rz(1.7858206) q[0];
rz(-pi) q[1];
rz(-3.0650862) q[2];
sx q[2];
rz(-1.6311809) q[2];
sx q[2];
rz(-3.1186207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0075078) q[1];
sx q[1];
rz(-0.99564394) q[1];
sx q[1];
rz(-1.6624727) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78531475) q[3];
sx q[3];
rz(-1.1668596) q[3];
sx q[3];
rz(-0.31111003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2498563) q[2];
sx q[2];
rz(-0.24954924) q[2];
sx q[2];
rz(-0.88805324) q[2];
rz(-0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(-2.1685062) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230098) q[0];
sx q[0];
rz(-0.50964481) q[0];
sx q[0];
rz(-0.50073671) q[0];
rz(1.1205193) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(-0.80148554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5699751) q[0];
sx q[0];
rz(-1.5727676) q[0];
sx q[0];
rz(0.070974102) q[0];
rz(2.0960484) q[2];
sx q[2];
rz(-2.0727951) q[2];
sx q[2];
rz(-1.0848622) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3858955) q[1];
sx q[1];
rz(-2.3508126) q[1];
sx q[1];
rz(-2.5189931) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3128223) q[3];
sx q[3];
rz(-1.7593062) q[3];
sx q[3];
rz(-2.9303355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72006172) q[2];
sx q[2];
rz(-1.820182) q[2];
sx q[2];
rz(-2.9113972) q[2];
rz(0.0097533334) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(0.15596381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941876) q[0];
sx q[0];
rz(-1.7740086) q[0];
sx q[0];
rz(0.76960027) q[0];
rz(-0.96254483) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(2.1544971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1302764) q[0];
sx q[0];
rz(-0.53142953) q[0];
sx q[0];
rz(-1.7087144) q[0];
rz(-pi) q[1];
rz(-1.8136421) q[2];
sx q[2];
rz(-1.7195903) q[2];
sx q[2];
rz(-0.75243581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7915114) q[1];
sx q[1];
rz(-1.529585) q[1];
sx q[1];
rz(0.033153127) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0760097) q[3];
sx q[3];
rz(-0.85799137) q[3];
sx q[3];
rz(-3.0312169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4276245) q[2];
sx q[2];
rz(-0.93728137) q[2];
sx q[2];
rz(0.06812185) q[2];
rz(2.6918329) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(1.4382039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.149067) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(-2.8515011) q[1];
sx q[1];
rz(-0.42872226) q[1];
sx q[1];
rz(0.52679481) q[1];
rz(2.6424825) q[2];
sx q[2];
rz(-1.2600095) q[2];
sx q[2];
rz(2.5799203) q[2];
rz(-0.3984821) q[3];
sx q[3];
rz(-1.553411) q[3];
sx q[3];
rz(-0.51391272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
