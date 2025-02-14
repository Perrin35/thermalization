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
rz(-1.8347972) q[0];
sx q[0];
rz(-2.2936294) q[0];
sx q[0];
rz(0.037394878) q[0];
rz(2.1283863) q[1];
sx q[1];
rz(-2.6599045) q[1];
sx q[1];
rz(-1.4451292) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8593711) q[0];
sx q[0];
rz(-1.9992775) q[0];
sx q[0];
rz(1.7707878) q[0];
x q[1];
rz(1.4679392) q[2];
sx q[2];
rz(-1.5005996) q[2];
sx q[2];
rz(1.4597033) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8822669) q[1];
sx q[1];
rz(-2.1082467) q[1];
sx q[1];
rz(0.56973704) q[1];
x q[2];
rz(1.7496893) q[3];
sx q[3];
rz(-2.5692392) q[3];
sx q[3];
rz(1.6088886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41363132) q[2];
sx q[2];
rz(-0.093158826) q[2];
sx q[2];
rz(1.6348582) q[2];
rz(0.19351752) q[3];
sx q[3];
rz(-0.7842803) q[3];
sx q[3];
rz(-0.94582742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042260878) q[0];
sx q[0];
rz(-1.3657382) q[0];
sx q[0];
rz(0.17901626) q[0];
rz(1.6049113) q[1];
sx q[1];
rz(-2.8004526) q[1];
sx q[1];
rz(1.590033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59373271) q[0];
sx q[0];
rz(-1.4734125) q[0];
sx q[0];
rz(-1.2800526) q[0];
x q[1];
rz(1.8297029) q[2];
sx q[2];
rz(-2.2635411) q[2];
sx q[2];
rz(1.9996114) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68374709) q[1];
sx q[1];
rz(-2.6553255) q[1];
sx q[1];
rz(-2.1621428) q[1];
x q[2];
rz(2.670774) q[3];
sx q[3];
rz(-2.0662466) q[3];
sx q[3];
rz(-1.9268056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67148525) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(3.1191077) q[2];
rz(2.2454028) q[3];
sx q[3];
rz(-2.360207) q[3];
sx q[3];
rz(-1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3317868) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(-1.5326387) q[0];
rz(1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(0.87361139) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667359) q[0];
sx q[0];
rz(-2.2169161) q[0];
sx q[0];
rz(-2.2227051) q[0];
rz(-pi) q[1];
x q[1];
rz(2.457946) q[2];
sx q[2];
rz(-2.3693759) q[2];
sx q[2];
rz(-1.5815841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.552678) q[1];
sx q[1];
rz(-0.50929087) q[1];
sx q[1];
rz(2.8224148) q[1];
rz(-pi) q[2];
rz(0.27287896) q[3];
sx q[3];
rz(-1.498184) q[3];
sx q[3];
rz(-0.85736322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56796271) q[2];
sx q[2];
rz(-2.6831388) q[2];
sx q[2];
rz(-2.6652794) q[2];
rz(-1.025398) q[3];
sx q[3];
rz(-1.7849331) q[3];
sx q[3];
rz(-2.8477113) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77182257) q[0];
sx q[0];
rz(-2.2875468) q[0];
sx q[0];
rz(-0.59637946) q[0];
rz(2.3884933) q[1];
sx q[1];
rz(-1.3558847) q[1];
sx q[1];
rz(-2.7900043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6839978) q[0];
sx q[0];
rz(-2.5980934) q[0];
sx q[0];
rz(0.25590055) q[0];
x q[1];
rz(0.48248191) q[2];
sx q[2];
rz(-2.1957779) q[2];
sx q[2];
rz(0.04908726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90801566) q[1];
sx q[1];
rz(-1.1777517) q[1];
sx q[1];
rz(0.24891757) q[1];
rz(-pi) q[2];
rz(2.5078689) q[3];
sx q[3];
rz(-1.0175287) q[3];
sx q[3];
rz(-2.3330757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72254649) q[2];
sx q[2];
rz(-1.498035) q[2];
sx q[2];
rz(0.92918116) q[2];
rz(1.2292713) q[3];
sx q[3];
rz(-0.036673948) q[3];
sx q[3];
rz(-1.8721972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96255985) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(-0.52781934) q[0];
rz(0.57012308) q[1];
sx q[1];
rz(-2.5716883) q[1];
sx q[1];
rz(-0.39400563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45371374) q[0];
sx q[0];
rz(-1.3766995) q[0];
sx q[0];
rz(-2.6537623) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6232004) q[2];
sx q[2];
rz(-2.4817991) q[2];
sx q[2];
rz(-2.348748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2886823) q[1];
sx q[1];
rz(-2.8894419) q[1];
sx q[1];
rz(2.0618477) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3287613) q[3];
sx q[3];
rz(-2.1413229) q[3];
sx q[3];
rz(1.3798483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3923308) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(-2.0070845) q[2];
rz(0.025731651) q[3];
sx q[3];
rz(-2.0870356) q[3];
sx q[3];
rz(-2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5662017) q[0];
sx q[0];
rz(-1.3302777) q[0];
sx q[0];
rz(-2.6824685) q[0];
rz(-2.3205914) q[1];
sx q[1];
rz(-2.2586925) q[1];
sx q[1];
rz(-0.51148907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82152459) q[0];
sx q[0];
rz(-1.8628768) q[0];
sx q[0];
rz(1.0758557) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3364001) q[2];
sx q[2];
rz(-2.2210741) q[2];
sx q[2];
rz(2.5087207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74384159) q[1];
sx q[1];
rz(-1.4201179) q[1];
sx q[1];
rz(2.56714) q[1];
x q[2];
rz(-0.13798519) q[3];
sx q[3];
rz(-2.249794) q[3];
sx q[3];
rz(0.0019055923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56200999) q[2];
sx q[2];
rz(-1.1082114) q[2];
sx q[2];
rz(-0.75135922) q[2];
rz(-2.5747418) q[3];
sx q[3];
rz(-1.5375429) q[3];
sx q[3];
rz(-1.8142627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71996561) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(0.0075465329) q[0];
rz(-0.94002062) q[1];
sx q[1];
rz(-0.66497856) q[1];
sx q[1];
rz(3.0090581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7203001) q[0];
sx q[0];
rz(-2.1425035) q[0];
sx q[0];
rz(2.6470143) q[0];
rz(-pi) q[1];
rz(2.7785702) q[2];
sx q[2];
rz(-1.6953812) q[2];
sx q[2];
rz(1.7302135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79515565) q[1];
sx q[1];
rz(-2.7043685) q[1];
sx q[1];
rz(2.6889685) q[1];
x q[2];
rz(-2.0293983) q[3];
sx q[3];
rz(-0.37887869) q[3];
sx q[3];
rz(-2.1357812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14564766) q[2];
sx q[2];
rz(-1.6463248) q[2];
sx q[2];
rz(3.0403467) q[2];
rz(-2.4920987) q[3];
sx q[3];
rz(-0.5972623) q[3];
sx q[3];
rz(-0.35972843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7262909) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(1.8743961) q[0];
rz(1.2665117) q[1];
sx q[1];
rz(-0.15788618) q[1];
sx q[1];
rz(-0.74323851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7605054) q[0];
sx q[0];
rz(-2.532027) q[0];
sx q[0];
rz(-2.7838092) q[0];
x q[1];
rz(0.27979346) q[2];
sx q[2];
rz(-0.81088561) q[2];
sx q[2];
rz(-1.5423797) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2253902) q[1];
sx q[1];
rz(-1.7280726) q[1];
sx q[1];
rz(1.3316403) q[1];
rz(-pi) q[2];
rz(-2.2970639) q[3];
sx q[3];
rz(-0.8958848) q[3];
sx q[3];
rz(2.9425987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0172952) q[2];
sx q[2];
rz(-0.20169078) q[2];
sx q[2];
rz(-1.7151493) q[2];
rz(2.1258449) q[3];
sx q[3];
rz(-1.4371212) q[3];
sx q[3];
rz(-2.3118741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.5947241) q[0];
sx q[0];
rz(-0.65383738) q[0];
sx q[0];
rz(-0.015901707) q[0];
rz(-1.5025899) q[1];
sx q[1];
rz(-0.67633164) q[1];
sx q[1];
rz(-2.1471088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2833982) q[0];
sx q[0];
rz(-0.97351979) q[0];
sx q[0];
rz(-2.6676763) q[0];
x q[1];
rz(-0.15764641) q[2];
sx q[2];
rz(-2.7589643) q[2];
sx q[2];
rz(1.9668818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3164036) q[1];
sx q[1];
rz(-2.0540407) q[1];
sx q[1];
rz(-0.8512599) q[1];
rz(1.307147) q[3];
sx q[3];
rz(-2.6126475) q[3];
sx q[3];
rz(-2.6890713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9859621) q[2];
sx q[2];
rz(-1.4302284) q[2];
sx q[2];
rz(2.963781) q[2];
rz(-1.6832247) q[3];
sx q[3];
rz(-0.42176133) q[3];
sx q[3];
rz(1.2678857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47852248) q[0];
sx q[0];
rz(-1.8926184) q[0];
sx q[0];
rz(-1.857969) q[0];
rz(0.78370699) q[1];
sx q[1];
rz(-2.3678534) q[1];
sx q[1];
rz(-1.3073889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6790153) q[0];
sx q[0];
rz(-1.1333915) q[0];
sx q[0];
rz(1.4694197) q[0];
rz(1.4965335) q[2];
sx q[2];
rz(-2.7527202) q[2];
sx q[2];
rz(-2.996794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5897023) q[1];
sx q[1];
rz(-1.3968727) q[1];
sx q[1];
rz(-0.059373418) q[1];
rz(-0.94759192) q[3];
sx q[3];
rz(-2.4952808) q[3];
sx q[3];
rz(1.0912967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54138294) q[2];
sx q[2];
rz(-1.7317438) q[2];
sx q[2];
rz(-0.48995885) q[2];
rz(2.0269035) q[3];
sx q[3];
rz(-1.1763108) q[3];
sx q[3];
rz(-0.32976845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212696) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(1.2976788) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(1.8280468) q[2];
sx q[2];
rz(-1.3197109) q[2];
sx q[2];
rz(-0.33725658) q[2];
rz(0.86434563) q[3];
sx q[3];
rz(-1.723524) q[3];
sx q[3];
rz(1.9836457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
