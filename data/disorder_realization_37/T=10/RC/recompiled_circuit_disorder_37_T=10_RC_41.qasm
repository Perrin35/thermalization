OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(-1.1480968) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6205665) q[0];
sx q[0];
rz(-1.5728083) q[0];
sx q[0];
rz(-2.8319671) q[0];
rz(-0.78656466) q[2];
sx q[2];
rz(-2.2872891) q[2];
sx q[2];
rz(-2.5094945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72109556) q[1];
sx q[1];
rz(-1.3008529) q[1];
sx q[1];
rz(-2.6707778) q[1];
rz(-pi) q[2];
rz(2.8586219) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(-0.52373826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(-2.9623048) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(0.13277408) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(2.9002088) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86357147) q[0];
sx q[0];
rz(-1.2153421) q[0];
sx q[0];
rz(1.1583207) q[0];
rz(-pi) q[1];
rz(0.74626211) q[2];
sx q[2];
rz(-0.66821874) q[2];
sx q[2];
rz(-1.8909188) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0028210359) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(0.9072733) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4003795) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(-1.9433598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(-2.4011491) q[0];
rz(-0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(-2.0887451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2967865) q[0];
sx q[0];
rz(-2.1150981) q[0];
sx q[0];
rz(1.143572) q[0];
rz(-pi) q[1];
rz(-1.2028678) q[2];
sx q[2];
rz(-2.433948) q[2];
sx q[2];
rz(2.5840273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82008541) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(2.4732175) q[1];
rz(-pi) q[2];
rz(3.1316109) q[3];
sx q[3];
rz(-1.1407033) q[3];
sx q[3];
rz(-0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(0.012399013) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(0.23637493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0683595) q[0];
sx q[0];
rz(-1.7422361) q[0];
sx q[0];
rz(1.839848) q[0];
x q[1];
rz(-0.9985853) q[2];
sx q[2];
rz(-0.59362715) q[2];
sx q[2];
rz(-1.3015391) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.826556) q[1];
sx q[1];
rz(-2.9737925) q[1];
sx q[1];
rz(-1.6102953) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27487367) q[3];
sx q[3];
rz(-2.636424) q[3];
sx q[3];
rz(1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(-2.5881055) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(-2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-2.4374403) q[0];
rz(1.0559121) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0740972) q[0];
sx q[0];
rz(-1.6620381) q[0];
sx q[0];
rz(1.3542961) q[0];
rz(-1.709183) q[2];
sx q[2];
rz(-0.90183898) q[2];
sx q[2];
rz(1.5630747) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.039918598) q[1];
sx q[1];
rz(-0.44821757) q[1];
sx q[1];
rz(0.43804534) q[1];
x q[2];
rz(1.7549873) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(-2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(2.5904783) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(-1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(2.1898848) q[0];
rz(2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2518894) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(2.9910812) q[0];
rz(-pi) q[1];
rz(-1.4163383) q[2];
sx q[2];
rz(-1.8784349) q[2];
sx q[2];
rz(-2.6442106) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8413137) q[1];
sx q[1];
rz(-1.0611371) q[1];
sx q[1];
rz(-2.120244) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0814249) q[3];
sx q[3];
rz(-1.6276629) q[3];
sx q[3];
rz(2.1631654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.6850083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(2.899535) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(0.40329969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2780669) q[0];
sx q[0];
rz(-2.6289872) q[0];
sx q[0];
rz(0.82786064) q[0];
rz(-1.8554011) q[2];
sx q[2];
rz(-0.91214123) q[2];
sx q[2];
rz(-1.8919485) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62337263) q[1];
sx q[1];
rz(-1.5689335) q[1];
sx q[1];
rz(-1.3606447) q[1];
rz(-pi) q[2];
rz(-3.1163252) q[3];
sx q[3];
rz(-1.7606887) q[3];
sx q[3];
rz(-0.59786284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(2.1408391) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-2.5411434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6143215) q[0];
sx q[0];
rz(-1.9691756) q[0];
sx q[0];
rz(-0.50171731) q[0];
rz(-pi) q[1];
rz(3.1152994) q[2];
sx q[2];
rz(-1.0458938) q[2];
sx q[2];
rz(-2.3412447) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32158347) q[1];
sx q[1];
rz(-1.4516186) q[1];
sx q[1];
rz(-2.9585341) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59372254) q[3];
sx q[3];
rz(-2.7799118) q[3];
sx q[3];
rz(-0.42632521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0042469) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-3.0252769) q[2];
rz(0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(-2.2019745) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-2.7240662) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4672887) q[0];
sx q[0];
rz(-2.1018565) q[0];
sx q[0];
rz(-0.051961016) q[0];
rz(-pi) q[1];
rz(-1.19403) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(-0.49809581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36996499) q[1];
sx q[1];
rz(-1.7376448) q[1];
sx q[1];
rz(-0.19373993) q[1];
rz(-pi) q[2];
rz(1.1674676) q[3];
sx q[3];
rz(-2.6780193) q[3];
sx q[3];
rz(-1.4232672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4801165) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(2.647906) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(-0.68558145) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(2.0102274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81998527) q[0];
sx q[0];
rz(-1.2258343) q[0];
sx q[0];
rz(2.5489775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2106902) q[2];
sx q[2];
rz(-0.97133342) q[2];
sx q[2];
rz(-2.714389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2724534) q[1];
sx q[1];
rz(-1.9483882) q[1];
sx q[1];
rz(-1.3294405) q[1];
rz(-pi) q[2];
rz(0.64254909) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(-2.97314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3161105) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(0.99954188) q[2];
sx q[2];
rz(-1.0456325) q[2];
sx q[2];
rz(0.13620302) q[2];
rz(-2.2255696) q[3];
sx q[3];
rz(-0.81867221) q[3];
sx q[3];
rz(1.4233854) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
