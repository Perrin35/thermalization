OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(5.4260317) q[0];
sx q[0];
rz(9.5572588) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(-0.58499709) q[1];
sx q[1];
rz(2.4490228) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4953295) q[0];
sx q[0];
rz(-1.4709934) q[0];
sx q[0];
rz(1.050759) q[0];
rz(2.4970826) q[2];
sx q[2];
rz(-1.95382) q[2];
sx q[2];
rz(-2.2565947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8707667) q[1];
sx q[1];
rz(-1.4987136) q[1];
sx q[1];
rz(-1.6986398) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0971783) q[3];
sx q[3];
rz(-1.3485104) q[3];
sx q[3];
rz(-1.8412631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0341558) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.6050603) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(-2.5800887) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(0.5805648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8037572) q[0];
sx q[0];
rz(-0.32807402) q[0];
sx q[0];
rz(0.34123811) q[0];
rz(-pi) q[1];
rz(2.850769) q[2];
sx q[2];
rz(-0.65867701) q[2];
sx q[2];
rz(2.3386699) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3851871) q[1];
sx q[1];
rz(-2.2573973) q[1];
sx q[1];
rz(2.3910206) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92506261) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(-0.95228449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4124174) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(0.87835971) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(-2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.47675258) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(-3.1402918) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(0.032827854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2680227) q[0];
sx q[0];
rz(-0.39186726) q[0];
sx q[0];
rz(-2.5001532) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29067729) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(-2.3454587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82959158) q[1];
sx q[1];
rz(-0.81273505) q[1];
sx q[1];
rz(-0.22711046) q[1];
rz(-pi) q[2];
rz(-1.0333943) q[3];
sx q[3];
rz(-2.5279547) q[3];
sx q[3];
rz(-2.9463241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(-0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70401496) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-2.8787956) q[0];
rz(-0.2335877) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(2.3707726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3931657) q[0];
sx q[0];
rz(-1.5512726) q[0];
sx q[0];
rz(3.1251735) q[0];
rz(-pi) q[1];
rz(3.0371573) q[2];
sx q[2];
rz(-2.0586788) q[2];
sx q[2];
rz(1.3654396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96344906) q[1];
sx q[1];
rz(-1.932486) q[1];
sx q[1];
rz(2.5835035) q[1];
x q[2];
rz(0.869107) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(-1.7770191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(0.7129933) q[2];
rz(-2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.4404526) q[0];
rz(-0.094093181) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1124135) q[0];
sx q[0];
rz(-0.83054435) q[0];
sx q[0];
rz(-1.2519757) q[0];
rz(-1.1414358) q[2];
sx q[2];
rz(-0.51157198) q[2];
sx q[2];
rz(2.9432952) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.070378455) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(2.1945164) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2682297) q[3];
sx q[3];
rz(-1.0358827) q[3];
sx q[3];
rz(1.3267335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99469441) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4858522) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(-2.6830542) q[0];
rz(-0.25587747) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(0.68516723) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348649) q[0];
sx q[0];
rz(-2.2198338) q[0];
sx q[0];
rz(3.0543324) q[0];
rz(-0.14641996) q[2];
sx q[2];
rz(-2.358987) q[2];
sx q[2];
rz(0.27461068) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9381147) q[1];
sx q[1];
rz(-2.1739829) q[1];
sx q[1];
rz(2.5342062) q[1];
rz(-pi) q[2];
x q[2];
rz(1.54322) q[3];
sx q[3];
rz(-2.555489) q[3];
sx q[3];
rz(0.54857777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3926065) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(2.4122453) q[0];
rz(2.8485281) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(-1.1475295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80619752) q[0];
sx q[0];
rz(-1.9837539) q[0];
sx q[0];
rz(1.6523244) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2044264) q[2];
sx q[2];
rz(-1.2680149) q[2];
sx q[2];
rz(-2.5123951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7123375) q[1];
sx q[1];
rz(-2.0550248) q[1];
sx q[1];
rz(-1.5244563) q[1];
rz(-pi) q[2];
rz(-0.64671867) q[3];
sx q[3];
rz(-2.5463856) q[3];
sx q[3];
rz(-1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3532233) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-2.3925171) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(-2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99326837) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.3759026) q[1];
sx q[1];
rz(-2.3283236) q[1];
sx q[1];
rz(2.7430699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24005323) q[0];
sx q[0];
rz(-1.1050637) q[0];
sx q[0];
rz(-0.17850152) q[0];
rz(-pi) q[1];
rz(-0.60110809) q[2];
sx q[2];
rz(-1.6582487) q[2];
sx q[2];
rz(2.2488038) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1426705) q[1];
sx q[1];
rz(-2.167836) q[1];
sx q[1];
rz(-3.1139657) q[1];
x q[2];
rz(-0.51316349) q[3];
sx q[3];
rz(-2.0675821) q[3];
sx q[3];
rz(2.7548807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9514256) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(-1.8438967) q[2];
rz(2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(2.4979533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012638906) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(-1.5147491) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(2.0432037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84868816) q[0];
sx q[0];
rz(-2.0540037) q[0];
sx q[0];
rz(2.8781761) q[0];
rz(-pi) q[1];
rz(-0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(-0.84504499) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92022773) q[1];
sx q[1];
rz(-2.5623294) q[1];
sx q[1];
rz(0.00097640493) q[1];
rz(-pi) q[2];
rz(1.7421726) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(-0.51952726) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39524233) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(2.675132) q[0];
rz(-2.9699504) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(0.62896532) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91613149) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(-1.7037237) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2538818) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(-0.72310477) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2168857) q[1];
sx q[1];
rz(-0.44210426) q[1];
sx q[1];
rz(3.1367723) q[1];
x q[2];
rz(-2.7528796) q[3];
sx q[3];
rz(-2.3832088) q[3];
sx q[3];
rz(-0.10314108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8978867) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(-2.0824599) q[2];
rz(2.4641666) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28329904) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(1.614744) q[2];
sx q[2];
rz(-2.2938001) q[2];
sx q[2];
rz(-3.0707035) q[2];
rz(-1.4678636) q[3];
sx q[3];
rz(-2.5522305) q[3];
sx q[3];
rz(1.1141368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
