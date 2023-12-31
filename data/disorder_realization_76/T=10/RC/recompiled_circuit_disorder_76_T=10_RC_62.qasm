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
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(-2.4490228) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(-1.371944) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64451005) q[2];
sx q[2];
rz(-1.1877726) q[2];
sx q[2];
rz(0.884998) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9308656) q[1];
sx q[1];
rz(-0.14666808) q[1];
sx q[1];
rz(-2.0861097) q[1];
rz(-pi) q[2];
rz(-1.1107221) q[3];
sx q[3];
rz(-0.51957909) q[3];
sx q[3];
rz(-2.4649064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(1.5365323) q[2];
rz(1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(-2.5610279) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5843129) q[0];
sx q[0];
rz(-1.4627539) q[0];
sx q[0];
rz(-2.8312107) q[0];
x q[1];
rz(2.850769) q[2];
sx q[2];
rz(-0.65867701) q[2];
sx q[2];
rz(-0.80292279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3851871) q[1];
sx q[1];
rz(-0.88419534) q[1];
sx q[1];
rz(-0.7505721) q[1];
rz(-pi) q[2];
rz(0.92506261) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(-2.1893082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4124174) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(-0.39204028) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(-0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(-3.1087648) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531909) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(1.8131959) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0228531) q[2];
sx q[2];
rz(-0.60048044) q[2];
sx q[2];
rz(-0.26417363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82959158) q[1];
sx q[1];
rz(-0.81273505) q[1];
sx q[1];
rz(0.22711046) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34605108) q[3];
sx q[3];
rz(-1.0533353) q[3];
sx q[3];
rz(-0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(0.39595655) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-0.26279703) q[0];
rz(-2.908005) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(0.77082005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3931657) q[0];
sx q[0];
rz(-1.5512726) q[0];
sx q[0];
rz(-0.016419134) q[0];
rz(-3.0371573) q[2];
sx q[2];
rz(-2.0586788) q[2];
sx q[2];
rz(-1.3654396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3168287) q[1];
sx q[1];
rz(-2.0889805) q[1];
sx q[1];
rz(1.9903241) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2724857) q[3];
sx q[3];
rz(-0.29033717) q[3];
sx q[3];
rz(1.7770191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16584855) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(2.4285994) q[2];
rz(2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825496) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(0.094093181) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(2.9715911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029179137) q[0];
sx q[0];
rz(-0.83054435) q[0];
sx q[0];
rz(1.2519757) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0988118) q[2];
sx q[2];
rz(-1.7760279) q[2];
sx q[2];
rz(2.1489378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4449094) q[1];
sx q[1];
rz(-1.0001567) q[1];
sx q[1];
rz(2.6615104) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3159536) q[3];
sx q[3];
rz(-0.85074433) q[3];
sx q[3];
rz(2.839098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(1.7374932) q[2];
rz(-1.5480301) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(-1.4353969) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4858522) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(2.4564254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72502575) q[0];
sx q[0];
rz(-1.6402813) q[0];
sx q[0];
rz(2.2216703) q[0];
x q[1];
rz(1.4267169) q[2];
sx q[2];
rz(-2.3428168) q[2];
sx q[2];
rz(3.0720403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.398715) q[1];
sx q[1];
rz(-1.0814953) q[1];
sx q[1];
rz(2.2687885) q[1];
rz(-pi) q[2];
rz(1.54322) q[3];
sx q[3];
rz(-2.555489) q[3];
sx q[3];
rz(-2.5930149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7489862) q[2];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(-0.72934735) q[0];
rz(2.8485281) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(1.1475295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3353951) q[0];
sx q[0];
rz(-1.1578387) q[0];
sx q[0];
rz(1.6523244) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93716623) q[2];
sx q[2];
rz(-1.2680149) q[2];
sx q[2];
rz(-2.5123951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16312576) q[1];
sx q[1];
rz(-1.5297869) q[1];
sx q[1];
rz(2.6569215) q[1];
rz(-pi) q[2];
rz(2.494874) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3532233) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(2.3925171) q[2];
rz(0.64368147) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-2.3102982) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-2.3283236) q[1];
sx q[1];
rz(2.7430699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4115899) q[0];
sx q[0];
rz(-1.7301136) q[0];
sx q[0];
rz(2.0429862) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60110809) q[2];
sx q[2];
rz(-1.6582487) q[2];
sx q[2];
rz(-2.2488038) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7290013) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(0.97357915) q[1];
x q[2];
rz(-2.6284292) q[3];
sx q[3];
rz(-2.0675821) q[3];
sx q[3];
rz(-2.7548807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9514256) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(-0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(1.7522316) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(-2.0432037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84868816) q[0];
sx q[0];
rz(-1.0875889) q[0];
sx q[0];
rz(2.8781761) q[0];
x q[1];
rz(-2.3214925) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(-2.1095914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2201982) q[1];
sx q[1];
rz(-0.99153334) q[1];
sx q[1];
rz(1.571435) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8875651) q[3];
sx q[3];
rz(-1.4048178) q[3];
sx q[3];
rz(-2.0167054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2417458) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(-2.6220654) q[2];
rz(-1.8064921) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-2.675132) q[0];
rz(2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254612) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(-1.4378689) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9009027) q[2];
sx q[2];
rz(-0.94229892) q[2];
sx q[2];
rz(-1.1185874) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91937376) q[1];
sx q[1];
rz(-1.1286976) q[1];
sx q[1];
rz(-1.5685146) q[1];
rz(-2.4217989) q[3];
sx q[3];
rz(-1.3070953) q[3];
sx q[3];
rz(1.3849474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8978867) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(-2.0824599) q[2];
rz(-2.4641666) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(2.8876866) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(-1.614744) q[2];
sx q[2];
rz(-0.84779253) q[2];
sx q[2];
rz(0.070889125) q[2];
rz(3.0729978) q[3];
sx q[3];
rz(-2.1566236) q[3];
sx q[3];
rz(-1.9038283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
