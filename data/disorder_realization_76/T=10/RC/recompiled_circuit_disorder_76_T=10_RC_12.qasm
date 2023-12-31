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
rz(5.6981882) q[1];
sx q[1];
rz(11.873801) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4953295) q[0];
sx q[0];
rz(-1.4709934) q[0];
sx q[0];
rz(1.050759) q[0];
x q[1];
rz(2.4970826) q[2];
sx q[2];
rz(-1.1877726) q[2];
sx q[2];
rz(2.2565947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9308656) q[1];
sx q[1];
rz(-0.14666808) q[1];
sx q[1];
rz(-2.0861097) q[1];
rz(-pi) q[2];
rz(-1.0971783) q[3];
sx q[3];
rz(-1.7930822) q[3];
sx q[3];
rz(1.3003295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0341558) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261616) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(1.2492299) q[0];
rz(-2.5800887) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(2.5610279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55727977) q[0];
sx q[0];
rz(-1.4627539) q[0];
sx q[0];
rz(0.31038196) q[0];
x q[1];
rz(-1.3524019) q[2];
sx q[2];
rz(-2.1973655) q[2];
sx q[2];
rz(0.44109694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9241582) q[1];
sx q[1];
rz(-0.96964753) q[1];
sx q[1];
rz(2.2648328) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-2.2632329) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(-0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(-0.032827854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093124495) q[0];
sx q[0];
rz(-1.3402407) q[0];
sx q[0];
rz(2.8218517) q[0];
x q[1];
rz(-0.29067729) q[2];
sx q[2];
rz(-1.0376087) q[2];
sx q[2];
rz(2.3454587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.153861) q[1];
sx q[1];
rz(-0.78480936) q[1];
sx q[1];
rz(1.8042817) q[1];
rz(2.7955416) q[3];
sx q[3];
rz(-2.0882574) q[3];
sx q[3];
rz(0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70401496) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3931657) q[0];
sx q[0];
rz(-1.5512726) q[0];
sx q[0];
rz(-0.016419134) q[0];
x q[1];
rz(1.7647676) q[2];
sx q[2];
rz(-0.49805222) q[2];
sx q[2];
rz(1.9961403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3168287) q[1];
sx q[1];
rz(-1.0526122) q[1];
sx q[1];
rz(1.9903241) q[1];
x q[2];
rz(-0.869107) q[3];
sx q[3];
rz(-0.29033717) q[3];
sx q[3];
rz(-1.7770191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9757441) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(-2.4285994) q[2];
rz(-2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(-0.094093181) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(2.9715911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3225587) q[0];
sx q[0];
rz(-1.337262) q[0];
sx q[0];
rz(0.76604953) q[0];
rz(-pi) q[1];
rz(-0.22959392) q[2];
sx q[2];
rz(-2.0320971) q[2];
sx q[2];
rz(0.68179828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.070378455) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(-0.94707625) q[1];
rz(2.2682297) q[3];
sx q[3];
rz(-2.1057099) q[3];
sx q[3];
rz(1.8148592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1468982) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(1.5935625) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(-1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4858522) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(2.8857152) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(2.4564254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79294369) q[0];
sx q[0];
rz(-2.2198338) q[0];
sx q[0];
rz(-0.087260212) q[0];
rz(2.3643656) q[2];
sx q[2];
rz(-1.6738552) q[2];
sx q[2];
rz(-1.4003786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31720686) q[1];
sx q[1];
rz(-0.82815352) q[1];
sx q[1];
rz(2.2627027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98486793) q[3];
sx q[3];
rz(-1.586048) q[3];
sx q[3];
rz(-1.0451942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3926065) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(1.7123327) q[2];
rz(-2.0424992) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(2.8485281) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(-1.9940631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3353951) q[0];
sx q[0];
rz(-1.1578387) q[0];
sx q[0];
rz(-1.6523244) q[0];
rz(1.0853026) q[2];
sx q[2];
rz(-2.448423) q[2];
sx q[2];
rz(1.327141) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8116248) q[1];
sx q[1];
rz(-2.6553272) q[1];
sx q[1];
rz(0.087841308) q[1];
rz(2.6461584) q[3];
sx q[3];
rz(-1.9154275) q[3];
sx q[3];
rz(-0.37638327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3532233) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7300028) q[0];
sx q[0];
rz(-1.7301136) q[0];
sx q[0];
rz(-2.0429862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4648828) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(2.523409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41259137) q[1];
sx q[1];
rz(-1.5479497) q[1];
sx q[1];
rz(2.1680135) q[1];
rz(-pi) q[2];
rz(0.51316349) q[3];
sx q[3];
rz(-2.0675821) q[3];
sx q[3];
rz(0.38671198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9514256) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.297696) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(-0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1289537) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(1.7522316) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(1.0983889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5441355) q[0];
sx q[0];
rz(-1.3381334) q[0];
sx q[0];
rz(-1.0730037) q[0];
rz(-pi) q[1];
rz(-0.82010014) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(2.1095914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92139447) q[1];
sx q[1];
rz(-2.1500593) q[1];
sx q[1];
rz(1.5701576) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8875651) q[3];
sx q[3];
rz(-1.7367749) q[3];
sx q[3];
rz(-1.1248873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(-2.6220654) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-0.46646068) q[0];
rz(-0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(0.62896532) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0590203) q[0];
sx q[0];
rz(-1.9452381) q[0];
sx q[0];
rz(3.0890205) q[0];
rz(1.2538818) q[2];
sx q[2];
rz(-0.66714087) q[2];
sx q[2];
rz(-2.4184879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92470691) q[1];
sx q[1];
rz(-2.6994884) q[1];
sx q[1];
rz(0.0048203992) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7528796) q[3];
sx q[3];
rz(-2.3832088) q[3];
sx q[3];
rz(3.0384516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(-1.0591327) q[2];
rz(-0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.8582936) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(-0.25390608) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(0.049747808) q[2];
sx q[2];
rz(-0.72409734) q[2];
sx q[2];
rz(-3.004336) q[2];
rz(1.6737291) q[3];
sx q[3];
rz(-2.5522305) q[3];
sx q[3];
rz(1.1141368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
