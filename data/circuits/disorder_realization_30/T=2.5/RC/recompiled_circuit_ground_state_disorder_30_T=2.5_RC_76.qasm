OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(2.365132) q[0];
sx q[0];
rz(9.790701) q[0];
rz(0.34691063) q[1];
sx q[1];
rz(-0.28471714) q[1];
sx q[1];
rz(-1.3887583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2251988) q[0];
sx q[0];
rz(-1.3687736) q[0];
sx q[0];
rz(-1.6388157) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2728782) q[2];
sx q[2];
rz(-2.2414444) q[2];
sx q[2];
rz(-0.16042319) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5947111) q[1];
sx q[1];
rz(-1.1818019) q[1];
sx q[1];
rz(-1.0884029) q[1];
x q[2];
rz(-0.71309375) q[3];
sx q[3];
rz(-2.4131107) q[3];
sx q[3];
rz(2.147004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1846788) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(2.206395) q[2];
rz(-0.30185559) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(-2.7479318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1644156) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(2.6432977) q[0];
rz(1.7138819) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(-1.0867585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1234498) q[0];
sx q[0];
rz(-0.88844127) q[0];
sx q[0];
rz(-1.8773954) q[0];
rz(-pi) q[1];
rz(0.94765122) q[2];
sx q[2];
rz(-2.8171879) q[2];
sx q[2];
rz(2.0632921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77176911) q[1];
sx q[1];
rz(-0.46716094) q[1];
sx q[1];
rz(-0.36231706) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5828176) q[3];
sx q[3];
rz(-1.88899) q[3];
sx q[3];
rz(-0.98230045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9429417) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(-1.9940935) q[3];
sx q[3];
rz(-1.5787326) q[3];
sx q[3];
rz(1.8089627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9320817) q[0];
sx q[0];
rz(-1.7387583) q[0];
sx q[0];
rz(-2.984356) q[0];
rz(1.2278185) q[1];
sx q[1];
rz(-2.9324052) q[1];
sx q[1];
rz(2.7195209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3738689) q[0];
sx q[0];
rz(-1.6607303) q[0];
sx q[0];
rz(1.7492245) q[0];
rz(-0.55251212) q[2];
sx q[2];
rz(-2.8185511) q[2];
sx q[2];
rz(-0.99977899) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6032519) q[1];
sx q[1];
rz(-1.9515688) q[1];
sx q[1];
rz(-0.095007665) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8019358) q[3];
sx q[3];
rz(-1.3462596) q[3];
sx q[3];
rz(-0.79659407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.813039) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(1.0775393) q[3];
sx q[3];
rz(-0.68818337) q[3];
sx q[3];
rz(0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828736) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(-3.1251113) q[0];
rz(3.0670498) q[1];
sx q[1];
rz(-0.45769474) q[1];
sx q[1];
rz(2.8864536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1265565) q[0];
sx q[0];
rz(-1.4269967) q[0];
sx q[0];
rz(-0.48499996) q[0];
rz(-pi) q[1];
x q[1];
rz(2.071923) q[2];
sx q[2];
rz(-0.6269031) q[2];
sx q[2];
rz(1.7834341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94261547) q[1];
sx q[1];
rz(-1.9991373) q[1];
sx q[1];
rz(2.4498536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.221352) q[3];
sx q[3];
rz(-1.4636875) q[3];
sx q[3];
rz(-2.8529594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61465803) q[2];
sx q[2];
rz(-0.69710985) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(-0.091863306) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(0.45472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62234539) q[0];
sx q[0];
rz(-2.3212101) q[0];
sx q[0];
rz(-1.1118332) q[0];
rz(2.9513997) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(1.9817188) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96944189) q[0];
sx q[0];
rz(-1.6052725) q[0];
sx q[0];
rz(-1.9146754) q[0];
rz(-pi) q[1];
rz(-2.2091523) q[2];
sx q[2];
rz(-1.4529422) q[2];
sx q[2];
rz(1.3502075) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3768757) q[1];
sx q[1];
rz(-2.5333166) q[1];
sx q[1];
rz(2.6981955) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0489063) q[3];
sx q[3];
rz(-1.5318222) q[3];
sx q[3];
rz(-3.1119973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5707034) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(-2.0795836) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45288169) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(-1.3558615) q[0];
rz(2.611825) q[1];
sx q[1];
rz(-0.7822839) q[1];
sx q[1];
rz(-1.9897602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36042155) q[0];
sx q[0];
rz(-0.36995212) q[0];
sx q[0];
rz(-0.2759224) q[0];
x q[1];
rz(2.1124338) q[2];
sx q[2];
rz(-1.6525606) q[2];
sx q[2];
rz(-2.0643016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.85688389) q[1];
sx q[1];
rz(-1.9344866) q[1];
sx q[1];
rz(-1.7784575) q[1];
rz(-pi) q[2];
rz(-1.4532667) q[3];
sx q[3];
rz(-1.3580048) q[3];
sx q[3];
rz(-2.8778587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6666224) q[2];
sx q[2];
rz(-2.5002561) q[2];
sx q[2];
rz(2.6744911) q[2];
rz(-2.1304255) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(-1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8491299) q[0];
sx q[0];
rz(-0.15501538) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(-2.977773) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(-0.64635578) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8317597) q[0];
sx q[0];
rz(-2.0947959) q[0];
sx q[0];
rz(-0.28964596) q[0];
rz(-pi) q[1];
rz(-2.9235513) q[2];
sx q[2];
rz(-1.8865117) q[2];
sx q[2];
rz(-1.9670847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19749459) q[1];
sx q[1];
rz(-1.0284327) q[1];
sx q[1];
rz(-1.4828862) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38742806) q[3];
sx q[3];
rz(-1.4447851) q[3];
sx q[3];
rz(-2.2962017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.15567638) q[2];
sx q[2];
rz(-1.6935657) q[2];
sx q[2];
rz(0.4099561) q[2];
rz(-0.83034596) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(1.693694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.603867) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(2.895288) q[0];
rz(2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(0.26396096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1896828) q[0];
sx q[0];
rz(-0.095746843) q[0];
sx q[0];
rz(1.4950947) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0954148) q[2];
sx q[2];
rz(-2.2959297) q[2];
sx q[2];
rz(-0.96607581) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3166271) q[1];
sx q[1];
rz(-2.3967603) q[1];
sx q[1];
rz(1.1492351) q[1];
x q[2];
rz(-2.554515) q[3];
sx q[3];
rz(-1.2012606) q[3];
sx q[3];
rz(-1.0472421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0366514) q[2];
sx q[2];
rz(-1.8737326) q[2];
sx q[2];
rz(0.71511739) q[2];
rz(-0.07130833) q[3];
sx q[3];
rz(-1.5846059) q[3];
sx q[3];
rz(-2.3947072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0095373) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(-0.43553964) q[0];
rz(-2.9938193) q[1];
sx q[1];
rz(-1.4316033) q[1];
sx q[1];
rz(-1.4685644) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2117812) q[0];
sx q[0];
rz(-0.47168293) q[0];
sx q[0];
rz(-2.7372772) q[0];
rz(-pi) q[1];
rz(-1.4917489) q[2];
sx q[2];
rz(-2.2666396) q[2];
sx q[2];
rz(-2.8082113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10994153) q[1];
sx q[1];
rz(-1.1835956) q[1];
sx q[1];
rz(2.9613858) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47504039) q[3];
sx q[3];
rz(-2.6938525) q[3];
sx q[3];
rz(-1.0645176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74464166) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(2.2486539) q[2];
rz(-1.8630155) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(-1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29499149) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(-0.60605979) q[0];
rz(-0.010295708) q[1];
sx q[1];
rz(-2.1538487) q[1];
sx q[1];
rz(3.0616679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.815926) q[0];
sx q[0];
rz(-1.0116966) q[0];
sx q[0];
rz(-2.8432106) q[0];
x q[1];
rz(2.299231) q[2];
sx q[2];
rz(-1.3536705) q[2];
sx q[2];
rz(-2.6922243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2458107) q[1];
sx q[1];
rz(-1.5930098) q[1];
sx q[1];
rz(2.7345045) q[1];
rz(-pi) q[2];
rz(0.84385033) q[3];
sx q[3];
rz(-0.30442134) q[3];
sx q[3];
rz(-0.756504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0004878) q[2];
sx q[2];
rz(-0.92685574) q[2];
sx q[2];
rz(-1.0151939) q[2];
rz(2.5403533) q[3];
sx q[3];
rz(-2.4398118) q[3];
sx q[3];
rz(1.0465485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.64423185) q[0];
sx q[0];
rz(-1.624122) q[0];
sx q[0];
rz(-2.4540785) q[0];
rz(-0.58912206) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(-2.9420058) q[2];
sx q[2];
rz(-0.45303405) q[2];
sx q[2];
rz(-0.41336811) q[2];
rz(2.3169869) q[3];
sx q[3];
rz(-1.6578703) q[3];
sx q[3];
rz(-0.35005611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
