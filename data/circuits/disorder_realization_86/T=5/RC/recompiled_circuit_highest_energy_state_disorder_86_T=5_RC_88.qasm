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
rz(1.6744094) q[0];
sx q[0];
rz(1.6011342) q[0];
sx q[0];
rz(8.1327333) q[0];
rz(2.0436824) q[1];
sx q[1];
rz(-0.75610375) q[1];
sx q[1];
rz(2.4296711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5626418) q[0];
sx q[0];
rz(-1.535784) q[0];
sx q[0];
rz(2.9981311) q[0];
x q[1];
rz(2.8991382) q[2];
sx q[2];
rz(-2.3676845) q[2];
sx q[2];
rz(2.7737308) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8375875) q[1];
sx q[1];
rz(-1.3437573) q[1];
sx q[1];
rz(2.5035437) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6308832) q[3];
sx q[3];
rz(-1.0282955) q[3];
sx q[3];
rz(1.6245019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11382515) q[2];
sx q[2];
rz(-1.0372459) q[2];
sx q[2];
rz(0.85565957) q[2];
rz(-1.8850108) q[3];
sx q[3];
rz(-0.62658566) q[3];
sx q[3];
rz(1.9470866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282495) q[0];
sx q[0];
rz(-0.11390587) q[0];
sx q[0];
rz(-2.2746427) q[0];
rz(-1.4504245) q[1];
sx q[1];
rz(-1.1879286) q[1];
sx q[1];
rz(-0.11122045) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1010872) q[0];
sx q[0];
rz(-1.936628) q[0];
sx q[0];
rz(0.58617298) q[0];
rz(-0.48008474) q[2];
sx q[2];
rz(-1.089555) q[2];
sx q[2];
rz(-0.54779886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4091089) q[1];
sx q[1];
rz(-0.95327158) q[1];
sx q[1];
rz(-0.78583048) q[1];
x q[2];
rz(-1.3830967) q[3];
sx q[3];
rz(-2.4259013) q[3];
sx q[3];
rz(-0.3970428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0337246) q[2];
sx q[2];
rz(-2.3794231) q[2];
sx q[2];
rz(-2.6503358) q[2];
rz(1.8852425) q[3];
sx q[3];
rz(-1.7137073) q[3];
sx q[3];
rz(2.6826732) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1678109) q[0];
sx q[0];
rz(-1.7901006) q[0];
sx q[0];
rz(3.1381881) q[0];
rz(-2.7963474) q[1];
sx q[1];
rz(-2.7216941) q[1];
sx q[1];
rz(-0.87510625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3728947) q[0];
sx q[0];
rz(-0.56124306) q[0];
sx q[0];
rz(3.0054557) q[0];
rz(-0.98710635) q[2];
sx q[2];
rz(-2.9265227) q[2];
sx q[2];
rz(1.2925066) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3899843) q[1];
sx q[1];
rz(-2.4953003) q[1];
sx q[1];
rz(-0.81089748) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3703824) q[3];
sx q[3];
rz(-1.589984) q[3];
sx q[3];
rz(-2.8018708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7615937) q[2];
sx q[2];
rz(-2.9271409) q[2];
sx q[2];
rz(-0.3271412) q[2];
rz(1.7548148) q[3];
sx q[3];
rz(-1.944647) q[3];
sx q[3];
rz(3.0676214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(0.73579329) q[0];
sx q[0];
rz(-2.1402335) q[0];
sx q[0];
rz(-1.5843947) q[0];
rz(2.7994075) q[1];
sx q[1];
rz(-1.0289959) q[1];
sx q[1];
rz(-2.7142966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9153684) q[0];
sx q[0];
rz(-1.5378204) q[0];
sx q[0];
rz(-0.054126496) q[0];
x q[1];
rz(-0.37125606) q[2];
sx q[2];
rz(-2.0780188) q[2];
sx q[2];
rz(-1.4812507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4327185) q[1];
sx q[1];
rz(-2.1418835) q[1];
sx q[1];
rz(2.4576748) q[1];
rz(-pi) q[2];
rz(-1.2702291) q[3];
sx q[3];
rz(-1.0247314) q[3];
sx q[3];
rz(-2.6566243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4348609) q[2];
sx q[2];
rz(-0.92919436) q[2];
sx q[2];
rz(0.84563196) q[2];
rz(-3.1137858) q[3];
sx q[3];
rz(-1.3549201) q[3];
sx q[3];
rz(-1.6644679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61331767) q[0];
sx q[0];
rz(-2.1149825) q[0];
sx q[0];
rz(0.051830526) q[0];
rz(-1.9822281) q[1];
sx q[1];
rz(-2.4720188) q[1];
sx q[1];
rz(2.5313098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51175115) q[0];
sx q[0];
rz(-1.7423706) q[0];
sx q[0];
rz(-2.9338994) q[0];
x q[1];
rz(-2.8437623) q[2];
sx q[2];
rz(-1.1874842) q[2];
sx q[2];
rz(0.9155067) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7597639) q[1];
sx q[1];
rz(-2.0501158) q[1];
sx q[1];
rz(1.3717531) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5984437) q[3];
sx q[3];
rz(-2.1787454) q[3];
sx q[3];
rz(-2.3635918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74756885) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(0.58689153) q[2];
rz(-1.3823357) q[3];
sx q[3];
rz(-2.0943677) q[3];
sx q[3];
rz(0.97572485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0368283) q[0];
sx q[0];
rz(-0.47127518) q[0];
sx q[0];
rz(0.18071827) q[0];
rz(-1.3982999) q[1];
sx q[1];
rz(-1.2340052) q[1];
sx q[1];
rz(-1.3999375) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619059) q[0];
sx q[0];
rz(-0.8041412) q[0];
sx q[0];
rz(1.8820394) q[0];
x q[1];
rz(-2.3467881) q[2];
sx q[2];
rz(-0.53404885) q[2];
sx q[2];
rz(0.42511031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.789098) q[1];
sx q[1];
rz(-2.2667655) q[1];
sx q[1];
rz(0.11274479) q[1];
rz(-pi) q[2];
rz(2.4816957) q[3];
sx q[3];
rz(-1.1142892) q[3];
sx q[3];
rz(-0.61566258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7515298) q[2];
sx q[2];
rz(-0.84948245) q[2];
sx q[2];
rz(-2.0709822) q[2];
rz(-1.4812329) q[3];
sx q[3];
rz(-0.79602066) q[3];
sx q[3];
rz(1.3720366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2299131) q[0];
sx q[0];
rz(-2.3340618) q[0];
sx q[0];
rz(-0.56021488) q[0];
rz(2.920976) q[1];
sx q[1];
rz(-0.95181528) q[1];
sx q[1];
rz(-0.87699786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91966832) q[0];
sx q[0];
rz(-3.1238365) q[0];
sx q[0];
rz(1.4229126) q[0];
rz(-pi) q[1];
rz(0.39773293) q[2];
sx q[2];
rz(-1.3215725) q[2];
sx q[2];
rz(-0.13532369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24755281) q[1];
sx q[1];
rz(-2.8553989) q[1];
sx q[1];
rz(1.9412089) q[1];
rz(1.6991922) q[3];
sx q[3];
rz(-2.4511325) q[3];
sx q[3];
rz(0.70790509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47473869) q[2];
sx q[2];
rz(-0.89484221) q[2];
sx q[2];
rz(-1.3372927) q[2];
rz(-2.9232591) q[3];
sx q[3];
rz(-2.1209013) q[3];
sx q[3];
rz(1.7201503) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6364994) q[0];
sx q[0];
rz(-1.8890843) q[0];
sx q[0];
rz(-0.1388347) q[0];
rz(0.88169634) q[1];
sx q[1];
rz(-1.1602743) q[1];
sx q[1];
rz(-0.82463157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449681) q[0];
sx q[0];
rz(-1.5486778) q[0];
sx q[0];
rz(1.478594) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6753372) q[2];
sx q[2];
rz(-1.6674768) q[2];
sx q[2];
rz(-0.51890512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3089203) q[1];
sx q[1];
rz(-1.7270346) q[1];
sx q[1];
rz(2.6392722) q[1];
x q[2];
rz(0.10004754) q[3];
sx q[3];
rz(-2.2216019) q[3];
sx q[3];
rz(-1.5318961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3489939) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(0.83354956) q[2];
rz(-0.27093568) q[3];
sx q[3];
rz(-2.0201594) q[3];
sx q[3];
rz(-0.84735316) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18868294) q[0];
sx q[0];
rz(-1.8930577) q[0];
sx q[0];
rz(0.46440014) q[0];
rz(0.69350997) q[1];
sx q[1];
rz(-0.92718569) q[1];
sx q[1];
rz(0.69673353) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2385134) q[0];
sx q[0];
rz(-0.12426025) q[0];
sx q[0];
rz(2.0342229) q[0];
rz(-1.3515662) q[2];
sx q[2];
rz(-0.18454177) q[2];
sx q[2];
rz(2.3957555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33524698) q[1];
sx q[1];
rz(-1.5981403) q[1];
sx q[1];
rz(1.3217447) q[1];
x q[2];
rz(0.018432004) q[3];
sx q[3];
rz(-0.64041172) q[3];
sx q[3];
rz(0.7225534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36737475) q[2];
sx q[2];
rz(-3.0493224) q[2];
sx q[2];
rz(-0.83887678) q[2];
rz(-1.687626) q[3];
sx q[3];
rz(-1.5980915) q[3];
sx q[3];
rz(1.4723697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0812747) q[0];
sx q[0];
rz(-1.769279) q[0];
sx q[0];
rz(-1.5244315) q[0];
rz(2.7853107) q[1];
sx q[1];
rz(-1.7226912) q[1];
sx q[1];
rz(1.8024532) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2196428) q[0];
sx q[0];
rz(-0.60254002) q[0];
sx q[0];
rz(1.6842808) q[0];
rz(-1.4993323) q[2];
sx q[2];
rz(-1.2515837) q[2];
sx q[2];
rz(2.6885374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4831743) q[1];
sx q[1];
rz(-0.69026679) q[1];
sx q[1];
rz(-0.19259318) q[1];
x q[2];
rz(1.8175199) q[3];
sx q[3];
rz(-2.9949247) q[3];
sx q[3];
rz(-0.038692729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39849207) q[2];
sx q[2];
rz(-0.83245459) q[2];
sx q[2];
rz(1.2644281) q[2];
rz(0.33665952) q[3];
sx q[3];
rz(-0.94069898) q[3];
sx q[3];
rz(1.8077067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3867415) q[0];
sx q[0];
rz(-1.516153) q[0];
sx q[0];
rz(-1.0990912) q[0];
rz(-0.59973888) q[1];
sx q[1];
rz(-2.6987684) q[1];
sx q[1];
rz(-0.59785688) q[1];
rz(-0.91942922) q[2];
sx q[2];
rz(-0.50242699) q[2];
sx q[2];
rz(-2.4519359) q[2];
rz(-2.609008) q[3];
sx q[3];
rz(-1.5859133) q[3];
sx q[3];
rz(1.5487158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
