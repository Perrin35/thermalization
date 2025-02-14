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
rz(1.037984) q[0];
sx q[0];
rz(-1.4572287) q[0];
sx q[0];
rz(1.7888223) q[0];
rz(1.1846722) q[1];
sx q[1];
rz(-0.67263043) q[1];
sx q[1];
rz(1.5157359) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11213451) q[0];
sx q[0];
rz(-1.7687651) q[0];
sx q[0];
rz(2.9517118) q[0];
rz(-2.6625161) q[2];
sx q[2];
rz(-2.276439) q[2];
sx q[2];
rz(2.3912663) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75396155) q[1];
sx q[1];
rz(-1.921614) q[1];
sx q[1];
rz(-2.6430348) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58324849) q[3];
sx q[3];
rz(-0.49213675) q[3];
sx q[3];
rz(-2.0898835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0036014) q[2];
sx q[2];
rz(-0.10819745) q[2];
sx q[2];
rz(0.18599621) q[2];
rz(0.018639175) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(-2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828534) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(-2.8232316) q[0];
rz(-0.63703498) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(-2.6723518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2923986) q[0];
sx q[0];
rz(-2.877826) q[0];
sx q[0];
rz(0.7103814) q[0];
rz(-1.9465916) q[2];
sx q[2];
rz(-1.4926806) q[2];
sx q[2];
rz(-2.5004435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0035867) q[1];
sx q[1];
rz(-0.39630656) q[1];
sx q[1];
rz(-2.519964) q[1];
rz(-0.0057159609) q[3];
sx q[3];
rz(-2.2607231) q[3];
sx q[3];
rz(1.5100513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8239173) q[2];
sx q[2];
rz(-0.21024148) q[2];
sx q[2];
rz(-0.69327411) q[2];
rz(1.3200101) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(-2.0461953) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005118) q[0];
sx q[0];
rz(-2.5040099) q[0];
sx q[0];
rz(0.47955036) q[0];
rz(-0.50411049) q[1];
sx q[1];
rz(-1.1481552) q[1];
sx q[1];
rz(-3.0832916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4625464) q[0];
sx q[0];
rz(-0.57314593) q[0];
sx q[0];
rz(-0.48991211) q[0];
rz(-pi) q[1];
rz(2.0928755) q[2];
sx q[2];
rz(-2.1914542) q[2];
sx q[2];
rz(-1.4246556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7229798) q[1];
sx q[1];
rz(-1.449391) q[1];
sx q[1];
rz(-1.0176119) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4543484) q[3];
sx q[3];
rz(-2.589041) q[3];
sx q[3];
rz(2.1228028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38640675) q[2];
sx q[2];
rz(-1.1815973) q[2];
sx q[2];
rz(0.03726658) q[2];
rz(1.5197598) q[3];
sx q[3];
rz(-0.45791364) q[3];
sx q[3];
rz(-0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.94654361) q[0];
sx q[0];
rz(-0.46849546) q[0];
sx q[0];
rz(3.0750437) q[0];
rz(1.5244124) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(0.7349416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733058) q[0];
sx q[0];
rz(-0.6319918) q[0];
sx q[0];
rz(-0.39936622) q[0];
rz(-pi) q[1];
rz(-1.8040015) q[2];
sx q[2];
rz(-0.8542904) q[2];
sx q[2];
rz(-0.75762123) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0385602) q[1];
sx q[1];
rz(-1.8502697) q[1];
sx q[1];
rz(-0.36906645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9959277) q[3];
sx q[3];
rz(-2.2861135) q[3];
sx q[3];
rz(0.13038929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28896004) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(-0.10870474) q[2];
rz(2.2147801) q[3];
sx q[3];
rz(-0.97065297) q[3];
sx q[3];
rz(1.5638117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0842593) q[0];
sx q[0];
rz(-0.63705343) q[0];
sx q[0];
rz(-1.0908701) q[0];
rz(-2.5690761) q[1];
sx q[1];
rz(-1.1895836) q[1];
sx q[1];
rz(-2.3172839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.123468) q[0];
sx q[0];
rz(-1.195692) q[0];
sx q[0];
rz(0.25183046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53218158) q[2];
sx q[2];
rz(-0.99747466) q[2];
sx q[2];
rz(0.5328726) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2870737) q[1];
sx q[1];
rz(-1.9654915) q[1];
sx q[1];
rz(-2.1629107) q[1];
rz(-pi) q[2];
rz(-2.3184954) q[3];
sx q[3];
rz(-0.56328008) q[3];
sx q[3];
rz(0.34264229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6606286) q[2];
sx q[2];
rz(-0.97233665) q[2];
sx q[2];
rz(0.29120293) q[2];
rz(0.63699841) q[3];
sx q[3];
rz(-1.4642508) q[3];
sx q[3];
rz(-1.8891107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3785192) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(0.13105233) q[0];
rz(-1.1669) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(-0.022620591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1194099) q[0];
sx q[0];
rz(-1.9227322) q[0];
sx q[0];
rz(0.58076136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6940704) q[2];
sx q[2];
rz(-1.7825812) q[2];
sx q[2];
rz(2.663954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5675332) q[1];
sx q[1];
rz(-2.1173734) q[1];
sx q[1];
rz(-0.36336561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7637445) q[3];
sx q[3];
rz(-1.4111327) q[3];
sx q[3];
rz(-0.99332383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8733946) q[2];
sx q[2];
rz(-0.70316535) q[2];
sx q[2];
rz(2.0568636) q[2];
rz(-1.1449413) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(2.1196608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5893843) q[0];
sx q[0];
rz(-1.6227868) q[0];
sx q[0];
rz(-0.4735221) q[0];
rz(-1.9206958) q[1];
sx q[1];
rz(-0.41243204) q[1];
sx q[1];
rz(-1.8730877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80618033) q[0];
sx q[0];
rz(-2.0783193) q[0];
sx q[0];
rz(-0.96854676) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94077295) q[2];
sx q[2];
rz(-0.93859276) q[2];
sx q[2];
rz(-0.8559627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.699325) q[1];
sx q[1];
rz(-1.066477) q[1];
sx q[1];
rz(2.0513351) q[1];
rz(-pi) q[2];
x q[2];
rz(1.394519) q[3];
sx q[3];
rz(-0.44197861) q[3];
sx q[3];
rz(0.37261841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97741693) q[2];
sx q[2];
rz(-1.5693393) q[2];
sx q[2];
rz(0.24641985) q[2];
rz(2.1988846) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(-2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3281658) q[0];
sx q[0];
rz(-1.8712217) q[0];
sx q[0];
rz(-3.0810007) q[0];
rz(-1.6391485) q[1];
sx q[1];
rz(-1.363058) q[1];
sx q[1];
rz(-1.5477808) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0133724) q[0];
sx q[0];
rz(-0.0039847535) q[0];
sx q[0];
rz(2.8811923) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8558989) q[2];
sx q[2];
rz(-2.5473352) q[2];
sx q[2];
rz(-0.93462925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1051576) q[1];
sx q[1];
rz(-0.97291017) q[1];
sx q[1];
rz(0.5374686) q[1];
rz(2.0261062) q[3];
sx q[3];
rz(-2.2336965) q[3];
sx q[3];
rz(-2.2772706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8760406) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(-1.457224) q[2];
rz(-0.11988104) q[3];
sx q[3];
rz(-0.20457743) q[3];
sx q[3];
rz(-2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(-2.0062398) q[0];
rz(-0.10336939) q[1];
sx q[1];
rz(-1.4048978) q[1];
sx q[1];
rz(-2.1939383) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1670939) q[0];
sx q[0];
rz(-1.5038953) q[0];
sx q[0];
rz(-0.6386021) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.060076) q[2];
sx q[2];
rz(-0.7502509) q[2];
sx q[2];
rz(0.81847755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7458742) q[1];
sx q[1];
rz(-1.304879) q[1];
sx q[1];
rz(-0.87541343) q[1];
x q[2];
rz(2.3941764) q[3];
sx q[3];
rz(-2.1028165) q[3];
sx q[3];
rz(-2.4126787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4244298) q[2];
sx q[2];
rz(-2.5134176) q[2];
sx q[2];
rz(-0.90432811) q[2];
rz(-1.2633911) q[3];
sx q[3];
rz(-1.9046013) q[3];
sx q[3];
rz(0.5164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.271027) q[0];
sx q[0];
rz(-0.75815433) q[0];
sx q[0];
rz(-0.98536056) q[0];
rz(0.35161463) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(2.7947289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064929334) q[0];
sx q[0];
rz(-2.4682625) q[0];
sx q[0];
rz(1.2162186) q[0];
x q[1];
rz(1.8029593) q[2];
sx q[2];
rz(-2.684991) q[2];
sx q[2];
rz(2.5176164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0488009) q[1];
sx q[1];
rz(-1.6597224) q[1];
sx q[1];
rz(-2.4503178) q[1];
x q[2];
rz(-1.204416) q[3];
sx q[3];
rz(-1.1561596) q[3];
sx q[3];
rz(0.36568991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88006192) q[2];
sx q[2];
rz(-2.4226649) q[2];
sx q[2];
rz(-0.13492179) q[2];
rz(0.12923446) q[3];
sx q[3];
rz(-0.40004572) q[3];
sx q[3];
rz(2.1028886) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226817) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(-0.46161721) q[1];
sx q[1];
rz(-1.7414265) q[1];
sx q[1];
rz(0.19207676) q[1];
rz(2.2444637) q[2];
sx q[2];
rz(-1.2012325) q[2];
sx q[2];
rz(2.7129632) q[2];
rz(0.57030525) q[3];
sx q[3];
rz(-2.8371596) q[3];
sx q[3];
rz(0.53573487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
