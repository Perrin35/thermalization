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
rz(-2.3075624) q[0];
sx q[0];
rz(2.7355255) q[0];
sx q[0];
rz(14.526215) q[0];
rz(1.740068) q[1];
sx q[1];
rz(5.8086173) q[1];
sx q[1];
rz(9.557815) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2849738) q[0];
sx q[0];
rz(-2.4388038) q[0];
sx q[0];
rz(-1.1590336) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1263761) q[2];
sx q[2];
rz(-2.0410178) q[2];
sx q[2];
rz(-2.0951867) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.592748) q[1];
sx q[1];
rz(-1.8858375) q[1];
sx q[1];
rz(-0.23579072) q[1];
rz(-pi) q[2];
rz(2.7772589) q[3];
sx q[3];
rz(-1.0291463) q[3];
sx q[3];
rz(-2.2026103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1349606) q[2];
sx q[2];
rz(-2.1217608) q[2];
sx q[2];
rz(0.68438619) q[2];
rz(1.1413752) q[3];
sx q[3];
rz(-2.8751774) q[3];
sx q[3];
rz(3.0548837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619277) q[0];
sx q[0];
rz(-0.68988887) q[0];
sx q[0];
rz(0.46098125) q[0];
rz(1.1191248) q[1];
sx q[1];
rz(-0.42367595) q[1];
sx q[1];
rz(0.31203312) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0658995) q[0];
sx q[0];
rz(-1.5830432) q[0];
sx q[0];
rz(3.1135983) q[0];
rz(-pi) q[1];
rz(0.076032555) q[2];
sx q[2];
rz(-2.3302148) q[2];
sx q[2];
rz(1.8234314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2300517) q[1];
sx q[1];
rz(-1.5649867) q[1];
sx q[1];
rz(-1.3951541) q[1];
x q[2];
rz(2.8617763) q[3];
sx q[3];
rz(-2.2035965) q[3];
sx q[3];
rz(-1.0029083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1877039) q[2];
sx q[2];
rz(-2.0362594) q[2];
sx q[2];
rz(-2.7375431) q[2];
rz(-0.36014253) q[3];
sx q[3];
rz(-1.6481684) q[3];
sx q[3];
rz(-2.4586316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3563511) q[0];
sx q[0];
rz(-2.1389565) q[0];
sx q[0];
rz(0.32401618) q[0];
rz(2.0815381) q[1];
sx q[1];
rz(-0.29393229) q[1];
sx q[1];
rz(0.70045984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7054847) q[0];
sx q[0];
rz(-1.7388845) q[0];
sx q[0];
rz(-1.4072493) q[0];
rz(-pi) q[1];
rz(-1.0251787) q[2];
sx q[2];
rz(-2.7433925) q[2];
sx q[2];
rz(0.37665483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30875242) q[1];
sx q[1];
rz(-0.82505095) q[1];
sx q[1];
rz(1.8925335) q[1];
x q[2];
rz(1.992393) q[3];
sx q[3];
rz(-1.3918878) q[3];
sx q[3];
rz(2.1960771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6274274) q[2];
sx q[2];
rz(-2.8943987) q[2];
sx q[2];
rz(0.79666454) q[2];
rz(1.1967777) q[3];
sx q[3];
rz(-1.5408885) q[3];
sx q[3];
rz(-0.59320199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6881123) q[0];
sx q[0];
rz(-1.0013094) q[0];
sx q[0];
rz(1.8185115) q[0];
rz(0.46609136) q[1];
sx q[1];
rz(-2.2345462) q[1];
sx q[1];
rz(1.1493433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1150165) q[0];
sx q[0];
rz(-1.1372024) q[0];
sx q[0];
rz(-2.114891) q[0];
x q[1];
rz(1.0428814) q[2];
sx q[2];
rz(-2.4529874) q[2];
sx q[2];
rz(-1.7797888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9117078) q[1];
sx q[1];
rz(-2.9344892) q[1];
sx q[1];
rz(2.9207499) q[1];
rz(-pi) q[2];
rz(-1.9051311) q[3];
sx q[3];
rz(-1.4282246) q[3];
sx q[3];
rz(2.4041374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.484802) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(-2.316324) q[2];
rz(1.2800062) q[3];
sx q[3];
rz(-1.620159) q[3];
sx q[3];
rz(-2.4204204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.698302) q[0];
sx q[0];
rz(-1.2606786) q[0];
sx q[0];
rz(1.3662421) q[0];
rz(1.0900991) q[1];
sx q[1];
rz(-1.5049728) q[1];
sx q[1];
rz(-1.3390138) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3109717) q[0];
sx q[0];
rz(-1.3960103) q[0];
sx q[0];
rz(2.4808148) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5191804) q[2];
sx q[2];
rz(-0.96660766) q[2];
sx q[2];
rz(-2.6780918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8112534) q[1];
sx q[1];
rz(-0.24345466) q[1];
sx q[1];
rz(-2.6324243) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3424266) q[3];
sx q[3];
rz(-0.46119565) q[3];
sx q[3];
rz(0.14652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1947386) q[2];
sx q[2];
rz(-1.7297435) q[2];
sx q[2];
rz(2.8080688) q[2];
rz(-2.5528095) q[3];
sx q[3];
rz(-2.6684561) q[3];
sx q[3];
rz(-2.3401882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81353417) q[0];
sx q[0];
rz(-1.3575587) q[0];
sx q[0];
rz(-1.8406156) q[0];
rz(-1.0282372) q[1];
sx q[1];
rz(-0.49003777) q[1];
sx q[1];
rz(0.44233826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8632354) q[0];
sx q[0];
rz(-1.8148737) q[0];
sx q[0];
rz(-0.4841448) q[0];
x q[1];
rz(1.0391159) q[2];
sx q[2];
rz(-1.9900512) q[2];
sx q[2];
rz(-0.13612572) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84917268) q[1];
sx q[1];
rz(-1.1467117) q[1];
sx q[1];
rz(1.7522274) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30486561) q[3];
sx q[3];
rz(-2.8072522) q[3];
sx q[3];
rz(1.5078733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77171317) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(-2.6698574) q[2];
rz(1.4419904) q[3];
sx q[3];
rz(-0.56157464) q[3];
sx q[3];
rz(0.26427463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1137375) q[0];
sx q[0];
rz(-1.4677445) q[0];
sx q[0];
rz(-0.16673985) q[0];
rz(-2.3475504) q[1];
sx q[1];
rz(-1.9414976) q[1];
sx q[1];
rz(0.099743191) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8309495) q[0];
sx q[0];
rz(-1.2590053) q[0];
sx q[0];
rz(-0.98592178) q[0];
rz(-2.0774309) q[2];
sx q[2];
rz(-2.4598413) q[2];
sx q[2];
rz(-0.1296986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9236026) q[1];
sx q[1];
rz(-2.1384708) q[1];
sx q[1];
rz(-1.4211378) q[1];
x q[2];
rz(-0.90899701) q[3];
sx q[3];
rz(-1.2332877) q[3];
sx q[3];
rz(0.26503978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2793067) q[2];
sx q[2];
rz(-1.0655094) q[2];
sx q[2];
rz(-0.12969895) q[2];
rz(-1.8779523) q[3];
sx q[3];
rz(-1.2010776) q[3];
sx q[3];
rz(-0.14550801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.170914) q[0];
sx q[0];
rz(-0.92249528) q[0];
sx q[0];
rz(-2.7150735) q[0];
rz(-2.049394) q[1];
sx q[1];
rz(-1.1323606) q[1];
sx q[1];
rz(-1.0027286) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195983) q[0];
sx q[0];
rz(-1.8289803) q[0];
sx q[0];
rz(0.95457558) q[0];
rz(-pi) q[1];
rz(-3.0363704) q[2];
sx q[2];
rz(-2.4159523) q[2];
sx q[2];
rz(1.9228455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1341237) q[1];
sx q[1];
rz(-1.136131) q[1];
sx q[1];
rz(2.77209) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0277983) q[3];
sx q[3];
rz(-2.6425231) q[3];
sx q[3];
rz(2.6036865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87454522) q[2];
sx q[2];
rz(-0.55570498) q[2];
sx q[2];
rz(0.12346539) q[2];
rz(0.63001436) q[3];
sx q[3];
rz(-1.8450626) q[3];
sx q[3];
rz(-2.4251078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89710871) q[0];
sx q[0];
rz(-0.69863313) q[0];
sx q[0];
rz(-0.83936349) q[0];
rz(-2.9821303) q[1];
sx q[1];
rz(-0.99226743) q[1];
sx q[1];
rz(-0.92963591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146067) q[0];
sx q[0];
rz(-0.85171284) q[0];
sx q[0];
rz(1.2308464) q[0];
x q[1];
rz(1.4484181) q[2];
sx q[2];
rz(-0.98320962) q[2];
sx q[2];
rz(1.1002968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1222555) q[1];
sx q[1];
rz(-2.4474065) q[1];
sx q[1];
rz(-1.7671142) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2655018) q[3];
sx q[3];
rz(-2.4809894) q[3];
sx q[3];
rz(-1.5098287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1975501) q[2];
sx q[2];
rz(-1.1669179) q[2];
sx q[2];
rz(0.76773947) q[2];
rz(2.7770216) q[3];
sx q[3];
rz(-1.3837827) q[3];
sx q[3];
rz(3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1345054) q[0];
sx q[0];
rz(-2.5275079) q[0];
sx q[0];
rz(-0.8836723) q[0];
rz(2.9332352) q[1];
sx q[1];
rz(-2.6388984) q[1];
sx q[1];
rz(1.3349894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55157298) q[0];
sx q[0];
rz(-3.0544222) q[0];
sx q[0];
rz(-2.4159191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6140974) q[2];
sx q[2];
rz(-1.732382) q[2];
sx q[2];
rz(1.8277663) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22182759) q[1];
sx q[1];
rz(-2.1073256) q[1];
sx q[1];
rz(-0.51862688) q[1];
x q[2];
rz(0.6649762) q[3];
sx q[3];
rz(-2.1998365) q[3];
sx q[3];
rz(1.4339989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9779382) q[2];
sx q[2];
rz(-1.5074707) q[2];
sx q[2];
rz(-0.046023544) q[2];
rz(-2.9764791) q[3];
sx q[3];
rz(-2.8486227) q[3];
sx q[3];
rz(1.2142115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5788427) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(2.9764755) q[1];
sx q[1];
rz(-1.4585635) q[1];
sx q[1];
rz(-0.30452902) q[1];
rz(-1.8971741) q[2];
sx q[2];
rz(-1.1481297) q[2];
sx q[2];
rz(1.4884244) q[2];
rz(-2.3143309) q[3];
sx q[3];
rz(-2.485459) q[3];
sx q[3];
rz(2.9970959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
