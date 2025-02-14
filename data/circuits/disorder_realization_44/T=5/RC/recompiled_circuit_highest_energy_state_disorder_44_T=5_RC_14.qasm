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
rz(2.9889838) q[0];
sx q[0];
rz(-1.9419365) q[0];
sx q[0];
rz(-2.8359523) q[0];
rz(-2.8513554) q[1];
sx q[1];
rz(-1.4479535) q[1];
sx q[1];
rz(-1.0040959) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1342524) q[0];
sx q[0];
rz(-1.5135845) q[0];
sx q[0];
rz(1.5438616) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3216686) q[2];
sx q[2];
rz(-1.3835356) q[2];
sx q[2];
rz(1.934777) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2724077) q[1];
sx q[1];
rz(-1.1704374) q[1];
sx q[1];
rz(0.34779208) q[1];
rz(-pi) q[2];
rz(0.53880764) q[3];
sx q[3];
rz(-1.5344005) q[3];
sx q[3];
rz(2.7766492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4940138) q[2];
sx q[2];
rz(-2.5250489) q[2];
sx q[2];
rz(-0.54440633) q[2];
rz(-2.7529649) q[3];
sx q[3];
rz(-1.2231239) q[3];
sx q[3];
rz(2.2012034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87043864) q[0];
sx q[0];
rz(-1.0259904) q[0];
sx q[0];
rz(2.0197268) q[0];
rz(2.0815966) q[1];
sx q[1];
rz(-0.50299877) q[1];
sx q[1];
rz(0.88567919) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6933357) q[0];
sx q[0];
rz(-1.8739032) q[0];
sx q[0];
rz(-1.8799549) q[0];
x q[1];
rz(0.76705025) q[2];
sx q[2];
rz(-0.97958744) q[2];
sx q[2];
rz(-1.5754099) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0421669) q[1];
sx q[1];
rz(-1.6337218) q[1];
sx q[1];
rz(-1.4896449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84195413) q[3];
sx q[3];
rz(-0.15781395) q[3];
sx q[3];
rz(2.006881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9649967) q[2];
sx q[2];
rz(-1.4623888) q[2];
sx q[2];
rz(-3.1122567) q[2];
rz(-2.7756179) q[3];
sx q[3];
rz(-2.6656606) q[3];
sx q[3];
rz(0.16987814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.733424) q[0];
sx q[0];
rz(-1.3329196) q[0];
sx q[0];
rz(-0.15342203) q[0];
rz(2.9452005) q[1];
sx q[1];
rz(-1.4686262) q[1];
sx q[1];
rz(0.82082716) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0401413) q[0];
sx q[0];
rz(-0.87803221) q[0];
sx q[0];
rz(0.67130868) q[0];
rz(-pi) q[1];
rz(0.94541855) q[2];
sx q[2];
rz(-2.4074005) q[2];
sx q[2];
rz(1.1933094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.046593) q[1];
sx q[1];
rz(-0.41727704) q[1];
sx q[1];
rz(-0.098578171) q[1];
x q[2];
rz(0.12755022) q[3];
sx q[3];
rz(-0.3693499) q[3];
sx q[3];
rz(2.7238059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11744943) q[2];
sx q[2];
rz(-1.8154181) q[2];
sx q[2];
rz(-0.93309012) q[2];
rz(0.34705958) q[3];
sx q[3];
rz(-2.0324028) q[3];
sx q[3];
rz(0.6412653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1953122) q[0];
sx q[0];
rz(-1.1271789) q[0];
sx q[0];
rz(2.0830925) q[0];
rz(0.095254101) q[1];
sx q[1];
rz(-1.9922099) q[1];
sx q[1];
rz(0.20406318) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6263555) q[0];
sx q[0];
rz(-1.7254819) q[0];
sx q[0];
rz(-1.3503287) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80486678) q[2];
sx q[2];
rz(-2.521732) q[2];
sx q[2];
rz(1.0484753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5285572) q[1];
sx q[1];
rz(-1.7712466) q[1];
sx q[1];
rz(2.5557808) q[1];
x q[2];
rz(1.3647365) q[3];
sx q[3];
rz(-2.2584174) q[3];
sx q[3];
rz(1.4233518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0846587) q[2];
sx q[2];
rz(-1.6997507) q[2];
sx q[2];
rz(-1.320768) q[2];
rz(-1.002958) q[3];
sx q[3];
rz(-1.2489677) q[3];
sx q[3];
rz(-2.5761719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.716662) q[0];
sx q[0];
rz(-1.2053763) q[0];
sx q[0];
rz(2.5004814) q[0];
rz(-2.5504316) q[1];
sx q[1];
rz(-1.4407651) q[1];
sx q[1];
rz(1.444918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49316367) q[0];
sx q[0];
rz(-0.38721353) q[0];
sx q[0];
rz(-3.0680543) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9282753) q[2];
sx q[2];
rz(-0.89175311) q[2];
sx q[2];
rz(-2.9625593) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7406297) q[1];
sx q[1];
rz(-2.0409806) q[1];
sx q[1];
rz(1.2931254) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0672449) q[3];
sx q[3];
rz(-2.5592521) q[3];
sx q[3];
rz(-2.989188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5833907) q[2];
sx q[2];
rz(-2.148197) q[2];
sx q[2];
rz(-0.85303419) q[2];
rz(2.9663626) q[3];
sx q[3];
rz(-1.8195189) q[3];
sx q[3];
rz(-0.62293735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22195062) q[0];
sx q[0];
rz(-2.3355244) q[0];
sx q[0];
rz(2.6958418) q[0];
rz(-0.99459612) q[1];
sx q[1];
rz(-1.0044731) q[1];
sx q[1];
rz(-1.4882784) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19916473) q[0];
sx q[0];
rz(-0.86948778) q[0];
sx q[0];
rz(-1.2851932) q[0];
x q[1];
rz(-3.0970552) q[2];
sx q[2];
rz(-1.6448529) q[2];
sx q[2];
rz(-3.0527405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46785746) q[1];
sx q[1];
rz(-0.36844353) q[1];
sx q[1];
rz(-0.96576565) q[1];
rz(-pi) q[2];
rz(3.0631243) q[3];
sx q[3];
rz(-1.244831) q[3];
sx q[3];
rz(-0.36074829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17066869) q[2];
sx q[2];
rz(-1.2430151) q[2];
sx q[2];
rz(2.1155913) q[2];
rz(2.3435727) q[3];
sx q[3];
rz(-0.84037104) q[3];
sx q[3];
rz(-2.87288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949424) q[0];
sx q[0];
rz(-1.3947399) q[0];
sx q[0];
rz(-2.7251439) q[0];
rz(1.7448447) q[1];
sx q[1];
rz(-1.5301219) q[1];
sx q[1];
rz(-1.6646615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48775169) q[0];
sx q[0];
rz(-1.6076822) q[0];
sx q[0];
rz(-1.0168309) q[0];
x q[1];
rz(0.32768048) q[2];
sx q[2];
rz(-1.9679356) q[2];
sx q[2];
rz(-2.5458592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.021183658) q[1];
sx q[1];
rz(-0.87090809) q[1];
sx q[1];
rz(-2.1645249) q[1];
x q[2];
rz(-2.9101679) q[3];
sx q[3];
rz(-0.34734472) q[3];
sx q[3];
rz(0.3350814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2485409) q[2];
sx q[2];
rz(-2.7613642) q[2];
sx q[2];
rz(0.60379544) q[2];
rz(-2.2327312) q[3];
sx q[3];
rz(-1.4330319) q[3];
sx q[3];
rz(1.0994256) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73868442) q[0];
sx q[0];
rz(-1.383902) q[0];
sx q[0];
rz(-2.4201194) q[0];
rz(-1.8062704) q[1];
sx q[1];
rz(-1.1386917) q[1];
sx q[1];
rz(-1.8849751) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5929991) q[0];
sx q[0];
rz(-1.8701435) q[0];
sx q[0];
rz(0.74798287) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0851791) q[2];
sx q[2];
rz(-2.0943421) q[2];
sx q[2];
rz(2.8741037) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80714204) q[1];
sx q[1];
rz(-2.6605407) q[1];
sx q[1];
rz(1.0472286) q[1];
x q[2];
rz(0.04560915) q[3];
sx q[3];
rz(-1.541434) q[3];
sx q[3];
rz(1.2330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98715544) q[2];
sx q[2];
rz(-2.4989231) q[2];
sx q[2];
rz(1.2745693) q[2];
rz(-0.68862033) q[3];
sx q[3];
rz(-1.6504811) q[3];
sx q[3];
rz(0.0037746599) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1275198) q[0];
sx q[0];
rz(-1.8919683) q[0];
sx q[0];
rz(2.4697812) q[0];
rz(1.6067243) q[1];
sx q[1];
rz(-0.22767362) q[1];
sx q[1];
rz(0.52275503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5005258) q[0];
sx q[0];
rz(-2.4198774) q[0];
sx q[0];
rz(-0.39875491) q[0];
rz(-0.050809697) q[2];
sx q[2];
rz(-1.9138971) q[2];
sx q[2];
rz(-1.9947413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.213428) q[1];
sx q[1];
rz(-1.1901179) q[1];
sx q[1];
rz(0.45427096) q[1];
rz(2.5795311) q[3];
sx q[3];
rz(-2.76427) q[3];
sx q[3];
rz(-2.5320208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25040024) q[2];
sx q[2];
rz(-0.6468536) q[2];
sx q[2];
rz(-2.7817173) q[2];
rz(1.322809) q[3];
sx q[3];
rz(-1.6836932) q[3];
sx q[3];
rz(-0.047164269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5322402) q[0];
sx q[0];
rz(-2.1873964) q[0];
sx q[0];
rz(0.94616079) q[0];
rz(-0.040945176) q[1];
sx q[1];
rz(-1.8649201) q[1];
sx q[1];
rz(-1.2650222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32308043) q[0];
sx q[0];
rz(-1.6744594) q[0];
sx q[0];
rz(-0.47804272) q[0];
x q[1];
rz(0.47750116) q[2];
sx q[2];
rz(-2.6386119) q[2];
sx q[2];
rz(-1.2625842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62502669) q[1];
sx q[1];
rz(-1.0426765) q[1];
sx q[1];
rz(-2.3305446) q[1];
x q[2];
rz(0.048790292) q[3];
sx q[3];
rz(-1.1975761) q[3];
sx q[3];
rz(-0.1540972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9762743) q[2];
sx q[2];
rz(-1.9860622) q[2];
sx q[2];
rz(-2.9602642) q[2];
rz(-0.85339648) q[3];
sx q[3];
rz(-2.3386164) q[3];
sx q[3];
rz(1.3388504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778018) q[0];
sx q[0];
rz(-1.3297357) q[0];
sx q[0];
rz(1.7672675) q[0];
rz(1.455066) q[1];
sx q[1];
rz(-1.4604026) q[1];
sx q[1];
rz(-0.30053465) q[1];
rz(0.65291885) q[2];
sx q[2];
rz(-0.75426415) q[2];
sx q[2];
rz(1.33367) q[2];
rz(2.4134514) q[3];
sx q[3];
rz(-2.7645088) q[3];
sx q[3];
rz(-1.0704452) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
