OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(7.6927778) q[0];
sx q[0];
rz(11.132244) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.943676) q[0];
sx q[0];
rz(-2.0352053) q[0];
sx q[0];
rz(2.8465413) q[0];
rz(-pi) q[1];
rz(1.7069874) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.9905123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0093706) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(-2.0396712) q[1];
rz(-2.9207346) q[3];
sx q[3];
rz(-2.0005895) q[3];
sx q[3];
rz(-2.5151099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(0.49301246) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(-3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74719602) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(-1.2778506) q[0];
rz(-0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(0.4321672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2910599) q[0];
sx q[0];
rz(-1.8467554) q[0];
sx q[0];
rz(0.54131298) q[0];
x q[1];
rz(-0.82168174) q[2];
sx q[2];
rz(-2.5857918) q[2];
sx q[2];
rz(2.1205714) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13465263) q[1];
sx q[1];
rz(-1.68774) q[1];
sx q[1];
rz(-1.1275396) q[1];
rz(-pi) q[2];
x q[2];
rz(0.024859419) q[3];
sx q[3];
rz(-2.1356574) q[3];
sx q[3];
rz(2.3853175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(-2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040722672) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(2.7242463) q[0];
rz(-1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(0.506385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2877809) q[0];
sx q[0];
rz(-1.7475442) q[0];
sx q[0];
rz(1.214083) q[0];
rz(-0.26489139) q[2];
sx q[2];
rz(-0.75220097) q[2];
sx q[2];
rz(-1.0571935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3370034) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(-0.6188436) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9403463) q[3];
sx q[3];
rz(-1.4003716) q[3];
sx q[3];
rz(0.91526645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69904077) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(2.3024094) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5485839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2617333) q[0];
sx q[0];
rz(-1.121959) q[0];
sx q[0];
rz(0.10106048) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4315005) q[2];
sx q[2];
rz(-2.2033764) q[2];
sx q[2];
rz(1.3131504) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81739391) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(2.6170931) q[1];
x q[2];
rz(2.7110093) q[3];
sx q[3];
rz(-1.8292556) q[3];
sx q[3];
rz(-2.5653198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(-2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(0.25594041) q[0];
rz(2.6804965) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(-0.76006132) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1003935) q[0];
sx q[0];
rz(-1.6738335) q[0];
sx q[0];
rz(-0.11739199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.063644479) q[2];
sx q[2];
rz(-1.9846989) q[2];
sx q[2];
rz(0.49583437) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45436817) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(1.441799) q[1];
rz(-2.2171668) q[3];
sx q[3];
rz(-1.8674208) q[3];
sx q[3];
rz(-0.38277205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(0.6742397) q[2];
rz(-2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-0.017344346) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.1791139) q[0];
rz(2.9367661) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(1.0669605) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63283352) q[0];
sx q[0];
rz(-1.5501223) q[0];
sx q[0];
rz(1.585292) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1322137) q[2];
sx q[2];
rz(-1.1819934) q[2];
sx q[2];
rz(0.68225451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8384335) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(0.82364239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9818929) q[3];
sx q[3];
rz(-0.80280639) q[3];
sx q[3];
rz(1.476895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.431488) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.7204826) q[0];
rz(2.9395318) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(2.2834159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1976397) q[0];
sx q[0];
rz(-1.7174935) q[0];
sx q[0];
rz(-0.12978817) q[0];
rz(-pi) q[1];
x q[1];
rz(0.045072149) q[2];
sx q[2];
rz(-1.1606693) q[2];
sx q[2];
rz(-0.28442597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.757829) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(1.5555192) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5704727) q[3];
sx q[3];
rz(-1.2763599) q[3];
sx q[3];
rz(1.7366228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28785607) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.3254335) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.3495061) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9106306) q[0];
sx q[0];
rz(-2.1719143) q[0];
sx q[0];
rz(-2.0448951) q[0];
x q[1];
rz(0.71808727) q[2];
sx q[2];
rz(-1.5880843) q[2];
sx q[2];
rz(-1.7503439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0950349) q[1];
sx q[1];
rz(-1.3890424) q[1];
sx q[1];
rz(-3.0912116) q[1];
rz(-pi) q[2];
rz(0.41534822) q[3];
sx q[3];
rz(-0.58598622) q[3];
sx q[3];
rz(-1.2522445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(1.9474585) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(-0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.5375686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81048548) q[0];
sx q[0];
rz(-1.2665505) q[0];
sx q[0];
rz(2.9097793) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6749009) q[2];
sx q[2];
rz(-1.9874007) q[2];
sx q[2];
rz(-0.57550752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2214339) q[1];
sx q[1];
rz(-0.71343525) q[1];
sx q[1];
rz(-0.17381298) q[1];
rz(2.0314625) q[3];
sx q[3];
rz(-1.8970282) q[3];
sx q[3];
rz(-0.5512475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1016772) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(0.90041655) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-2.9246869) q[0];
rz(0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(-2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2749598) q[0];
sx q[0];
rz(-1.8909591) q[0];
sx q[0];
rz(0.76036705) q[0];
rz(-0.54951349) q[2];
sx q[2];
rz(-1.9351442) q[2];
sx q[2];
rz(3.0669616) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.875647) q[1];
sx q[1];
rz(-1.4946283) q[1];
sx q[1];
rz(0.45653685) q[1];
x q[2];
rz(-2.1438164) q[3];
sx q[3];
rz(-0.7191092) q[3];
sx q[3];
rz(2.5620808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(-2.1968502) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(-2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.5007297) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-2.5169218) q[2];
sx q[2];
rz(-2.2041337) q[2];
sx q[2];
rz(0.49908257) q[2];
rz(-1.1576204) q[3];
sx q[3];
rz(-0.4784085) q[3];
sx q[3];
rz(2.9984409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
