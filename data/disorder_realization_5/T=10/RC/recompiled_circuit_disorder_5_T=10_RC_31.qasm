OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514483) q[0];
sx q[0];
rz(-1.7062618) q[0];
sx q[0];
rz(-2.9209903) q[0];
rz(-0.692042) q[2];
sx q[2];
rz(-3.0243052) q[2];
sx q[2];
rz(-0.28809822) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4133271) q[1];
sx q[1];
rz(-1.773355) q[1];
sx q[1];
rz(-0.061518016) q[1];
rz(-pi) q[2];
rz(-1.0702707) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4883603) q[0];
sx q[0];
rz(-0.45408861) q[0];
sx q[0];
rz(-2.1064208) q[0];
x q[1];
rz(-2.0741834) q[2];
sx q[2];
rz(-1.6620087) q[2];
sx q[2];
rz(-2.3561321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2903366) q[1];
sx q[1];
rz(-2.1888121) q[1];
sx q[1];
rz(3.0783155) q[1];
rz(3.0117399) q[3];
sx q[3];
rz(-2.2359071) q[3];
sx q[3];
rz(-0.065404281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(-2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547566) q[0];
sx q[0];
rz(-1.1491547) q[0];
sx q[0];
rz(0.7936759) q[0];
rz(-pi) q[1];
rz(-1.2013024) q[2];
sx q[2];
rz(-0.72417799) q[2];
sx q[2];
rz(-2.5232814) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(1.3819441) q[1];
rz(-pi) q[2];
rz(2.6684746) q[3];
sx q[3];
rz(-2.3948673) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.408067) q[2];
rz(0.18313289) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(2.4750211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4857805) q[0];
sx q[0];
rz(-0.27944767) q[0];
sx q[0];
rz(0.73894545) q[0];
rz(-1.2810983) q[2];
sx q[2];
rz(-2.2640267) q[2];
sx q[2];
rz(1.9300269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8909292) q[1];
sx q[1];
rz(-2.2729985) q[1];
sx q[1];
rz(2.3198747) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7188844) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8551222) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(0.83706013) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(-1.9821092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5975208) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(1.9835299) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(-1.3364524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3573208) q[1];
sx q[1];
rz(-1.9958152) q[1];
sx q[1];
rz(2.1678863) q[1];
rz(-pi) q[2];
rz(-0.63727832) q[3];
sx q[3];
rz(-1.1253469) q[3];
sx q[3];
rz(-0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7065457) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(1.2716028) q[2];
rz(1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5746675) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(3.040722) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6504165) q[0];
sx q[0];
rz(-2.4917779) q[0];
sx q[0];
rz(-0.25214904) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51668824) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(1.237243) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-0.15921758) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50562596) q[3];
sx q[3];
rz(-0.45590948) q[3];
sx q[3];
rz(2.0924007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(-2.4538453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6531658) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(1.2290918) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8887927) q[2];
sx q[2];
rz(-1.3557634) q[2];
sx q[2];
rz(0.064934864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3878165) q[1];
sx q[1];
rz(-0.13145914) q[1];
sx q[1];
rz(0.17863518) q[1];
x q[2];
rz(-2.1391159) q[3];
sx q[3];
rz(-1.5573313) q[3];
sx q[3];
rz(2.3426567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(1.942873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872902) q[0];
sx q[0];
rz(-2.0044921) q[0];
sx q[0];
rz(1.646423) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0555658) q[2];
sx q[2];
rz(-1.9326397) q[2];
sx q[2];
rz(1.8334243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3938013) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(1.6277114) q[1];
rz(-pi) q[2];
rz(1.8317322) q[3];
sx q[3];
rz(-2.0203777) q[3];
sx q[3];
rz(0.98228067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.566074) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(0.9334329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8342455) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(3.0349915) q[0];
rz(-2.2869296) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(1.9180627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6586106) q[1];
sx q[1];
rz(-0.88739363) q[1];
sx q[1];
rz(1.5233836) q[1];
rz(2.6870319) q[3];
sx q[3];
rz(-2.436736) q[3];
sx q[3];
rz(-0.2017894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(0.36994568) q[2];
rz(-0.89921078) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-0.95329154) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1870435) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(2.488234) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606124) q[0];
sx q[0];
rz(-1.6913337) q[0];
sx q[0];
rz(-1.6880077) q[0];
rz(0.75283639) q[2];
sx q[2];
rz(-0.92217731) q[2];
sx q[2];
rz(2.2019049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(0.88840719) q[1];
rz(-3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-0.72485483) q[2];
sx q[2];
rz(-2.5584494) q[2];
sx q[2];
rz(0.7128788) q[2];
rz(-1.8300874) q[3];
sx q[3];
rz(-2.5355831) q[3];
sx q[3];
rz(-2.2021962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
