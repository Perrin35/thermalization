OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8074789) q[0];
sx q[0];
rz(-0.69598323) q[0];
sx q[0];
rz(-0.88710436) q[0];
rz(-0.71169418) q[1];
sx q[1];
rz(-1.0625755) q[1];
sx q[1];
rz(-0.33023155) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7831551) q[0];
sx q[0];
rz(-2.3312954) q[0];
sx q[0];
rz(0.94457074) q[0];
rz(2.5862972) q[2];
sx q[2];
rz(-1.2568297) q[2];
sx q[2];
rz(1.2295251) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6373693) q[1];
sx q[1];
rz(-2.0614701) q[1];
sx q[1];
rz(1.4549903) q[1];
rz(-1.4012778) q[3];
sx q[3];
rz(-2.7054686) q[3];
sx q[3];
rz(0.90986246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3868788) q[2];
sx q[2];
rz(-0.42749307) q[2];
sx q[2];
rz(-1.0270366) q[2];
rz(-2.2586281) q[3];
sx q[3];
rz(-1.797978) q[3];
sx q[3];
rz(0.55356717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.7850194) q[0];
sx q[0];
rz(-0.078332575) q[0];
sx q[0];
rz(0.26764348) q[0];
rz(1.74125) q[1];
sx q[1];
rz(-0.61336556) q[1];
sx q[1];
rz(-2.5852481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436377) q[0];
sx q[0];
rz(-2.58444) q[0];
sx q[0];
rz(2.055026) q[0];
rz(-pi) q[1];
rz(-2.7423233) q[2];
sx q[2];
rz(-0.74490163) q[2];
sx q[2];
rz(-2.5585554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2399064) q[1];
sx q[1];
rz(-1.0178029) q[1];
sx q[1];
rz(-2.6917975) q[1];
rz(2.56145) q[3];
sx q[3];
rz(-2.1506967) q[3];
sx q[3];
rz(0.69441352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23655015) q[2];
sx q[2];
rz(-1.6908129) q[2];
sx q[2];
rz(-1.2518008) q[2];
rz(-1.8768138) q[3];
sx q[3];
rz(-2.4301961) q[3];
sx q[3];
rz(-1.0714162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.4816137) q[0];
sx q[0];
rz(-2.5857506) q[0];
sx q[0];
rz(-0.83928338) q[0];
rz(-1.6370157) q[1];
sx q[1];
rz(-1.4312276) q[1];
sx q[1];
rz(-1.379871) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312442) q[0];
sx q[0];
rz(-1.4891242) q[0];
sx q[0];
rz(-1.6947338) q[0];
x q[1];
rz(-2.9271651) q[2];
sx q[2];
rz(-2.1058208) q[2];
sx q[2];
rz(2.9012539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31095895) q[1];
sx q[1];
rz(-1.1006025) q[1];
sx q[1];
rz(-0.36010919) q[1];
rz(-pi) q[2];
rz(2.1228509) q[3];
sx q[3];
rz(-2.0716487) q[3];
sx q[3];
rz(-0.93397442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.69924361) q[2];
sx q[2];
rz(-1.2171429) q[2];
sx q[2];
rz(-2.4231363) q[2];
rz(-2.5213304) q[3];
sx q[3];
rz(-0.93554997) q[3];
sx q[3];
rz(-0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45850596) q[0];
sx q[0];
rz(-1.9106671) q[0];
sx q[0];
rz(-1.7339535) q[0];
rz(0.55157026) q[1];
sx q[1];
rz(-3.0307814) q[1];
sx q[1];
rz(1.8254021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5410446) q[0];
sx q[0];
rz(-1.2886314) q[0];
sx q[0];
rz(0.94775424) q[0];
rz(-pi) q[1];
rz(3.1291601) q[2];
sx q[2];
rz(-1.2913449) q[2];
sx q[2];
rz(0.81571992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1626959) q[1];
sx q[1];
rz(-2.7172019) q[1];
sx q[1];
rz(-0.32482206) q[1];
x q[2];
rz(0.67263453) q[3];
sx q[3];
rz(-2.1676873) q[3];
sx q[3];
rz(-1.2850645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59336415) q[2];
sx q[2];
rz(-2.7061988) q[2];
sx q[2];
rz(-1.0080053) q[2];
rz(-0.27901444) q[3];
sx q[3];
rz(-0.81441003) q[3];
sx q[3];
rz(0.45849714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27610436) q[0];
sx q[0];
rz(-1.3300329) q[0];
sx q[0];
rz(-0.27648595) q[0];
rz(2.8522988) q[1];
sx q[1];
rz(-1.1951059) q[1];
sx q[1];
rz(-0.2624661) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55928265) q[0];
sx q[0];
rz(-1.773343) q[0];
sx q[0];
rz(1.7515394) q[0];
rz(-pi) q[1];
rz(0.010774537) q[2];
sx q[2];
rz(-0.97948217) q[2];
sx q[2];
rz(-0.39515218) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20690675) q[1];
sx q[1];
rz(-2.1175724) q[1];
sx q[1];
rz(1.9711167) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9093303) q[3];
sx q[3];
rz(-2.5506936) q[3];
sx q[3];
rz(0.10400203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7485973) q[2];
sx q[2];
rz(-2.5454919) q[2];
sx q[2];
rz(-1.3470915) q[2];
rz(0.10284452) q[3];
sx q[3];
rz(-1.9072073) q[3];
sx q[3];
rz(-1.6245406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1384001) q[0];
sx q[0];
rz(-1.5035368) q[0];
sx q[0];
rz(0.060977161) q[0];
rz(-1.2334476) q[1];
sx q[1];
rz(-1.7465218) q[1];
sx q[1];
rz(1.5256418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4075482) q[0];
sx q[0];
rz(-1.6400695) q[0];
sx q[0];
rz(2.6114527) q[0];
rz(-pi) q[1];
rz(1.5360918) q[2];
sx q[2];
rz(-2.3922709) q[2];
sx q[2];
rz(2.9108436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96887302) q[1];
sx q[1];
rz(-0.27471009) q[1];
sx q[1];
rz(1.9897047) q[1];
rz(1.3539702) q[3];
sx q[3];
rz(-2.4783756) q[3];
sx q[3];
rz(1.7182223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47152758) q[2];
sx q[2];
rz(-1.7933041) q[2];
sx q[2];
rz(1.8111551) q[2];
rz(2.0165675) q[3];
sx q[3];
rz(-2.2671813) q[3];
sx q[3];
rz(0.037954656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838945) q[0];
sx q[0];
rz(-1.6989468) q[0];
sx q[0];
rz(-2.884602) q[0];
rz(-0.43486241) q[1];
sx q[1];
rz(-0.75006524) q[1];
sx q[1];
rz(0.90352568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57243915) q[0];
sx q[0];
rz(-1.4785197) q[0];
sx q[0];
rz(3.1374428) q[0];
rz(1.4353739) q[2];
sx q[2];
rz(-0.82937448) q[2];
sx q[2];
rz(-3.137085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7419646) q[1];
sx q[1];
rz(-3.0638712) q[1];
sx q[1];
rz(-0.6570973) q[1];
rz(-2.2478476) q[3];
sx q[3];
rz(-2.4747958) q[3];
sx q[3];
rz(2.9725985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32288185) q[2];
sx q[2];
rz(-1.633753) q[2];
sx q[2];
rz(-1.7203077) q[2];
rz(-0.70982248) q[3];
sx q[3];
rz(-1.1387419) q[3];
sx q[3];
rz(0.53409725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1555971) q[0];
sx q[0];
rz(-0.43102145) q[0];
sx q[0];
rz(1.0850329) q[0];
rz(2.494508) q[1];
sx q[1];
rz(-1.9661463) q[1];
sx q[1];
rz(3.0771902) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2213106) q[0];
sx q[0];
rz(-2.1951809) q[0];
sx q[0];
rz(1.1657752) q[0];
rz(-2.4738381) q[2];
sx q[2];
rz(-2.1760094) q[2];
sx q[2];
rz(-0.79082327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6820087) q[1];
sx q[1];
rz(-1.1136076) q[1];
sx q[1];
rz(2.3383635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9265979) q[3];
sx q[3];
rz(-1.8321773) q[3];
sx q[3];
rz(1.3047117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.881968) q[2];
sx q[2];
rz(-2.7590064) q[2];
sx q[2];
rz(3.0982223) q[2];
rz(-2.4671593) q[3];
sx q[3];
rz(-1.3573656) q[3];
sx q[3];
rz(-1.1446713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6944273) q[0];
sx q[0];
rz(-2.5612216) q[0];
sx q[0];
rz(-2.8149783) q[0];
rz(-2.0475552) q[1];
sx q[1];
rz(-0.84427619) q[1];
sx q[1];
rz(0.48386595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1492067) q[0];
sx q[0];
rz(-2.0786374) q[0];
sx q[0];
rz(0.81224982) q[0];
rz(2.6695374) q[2];
sx q[2];
rz(-0.46923551) q[2];
sx q[2];
rz(0.068048565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30712515) q[1];
sx q[1];
rz(-0.79453429) q[1];
sx q[1];
rz(-1.267872) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3392162) q[3];
sx q[3];
rz(-2.2001523) q[3];
sx q[3];
rz(1.6762102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6924374) q[2];
sx q[2];
rz(-2.1026976) q[2];
sx q[2];
rz(3.0628487) q[2];
rz(-0.76830831) q[3];
sx q[3];
rz(-1.2945622) q[3];
sx q[3];
rz(-2.9505742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67097265) q[0];
sx q[0];
rz(-2.5363531) q[0];
sx q[0];
rz(-3.0249) q[0];
rz(-0.8600046) q[1];
sx q[1];
rz(-1.3415033) q[1];
sx q[1];
rz(0.87342993) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610342) q[0];
sx q[0];
rz(-1.8115028) q[0];
sx q[0];
rz(0.81674256) q[0];
rz(-pi) q[1];
rz(2.6038613) q[2];
sx q[2];
rz(-0.5774379) q[2];
sx q[2];
rz(2.2619907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5705986) q[1];
sx q[1];
rz(-2.6254404) q[1];
sx q[1];
rz(-2.8883449) q[1];
rz(-1.5029346) q[3];
sx q[3];
rz(-0.96047365) q[3];
sx q[3];
rz(0.44059702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42150911) q[2];
sx q[2];
rz(-2.8953711) q[2];
sx q[2];
rz(-2.1923547) q[2];
rz(1.2468437) q[3];
sx q[3];
rz(-1.4693762) q[3];
sx q[3];
rz(-2.5505572) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0134037) q[0];
sx q[0];
rz(-2.1407776) q[0];
sx q[0];
rz(1.5110973) q[0];
rz(2.7746101) q[1];
sx q[1];
rz(-1.783168) q[1];
sx q[1];
rz(2.5617243) q[1];
rz(-2.1268216) q[2];
sx q[2];
rz(-2.6064739) q[2];
sx q[2];
rz(0.50630416) q[2];
rz(-1.1506769) q[3];
sx q[3];
rz(-0.92460604) q[3];
sx q[3];
rz(-2.9210319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
