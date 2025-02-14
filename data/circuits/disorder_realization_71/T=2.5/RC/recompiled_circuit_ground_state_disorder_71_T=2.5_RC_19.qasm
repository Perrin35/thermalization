OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3040721) q[0];
sx q[0];
rz(3.8068258) q[0];
sx q[0];
rz(8.1991631) q[0];
rz(2.4576814) q[1];
sx q[1];
rz(-0.39165762) q[1];
sx q[1];
rz(0.70206967) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0618781) q[0];
sx q[0];
rz(-1.4532545) q[0];
sx q[0];
rz(-2.4410309) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4058787) q[2];
sx q[2];
rz(-1.5497297) q[2];
sx q[2];
rz(-1.5738459) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7309642) q[1];
sx q[1];
rz(-2.0515039) q[1];
sx q[1];
rz(2.5334873) q[1];
x q[2];
rz(1.8075712) q[3];
sx q[3];
rz(-0.70565685) q[3];
sx q[3];
rz(1.9183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94032732) q[2];
sx q[2];
rz(-1.5387115) q[2];
sx q[2];
rz(-2.8931457) q[2];
rz(-1.1343608) q[3];
sx q[3];
rz(-0.13392197) q[3];
sx q[3];
rz(-2.0094481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(3.0487173) q[0];
sx q[0];
rz(-0.55260125) q[0];
sx q[0];
rz(1.5485113) q[0];
rz(0.50367194) q[1];
sx q[1];
rz(-1.9182938) q[1];
sx q[1];
rz(-2.7647387) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5708964) q[0];
sx q[0];
rz(-2.6691169) q[0];
sx q[0];
rz(3.0993479) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69825577) q[2];
sx q[2];
rz(-0.94329903) q[2];
sx q[2];
rz(-2.1098532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7347757) q[1];
sx q[1];
rz(-1.3208913) q[1];
sx q[1];
rz(2.6491449) q[1];
x q[2];
rz(1.6474071) q[3];
sx q[3];
rz(-1.239068) q[3];
sx q[3];
rz(-2.2463617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5571931) q[2];
sx q[2];
rz(-2.445745) q[2];
sx q[2];
rz(-0.86563555) q[2];
rz(-1.4731167) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(-2.1540811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(0.57403785) q[0];
sx q[0];
rz(-1.8542629) q[0];
sx q[0];
rz(-1.0512742) q[0];
rz(1.6846664) q[1];
sx q[1];
rz(-1.3070062) q[1];
sx q[1];
rz(1.6251224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.595436) q[0];
sx q[0];
rz(-1.8519028) q[0];
sx q[0];
rz(0.77451046) q[0];
rz(-pi) q[1];
rz(-2.3651354) q[2];
sx q[2];
rz(-0.8469905) q[2];
sx q[2];
rz(-0.077082917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1976002) q[1];
sx q[1];
rz(-1.9514582) q[1];
sx q[1];
rz(2.073111) q[1];
x q[2];
rz(-0.061013075) q[3];
sx q[3];
rz(-1.5170238) q[3];
sx q[3];
rz(-2.6110569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39525825) q[2];
sx q[2];
rz(-2.0863775) q[2];
sx q[2];
rz(-2.9193817) q[2];
rz(2.639751) q[3];
sx q[3];
rz(-1.4072489) q[3];
sx q[3];
rz(0.80880222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489081) q[0];
sx q[0];
rz(-0.48766708) q[0];
sx q[0];
rz(1.9768313) q[0];
rz(-0.89272967) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(2.8597615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4455216) q[0];
sx q[0];
rz(-2.0508678) q[0];
sx q[0];
rz(-1.4351373) q[0];
x q[1];
rz(-0.71433432) q[2];
sx q[2];
rz(-2.0296944) q[2];
sx q[2];
rz(2.4752576) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52543312) q[1];
sx q[1];
rz(-1.6405917) q[1];
sx q[1];
rz(-1.4074123) q[1];
rz(-pi) q[2];
x q[2];
rz(1.124566) q[3];
sx q[3];
rz(-1.5939624) q[3];
sx q[3];
rz(1.5737826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0622327) q[2];
sx q[2];
rz(-2.5793109) q[2];
sx q[2];
rz(-2.6083561) q[2];
rz(2.0594635) q[3];
sx q[3];
rz(-2.5366668) q[3];
sx q[3];
rz(-0.81348872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37394062) q[0];
sx q[0];
rz(-2.7426608) q[0];
sx q[0];
rz(1.8178513) q[0];
rz(-0.81870493) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(2.0735819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477953) q[0];
sx q[0];
rz(-1.6500705) q[0];
sx q[0];
rz(-2.981841) q[0];
rz(0.85535149) q[2];
sx q[2];
rz(-0.87424874) q[2];
sx q[2];
rz(1.0698505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9449759) q[1];
sx q[1];
rz(-2.2240337) q[1];
sx q[1];
rz(1.7004299) q[1];
rz(1.8858484) q[3];
sx q[3];
rz(-0.49335155) q[3];
sx q[3];
rz(-0.82786614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9615122) q[2];
sx q[2];
rz(-2.027812) q[2];
sx q[2];
rz(0.86165825) q[2];
rz(-1.0263475) q[3];
sx q[3];
rz(-1.15871) q[3];
sx q[3];
rz(-2.0238743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79477972) q[0];
sx q[0];
rz(-2.3220799) q[0];
sx q[0];
rz(2.2487707) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-0.67829689) q[1];
sx q[1];
rz(0.13042626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3149726) q[0];
sx q[0];
rz(-2.9493647) q[0];
sx q[0];
rz(-0.97127025) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4685235) q[2];
sx q[2];
rz(-2.4166738) q[2];
sx q[2];
rz(-2.6672305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6857026) q[1];
sx q[1];
rz(-1.3336542) q[1];
sx q[1];
rz(-3.1115467) q[1];
rz(-pi) q[2];
rz(-2.942286) q[3];
sx q[3];
rz(-0.79108566) q[3];
sx q[3];
rz(-1.4996604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22455198) q[2];
sx q[2];
rz(-1.2265393) q[2];
sx q[2];
rz(2.2652333) q[2];
rz(-1.2419491) q[3];
sx q[3];
rz(-0.94049898) q[3];
sx q[3];
rz(-1.0103286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7414311) q[0];
sx q[0];
rz(-3.042996) q[0];
sx q[0];
rz(0.36369351) q[0];
rz(-2.3997276) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(-2.738764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2296933) q[0];
sx q[0];
rz(-0.48297627) q[0];
sx q[0];
rz(-2.1080024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9697519) q[2];
sx q[2];
rz(-1.2201933) q[2];
sx q[2];
rz(-0.95870852) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2022823) q[1];
sx q[1];
rz(-2.99761) q[1];
sx q[1];
rz(1.0054863) q[1];
rz(-pi) q[2];
rz(-0.36084036) q[3];
sx q[3];
rz(-2.4746102) q[3];
sx q[3];
rz(-2.9426394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3591298) q[2];
sx q[2];
rz(-0.368258) q[2];
sx q[2];
rz(-1.5472319) q[2];
rz(0.83526978) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(-0.48867759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816724) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(-0.51505995) q[0];
rz(1.7022279) q[1];
sx q[1];
rz(-1.2934338) q[1];
sx q[1];
rz(2.7424367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24267749) q[0];
sx q[0];
rz(-1.5689881) q[0];
sx q[0];
rz(-0.13092069) q[0];
rz(-pi) q[1];
rz(1.8410583) q[2];
sx q[2];
rz(-2.3944602) q[2];
sx q[2];
rz(3.0556553) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.040144939) q[1];
sx q[1];
rz(-1.8310837) q[1];
sx q[1];
rz(0.066802967) q[1];
x q[2];
rz(0.065654556) q[3];
sx q[3];
rz(-1.5720417) q[3];
sx q[3];
rz(2.7630382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9814375) q[2];
sx q[2];
rz(-1.8835386) q[2];
sx q[2];
rz(-2.0951648) q[2];
rz(-2.4753172) q[3];
sx q[3];
rz(-0.25686887) q[3];
sx q[3];
rz(2.090914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.94205725) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(-0.60607213) q[0];
rz(2.3145158) q[1];
sx q[1];
rz(-1.8111633) q[1];
sx q[1];
rz(-1.3333295) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2829958) q[0];
sx q[0];
rz(-1.9278316) q[0];
sx q[0];
rz(1.6813361) q[0];
rz(-1.1252046) q[2];
sx q[2];
rz(-2.2632709) q[2];
sx q[2];
rz(-2.8852579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3279475) q[1];
sx q[1];
rz(-2.4634117) q[1];
sx q[1];
rz(-1.3840883) q[1];
rz(0.17154947) q[3];
sx q[3];
rz(-2.6772873) q[3];
sx q[3];
rz(-2.7363079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0102319) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(-2.7566747) q[2];
rz(-0.63079232) q[3];
sx q[3];
rz(-1.8600978) q[3];
sx q[3];
rz(1.7990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.25403062) q[0];
sx q[0];
rz(-1.7429202) q[0];
sx q[0];
rz(1.4400462) q[0];
rz(0.74132672) q[1];
sx q[1];
rz(-1.7404375) q[1];
sx q[1];
rz(-0.25873605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1848768) q[0];
sx q[0];
rz(-1.3248845) q[0];
sx q[0];
rz(3.1374404) q[0];
rz(-pi) q[1];
rz(0.3157987) q[2];
sx q[2];
rz(-1.074203) q[2];
sx q[2];
rz(1.1441355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1867793) q[1];
sx q[1];
rz(-0.80866058) q[1];
sx q[1];
rz(1.0891799) q[1];
x q[2];
rz(-1.9680989) q[3];
sx q[3];
rz(-2.3081995) q[3];
sx q[3];
rz(-1.5381588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72186333) q[2];
sx q[2];
rz(-2.0294919) q[2];
sx q[2];
rz(-1.8224243) q[2];
rz(2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(0.54660249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.4010314) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(-1.7026547) q[1];
sx q[1];
rz(-0.66822744) q[1];
sx q[1];
rz(-2.6034036) q[1];
rz(1.5580675) q[2];
sx q[2];
rz(-2.4314778) q[2];
sx q[2];
rz(1.1927803) q[2];
rz(-1.9520252) q[3];
sx q[3];
rz(-2.5459131) q[3];
sx q[3];
rz(-0.66543647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
