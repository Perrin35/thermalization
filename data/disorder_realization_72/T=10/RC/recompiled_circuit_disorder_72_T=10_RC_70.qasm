OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(-0.49931061) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0190174) q[0];
sx q[0];
rz(-0.28261533) q[0];
sx q[0];
rz(-2.9119888) q[0];
x q[1];
rz(0.17272858) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(1.8634862) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0608622) q[1];
sx q[1];
rz(-2.3002491) q[1];
sx q[1];
rz(1.7727477) q[1];
rz(1.6742168) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(-0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(-0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(-1.0027592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28264499) q[0];
sx q[0];
rz(-0.98416057) q[0];
sx q[0];
rz(2.100201) q[0];
rz(-pi) q[1];
rz(-1.3833212) q[2];
sx q[2];
rz(-1.7406929) q[2];
sx q[2];
rz(-2.1263009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7748389) q[1];
sx q[1];
rz(-0.60791053) q[1];
sx q[1];
rz(-1.2971836) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6456804) q[3];
sx q[3];
rz(-0.65973385) q[3];
sx q[3];
rz(2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(2.3625968) q[0];
rz(-2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4602063) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(1.7494739) q[0];
rz(-1.1728889) q[2];
sx q[2];
rz(-1.298561) q[2];
sx q[2];
rz(2.9583601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1556342) q[1];
sx q[1];
rz(-1.2859207) q[1];
sx q[1];
rz(0.023521544) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65914764) q[3];
sx q[3];
rz(-1.6640321) q[3];
sx q[3];
rz(0.4841218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-0.86301962) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097968) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(-0.61169949) q[0];
rz(-1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-2.591419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69458354) q[0];
sx q[0];
rz(-2.3059079) q[0];
sx q[0];
rz(-1.4674594) q[0];
rz(-pi) q[1];
rz(1.6646531) q[2];
sx q[2];
rz(-1.8641194) q[2];
sx q[2];
rz(-0.23767995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1830814) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(2.1469841) q[1];
x q[2];
rz(-0.21869603) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(-0.72363879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-0.91526389) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4782151) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-0.94232725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.077276262) q[2];
sx q[2];
rz(-2.6151267) q[2];
sx q[2];
rz(-0.72409814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.535714) q[1];
sx q[1];
rz(-1.6026346) q[1];
sx q[1];
rz(-1.72623) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81178938) q[0];
sx q[0];
rz(-1.3178696) q[0];
sx q[0];
rz(-1.1008218) q[0];
x q[1];
rz(-1.8814538) q[2];
sx q[2];
rz(-2.0271745) q[2];
sx q[2];
rz(-3.0820456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.10627667) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(-1.4186526) q[1];
rz(-2.5572705) q[3];
sx q[3];
rz(-1.366368) q[3];
sx q[3];
rz(-1.0014597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5863824) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267245) q[0];
sx q[0];
rz(-1.4818026) q[0];
sx q[0];
rz(-0.29863775) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(-2.8587647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3583402) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(1.7794442) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2231636) q[3];
sx q[3];
rz(-1.6602483) q[3];
sx q[3];
rz(1.7355433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(3.0287108) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(0.032756068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54887933) q[0];
sx q[0];
rz(-1.2476876) q[0];
sx q[0];
rz(-0.26922853) q[0];
x q[1];
rz(-0.92808) q[2];
sx q[2];
rz(-2.4734481) q[2];
sx q[2];
rz(1.1754787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36383648) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(-1.1532409) q[1];
x q[2];
rz(1.1912187) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(0.4656725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(1.9744251) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.099982925) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(2.656235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55886666) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(-2.7607714) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.543407) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(2.8796632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34708729) q[1];
sx q[1];
rz(-0.80650389) q[1];
sx q[1];
rz(-2.0054714) q[1];
rz(-0.58166196) q[3];
sx q[3];
rz(-1.9546486) q[3];
sx q[3];
rz(1.7285085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.52206802) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-0.5982582) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
x q[1];
rz(1.422545) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(-1.6809747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.6655386) q[1];
sx q[1];
rz(1.2050864) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9833097) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(2.3815001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(2.9700206) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(2.2451154) q[3];
sx q[3];
rz(-1.258068) q[3];
sx q[3];
rz(0.51545959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
