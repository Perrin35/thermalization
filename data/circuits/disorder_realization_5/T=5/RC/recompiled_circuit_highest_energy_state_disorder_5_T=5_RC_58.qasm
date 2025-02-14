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
rz(0.24231237) q[0];
sx q[0];
rz(-2.2678312) q[0];
sx q[0];
rz(1.4886966) q[0];
rz(-2.0714662) q[1];
sx q[1];
rz(-0.5572328) q[1];
sx q[1];
rz(0.4692404) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0996286) q[0];
sx q[0];
rz(-1.6828186) q[0];
sx q[0];
rz(-2.4449181) q[0];
rz(-2.4503684) q[2];
sx q[2];
rz(-0.35688734) q[2];
sx q[2];
rz(1.7920902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1716071) q[1];
sx q[1];
rz(-2.0678221) q[1];
sx q[1];
rz(-2.7812048) q[1];
x q[2];
rz(-1.038299) q[3];
sx q[3];
rz(-0.91404741) q[3];
sx q[3];
rz(2.9503502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0327518) q[2];
sx q[2];
rz(-2.8665172) q[2];
sx q[2];
rz(1.8185505) q[2];
rz(-2.1950586) q[3];
sx q[3];
rz(-0.60775477) q[3];
sx q[3];
rz(0.058145903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7846947) q[0];
sx q[0];
rz(-0.31814831) q[0];
sx q[0];
rz(-1.1784877) q[0];
rz(-2.5937041) q[1];
sx q[1];
rz(-1.9375487) q[1];
sx q[1];
rz(-1.2335802) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6797076) q[0];
sx q[0];
rz(-1.8543482) q[0];
sx q[0];
rz(-0.028282982) q[0];
x q[1];
rz(1.7370102) q[2];
sx q[2];
rz(-1.4912327) q[2];
sx q[2];
rz(-0.97841371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0280721) q[1];
sx q[1];
rz(-1.6238316) q[1];
sx q[1];
rz(3.0854532) q[1];
x q[2];
rz(0.61290755) q[3];
sx q[3];
rz(-2.7281085) q[3];
sx q[3];
rz(-0.72847937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2445688) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(1.7986253) q[2];
rz(-0.025394214) q[3];
sx q[3];
rz(-1.3586905) q[3];
sx q[3];
rz(3.0620745) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19567604) q[0];
sx q[0];
rz(-1.8657277) q[0];
sx q[0];
rz(-0.50862408) q[0];
rz(1.7018082) q[1];
sx q[1];
rz(-1.6248117) q[1];
sx q[1];
rz(-1.5135099) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.968983) q[0];
sx q[0];
rz(-1.727761) q[0];
sx q[0];
rz(-0.60625546) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3966826) q[2];
sx q[2];
rz(-1.422494) q[2];
sx q[2];
rz(2.1789511) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.73838961) q[1];
sx q[1];
rz(-1.566267) q[1];
sx q[1];
rz(0.036935115) q[1];
x q[2];
rz(-1.1636038) q[3];
sx q[3];
rz(-1.9137205) q[3];
sx q[3];
rz(-0.69780998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9958008) q[2];
sx q[2];
rz(-2.3172816) q[2];
sx q[2];
rz(-2.9902747) q[2];
rz(0.96210903) q[3];
sx q[3];
rz(-1.5113219) q[3];
sx q[3];
rz(0.37402737) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91604084) q[0];
sx q[0];
rz(-0.89260888) q[0];
sx q[0];
rz(1.4501866) q[0];
rz(-1.0610896) q[1];
sx q[1];
rz(-2.7174945) q[1];
sx q[1];
rz(-2.012097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0267222) q[0];
sx q[0];
rz(-0.78512275) q[0];
sx q[0];
rz(-0.27717395) q[0];
x q[1];
rz(1.1808632) q[2];
sx q[2];
rz(-2.0695114) q[2];
sx q[2];
rz(2.5759144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7348014) q[1];
sx q[1];
rz(-0.86002058) q[1];
sx q[1];
rz(-2.9578277) q[1];
rz(-pi) q[2];
rz(-1.8742423) q[3];
sx q[3];
rz(-3.0253694) q[3];
sx q[3];
rz(-0.35467142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9811161) q[2];
sx q[2];
rz(-0.99157292) q[2];
sx q[2];
rz(2.0151095) q[2];
rz(0.22731656) q[3];
sx q[3];
rz(-1.428363) q[3];
sx q[3];
rz(0.050392438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324976) q[0];
sx q[0];
rz(-0.80642527) q[0];
sx q[0];
rz(-2.4438013) q[0];
rz(1.5296439) q[1];
sx q[1];
rz(-0.35747129) q[1];
sx q[1];
rz(-0.43561092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39640204) q[0];
sx q[0];
rz(-1.309309) q[0];
sx q[0];
rz(-3.0130062) q[0];
rz(-pi) q[1];
rz(-2.7694791) q[2];
sx q[2];
rz(-0.96192322) q[2];
sx q[2];
rz(3.0022668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.10098329) q[1];
sx q[1];
rz(-2.6517649) q[1];
sx q[1];
rz(0.71160729) q[1];
rz(-pi) q[2];
rz(1.645363) q[3];
sx q[3];
rz(-1.8454937) q[3];
sx q[3];
rz(-2.4303586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39516732) q[2];
sx q[2];
rz(-2.1992511) q[2];
sx q[2];
rz(0.64685267) q[2];
rz(0.70710373) q[3];
sx q[3];
rz(-3.0349858) q[3];
sx q[3];
rz(2.148237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7917204) q[0];
sx q[0];
rz(-0.40769044) q[0];
sx q[0];
rz(0.83734018) q[0];
rz(2.1912241) q[1];
sx q[1];
rz(-1.252754) q[1];
sx q[1];
rz(-0.18384917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35826916) q[0];
sx q[0];
rz(-1.2766495) q[0];
sx q[0];
rz(-1.8479895) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2410795) q[2];
sx q[2];
rz(-1.7779254) q[2];
sx q[2];
rz(1.1348534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5317418) q[1];
sx q[1];
rz(-0.86656266) q[1];
sx q[1];
rz(2.0805712) q[1];
rz(1.4687187) q[3];
sx q[3];
rz(-1.3390171) q[3];
sx q[3];
rz(-1.8382324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4652555) q[2];
sx q[2];
rz(-2.4691171) q[2];
sx q[2];
rz(2.903751) q[2];
rz(3.0658718) q[3];
sx q[3];
rz(-3.0140311) q[3];
sx q[3];
rz(-2.7231176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38943648) q[0];
sx q[0];
rz(-0.51920033) q[0];
sx q[0];
rz(0.12553781) q[0];
rz(2.708066) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(-1.7009521) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020156064) q[0];
sx q[0];
rz(-1.5227063) q[0];
sx q[0];
rz(-1.5096971) q[0];
x q[1];
rz(2.6532778) q[2];
sx q[2];
rz(-1.7124694) q[2];
sx q[2];
rz(-0.087208793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6969917) q[1];
sx q[1];
rz(-1.9755413) q[1];
sx q[1];
rz(0.54411097) q[1];
x q[2];
rz(2.3818199) q[3];
sx q[3];
rz(-2.2380412) q[3];
sx q[3];
rz(-1.5189288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0354249) q[2];
sx q[2];
rz(-1.9751534) q[2];
sx q[2];
rz(0.20797813) q[2];
rz(-2.4584127) q[3];
sx q[3];
rz(-0.71931374) q[3];
sx q[3];
rz(1.6600367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55027562) q[0];
sx q[0];
rz(-0.90917259) q[0];
sx q[0];
rz(-2.7081178) q[0];
rz(-2.6705961) q[1];
sx q[1];
rz(-0.29312557) q[1];
sx q[1];
rz(2.8004144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29487404) q[0];
sx q[0];
rz(-0.12285422) q[0];
sx q[0];
rz(0.8126749) q[0];
rz(2.3287235) q[2];
sx q[2];
rz(-1.4238266) q[2];
sx q[2];
rz(-0.16260553) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54509974) q[1];
sx q[1];
rz(-2.5148646) q[1];
sx q[1];
rz(-1.2709446) q[1];
rz(-pi) q[2];
rz(1.6834931) q[3];
sx q[3];
rz(-2.1612565) q[3];
sx q[3];
rz(-2.0463508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6097673) q[2];
sx q[2];
rz(-0.34588459) q[2];
sx q[2];
rz(2.7609265) q[2];
rz(-0.62764132) q[3];
sx q[3];
rz(-2.6142879) q[3];
sx q[3];
rz(0.92831534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31955433) q[0];
sx q[0];
rz(-1.1975937) q[0];
sx q[0];
rz(-2.8354216) q[0];
rz(2.6272558) q[1];
sx q[1];
rz(-2.3898333) q[1];
sx q[1];
rz(-2.9591282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7577359) q[0];
sx q[0];
rz(-2.3691874) q[0];
sx q[0];
rz(-0.17018233) q[0];
rz(-pi) q[1];
rz(0.20240558) q[2];
sx q[2];
rz(-0.53779624) q[2];
sx q[2];
rz(-1.1059731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.565158) q[1];
sx q[1];
rz(-1.7243313) q[1];
sx q[1];
rz(1.1576732) q[1];
x q[2];
rz(1.9133041) q[3];
sx q[3];
rz(-1.1376024) q[3];
sx q[3];
rz(-1.1802281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3719486) q[2];
sx q[2];
rz(-0.19522788) q[2];
sx q[2];
rz(0.39607421) q[2];
rz(1.1560446) q[3];
sx q[3];
rz(-2.5542185) q[3];
sx q[3];
rz(-2.7189232) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7523772) q[0];
sx q[0];
rz(-0.14425819) q[0];
sx q[0];
rz(-0.17917646) q[0];
rz(-1.3009118) q[1];
sx q[1];
rz(-0.4549028) q[1];
sx q[1];
rz(0.5295583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048470174) q[0];
sx q[0];
rz(-2.7137965) q[0];
sx q[0];
rz(-0.64944511) q[0];
rz(-pi) q[1];
rz(1.851469) q[2];
sx q[2];
rz(-1.5641128) q[2];
sx q[2];
rz(1.815442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22687616) q[1];
sx q[1];
rz(-0.90639948) q[1];
sx q[1];
rz(-1.9383345) q[1];
x q[2];
rz(-0.29278585) q[3];
sx q[3];
rz(-2.5857627) q[3];
sx q[3];
rz(2.0563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.788488) q[2];
sx q[2];
rz(-2.4994734) q[2];
sx q[2];
rz(-2.8118964) q[2];
rz(-1.6126136) q[3];
sx q[3];
rz(-1.8776882) q[3];
sx q[3];
rz(0.32040709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0223087) q[0];
sx q[0];
rz(-1.4531463) q[0];
sx q[0];
rz(-1.4896738) q[0];
rz(0.48619167) q[1];
sx q[1];
rz(-1.9959027) q[1];
sx q[1];
rz(-2.1015658) q[1];
rz(0.1472856) q[2];
sx q[2];
rz(-2.4596557) q[2];
sx q[2];
rz(-0.13849592) q[2];
rz(-0.43735237) q[3];
sx q[3];
rz(-1.4244546) q[3];
sx q[3];
rz(-0.33311346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
