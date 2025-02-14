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
rz(-1.9395186) q[0];
sx q[0];
rz(-1.0459067) q[0];
sx q[0];
rz(0.75049555) q[0];
rz(6.5039492) q[1];
sx q[1];
rz(1.060744) q[1];
sx q[1];
rz(4.3252601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25607294) q[0];
sx q[0];
rz(-1.1587496) q[0];
sx q[0];
rz(2.4072007) q[0];
rz(-0.92249845) q[2];
sx q[2];
rz(-2.9335409) q[2];
sx q[2];
rz(1.7788943) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.466063) q[1];
sx q[1];
rz(-2.4616082) q[1];
sx q[1];
rz(-3.1230981) q[1];
rz(-pi) q[2];
rz(-2.9394807) q[3];
sx q[3];
rz(-1.8732866) q[3];
sx q[3];
rz(-1.9155531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.815328) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(-3.0495138) q[2];
rz(1.1312671) q[3];
sx q[3];
rz(-1.0695894) q[3];
sx q[3];
rz(2.2856975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6446514) q[0];
sx q[0];
rz(-2.6583789) q[0];
sx q[0];
rz(0.53571969) q[0];
rz(3.1247395) q[1];
sx q[1];
rz(-2.2622175) q[1];
sx q[1];
rz(-0.0171612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0417174) q[0];
sx q[0];
rz(-1.2493396) q[0];
sx q[0];
rz(-2.0896555) q[0];
rz(0.13249915) q[2];
sx q[2];
rz(-0.4702869) q[2];
sx q[2];
rz(-2.7072226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29083272) q[1];
sx q[1];
rz(-0.28181812) q[1];
sx q[1];
rz(-2.5435104) q[1];
rz(-pi) q[2];
rz(-1.4799177) q[3];
sx q[3];
rz(-1.7044997) q[3];
sx q[3];
rz(2.9667678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0422334) q[2];
sx q[2];
rz(-1.4734273) q[2];
sx q[2];
rz(-2.8933706) q[2];
rz(0.92784268) q[3];
sx q[3];
rz(-0.78301269) q[3];
sx q[3];
rz(-0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8645653) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(0.16265854) q[0];
rz(0.76205572) q[1];
sx q[1];
rz(-2.9167852) q[1];
sx q[1];
rz(2.3077097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974166) q[0];
sx q[0];
rz(-0.58338291) q[0];
sx q[0];
rz(0.22384217) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35807407) q[2];
sx q[2];
rz(-0.34293567) q[2];
sx q[2];
rz(-1.5141405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1385041) q[1];
sx q[1];
rz(-0.89996594) q[1];
sx q[1];
rz(2.6062232) q[1];
x q[2];
rz(-1.6271851) q[3];
sx q[3];
rz(-1.1306376) q[3];
sx q[3];
rz(1.7007916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0791846) q[2];
sx q[2];
rz(-1.1391613) q[2];
sx q[2];
rz(-2.8774234) q[2];
rz(1.0583813) q[3];
sx q[3];
rz(-2.1975785) q[3];
sx q[3];
rz(-0.03037608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3707378) q[0];
sx q[0];
rz(-0.83468947) q[0];
sx q[0];
rz(-1.7133065) q[0];
rz(0.76438534) q[1];
sx q[1];
rz(-1.1420219) q[1];
sx q[1];
rz(-1.3216602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.747731) q[0];
sx q[0];
rz(-1.9698304) q[0];
sx q[0];
rz(-1.948843) q[0];
rz(-1.9025397) q[2];
sx q[2];
rz(-1.52196) q[2];
sx q[2];
rz(1.2570753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0881867) q[1];
sx q[1];
rz(-1.2866396) q[1];
sx q[1];
rz(2.8233145) q[1];
rz(-pi) q[2];
rz(2.2854487) q[3];
sx q[3];
rz(-1.9106746) q[3];
sx q[3];
rz(-1.4902756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6617714) q[2];
sx q[2];
rz(-1.4134553) q[2];
sx q[2];
rz(0.28523764) q[2];
rz(-1.8789004) q[3];
sx q[3];
rz(-1.217239) q[3];
sx q[3];
rz(-1.4234022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406463) q[0];
sx q[0];
rz(-0.74746376) q[0];
sx q[0];
rz(-2.2593011) q[0];
rz(2.4709985) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(2.1569599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1272972) q[0];
sx q[0];
rz(-1.5080351) q[0];
sx q[0];
rz(1.602229) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6419471) q[2];
sx q[2];
rz(-2.5987968) q[2];
sx q[2];
rz(1.1384447) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9198926) q[1];
sx q[1];
rz(-1.0439149) q[1];
sx q[1];
rz(2.5008898) q[1];
x q[2];
rz(0.71159848) q[3];
sx q[3];
rz(-2.4715373) q[3];
sx q[3];
rz(-1.3942476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8661477) q[2];
sx q[2];
rz(-0.53469849) q[2];
sx q[2];
rz(0.61694413) q[2];
rz(2.3894737) q[3];
sx q[3];
rz(-1.3277466) q[3];
sx q[3];
rz(2.6822958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1521456) q[0];
sx q[0];
rz(-0.716827) q[0];
sx q[0];
rz(-2.3710807) q[0];
rz(-0.72692263) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(2.4368584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7677697) q[0];
sx q[0];
rz(-2.2859068) q[0];
sx q[0];
rz(2.3138347) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36811916) q[2];
sx q[2];
rz(-1.0794753) q[2];
sx q[2];
rz(1.7047326) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.875346) q[1];
sx q[1];
rz(-0.35142144) q[1];
sx q[1];
rz(-3.0446192) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25382692) q[3];
sx q[3];
rz(-1.0146662) q[3];
sx q[3];
rz(-1.1185916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7695693) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(2.7083569) q[2];
rz(-2.9406934) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(-2.6065629) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3055426) q[0];
sx q[0];
rz(-1.9925646) q[0];
sx q[0];
rz(-2.7591163) q[0];
rz(2.1005232) q[1];
sx q[1];
rz(-1.6761227) q[1];
sx q[1];
rz(2.6952851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6093369) q[0];
sx q[0];
rz(-2.8777538) q[0];
sx q[0];
rz(0.20045073) q[0];
rz(0.52526229) q[2];
sx q[2];
rz(-2.7064311) q[2];
sx q[2];
rz(2.9829142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.32532) q[1];
sx q[1];
rz(-2.153646) q[1];
sx q[1];
rz(0.42709777) q[1];
rz(-2.2931523) q[3];
sx q[3];
rz(-2.1506243) q[3];
sx q[3];
rz(0.27001303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.81991601) q[2];
sx q[2];
rz(-2.0240462) q[2];
sx q[2];
rz(2.4134911) q[2];
rz(-0.32650945) q[3];
sx q[3];
rz(-1.510334) q[3];
sx q[3];
rz(-0.85730332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.2364872) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(1.2976728) q[0];
rz(-2.1406651) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(-0.40850684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1649969) q[0];
sx q[0];
rz(-1.5566096) q[0];
sx q[0];
rz(-3.1407457) q[0];
rz(2.0493919) q[2];
sx q[2];
rz(-1.9754306) q[2];
sx q[2];
rz(1.7737845) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4104011) q[1];
sx q[1];
rz(-2.9653314) q[1];
sx q[1];
rz(-0.86430637) q[1];
rz(0.85902416) q[3];
sx q[3];
rz(-2.6005496) q[3];
sx q[3];
rz(-2.4699901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11477509) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(2.3180023) q[2];
rz(2.1987727) q[3];
sx q[3];
rz(-1.4225682) q[3];
sx q[3];
rz(-1.5028809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59458643) q[0];
sx q[0];
rz(-0.44438812) q[0];
sx q[0];
rz(1.8894926) q[0];
rz(1.7777255) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(-1.8069256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20864381) q[0];
sx q[0];
rz(-1.874039) q[0];
sx q[0];
rz(-1.6088845) q[0];
rz(-pi) q[1];
rz(-1.3965071) q[2];
sx q[2];
rz(-2.2054923) q[2];
sx q[2];
rz(1.8746992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5022565) q[1];
sx q[1];
rz(-0.4570804) q[1];
sx q[1];
rz(0.16641081) q[1];
x q[2];
rz(2.8765479) q[3];
sx q[3];
rz(-1.1172034) q[3];
sx q[3];
rz(2.3475304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1198662) q[2];
sx q[2];
rz(-2.848337) q[2];
sx q[2];
rz(-2.8450656) q[2];
rz(1.4793388) q[3];
sx q[3];
rz(-1.5115503) q[3];
sx q[3];
rz(0.1571981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6327561) q[0];
sx q[0];
rz(-2.4907676) q[0];
sx q[0];
rz(0.79488361) q[0];
rz(2.5859313) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(1.6937675) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7894365) q[0];
sx q[0];
rz(-1.1934501) q[0];
sx q[0];
rz(-2.3579602) q[0];
rz(0.49782355) q[2];
sx q[2];
rz(-1.6256672) q[2];
sx q[2];
rz(-0.16615088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9196133) q[1];
sx q[1];
rz(-1.3718425) q[1];
sx q[1];
rz(1.0859906) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85526222) q[3];
sx q[3];
rz(-0.73270933) q[3];
sx q[3];
rz(0.39855188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8453703) q[2];
sx q[2];
rz(-1.5105379) q[2];
sx q[2];
rz(-1.7108062) q[2];
rz(0.75302643) q[3];
sx q[3];
rz(-2.3236661) q[3];
sx q[3];
rz(0.16451612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13567781) q[0];
sx q[0];
rz(-0.68886859) q[0];
sx q[0];
rz(-1.9266358) q[0];
rz(-1.4665435) q[1];
sx q[1];
rz(-2.2129682) q[1];
sx q[1];
rz(0.60563544) q[1];
rz(-0.56500021) q[2];
sx q[2];
rz(-2.0774659) q[2];
sx q[2];
rz(-2.5027711) q[2];
rz(-0.10726992) q[3];
sx q[3];
rz(-1.4013877) q[3];
sx q[3];
rz(0.79584484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
