OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27788568) q[0];
sx q[0];
rz(3.570896) q[0];
sx q[0];
rz(12.255393) q[0];
rz(-1.3485981) q[1];
sx q[1];
rz(4.2378874) q[1];
sx q[1];
rz(8.8506946) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35168655) q[0];
sx q[0];
rz(-1.154019) q[0];
sx q[0];
rz(-0.4051286) q[0];
rz(1.5181731) q[2];
sx q[2];
rz(-2.5580018) q[2];
sx q[2];
rz(-0.4837732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52770268) q[1];
sx q[1];
rz(-0.37606323) q[1];
sx q[1];
rz(1.3251873) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0776377) q[3];
sx q[3];
rz(-1.9465145) q[3];
sx q[3];
rz(0.8448173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2744039) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(-1.8223507) q[2];
rz(-0.073171767) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0480334) q[0];
sx q[0];
rz(-1.1807384) q[0];
sx q[0];
rz(2.4617885) q[0];
rz(-0.80348429) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(0.63978535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2701975) q[0];
sx q[0];
rz(-1.2551487) q[0];
sx q[0];
rz(2.2391566) q[0];
x q[1];
rz(-2.4117208) q[2];
sx q[2];
rz(-1.6449071) q[2];
sx q[2];
rz(-2.6593936) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82347875) q[1];
sx q[1];
rz(-2.1770855) q[1];
sx q[1];
rz(-0.73474523) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7946258) q[3];
sx q[3];
rz(-1.2696243) q[3];
sx q[3];
rz(-2.6909242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60376072) q[2];
sx q[2];
rz(-0.49850264) q[2];
sx q[2];
rz(2.4375088) q[2];
rz(3.0294561) q[3];
sx q[3];
rz(-1.4487368) q[3];
sx q[3];
rz(-1.0224379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10244399) q[0];
sx q[0];
rz(-1.1071438) q[0];
sx q[0];
rz(2.802134) q[0];
rz(1.1184232) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(-3.1059473) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064020412) q[0];
sx q[0];
rz(-1.3638931) q[0];
sx q[0];
rz(1.4830566) q[0];
rz(-0.33907922) q[2];
sx q[2];
rz(-1.2847752) q[2];
sx q[2];
rz(-1.7947444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34295568) q[1];
sx q[1];
rz(-1.3561212) q[1];
sx q[1];
rz(-0.66185419) q[1];
rz(2.6958353) q[3];
sx q[3];
rz(-0.66346079) q[3];
sx q[3];
rz(1.9395246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3493335) q[2];
sx q[2];
rz(-2.4534241) q[2];
sx q[2];
rz(-2.2815857) q[2];
rz(-2.3490014) q[3];
sx q[3];
rz(-0.20321295) q[3];
sx q[3];
rz(-1.0205166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0122796) q[0];
sx q[0];
rz(-1.7565933) q[0];
sx q[0];
rz(-2.5087575) q[0];
rz(-1.9889779) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(1.4702183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2180207) q[0];
sx q[0];
rz(-1.5578736) q[0];
sx q[0];
rz(-2.2201029) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15369065) q[2];
sx q[2];
rz(-2.3665303) q[2];
sx q[2];
rz(2.3714921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51293514) q[1];
sx q[1];
rz(-0.75408616) q[1];
sx q[1];
rz(-1.1910466) q[1];
rz(2.4864462) q[3];
sx q[3];
rz(-2.0787424) q[3];
sx q[3];
rz(-2.9367713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4319438) q[2];
sx q[2];
rz(-0.96082965) q[2];
sx q[2];
rz(0.43080899) q[2];
rz(-2.4591947) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(0.94990134) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34669852) q[0];
sx q[0];
rz(-1.2656724) q[0];
sx q[0];
rz(-1.9516113) q[0];
rz(-0.31040141) q[1];
sx q[1];
rz(-2.0605395) q[1];
sx q[1];
rz(0.50500542) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81722084) q[0];
sx q[0];
rz(-2.007487) q[0];
sx q[0];
rz(0.94193108) q[0];
rz(-pi) q[1];
rz(0.1311028) q[2];
sx q[2];
rz(-1.1590247) q[2];
sx q[2];
rz(0.51136651) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6663899) q[1];
sx q[1];
rz(-2.0670927) q[1];
sx q[1];
rz(1.6700891) q[1];
rz(-1.4851863) q[3];
sx q[3];
rz(-0.84692779) q[3];
sx q[3];
rz(1.6023265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7348822) q[2];
sx q[2];
rz(-2.4411185) q[2];
sx q[2];
rz(-0.90275466) q[2];
rz(-2.086575) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(-0.46190754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91586739) q[0];
sx q[0];
rz(-0.71284717) q[0];
sx q[0];
rz(-2.0676887) q[0];
rz(1.7621) q[1];
sx q[1];
rz(-2.4122489) q[1];
sx q[1];
rz(3.0440547) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.486575) q[0];
sx q[0];
rz(-1.3897822) q[0];
sx q[0];
rz(0.54365309) q[0];
x q[1];
rz(1.9768841) q[2];
sx q[2];
rz(-2.7345719) q[2];
sx q[2];
rz(0.87813745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40332023) q[1];
sx q[1];
rz(-0.077906713) q[1];
sx q[1];
rz(2.2408443) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0076526) q[3];
sx q[3];
rz(-2.030367) q[3];
sx q[3];
rz(-2.8569222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0292616) q[2];
sx q[2];
rz(-1.6882201) q[2];
sx q[2];
rz(2.7304999) q[2];
rz(1.4136275) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(-1.5467862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(0.77808648) q[0];
sx q[0];
rz(-3.111105) q[0];
sx q[0];
rz(-1.4084858) q[0];
rz(2.1259437) q[1];
sx q[1];
rz(-1.3280832) q[1];
sx q[1];
rz(-1.4782864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9736098) q[0];
sx q[0];
rz(-2.4496299) q[0];
sx q[0];
rz(0.47398098) q[0];
rz(-0.99009606) q[2];
sx q[2];
rz(-1.0063356) q[2];
sx q[2];
rz(-1.9866458) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95132213) q[1];
sx q[1];
rz(-0.45004216) q[1];
sx q[1];
rz(-1.3954891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0706581) q[3];
sx q[3];
rz(-1.4551509) q[3];
sx q[3];
rz(2.107419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1136855) q[2];
sx q[2];
rz(-2.5137641) q[2];
sx q[2];
rz(-2.9713463) q[2];
rz(1.2804735) q[3];
sx q[3];
rz(-1.4743285) q[3];
sx q[3];
rz(-1.1580275) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.668642) q[0];
sx q[0];
rz(-2.1773715) q[0];
sx q[0];
rz(0.52841312) q[0];
rz(-2.2083185) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(0.55567137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816102) q[0];
sx q[0];
rz(-2.1851808) q[0];
sx q[0];
rz(2.9016764) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3494928) q[2];
sx q[2];
rz(-1.6413771) q[2];
sx q[2];
rz(0.5943228) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1126777) q[1];
sx q[1];
rz(-2.2820447) q[1];
sx q[1];
rz(-1.6131141) q[1];
rz(-pi) q[2];
rz(-0.37013388) q[3];
sx q[3];
rz(-0.72648337) q[3];
sx q[3];
rz(-0.47011061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3282503) q[2];
sx q[2];
rz(-0.80208653) q[2];
sx q[2];
rz(-0.30725202) q[2];
rz(0.79636374) q[3];
sx q[3];
rz(-0.73085228) q[3];
sx q[3];
rz(-0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873213) q[0];
sx q[0];
rz(-1.048943) q[0];
sx q[0];
rz(1.3294504) q[0];
rz(0.21996552) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(-1.3185917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057379) q[0];
sx q[0];
rz(-1.8377343) q[0];
sx q[0];
rz(0.22375317) q[0];
rz(-2.6418453) q[2];
sx q[2];
rz(-1.6732273) q[2];
sx q[2];
rz(-0.62291716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3941806) q[1];
sx q[1];
rz(-2.5492382) q[1];
sx q[1];
rz(-0.92983957) q[1];
x q[2];
rz(2.044146) q[3];
sx q[3];
rz(-2.2339377) q[3];
sx q[3];
rz(1.1569661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69507504) q[2];
sx q[2];
rz(-1.4277642) q[2];
sx q[2];
rz(1.2850777) q[2];
rz(-2.255693) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.5281965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000243) q[0];
sx q[0];
rz(-0.73606857) q[0];
sx q[0];
rz(0.31684434) q[0];
rz(-0.44081229) q[1];
sx q[1];
rz(-0.79439729) q[1];
sx q[1];
rz(-2.7712834) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.142684) q[0];
sx q[0];
rz(-0.38247358) q[0];
sx q[0];
rz(-1.6539025) q[0];
x q[1];
rz(0.84443386) q[2];
sx q[2];
rz(-2.0877247) q[2];
sx q[2];
rz(-0.58973613) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2930544) q[1];
sx q[1];
rz(-1.8369743) q[1];
sx q[1];
rz(1.2598424) q[1];
rz(-pi) q[2];
rz(-3.0906865) q[3];
sx q[3];
rz(-1.6744782) q[3];
sx q[3];
rz(3.0391676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.73654282) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(1.6323818) q[2];
rz(0.80471188) q[3];
sx q[3];
rz(-1.5038306) q[3];
sx q[3];
rz(-0.10778431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1152773) q[0];
sx q[0];
rz(-1.7290709) q[0];
sx q[0];
rz(1.2363634) q[0];
rz(0.23030494) q[1];
sx q[1];
rz(-0.31619148) q[1];
sx q[1];
rz(-2.2264623) q[1];
rz(-0.7028107) q[2];
sx q[2];
rz(-2.3431449) q[2];
sx q[2];
rz(0.14771067) q[2];
rz(0.32634278) q[3];
sx q[3];
rz(-0.97728609) q[3];
sx q[3];
rz(1.4090007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
