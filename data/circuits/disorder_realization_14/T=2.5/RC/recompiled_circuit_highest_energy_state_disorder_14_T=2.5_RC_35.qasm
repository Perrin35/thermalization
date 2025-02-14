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
rz(-0.7402339) q[0];
sx q[0];
rz(4.8010173) q[0];
sx q[0];
rz(12.231449) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(-1.0022751) q[1];
sx q[1];
rz(0.60751539) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7142732) q[0];
sx q[0];
rz(-1.8022984) q[0];
sx q[0];
rz(-1.0697068) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11694853) q[2];
sx q[2];
rz(-2.0772572) q[2];
sx q[2];
rz(0.036265515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40558896) q[1];
sx q[1];
rz(-2.0479255) q[1];
sx q[1];
rz(1.3911584) q[1];
rz(-pi) q[2];
rz(-2.7685952) q[3];
sx q[3];
rz(-1.1524876) q[3];
sx q[3];
rz(-1.7278863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-0.66317916) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(1.1801571) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8993768) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(-1.2902749) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(1.7346409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5840184) q[0];
sx q[0];
rz(-1.0913335) q[0];
sx q[0];
rz(-2.8842501) q[0];
x q[1];
rz(0.036769899) q[2];
sx q[2];
rz(-1.2388133) q[2];
sx q[2];
rz(-0.22184243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9216361) q[1];
sx q[1];
rz(-2.0825279) q[1];
sx q[1];
rz(0.10351609) q[1];
rz(-pi) q[2];
rz(1.155726) q[3];
sx q[3];
rz(-1.0336116) q[3];
sx q[3];
rz(-1.146917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0824288) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(-1.814369) q[2];
rz(2.0969157) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(-0.4253934) q[0];
rz(-1.3771903) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(-1.707071) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73983708) q[0];
sx q[0];
rz(-0.65261894) q[0];
sx q[0];
rz(-2.9463861) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1026406) q[2];
sx q[2];
rz(-0.70868451) q[2];
sx q[2];
rz(1.5295636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.19615281) q[1];
sx q[1];
rz(-1.9827794) q[1];
sx q[1];
rz(-0.71228551) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8355153) q[3];
sx q[3];
rz(-2.1669905) q[3];
sx q[3];
rz(-1.5706737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(1.8505992) q[2];
rz(-1.7839606) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.4972081) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62035471) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(2.3714016) q[0];
rz(0.99984804) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(-1.3410478) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65085852) q[0];
sx q[0];
rz(-2.5528209) q[0];
sx q[0];
rz(2.7522699) q[0];
rz(0.88444986) q[2];
sx q[2];
rz(-1.9362861) q[2];
sx q[2];
rz(2.3523503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51777202) q[1];
sx q[1];
rz(-1.2720577) q[1];
sx q[1];
rz(-1.8483398) q[1];
x q[2];
rz(-1.2433231) q[3];
sx q[3];
rz(-0.30042111) q[3];
sx q[3];
rz(-1.5113561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(-0.10041222) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95555821) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(2.9146063) q[0];
rz(-1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(-1.7452128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76043789) q[0];
sx q[0];
rz(-1.210307) q[0];
sx q[0];
rz(-1.9442149) q[0];
x q[1];
rz(1.7368083) q[2];
sx q[2];
rz(-1.1653656) q[2];
sx q[2];
rz(-2.977598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5310881) q[1];
sx q[1];
rz(-2.5042494) q[1];
sx q[1];
rz(-0.12048851) q[1];
x q[2];
rz(-0.8245411) q[3];
sx q[3];
rz(-1.5932114) q[3];
sx q[3];
rz(0.21460545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22623006) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(-0.076233141) q[2];
rz(1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(-1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48080322) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.3355108) q[0];
rz(-2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(2.9551771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.473523) q[0];
sx q[0];
rz(-2.0041564) q[0];
sx q[0];
rz(0.98287232) q[0];
x q[1];
rz(-1.3720296) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(-1.1793062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5552206) q[1];
sx q[1];
rz(-2.456291) q[1];
sx q[1];
rz(-1.1838811) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3535054) q[3];
sx q[3];
rz(-2.4895146) q[3];
sx q[3];
rz(1.2913845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4621801) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(1.3235486) q[2];
rz(-1.1421674) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(2.2511258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46397504) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(-2.5073945) q[0];
rz(1.0448666) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(0.96010906) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8909559) q[0];
sx q[0];
rz(-1.5433558) q[0];
sx q[0];
rz(1.1720042) q[0];
x q[1];
rz(-1.5788001) q[2];
sx q[2];
rz(-1.8718613) q[2];
sx q[2];
rz(2.0640399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81183103) q[1];
sx q[1];
rz(-2.0961746) q[1];
sx q[1];
rz(-3.1282022) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3784901) q[3];
sx q[3];
rz(-2.1402485) q[3];
sx q[3];
rz(2.4214937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94379696) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(1.724285) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(1.5203169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6087795) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(3.0335562) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.9727762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10442142) q[0];
sx q[0];
rz(-1.1544495) q[0];
sx q[0];
rz(-1.0486616) q[0];
rz(-1.5708099) q[2];
sx q[2];
rz(-2.332649) q[2];
sx q[2];
rz(-2.1894313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16202422) q[1];
sx q[1];
rz(-1.3378394) q[1];
sx q[1];
rz(0.090090171) q[1];
rz(-pi) q[2];
rz(-3.1049012) q[3];
sx q[3];
rz(-0.49202575) q[3];
sx q[3];
rz(1.2855315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10890266) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(0.6558134) q[2];
rz(3.0374895) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(-0.18812215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8769237) q[0];
sx q[0];
rz(-1.8380565) q[0];
sx q[0];
rz(-1.4917829) q[0];
rz(2.1022294) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(0.46844354) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2830848) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(1.5664738) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2861112) q[2];
sx q[2];
rz(-1.9980901) q[2];
sx q[2];
rz(-2.9197249) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2704308) q[1];
sx q[1];
rz(-2.3498145) q[1];
sx q[1];
rz(-1.3429848) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78448589) q[3];
sx q[3];
rz(-2.753559) q[3];
sx q[3];
rz(-1.3875543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.09482065) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(0.92791933) q[2];
rz(1.8038484) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-3.0874262) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(-1.077865) q[0];
rz(-2.7742591) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(1.8574538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8074984) q[0];
sx q[0];
rz(-2.0413412) q[0];
sx q[0];
rz(0.28326359) q[0];
x q[1];
rz(1.4201035) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(-0.53021741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2312647) q[1];
sx q[1];
rz(-0.18016768) q[1];
sx q[1];
rz(0.50039165) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23350164) q[3];
sx q[3];
rz(-2.3931008) q[3];
sx q[3];
rz(-3.0155663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(-2.7122811) q[2];
rz(2.9863827) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(1.8163053) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406381) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(-0.86391972) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(-0.39045329) q[2];
sx q[2];
rz(-0.57885546) q[2];
sx q[2];
rz(-0.59138966) q[2];
rz(0.94384296) q[3];
sx q[3];
rz(-1.8494434) q[3];
sx q[3];
rz(-0.2260126) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
