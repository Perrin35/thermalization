OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(-3.1080973) q[0];
sx q[0];
rz(-1.7749696) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(4.7935901) q[1];
sx q[1];
rz(10.556769) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7819408) q[0];
sx q[0];
rz(-1.9160761) q[0];
sx q[0];
rz(-2.3495673) q[0];
rz(-pi) q[1];
rz(2.9801324) q[2];
sx q[2];
rz(-1.7345288) q[2];
sx q[2];
rz(-1.8979567) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0238029) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(2.6610713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55478401) q[3];
sx q[3];
rz(-1.8150419) q[3];
sx q[3];
rz(-0.82575219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(-1.988391) q[2];
rz(0.48405805) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(2.6822065) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83950481) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(-2.6385345) q[0];
rz(1.5548276) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(0.15393004) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.670823) q[0];
sx q[0];
rz(-1.5851067) q[0];
sx q[0];
rz(-1.0679246) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8106861) q[2];
sx q[2];
rz(-0.94793301) q[2];
sx q[2];
rz(-1.1922342) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1679722) q[1];
sx q[1];
rz(-2.4929843) q[1];
sx q[1];
rz(-0.82778511) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50977317) q[3];
sx q[3];
rz(-1.5700649) q[3];
sx q[3];
rz(0.36611205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8872035) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(1.8187693) q[2];
rz(-1.4860738) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(-2.4310908) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1948497) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(-0.90993607) q[0];
rz(2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(2.1562703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5062149) q[0];
sx q[0];
rz(-0.92184508) q[0];
sx q[0];
rz(-0.20091591) q[0];
rz(1.9820205) q[2];
sx q[2];
rz(-0.69934884) q[2];
sx q[2];
rz(1.5027836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31919033) q[1];
sx q[1];
rz(-1.625758) q[1];
sx q[1];
rz(0.24939202) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3722234) q[3];
sx q[3];
rz(-2.3017075) q[3];
sx q[3];
rz(1.6826671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.862792) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(-1.2476236) q[2];
rz(-0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(0.32143337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9001532) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(1.1234269) q[0];
rz(-0.57888794) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.7480063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508586) q[0];
sx q[0];
rz(-2.3669741) q[0];
sx q[0];
rz(-2.5289092) q[0];
x q[1];
rz(-2.706714) q[2];
sx q[2];
rz(-1.7231427) q[2];
sx q[2];
rz(-0.16532126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77809282) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-2.5219445) q[1];
x q[2];
rz(-1.8279151) q[3];
sx q[3];
rz(-2.3895279) q[3];
sx q[3];
rz(-2.5079692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-3.139479) q[2];
rz(0.56143108) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(0.58825618) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.026022) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(1.2258688) q[0];
rz(-1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.7747169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5565235) q[0];
sx q[0];
rz(-1.1856623) q[0];
sx q[0];
rz(-1.3655846) q[0];
rz(-pi) q[1];
rz(-2.3943704) q[2];
sx q[2];
rz(-1.5486273) q[2];
sx q[2];
rz(1.6037461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99202427) q[1];
sx q[1];
rz(-0.7115041) q[1];
sx q[1];
rz(1.6929388) q[1];
rz(-pi) q[2];
rz(0.59646888) q[3];
sx q[3];
rz(-0.58613741) q[3];
sx q[3];
rz(-2.6479207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(0.40536353) q[2];
rz(2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(1.595114) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(-2.6382085) q[0];
rz(-0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(-0.65178451) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3394649) q[0];
sx q[0];
rz(-1.8101242) q[0];
sx q[0];
rz(-1.375074) q[0];
rz(-pi) q[1];
rz(0.82181828) q[2];
sx q[2];
rz(-1.6449252) q[2];
sx q[2];
rz(-0.2766343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5469049) q[1];
sx q[1];
rz(-1.6225796) q[1];
sx q[1];
rz(-1.7054218) q[1];
rz(-pi) q[2];
rz(2.7995293) q[3];
sx q[3];
rz(-2.8422589) q[3];
sx q[3];
rz(-0.39810668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51263222) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(0.39166489) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(0.96965924) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-3.0113509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0129583) q[0];
sx q[0];
rz(-0.35612125) q[0];
sx q[0];
rz(-2.998259) q[0];
rz(-pi) q[1];
rz(-2.4728751) q[2];
sx q[2];
rz(-2.0520743) q[2];
sx q[2];
rz(0.95578335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3324193) q[1];
sx q[1];
rz(-1.0671339) q[1];
sx q[1];
rz(-1.3219576) q[1];
rz(-pi) q[2];
rz(2.7534915) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(-2.2263118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-1.0858034) q[2];
rz(-0.052224934) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645585) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(-1.9751256) q[0];
rz(0.72558609) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(2.3988147) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60177875) q[0];
sx q[0];
rz(-0.83866461) q[0];
sx q[0];
rz(1.3108805) q[0];
x q[1];
rz(2.5119945) q[2];
sx q[2];
rz(-0.78231914) q[2];
sx q[2];
rz(-0.19537374) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.038866) q[1];
sx q[1];
rz(-0.83737774) q[1];
sx q[1];
rz(-2.1329761) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5524213) q[3];
sx q[3];
rz(-0.99654752) q[3];
sx q[3];
rz(-1.9598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(0.26088866) q[2];
rz(-2.0339113) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(-2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(0.93604952) q[0];
rz(-3.127457) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(-0.65151185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441898) q[0];
sx q[0];
rz(-1.5751) q[0];
sx q[0];
rz(-1.5944907) q[0];
rz(-pi) q[1];
rz(-0.94349761) q[2];
sx q[2];
rz(-2.128696) q[2];
sx q[2];
rz(-3.0859335) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9459322) q[1];
sx q[1];
rz(-1.6911427) q[1];
sx q[1];
rz(1.8912998) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7016919) q[3];
sx q[3];
rz(-1.2370046) q[3];
sx q[3];
rz(1.7037638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8462048) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(-0.52337581) q[2];
rz(-2.8335617) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(-1.403632) q[0];
rz(-2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(1.7636991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4353367) q[0];
sx q[0];
rz(-2.113111) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(-pi) q[1];
rz(-0.84656287) q[2];
sx q[2];
rz(-2.2228974) q[2];
sx q[2];
rz(-2.6400499) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5046895) q[1];
sx q[1];
rz(-2.6733589) q[1];
sx q[1];
rz(0.56028985) q[1];
x q[2];
rz(1.4078478) q[3];
sx q[3];
rz(-1.019422) q[3];
sx q[3];
rz(2.5577302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.750981) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(1.5314468) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(-2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703341) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(2.4709672) q[2];
sx q[2];
rz(-1.8319596) q[2];
sx q[2];
rz(-1.7419227) q[2];
rz(0.74606568) q[3];
sx q[3];
rz(-2.2545771) q[3];
sx q[3];
rz(1.547326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];