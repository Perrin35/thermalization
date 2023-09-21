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
rz(3.175088) q[0];
sx q[0];
rz(7.6498084) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(4.7935901) q[1];
sx q[1];
rz(10.556769) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53349797) q[0];
sx q[0];
rz(-2.29288) q[0];
sx q[0];
rz(0.46790926) q[0];
x q[1];
rz(-2.9801324) q[2];
sx q[2];
rz(-1.4070639) q[2];
sx q[2];
rz(-1.8979567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4181617) q[1];
sx q[1];
rz(-1.1685813) q[1];
sx q[1];
rz(2.1869786) q[1];
x q[2];
rz(-0.44191435) q[3];
sx q[3];
rz(-0.60097296) q[3];
sx q[3];
rz(2.0244983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91036096) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(-1.988391) q[2];
rz(-0.48405805) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(0.4593862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020878) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(-0.50305811) q[0];
rz(-1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(2.9876626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0921558) q[0];
sx q[0];
rz(-2.0736118) q[0];
sx q[0];
rz(-0.016331971) q[0];
rz(-2.2203127) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(0.18077476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1679722) q[1];
sx q[1];
rz(-2.4929843) q[1];
sx q[1];
rz(-0.82778511) q[1];
rz(-pi) q[2];
rz(-2.6318195) q[3];
sx q[3];
rz(-1.5715277) q[3];
sx q[3];
rz(-0.36611205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8872035) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(1.3228234) q[2];
rz(-1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-2.4310908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1948497) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(2.2316566) q[0];
rz(0.77727708) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(-2.1562703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9604208) q[0];
sx q[0];
rz(-0.67502484) q[0];
sx q[0];
rz(1.3135364) q[0];
rz(-pi) q[1];
rz(-0.3243685) q[2];
sx q[2];
rz(-2.2019221) q[2];
sx q[2];
rz(-0.98482519) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6775551) q[1];
sx q[1];
rz(-2.8863393) q[1];
sx q[1];
rz(0.21932253) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4661029) q[3];
sx q[3];
rz(-1.0263066) q[3];
sx q[3];
rz(-0.68577037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.862792) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(-1.2476236) q[2];
rz(0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-0.32143337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2414395) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.4588979) q[1];
sx q[1];
rz(1.7480063) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907341) q[0];
sx q[0];
rz(-0.77461857) q[0];
sx q[0];
rz(2.5289092) q[0];
rz(1.4030928) q[2];
sx q[2];
rz(-1.1412914) q[2];
sx q[2];
rz(1.6657366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77809282) q[1];
sx q[1];
rz(-0.85077319) q[1];
sx q[1];
rz(2.5219445) q[1];
rz(2.9080503) q[3];
sx q[3];
rz(-0.84905784) q[3];
sx q[3];
rz(-2.8535709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0397772) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(0.56143108) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11557065) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(1.9157238) q[0];
rz(-1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(-1.3668758) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5850692) q[0];
sx q[0];
rz(-1.1856623) q[0];
sx q[0];
rz(-1.7760081) q[0];
rz(-1.6010124) q[2];
sx q[2];
rz(-0.82380166) q[2];
sx q[2];
rz(0.012416427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83134507) q[1];
sx q[1];
rz(-2.275895) q[1];
sx q[1];
rz(0.10465937) q[1];
x q[2];
rz(2.5451238) q[3];
sx q[3];
rz(-2.5554552) q[3];
sx q[3];
rz(0.49367192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2772969) q[2];
sx q[2];
rz(-2.0531451) q[2];
sx q[2];
rz(0.40536353) q[2];
rz(2.6799485) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(-1.595114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(2.6382085) q[0];
rz(-2.9227496) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(2.4898081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9572243) q[0];
sx q[0];
rz(-1.7608709) q[0];
sx q[0];
rz(0.24380542) q[0];
x q[1];
rz(-0.10105614) q[2];
sx q[2];
rz(-2.3172242) q[2];
sx q[2];
rz(-1.2254168) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8004605) q[1];
sx q[1];
rz(-2.9974077) q[1];
sx q[1];
rz(1.202281) q[1];
rz(-pi) q[2];
rz(1.6739453) q[3];
sx q[3];
rz(-1.2892937) q[3];
sx q[3];
rz(-2.3867949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51263222) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(-2.7499278) q[2];
rz(-2.9351249) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7713292) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(0.96965924) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-3.0113509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0129583) q[0];
sx q[0];
rz(-0.35612125) q[0];
sx q[0];
rz(0.14333368) q[0];
rz(-2.1580556) q[2];
sx q[2];
rz(-0.98896719) q[2];
sx q[2];
rz(0.26435095) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3245557) q[1];
sx q[1];
rz(-0.55700028) q[1];
sx q[1];
rz(0.42028285) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98826615) q[3];
sx q[3];
rz(-1.2417955) q[3];
sx q[3];
rz(0.44579166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51817259) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(-1.0858034) q[2];
rz(0.052224934) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(-1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645585) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(0.72558609) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(-2.3988147) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9966492) q[0];
sx q[0];
rz(-1.7631233) q[0];
sx q[0];
rz(0.74914519) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1002662) q[2];
sx q[2];
rz(-2.1770182) q[2];
sx q[2];
rz(-0.60339061) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85779491) q[1];
sx q[1];
rz(-0.89110121) q[1];
sx q[1];
rz(0.53417511) q[1];
x q[2];
rz(2.2804184) q[3];
sx q[3];
rz(-0.79813254) q[3];
sx q[3];
rz(2.8482311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(2.880704) q[2];
rz(2.0339113) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35659197) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(-0.93604952) q[0];
rz(-3.127457) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(2.4900808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441898) q[0];
sx q[0];
rz(-1.5751) q[0];
sx q[0];
rz(-1.547102) q[0];
rz(0.94349761) q[2];
sx q[2];
rz(-2.128696) q[2];
sx q[2];
rz(-0.055659143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1134909) q[1];
sx q[1];
rz(-0.3416225) q[1];
sx q[1];
rz(-1.2042868) q[1];
rz(-pi) q[2];
rz(-2.4581962) q[3];
sx q[3];
rz(-2.5960687) q[3];
sx q[3];
rz(-0.47513902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(0.30803099) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(-1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734252) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(0.47406468) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(-1.3778936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70625593) q[0];
sx q[0];
rz(-2.113111) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(2.4268742) q[2];
sx q[2];
rz(-0.93313365) q[2];
sx q[2];
rz(-1.4710466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4442515) q[1];
sx q[1];
rz(-1.3285944) q[1];
sx q[1];
rz(0.40476207) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7337448) q[3];
sx q[3];
rz(-1.019422) q[3];
sx q[3];
rz(-0.58386246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.6837998) q[3];
sx q[3];
rz(-1.0554353) q[3];
sx q[3];
rz(2.2916601) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67125852) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(-2.4709672) q[2];
sx q[2];
rz(-1.309633) q[2];
sx q[2];
rz(1.39967) q[2];
rz(2.408151) q[3];
sx q[3];
rz(-1.0167132) q[3];
sx q[3];
rz(0.50501962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
