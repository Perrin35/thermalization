OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5151514) q[0];
sx q[0];
rz(-0.03349537) q[0];
sx q[0];
rz(-1.3666231) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3596518) q[0];
sx q[0];
rz(-1.2255166) q[0];
sx q[0];
rz(0.79202534) q[0];
x q[1];
rz(1.4049454) q[2];
sx q[2];
rz(-1.7300786) q[2];
sx q[2];
rz(0.35370358) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72343091) q[1];
sx q[1];
rz(-1.9730113) q[1];
sx q[1];
rz(2.1869786) q[1];
x q[2];
rz(2.5868086) q[3];
sx q[3];
rz(-1.8150419) q[3];
sx q[3];
rz(2.3158405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2312317) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.988391) q[2];
rz(-2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(2.6385345) q[0];
rz(1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(0.15393004) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.126037) q[0];
sx q[0];
rz(-2.6385348) q[0];
sx q[0];
rz(1.5411099) q[0];
rz(-pi) q[1];
rz(-1.1459848) q[2];
sx q[2];
rz(-0.69485352) q[2];
sx q[2];
rz(1.7244463) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9736204) q[1];
sx q[1];
rz(-2.4929843) q[1];
sx q[1];
rz(-0.82778511) q[1];
rz(-pi) q[2];
rz(3.1400938) q[3];
sx q[3];
rz(-0.50977364) q[3];
sx q[3];
rz(1.203376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2543891) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(-1.3228234) q[2];
rz(1.4860738) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(-0.71050182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(-2.2316566) q[0];
rz(-2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(0.98532239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6353778) q[0];
sx q[0];
rz(-2.2197476) q[0];
sx q[0];
rz(0.20091591) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2276332) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(0.39011207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9039771) q[1];
sx q[1];
rz(-1.8198038) q[1];
sx q[1];
rz(1.6275089) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91089532) q[3];
sx q[3];
rz(-1.0064555) q[3];
sx q[3];
rz(-0.49163715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2788006) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(2.6990081) q[3];
sx q[3];
rz(-1.7740039) q[3];
sx q[3];
rz(2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(-0.57888794) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.3935864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907341) q[0];
sx q[0];
rz(-0.77461857) q[0];
sx q[0];
rz(-0.61268341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7384999) q[2];
sx q[2];
rz(-2.0003013) q[2];
sx q[2];
rz(1.6657366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.046603831) q[1];
sx q[1];
rz(-2.2294083) q[1];
sx q[1];
rz(2.1556426) q[1];
rz(2.3062069) q[3];
sx q[3];
rz(-1.3961892) q[3];
sx q[3];
rz(-1.1268827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(0.0021136443) q[2];
rz(-2.5801616) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11557065) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(1.9157238) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.7747169) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5850692) q[0];
sx q[0];
rz(-1.9559304) q[0];
sx q[0];
rz(1.7760081) q[0];
rz(-pi) q[1];
rz(1.5405802) q[2];
sx q[2];
rz(-0.82380166) q[2];
sx q[2];
rz(-3.1291762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3102476) q[1];
sx q[1];
rz(-0.86569769) q[1];
sx q[1];
rz(3.0369333) q[1];
rz(2.5451238) q[3];
sx q[3];
rz(-0.58613741) q[3];
sx q[3];
rz(-0.49367192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2772969) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(2.7362291) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(-1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49155238) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(-0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(-0.65178451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0359905) q[0];
sx q[0];
rz(-0.30797568) q[0];
sx q[0];
rz(-0.67291798) q[0];
rz(-1.4621554) q[2];
sx q[2];
rz(-2.3896653) q[2];
sx q[2];
rz(1.7679364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.3411322) q[1];
sx q[1];
rz(-2.9974077) q[1];
sx q[1];
rz(1.9393117) q[1];
x q[2];
rz(-1.4676474) q[3];
sx q[3];
rz(-1.2892937) q[3];
sx q[3];
rz(-2.3867949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51263222) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(-0.39166489) q[2];
rz(0.20646778) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-0.68968836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(-0.96965924) q[0];
rz(-2.5580653) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(0.13024174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286344) q[0];
sx q[0];
rz(-0.35612125) q[0];
sx q[0];
rz(2.998259) q[0];
rz(-pi) q[1];
rz(-2.1580556) q[2];
sx q[2];
rz(-2.1526255) q[2];
sx q[2];
rz(-0.26435095) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3245557) q[1];
sx q[1];
rz(-2.5845924) q[1];
sx q[1];
rz(0.42028285) q[1];
rz(0.98826615) q[3];
sx q[3];
rz(-1.8997972) q[3];
sx q[3];
rz(-0.44579166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(3.0893677) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(-1.0857371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6770342) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(-1.1664671) q[0];
rz(-0.72558609) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(0.74277791) q[1];
rz(-pi) q[2];
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
rz(-2.4648415) q[2];
sx q[2];
rz(-1.1427715) q[2];
sx q[2];
rz(1.2892436) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2837977) q[1];
sx q[1];
rz(-2.2504914) q[1];
sx q[1];
rz(2.6074175) q[1];
x q[2];
rz(0.86117427) q[3];
sx q[3];
rz(-0.79813254) q[3];
sx q[3];
rz(-2.8482311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(2.880704) q[2];
rz(1.1076814) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(-1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7850007) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(2.2055431) q[0];
rz(-0.014135663) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(-2.4900808) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441898) q[0];
sx q[0];
rz(-1.5664926) q[0];
sx q[0];
rz(1.5944907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.198095) q[2];
sx q[2];
rz(-2.128696) q[2];
sx q[2];
rz(0.055659143) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7266255) q[1];
sx q[1];
rz(-1.2526928) q[1];
sx q[1];
rz(0.12673641) q[1];
x q[2];
rz(2.4581962) q[3];
sx q[3];
rz(-2.5960687) q[3];
sx q[3];
rz(-2.6664536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2953879) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(2.6182168) q[2];
rz(2.8335617) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4073407) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(-1.403632) q[0];
rz(0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(-1.7636991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5546075) q[0];
sx q[0];
rz(-0.73903144) q[0];
sx q[0];
rz(2.2686195) q[0];
x q[1];
rz(-2.4268742) q[2];
sx q[2];
rz(-0.93313365) q[2];
sx q[2];
rz(1.4710466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6369032) q[1];
sx q[1];
rz(-2.6733589) q[1];
sx q[1];
rz(-2.5813028) q[1];
x q[2];
rz(-1.7337448) q[3];
sx q[3];
rz(-1.019422) q[3];
sx q[3];
rz(2.5577302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.750981) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.6101458) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(0.40614265) q[2];
sx q[2];
rz(-0.71229013) q[2];
sx q[2];
rz(-2.9980414) q[2];
rz(-0.73344161) q[3];
sx q[3];
rz(-1.0167132) q[3];
sx q[3];
rz(0.50501962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
