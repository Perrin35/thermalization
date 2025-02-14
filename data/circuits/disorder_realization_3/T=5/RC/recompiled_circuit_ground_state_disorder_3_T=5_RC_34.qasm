OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(-0.29932061) q[0];
sx q[0];
rz(-2.646995) q[0];
rz(1.142113) q[1];
sx q[1];
rz(-1.0057058) q[1];
sx q[1];
rz(-2.0118654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0725192) q[0];
sx q[0];
rz(-1.7290218) q[0];
sx q[0];
rz(0.67969602) q[0];
x q[1];
rz(-0.45290516) q[2];
sx q[2];
rz(-1.1158841) q[2];
sx q[2];
rz(0.95547966) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3898824) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(2.8044279) q[1];
x q[2];
rz(0.42721984) q[3];
sx q[3];
rz(-0.96517206) q[3];
sx q[3];
rz(-0.96715121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8226681) q[2];
sx q[2];
rz(-1.8225887) q[2];
sx q[2];
rz(-2.2887716) q[2];
rz(-1.3301814) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3332719) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(1.464123) q[0];
rz(-1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(-1.2082072) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0776652) q[0];
sx q[0];
rz(-2.0443235) q[0];
sx q[0];
rz(-2.5350476) q[0];
x q[1];
rz(-2.9242427) q[2];
sx q[2];
rz(-2.0361379) q[2];
sx q[2];
rz(0.67721043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.959001) q[1];
sx q[1];
rz(-1.4495069) q[1];
sx q[1];
rz(1.3022997) q[1];
rz(-2.6868846) q[3];
sx q[3];
rz(-1.7268001) q[3];
sx q[3];
rz(-1.7059513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53604424) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(-1.6571244) q[2];
rz(-0.015497192) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(-1.9626455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.991796) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-2.8699744) q[0];
rz(2.244921) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(-1.2976049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0947726) q[0];
sx q[0];
rz(-1.5607139) q[0];
sx q[0];
rz(-2.0338414) q[0];
rz(1.9515368) q[2];
sx q[2];
rz(-1.0972766) q[2];
sx q[2];
rz(1.2981594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5813959) q[1];
sx q[1];
rz(-2.8688258) q[1];
sx q[1];
rz(-1.1352989) q[1];
rz(-pi) q[2];
rz(-3.1362757) q[3];
sx q[3];
rz(-2.379619) q[3];
sx q[3];
rz(-0.055082037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74035949) q[2];
sx q[2];
rz(-2.3637502) q[2];
sx q[2];
rz(3.0878301) q[2];
rz(1.2939804) q[3];
sx q[3];
rz(-2.3432178) q[3];
sx q[3];
rz(-2.9062041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.9409222) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(-1.6773552) q[0];
rz(1.1306521) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(-3.0272223) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0972041) q[0];
sx q[0];
rz(-2.1425793) q[0];
sx q[0];
rz(-1.0685789) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0691772) q[2];
sx q[2];
rz(-0.97626462) q[2];
sx q[2];
rz(-1.6785113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3721766) q[1];
sx q[1];
rz(-1.5443364) q[1];
sx q[1];
rz(-2.2464804) q[1];
rz(1.5642197) q[3];
sx q[3];
rz(-1.7688171) q[3];
sx q[3];
rz(-1.1634367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4734681) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(0.075695666) q[2];
rz(-0.43478742) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39744034) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(2.721526) q[0];
rz(-0.21332598) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(1.8992281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5899701) q[0];
sx q[0];
rz(-1.9931355) q[0];
sx q[0];
rz(2.3809654) q[0];
x q[1];
rz(-0.43102805) q[2];
sx q[2];
rz(-1.8371474) q[2];
sx q[2];
rz(2.9961627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8458357) q[1];
sx q[1];
rz(-1.6713665) q[1];
sx q[1];
rz(-0.015458903) q[1];
rz(-1.5228988) q[3];
sx q[3];
rz(-1.7918799) q[3];
sx q[3];
rz(0.97122279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8215948) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(0.37459174) q[2];
rz(1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(-0.64468002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11151611) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(-0.43055713) q[0];
rz(-0.14398362) q[1];
sx q[1];
rz(-2.0591683) q[1];
sx q[1];
rz(0.63124257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66723204) q[0];
sx q[0];
rz(-2.1788414) q[0];
sx q[0];
rz(-1.6365746) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14603931) q[2];
sx q[2];
rz(-1.1542873) q[2];
sx q[2];
rz(-3.0768968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5938607) q[1];
sx q[1];
rz(-1.0775078) q[1];
sx q[1];
rz(-2.4592722) q[1];
rz(-pi) q[2];
rz(2.1926375) q[3];
sx q[3];
rz(-2.6027711) q[3];
sx q[3];
rz(-0.99297374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97079078) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(1.4405174) q[2];
rz(1.7801646) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(-2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.454527) q[0];
sx q[0];
rz(-2.8519958) q[0];
sx q[0];
rz(-0.92426306) q[0];
rz(-0.29280064) q[1];
sx q[1];
rz(-1.9127138) q[1];
sx q[1];
rz(-0.80088314) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7445114) q[0];
sx q[0];
rz(-1.1985221) q[0];
sx q[0];
rz(-0.27639322) q[0];
rz(-pi) q[1];
rz(2.2325781) q[2];
sx q[2];
rz(-1.2940727) q[2];
sx q[2];
rz(-0.02034517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.427721) q[1];
sx q[1];
rz(-1.5414951) q[1];
sx q[1];
rz(-0.10515611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5436567) q[3];
sx q[3];
rz(-1.5487897) q[3];
sx q[3];
rz(-2.1834971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31876365) q[2];
sx q[2];
rz(-1.9209361) q[2];
sx q[2];
rz(1.8611543) q[2];
rz(2.473623) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(-2.7282696) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845344) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(-1.1827693) q[0];
rz(0.23712748) q[1];
sx q[1];
rz(-0.10219899) q[1];
sx q[1];
rz(3.0822486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5049669) q[0];
sx q[0];
rz(-1.4511746) q[0];
sx q[0];
rz(-1.2960394) q[0];
rz(2.0996712) q[2];
sx q[2];
rz(-2.2239074) q[2];
sx q[2];
rz(-2.3082341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0036111) q[1];
sx q[1];
rz(-1.4661769) q[1];
sx q[1];
rz(1.8358747) q[1];
x q[2];
rz(-1.343325) q[3];
sx q[3];
rz(-0.331628) q[3];
sx q[3];
rz(-2.22081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55988971) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(1.3524559) q[2];
rz(1.3672359) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(-0.8555612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.2818114) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(2.6819041) q[0];
rz(0.028060878) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(-1.9557767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494649) q[0];
sx q[0];
rz(-0.26927265) q[0];
sx q[0];
rz(2.9721391) q[0];
rz(-0.90699754) q[2];
sx q[2];
rz(-3.0131648) q[2];
sx q[2];
rz(-0.80824404) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2433127) q[1];
sx q[1];
rz(-2.4377258) q[1];
sx q[1];
rz(2.3614592) q[1];
rz(-pi) q[2];
rz(2.7544153) q[3];
sx q[3];
rz(-0.65473352) q[3];
sx q[3];
rz(-2.2019405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(0.15677491) q[2];
rz(-2.5458941) q[3];
sx q[3];
rz(-2.7955293) q[3];
sx q[3];
rz(-0.025207635) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51782411) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(0.18381707) q[0];
rz(-0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-2.877291) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1109836) q[0];
sx q[0];
rz(-1.302659) q[0];
sx q[0];
rz(-2.5546352) q[0];
rz(-pi) q[1];
rz(-1.043522) q[2];
sx q[2];
rz(-0.95607483) q[2];
sx q[2];
rz(1.6312903) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7122927) q[1];
sx q[1];
rz(-1.7279062) q[1];
sx q[1];
rz(0.18204851) q[1];
rz(2.6464858) q[3];
sx q[3];
rz(-1.6551842) q[3];
sx q[3];
rz(1.6364678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(-1.7372355) q[3];
sx q[3];
rz(-0.93969932) q[3];
sx q[3];
rz(-3.0586045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(0.71113853) q[1];
sx q[1];
rz(-1.8322721) q[1];
sx q[1];
rz(-2.4087404) q[1];
rz(1.0719094) q[2];
sx q[2];
rz(-1.6098534) q[2];
sx q[2];
rz(0.44865566) q[2];
rz(2.5855999) q[3];
sx q[3];
rz(-1.0992194) q[3];
sx q[3];
rz(1.9627375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
