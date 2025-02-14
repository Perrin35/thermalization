OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8102326) q[0];
sx q[0];
rz(4.3963764) q[0];
sx q[0];
rz(8.7788361) q[0];
rz(1.6826001) q[1];
sx q[1];
rz(3.269722) q[1];
sx q[1];
rz(10.092957) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2006045) q[0];
sx q[0];
rz(-1.683569) q[0];
sx q[0];
rz(2.8975945) q[0];
rz(-pi) q[1];
rz(-1.1485841) q[2];
sx q[2];
rz(-0.6919043) q[2];
sx q[2];
rz(0.9949323) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5234452) q[1];
sx q[1];
rz(-2.0904358) q[1];
sx q[1];
rz(-2.6818399) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1124437) q[3];
sx q[3];
rz(-2.3437269) q[3];
sx q[3];
rz(3.0276379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51840034) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(-1.0830967) q[2];
rz(0.25534233) q[3];
sx q[3];
rz(-2.2446003) q[3];
sx q[3];
rz(-1.9820836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.979368) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-3.1080143) q[0];
rz(1.1383188) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(2.9046362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19582307) q[0];
sx q[0];
rz(-1.3422215) q[0];
sx q[0];
rz(1.4287455) q[0];
rz(-pi) q[1];
x q[1];
rz(1.218538) q[2];
sx q[2];
rz(-0.50261231) q[2];
sx q[2];
rz(1.9723122) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8023876) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(1.7632381) q[1];
x q[2];
rz(2.6941766) q[3];
sx q[3];
rz(-1.0378812) q[3];
sx q[3];
rz(0.79870236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.529155) q[2];
sx q[2];
rz(-1.0470904) q[2];
sx q[2];
rz(0.1203514) q[2];
rz(-2.555661) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(-1.8732635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2633857) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(-2.1534488) q[0];
rz(1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(1.5879226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7415888) q[0];
sx q[0];
rz(-0.015946139) q[0];
sx q[0];
rz(0.069132968) q[0];
x q[1];
rz(1.5018671) q[2];
sx q[2];
rz(-1.7048827) q[2];
sx q[2];
rz(-2.9390735) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37582477) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(0.029726083) q[1];
rz(0.56769811) q[3];
sx q[3];
rz(-1.9318707) q[3];
sx q[3];
rz(-0.89601024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83739088) q[2];
sx q[2];
rz(-1.7981671) q[2];
sx q[2];
rz(-3.0710132) q[2];
rz(0.48921674) q[3];
sx q[3];
rz(-0.990812) q[3];
sx q[3];
rz(2.9941471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9551142) q[0];
sx q[0];
rz(-0.69545737) q[0];
sx q[0];
rz(-1.7406933) q[0];
rz(-2.9715624) q[1];
sx q[1];
rz(-1.731512) q[1];
sx q[1];
rz(1.9084575) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5197553) q[0];
sx q[0];
rz(-0.4404236) q[0];
sx q[0];
rz(1.2746524) q[0];
rz(2.4242086) q[2];
sx q[2];
rz(-1.1108494) q[2];
sx q[2];
rz(-1.9460974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4608232) q[1];
sx q[1];
rz(-1.555838) q[1];
sx q[1];
rz(-1.8504924) q[1];
x q[2];
rz(-1.5587658) q[3];
sx q[3];
rz(-2.0352484) q[3];
sx q[3];
rz(2.8248276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79973811) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(-2.5119761) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-2.5119731) q[3];
sx q[3];
rz(-1.4153882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.28302309) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(0.10230219) q[0];
rz(1.6434068) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(2.7844875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84182318) q[0];
sx q[0];
rz(-1.7938611) q[0];
sx q[0];
rz(-1.3510432) q[0];
x q[1];
rz(-1.5979366) q[2];
sx q[2];
rz(-1.4476336) q[2];
sx q[2];
rz(2.9371098) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6372924) q[1];
sx q[1];
rz(-1.8691113) q[1];
sx q[1];
rz(0.53108414) q[1];
x q[2];
rz(-2.6980347) q[3];
sx q[3];
rz(-2.0152115) q[3];
sx q[3];
rz(0.23469521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47306481) q[2];
sx q[2];
rz(-1.4184971) q[2];
sx q[2];
rz(-1.3225887) q[2];
rz(1.4332917) q[3];
sx q[3];
rz(-1.3168443) q[3];
sx q[3];
rz(-2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54356164) q[0];
sx q[0];
rz(-0.99995166) q[0];
sx q[0];
rz(-1.438197) q[0];
rz(1.9012798) q[1];
sx q[1];
rz(-2.4867609) q[1];
sx q[1];
rz(1.7887438) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677296) q[0];
sx q[0];
rz(-0.41711787) q[0];
sx q[0];
rz(-1.3365082) q[0];
rz(-pi) q[1];
rz(1.7393624) q[2];
sx q[2];
rz(-2.0907846) q[2];
sx q[2];
rz(-1.2181768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8554012) q[1];
sx q[1];
rz(-2.0778065) q[1];
sx q[1];
rz(-0.42393522) q[1];
rz(-pi) q[2];
rz(-0.39522533) q[3];
sx q[3];
rz(-1.4661643) q[3];
sx q[3];
rz(-0.44236174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.214434) q[2];
sx q[2];
rz(-2.7775601) q[2];
sx q[2];
rz(-0.39043179) q[2];
rz(-1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(2.3217679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696233) q[0];
sx q[0];
rz(-2.2177028) q[0];
sx q[0];
rz(-1.7286638) q[0];
rz(-1.4631118) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(0.63953343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4832925) q[0];
sx q[0];
rz(-2.4118703) q[0];
sx q[0];
rz(-0.99131063) q[0];
rz(-pi) q[1];
rz(-0.15070559) q[2];
sx q[2];
rz(-1.4960714) q[2];
sx q[2];
rz(-3.0091803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4762791) q[1];
sx q[1];
rz(-1.3060044) q[1];
sx q[1];
rz(2.4864343) q[1];
x q[2];
rz(0.11951167) q[3];
sx q[3];
rz(-1.9627213) q[3];
sx q[3];
rz(-2.3562795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9245236) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(-2.6893943) q[2];
rz(0.2229812) q[3];
sx q[3];
rz(-2.860234) q[3];
sx q[3];
rz(-2.9862459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13101354) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(0.16192326) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(-2.9071992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1680555) q[0];
sx q[0];
rz(-0.99322766) q[0];
sx q[0];
rz(2.4650991) q[0];
rz(-pi) q[1];
rz(-1.8816119) q[2];
sx q[2];
rz(-2.2665885) q[2];
sx q[2];
rz(-2.0134848) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.80658) q[1];
sx q[1];
rz(-1.7579105) q[1];
sx q[1];
rz(-0.92430618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27674562) q[3];
sx q[3];
rz(-1.5066772) q[3];
sx q[3];
rz(1.319805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-1.0711292) q[2];
sx q[2];
rz(2.5174649) q[2];
rz(2.3983119) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(-0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2766992) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(-1.0869166) q[0];
rz(-1.2326321) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(-1.3023652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75847699) q[0];
sx q[0];
rz(-1.7799151) q[0];
sx q[0];
rz(1.6864548) q[0];
x q[1];
rz(-2.7872269) q[2];
sx q[2];
rz(-2.1379037) q[2];
sx q[2];
rz(1.3256377) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2375396) q[1];
sx q[1];
rz(-1.1773279) q[1];
sx q[1];
rz(-2.5649125) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1152505) q[3];
sx q[3];
rz(-2.5183916) q[3];
sx q[3];
rz(-1.88456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46633259) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(-0.41435286) q[2];
rz(-2.3710592) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122124) q[0];
sx q[0];
rz(-0.82013622) q[0];
sx q[0];
rz(0.46863753) q[0];
rz(1.8679484) q[1];
sx q[1];
rz(-1.8654774) q[1];
sx q[1];
rz(-2.830107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0429223) q[0];
sx q[0];
rz(-2.2880974) q[0];
sx q[0];
rz(3.0649661) q[0];
rz(-pi) q[1];
rz(-0.63705541) q[2];
sx q[2];
rz(-2.9961176) q[2];
sx q[2];
rz(1.6033974) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1232417) q[1];
sx q[1];
rz(-2.6518378) q[1];
sx q[1];
rz(1.4714929) q[1];
x q[2];
rz(2.4021637) q[3];
sx q[3];
rz(-0.66777705) q[3];
sx q[3];
rz(1.8831933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0353388) q[2];
sx q[2];
rz(-1.9474578) q[2];
sx q[2];
rz(-1.3930456) q[2];
rz(2.2140908) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(2.2813796) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911693) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-2.9604079) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(-2.1045024) q[2];
sx q[2];
rz(-1.4363534) q[2];
sx q[2];
rz(-1.7912122) q[2];
rz(2.2374033) q[3];
sx q[3];
rz(-1.0720357) q[3];
sx q[3];
rz(-2.8513249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
