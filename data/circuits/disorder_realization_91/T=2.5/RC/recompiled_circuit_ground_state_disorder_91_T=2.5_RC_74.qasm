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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2006045) q[0];
sx q[0];
rz(-1.683569) q[0];
sx q[0];
rz(2.8975945) q[0];
rz(-2.8142848) q[2];
sx q[2];
rz(-0.94963726) q[2];
sx q[2];
rz(2.6747764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1937374) q[1];
sx q[1];
rz(-1.1754218) q[1];
sx q[1];
rz(-2.1389524) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0291489) q[3];
sx q[3];
rz(-2.3437269) q[3];
sx q[3];
rz(-3.0276379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51840034) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(-1.0830967) q[2];
rz(2.8862503) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(-1.9820836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.979368) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-0.033578385) q[0];
rz(-2.0032739) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(-0.23695645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9457696) q[0];
sx q[0];
rz(-1.3422215) q[0];
sx q[0];
rz(1.7128471) q[0];
rz(-2.9541624) q[2];
sx q[2];
rz(-2.0399562) q[2];
sx q[2];
rz(-2.3694866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8023876) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(1.7632381) q[1];
x q[2];
rz(-0.93795012) q[3];
sx q[3];
rz(-0.68162912) q[3];
sx q[3];
rz(-1.5860032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61243764) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(0.1203514) q[2];
rz(0.58593166) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(1.2683292) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2633857) q[0];
sx q[0];
rz(-2.1757941) q[0];
sx q[0];
rz(-0.98814386) q[0];
rz(-1.490961) q[1];
sx q[1];
rz(-0.71690503) q[1];
sx q[1];
rz(-1.5536701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7415888) q[0];
sx q[0];
rz(-0.015946139) q[0];
sx q[0];
rz(3.0724597) q[0];
x q[1];
rz(2.6695197) q[2];
sx q[2];
rz(-2.9909212) q[2];
sx q[2];
rz(-2.8674088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9389438) q[1];
sx q[1];
rz(-1.5420785) q[1];
sx q[1];
rz(-1.8320626) q[1];
x q[2];
rz(-2.5738945) q[3];
sx q[3];
rz(-1.2097219) q[3];
sx q[3];
rz(-2.2455824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3042018) q[2];
sx q[2];
rz(-1.3434255) q[2];
sx q[2];
rz(0.070579441) q[2];
rz(-0.48921674) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(-0.14744559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9551142) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(1.4008993) q[0];
rz(0.1700302) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(1.2331351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5197553) q[0];
sx q[0];
rz(-2.701169) q[0];
sx q[0];
rz(-1.8669403) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1523684) q[2];
sx q[2];
rz(-2.2007341) q[2];
sx q[2];
rz(3.1357855) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88573061) q[1];
sx q[1];
rz(-1.2911324) q[1];
sx q[1];
rz(0.015563029) q[1];
x q[2];
rz(3.1175851) q[3];
sx q[3];
rz(-2.6769961) q[3];
sx q[3];
rz(2.8516804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3418545) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(2.5119761) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(1.4153882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302309) q[0];
sx q[0];
rz(-2.6333599) q[0];
sx q[0];
rz(-0.10230219) q[0];
rz(1.6434068) q[1];
sx q[1];
rz(-1.3003277) q[1];
sx q[1];
rz(-2.7844875) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2997695) q[0];
sx q[0];
rz(-1.3477316) q[0];
sx q[0];
rz(-1.7905495) q[0];
rz(0.12320766) q[2];
sx q[2];
rz(-1.5438617) q[2];
sx q[2];
rz(1.3696485) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2374463) q[1];
sx q[1];
rz(-1.0654628) q[1];
sx q[1];
rz(1.9133486) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6980347) q[3];
sx q[3];
rz(-2.0152115) q[3];
sx q[3];
rz(2.9068974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47306481) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.3225887) q[2];
rz(-1.4332917) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(-2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54356164) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(-1.7033956) q[0];
rz(1.2403129) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(-1.3528489) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2738631) q[0];
sx q[0];
rz(-2.7244748) q[0];
sx q[0];
rz(1.8050844) q[0];
rz(-pi) q[1];
rz(2.8565494) q[2];
sx q[2];
rz(-0.54423287) q[2];
sx q[2];
rz(0.88819347) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50033113) q[1];
sx q[1];
rz(-1.2029543) q[1];
sx q[1];
rz(2.1180875) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8753136) q[3];
sx q[3];
rz(-2.7334573) q[3];
sx q[3];
rz(-0.88312393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9271586) q[2];
sx q[2];
rz(-2.7775601) q[2];
sx q[2];
rz(-2.7511609) q[2];
rz(1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(-2.3217679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696233) q[0];
sx q[0];
rz(-2.2177028) q[0];
sx q[0];
rz(-1.7286638) q[0];
rz(1.4631118) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(2.5020592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36632682) q[0];
sx q[0];
rz(-1.9444939) q[0];
sx q[0];
rz(-2.213272) q[0];
rz(-2.9908871) q[2];
sx q[2];
rz(-1.6455212) q[2];
sx q[2];
rz(0.13241235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8486316) q[1];
sx q[1];
rz(-2.1994414) q[1];
sx q[1];
rz(-1.2413003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9652548) q[3];
sx q[3];
rz(-1.4603851) q[3];
sx q[3];
rz(0.73964707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21706906) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(-0.45219839) q[2];
rz(2.9186115) q[3];
sx q[3];
rz(-2.860234) q[3];
sx q[3];
rz(-0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0105791) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(-0.16192326) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.2740998) q[1];
sx q[1];
rz(-0.23439342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1944626) q[0];
sx q[0];
rz(-0.85887733) q[0];
sx q[0];
rz(2.3361337) q[0];
x q[1];
rz(1.2599808) q[2];
sx q[2];
rz(-0.87500415) q[2];
sx q[2];
rz(-1.1281079) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.99410759) q[1];
sx q[1];
rz(-2.4723158) q[1];
sx q[1];
rz(-1.8753176) q[1];
x q[2];
rz(2.864847) q[3];
sx q[3];
rz(-1.5066772) q[3];
sx q[3];
rz(-1.8217877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28635412) q[2];
sx q[2];
rz(-1.0711292) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(-0.74328077) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(-2.4519517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2766992) q[0];
sx q[0];
rz(-1.6962637) q[0];
sx q[0];
rz(-1.0869166) q[0];
rz(-1.9089606) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(1.8392275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83643276) q[0];
sx q[0];
rz(-1.4576685) q[0];
sx q[0];
rz(-2.9311084) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97424284) q[2];
sx q[2];
rz(-1.2737717) q[2];
sx q[2];
rz(-3.0926306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2757146) q[1];
sx q[1];
rz(-0.68531407) q[1];
sx q[1];
rz(2.4908743) q[1];
x q[2];
rz(1.5897254) q[3];
sx q[3];
rz(-2.193748) q[3];
sx q[3];
rz(-1.2894693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46633259) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(-2.7272398) q[2];
rz(-0.77053344) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-2.2079302) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82938021) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(0.46863753) q[0];
rz(-1.8679484) q[1];
sx q[1];
rz(-1.8654774) q[1];
sx q[1];
rz(-0.31148568) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5632919) q[0];
sx q[0];
rz(-1.5130763) q[0];
sx q[0];
rz(2.2895534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6577254) q[2];
sx q[2];
rz(-1.6875899) q[2];
sx q[2];
rz(2.180336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35986082) q[1];
sx q[1];
rz(-1.5241429) q[1];
sx q[1];
rz(-2.0585039) q[1];
x q[2];
rz(-0.52759513) q[3];
sx q[3];
rz(-2.0012534) q[3];
sx q[3];
rz(-0.93387077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(1.3930456) q[2];
rz(0.92750183) q[3];
sx q[3];
rz(-1.564097) q[3];
sx q[3];
rz(2.2813796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504234) q[0];
sx q[0];
rz(-0.32857729) q[0];
sx q[0];
rz(1.487442) q[0];
rz(0.18118478) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(-1.8306611) q[2];
sx q[2];
rz(-0.54878546) q[2];
sx q[2];
rz(2.6981163) q[2];
rz(2.2929706) q[3];
sx q[3];
rz(-0.80905882) q[3];
sx q[3];
rz(-0.73425135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
