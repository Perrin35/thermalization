OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2663015) q[0];
sx q[0];
rz(-0.40691352) q[0];
sx q[0];
rz(-2.9391675) q[0];
rz(-1.3099194) q[1];
sx q[1];
rz(-1.1453495) q[1];
sx q[1];
rz(0.98869148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3770954) q[0];
sx q[0];
rz(-2.3125153) q[0];
sx q[0];
rz(-0.25804706) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0634544) q[2];
sx q[2];
rz(-1.4333908) q[2];
sx q[2];
rz(-1.5619141) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0571756) q[1];
sx q[1];
rz(-1.7965846) q[1];
sx q[1];
rz(0.94364072) q[1];
rz(-pi) q[2];
rz(-1.3019788) q[3];
sx q[3];
rz(-0.43987396) q[3];
sx q[3];
rz(-1.1209821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6053091) q[2];
sx q[2];
rz(-0.80159694) q[2];
sx q[2];
rz(1.0514222) q[2];
rz(1.6384404) q[3];
sx q[3];
rz(-1.7283311) q[3];
sx q[3];
rz(0.2352636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1992092) q[0];
sx q[0];
rz(-1.0301882) q[0];
sx q[0];
rz(-0.27401608) q[0];
rz(2.7577397) q[1];
sx q[1];
rz(-2.5089788) q[1];
sx q[1];
rz(-1.3207818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9250671) q[0];
sx q[0];
rz(-1.4122047) q[0];
sx q[0];
rz(-3.0761883) q[0];
rz(-pi) q[1];
rz(-1.8019559) q[2];
sx q[2];
rz(-2.0589239) q[2];
sx q[2];
rz(1.3991935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1283001) q[1];
sx q[1];
rz(-1.6183524) q[1];
sx q[1];
rz(2.2940728) q[1];
x q[2];
rz(0.78944309) q[3];
sx q[3];
rz(-2.200921) q[3];
sx q[3];
rz(0.9073782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5447834) q[2];
sx q[2];
rz(-1.7388672) q[2];
sx q[2];
rz(1.7421494) q[2];
rz(0.60747373) q[3];
sx q[3];
rz(-1.6739269) q[3];
sx q[3];
rz(0.39052159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52943277) q[0];
sx q[0];
rz(-2.2859892) q[0];
sx q[0];
rz(0.31923077) q[0];
rz(-3.0901129) q[1];
sx q[1];
rz(-1.9620644) q[1];
sx q[1];
rz(2.0850339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7565472) q[0];
sx q[0];
rz(-0.17160417) q[0];
sx q[0];
rz(1.8770939) q[0];
rz(-pi) q[1];
rz(-1.5139903) q[2];
sx q[2];
rz(-0.70630766) q[2];
sx q[2];
rz(2.9212223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33256862) q[1];
sx q[1];
rz(-1.3961375) q[1];
sx q[1];
rz(0.93775429) q[1];
rz(3.0491054) q[3];
sx q[3];
rz(-1.1690946) q[3];
sx q[3];
rz(0.71850151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9487379) q[2];
sx q[2];
rz(-1.4317908) q[2];
sx q[2];
rz(1.3445492) q[2];
rz(0.91790849) q[3];
sx q[3];
rz(-0.63799208) q[3];
sx q[3];
rz(-1.8857672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767839) q[0];
sx q[0];
rz(-0.5905686) q[0];
sx q[0];
rz(-1.6326686) q[0];
rz(-2.6347939) q[1];
sx q[1];
rz(-1.9281887) q[1];
sx q[1];
rz(-0.39055821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9386425) q[0];
sx q[0];
rz(-3.1331867) q[0];
sx q[0];
rz(2.6132843) q[0];
rz(-pi) q[1];
rz(2.7616169) q[2];
sx q[2];
rz(-1.7479428) q[2];
sx q[2];
rz(0.65922696) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4949172) q[1];
sx q[1];
rz(-2.5220089) q[1];
sx q[1];
rz(1.5219206) q[1];
x q[2];
rz(-2.6253013) q[3];
sx q[3];
rz(-0.18169475) q[3];
sx q[3];
rz(-2.5207375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.249923) q[2];
sx q[2];
rz(-2.2107783) q[2];
sx q[2];
rz(-0.8963975) q[2];
rz(-2.2253288) q[3];
sx q[3];
rz(-2.2057585) q[3];
sx q[3];
rz(-0.96074444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5498098) q[0];
sx q[0];
rz(-2.7857605) q[0];
sx q[0];
rz(0.47019666) q[0];
rz(1.7377986) q[1];
sx q[1];
rz(-2.4411968) q[1];
sx q[1];
rz(1.2276924) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.550023) q[0];
sx q[0];
rz(-1.8409022) q[0];
sx q[0];
rz(-0.0065992532) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9982905) q[2];
sx q[2];
rz(-1.059747) q[2];
sx q[2];
rz(-1.4882027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44547651) q[1];
sx q[1];
rz(-1.646478) q[1];
sx q[1];
rz(-2.2161525) q[1];
x q[2];
rz(0.49318243) q[3];
sx q[3];
rz(-2.5458284) q[3];
sx q[3];
rz(-0.34957928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1710743) q[2];
sx q[2];
rz(-0.71275622) q[2];
sx q[2];
rz(-1.7573382) q[2];
rz(-2.5303043) q[3];
sx q[3];
rz(-1.3795815) q[3];
sx q[3];
rz(-2.7454564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.1432081) q[0];
sx q[0];
rz(-1.0118326) q[0];
sx q[0];
rz(-2.5894453) q[0];
rz(2.4119008) q[1];
sx q[1];
rz(-2.153986) q[1];
sx q[1];
rz(2.139835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.998913) q[0];
sx q[0];
rz(-2.985654) q[0];
sx q[0];
rz(-2.1151311) q[0];
rz(-0.92879734) q[2];
sx q[2];
rz(-1.1603519) q[2];
sx q[2];
rz(2.9985119) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3730091) q[1];
sx q[1];
rz(-0.78861744) q[1];
sx q[1];
rz(-2.8810701) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4625307) q[3];
sx q[3];
rz(-1.9894532) q[3];
sx q[3];
rz(-2.1931894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45727649) q[2];
sx q[2];
rz(-1.7924954) q[2];
sx q[2];
rz(-1.7982193) q[2];
rz(2.4172879) q[3];
sx q[3];
rz(-1.2035921) q[3];
sx q[3];
rz(-0.30266416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2398719) q[0];
sx q[0];
rz(-1.8093103) q[0];
sx q[0];
rz(-2.4161762) q[0];
rz(1.2984917) q[1];
sx q[1];
rz(-0.45227554) q[1];
sx q[1];
rz(0.27870146) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50510079) q[0];
sx q[0];
rz(-1.0351702) q[0];
sx q[0];
rz(1.3938175) q[0];
rz(-2.2594035) q[2];
sx q[2];
rz(-0.49748245) q[2];
sx q[2];
rz(-2.62219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1126453) q[1];
sx q[1];
rz(-0.7794081) q[1];
sx q[1];
rz(-2.0112627) q[1];
x q[2];
rz(-1.0213357) q[3];
sx q[3];
rz(-0.82290339) q[3];
sx q[3];
rz(-1.2513892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7574888) q[2];
sx q[2];
rz(-2.8748547) q[2];
sx q[2];
rz(2.1534446) q[2];
rz(0.10495505) q[3];
sx q[3];
rz(-1.8433439) q[3];
sx q[3];
rz(2.9102563) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4526378) q[0];
sx q[0];
rz(-0.80113688) q[0];
sx q[0];
rz(-0.6193921) q[0];
rz(3.1267005) q[1];
sx q[1];
rz(-2.2512524) q[1];
sx q[1];
rz(2.4598918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3103763) q[0];
sx q[0];
rz(-1.8583603) q[0];
sx q[0];
rz(2.1307261) q[0];
x q[1];
rz(3.0613941) q[2];
sx q[2];
rz(-2.2343582) q[2];
sx q[2];
rz(-0.71127438) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8047844) q[1];
sx q[1];
rz(-0.35165916) q[1];
sx q[1];
rz(-3.1218887) q[1];
rz(-pi) q[2];
rz(-2.3084749) q[3];
sx q[3];
rz(-1.7135812) q[3];
sx q[3];
rz(-2.8974183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90998489) q[2];
sx q[2];
rz(-2.3384428) q[2];
sx q[2];
rz(-0.33965674) q[2];
rz(2.5275982) q[3];
sx q[3];
rz(-1.7620148) q[3];
sx q[3];
rz(-1.5617255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.9851819) q[0];
sx q[0];
rz(-2.7934533) q[0];
sx q[0];
rz(-0.85439318) q[0];
rz(-0.59335452) q[1];
sx q[1];
rz(-1.0320458) q[1];
sx q[1];
rz(-2.8462483) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4446553) q[0];
sx q[0];
rz(-1.5155223) q[0];
sx q[0];
rz(1.6037206) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6215084) q[2];
sx q[2];
rz(-0.4182294) q[2];
sx q[2];
rz(-0.39722463) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1258537) q[1];
sx q[1];
rz(-2.2496114) q[1];
sx q[1];
rz(0.62894459) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2611102) q[3];
sx q[3];
rz(-1.6777552) q[3];
sx q[3];
rz(0.84236077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16193834) q[2];
sx q[2];
rz(-1.6831968) q[2];
sx q[2];
rz(-2.4597607) q[2];
rz(-1.2006867) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(-0.43584263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0308663) q[0];
sx q[0];
rz(-1.8616572) q[0];
sx q[0];
rz(-2.6357292) q[0];
rz(2.1234296) q[1];
sx q[1];
rz(-1.4351427) q[1];
sx q[1];
rz(0.86046576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8296536) q[0];
sx q[0];
rz(-1.5235434) q[0];
sx q[0];
rz(-0.67619369) q[0];
rz(-pi) q[1];
rz(2.0138836) q[2];
sx q[2];
rz(-0.32226105) q[2];
sx q[2];
rz(-1.6784422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0588433) q[1];
sx q[1];
rz(-0.25139055) q[1];
sx q[1];
rz(-3.0091834) q[1];
rz(-2.5768336) q[3];
sx q[3];
rz(-1.0339206) q[3];
sx q[3];
rz(0.91954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0818417) q[2];
sx q[2];
rz(-1.1885208) q[2];
sx q[2];
rz(-0.58319485) q[2];
rz(0.0010631337) q[3];
sx q[3];
rz(-0.93102264) q[3];
sx q[3];
rz(-0.49334905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6625593) q[0];
sx q[0];
rz(-1.4190577) q[0];
sx q[0];
rz(0.90202913) q[0];
rz(-2.5262911) q[1];
sx q[1];
rz(-1.6533783) q[1];
sx q[1];
rz(1.3396214) q[1];
rz(2.6189005) q[2];
sx q[2];
rz(-1.3466071) q[2];
sx q[2];
rz(-1.9853332) q[2];
rz(-1.2801805) q[3];
sx q[3];
rz(-0.59528412) q[3];
sx q[3];
rz(-3.0045454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
