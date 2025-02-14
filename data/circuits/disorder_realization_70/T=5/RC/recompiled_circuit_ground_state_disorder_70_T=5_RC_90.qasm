OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8283451) q[0];
sx q[0];
rz(-0.77594835) q[0];
sx q[0];
rz(-1.1021855) q[0];
rz(1.6232396) q[1];
sx q[1];
rz(2.0145388) q[1];
sx q[1];
rz(9.1993499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878712) q[0];
sx q[0];
rz(-1.5086156) q[0];
sx q[0];
rz(2.4440701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61580832) q[2];
sx q[2];
rz(-1.4501418) q[2];
sx q[2];
rz(3.0349611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77535875) q[1];
sx q[1];
rz(-2.6817396) q[1];
sx q[1];
rz(-2.2103106) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59360154) q[3];
sx q[3];
rz(-2.5656325) q[3];
sx q[3];
rz(-2.0570448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48308358) q[2];
sx q[2];
rz(-0.67407125) q[2];
sx q[2];
rz(-3.1209514) q[2];
rz(2.9523197) q[3];
sx q[3];
rz(-2.7920189) q[3];
sx q[3];
rz(1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610158) q[0];
sx q[0];
rz(-0.78106946) q[0];
sx q[0];
rz(0.97907698) q[0];
rz(-3.0304404) q[1];
sx q[1];
rz(-2.4039098) q[1];
sx q[1];
rz(1.0609863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95288657) q[0];
sx q[0];
rz(-2.0659819) q[0];
sx q[0];
rz(1.0409036) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0041372) q[2];
sx q[2];
rz(-1.571621) q[2];
sx q[2];
rz(2.0440136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0909522) q[1];
sx q[1];
rz(-2.0555858) q[1];
sx q[1];
rz(-1.3371972) q[1];
rz(-2.8021566) q[3];
sx q[3];
rz(-0.41012479) q[3];
sx q[3];
rz(-2.8462178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7079118) q[2];
sx q[2];
rz(-1.3206864) q[2];
sx q[2];
rz(-2.743538) q[2];
rz(1.7178644) q[3];
sx q[3];
rz(-2.8545696) q[3];
sx q[3];
rz(2.6557693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240391) q[0];
sx q[0];
rz(-0.51152873) q[0];
sx q[0];
rz(-2.0892573) q[0];
rz(-2.6984093) q[1];
sx q[1];
rz(-2.1801703) q[1];
sx q[1];
rz(-2.4876432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0130413) q[0];
sx q[0];
rz(-1.6975901) q[0];
sx q[0];
rz(2.3953505) q[0];
rz(-pi) q[1];
rz(2.1939799) q[2];
sx q[2];
rz(-1.6296367) q[2];
sx q[2];
rz(0.72966444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7591651) q[1];
sx q[1];
rz(-0.8635206) q[1];
sx q[1];
rz(1.00818) q[1];
x q[2];
rz(1.3929358) q[3];
sx q[3];
rz(-1.7962958) q[3];
sx q[3];
rz(-2.5009843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7707278) q[2];
sx q[2];
rz(-0.76954049) q[2];
sx q[2];
rz(0.23180836) q[2];
rz(-1.9306463) q[3];
sx q[3];
rz(-1.959266) q[3];
sx q[3];
rz(0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4732707) q[0];
sx q[0];
rz(-2.7372161) q[0];
sx q[0];
rz(-2.7647198) q[0];
rz(-0.51141524) q[1];
sx q[1];
rz(-1.8835953) q[1];
sx q[1];
rz(0.84397856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91746861) q[0];
sx q[0];
rz(-0.11118764) q[0];
sx q[0];
rz(1.7286517) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2449712) q[2];
sx q[2];
rz(-2.5375536) q[2];
sx q[2];
rz(0.86261311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62629265) q[1];
sx q[1];
rz(-2.3808783) q[1];
sx q[1];
rz(-2.4079851) q[1];
rz(1.549349) q[3];
sx q[3];
rz(-2.0923695) q[3];
sx q[3];
rz(2.4873268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1692928) q[2];
sx q[2];
rz(-2.1268714) q[2];
sx q[2];
rz(2.443552) q[2];
rz(1.9418779) q[3];
sx q[3];
rz(-2.2359087) q[3];
sx q[3];
rz(2.5419295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19071628) q[0];
sx q[0];
rz(-0.27442014) q[0];
sx q[0];
rz(3.0158667) q[0];
rz(-1.0267286) q[1];
sx q[1];
rz(-1.3871565) q[1];
sx q[1];
rz(2.6153807) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77718098) q[0];
sx q[0];
rz(-0.48425774) q[0];
sx q[0];
rz(-0.15987349) q[0];
x q[1];
rz(-1.7650928) q[2];
sx q[2];
rz(-0.70617968) q[2];
sx q[2];
rz(-0.45748177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3359068) q[1];
sx q[1];
rz(-1.1912002) q[1];
sx q[1];
rz(1.6811045) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79546572) q[3];
sx q[3];
rz(-1.0364) q[3];
sx q[3];
rz(2.5837542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0937664) q[2];
sx q[2];
rz(-2.3036849) q[2];
sx q[2];
rz(-0.31025904) q[2];
rz(0.034916498) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(-2.3518899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652311) q[0];
sx q[0];
rz(-1.6595474) q[0];
sx q[0];
rz(0.37102997) q[0];
rz(3.0767483) q[1];
sx q[1];
rz(-2.3193181) q[1];
sx q[1];
rz(-1.8745905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3716301) q[0];
sx q[0];
rz(-0.93608674) q[0];
sx q[0];
rz(0.43031613) q[0];
x q[1];
rz(-2.7899988) q[2];
sx q[2];
rz(-0.97132713) q[2];
sx q[2];
rz(-2.4945187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6325841) q[1];
sx q[1];
rz(-1.9688089) q[1];
sx q[1];
rz(-1.5334315) q[1];
rz(-pi) q[2];
rz(-2.6755813) q[3];
sx q[3];
rz(-1.8309085) q[3];
sx q[3];
rz(2.8046908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0176257) q[2];
sx q[2];
rz(-1.4767246) q[2];
sx q[2];
rz(1.5383447) q[2];
rz(2.7247143) q[3];
sx q[3];
rz(-2.6346801) q[3];
sx q[3];
rz(0.1703593) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040319547) q[0];
sx q[0];
rz(-2.9077001) q[0];
sx q[0];
rz(2.7129569) q[0];
rz(-1.1242695) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(0.7695778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3316318) q[0];
sx q[0];
rz(-1.5726345) q[0];
sx q[0];
rz(1.5738945) q[0];
rz(-pi) q[1];
rz(2.6181302) q[2];
sx q[2];
rz(-2.0085196) q[2];
sx q[2];
rz(-1.4283534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10248549) q[1];
sx q[1];
rz(-1.277248) q[1];
sx q[1];
rz(2.5204646) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1355204) q[3];
sx q[3];
rz(-0.58593633) q[3];
sx q[3];
rz(-2.1028818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6577242) q[2];
sx q[2];
rz(-1.8572073) q[2];
sx q[2];
rz(-0.095495187) q[2];
rz(-0.5419845) q[3];
sx q[3];
rz(-0.46613765) q[3];
sx q[3];
rz(-0.12921648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1456881) q[0];
sx q[0];
rz(-2.2112084) q[0];
sx q[0];
rz(-3.0297739) q[0];
rz(-1.8226786) q[1];
sx q[1];
rz(-1.8967862) q[1];
sx q[1];
rz(1.8359312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625264) q[0];
sx q[0];
rz(-0.56076196) q[0];
sx q[0];
rz(1.3365585) q[0];
rz(1.4292923) q[2];
sx q[2];
rz(-1.8496528) q[2];
sx q[2];
rz(-0.57153801) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9247465) q[1];
sx q[1];
rz(-1.9439183) q[1];
sx q[1];
rz(1.2678964) q[1];
rz(2.5091845) q[3];
sx q[3];
rz(-2.494209) q[3];
sx q[3];
rz(1.6066911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32885113) q[2];
sx q[2];
rz(-1.988827) q[2];
sx q[2];
rz(1.1197155) q[2];
rz(-2.8730734) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(-1.1855116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28065228) q[0];
sx q[0];
rz(-2.379731) q[0];
sx q[0];
rz(-2.5326488) q[0];
rz(-2.2641585) q[1];
sx q[1];
rz(-1.0017706) q[1];
sx q[1];
rz(-0.61224365) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5800047) q[0];
sx q[0];
rz(-1.1682434) q[0];
sx q[0];
rz(2.3775565) q[0];
x q[1];
rz(1.2127146) q[2];
sx q[2];
rz(-1.0306669) q[2];
sx q[2];
rz(-3.0703406) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7944736) q[1];
sx q[1];
rz(-2.0277727) q[1];
sx q[1];
rz(0.73758482) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4928451) q[3];
sx q[3];
rz(-2.1105937) q[3];
sx q[3];
rz(-2.0920366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.073254243) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(0.18745984) q[2];
rz(-1.1120262) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(-2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336695) q[0];
sx q[0];
rz(-0.81403553) q[0];
sx q[0];
rz(-2.8975876) q[0];
rz(-1.4834652) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(2.3689178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940154) q[0];
sx q[0];
rz(-2.4227067) q[0];
sx q[0];
rz(0.80475828) q[0];
x q[1];
rz(-1.4336606) q[2];
sx q[2];
rz(-1.8923645) q[2];
sx q[2];
rz(-2.0069881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9909042) q[1];
sx q[1];
rz(-2.7953618) q[1];
sx q[1];
rz(-2.8105286) q[1];
x q[2];
rz(-2.3599125) q[3];
sx q[3];
rz(-2.1865851) q[3];
sx q[3];
rz(0.35889949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3385758) q[2];
sx q[2];
rz(-2.1803941) q[2];
sx q[2];
rz(0.72688603) q[2];
rz(2.0342942) q[3];
sx q[3];
rz(-0.50873435) q[3];
sx q[3];
rz(-1.0072964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0016639391) q[0];
sx q[0];
rz(-1.7353084) q[0];
sx q[0];
rz(-1.1575862) q[0];
rz(0.4460371) q[1];
sx q[1];
rz(-2.0560494) q[1];
sx q[1];
rz(2.5037419) q[1];
rz(1.4459417) q[2];
sx q[2];
rz(-1.1706252) q[2];
sx q[2];
rz(-1.3782383) q[2];
rz(2.7160326) q[3];
sx q[3];
rz(-1.4360518) q[3];
sx q[3];
rz(-0.19946972) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
