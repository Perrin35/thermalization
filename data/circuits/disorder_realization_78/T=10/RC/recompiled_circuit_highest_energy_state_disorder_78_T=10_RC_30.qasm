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
rz(1.3820833) q[0];
sx q[0];
rz(-0.74892646) q[0];
sx q[0];
rz(1.7480667) q[0];
rz(0.8028318) q[1];
sx q[1];
rz(-0.79610151) q[1];
sx q[1];
rz(0.5308477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323468) q[0];
sx q[0];
rz(-2.0589925) q[0];
sx q[0];
rz(2.7335579) q[0];
rz(3.0391998) q[2];
sx q[2];
rz(-2.9739485) q[2];
sx q[2];
rz(-1.8554115) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72076529) q[1];
sx q[1];
rz(-1.4958352) q[1];
sx q[1];
rz(-2.1032745) q[1];
x q[2];
rz(-2.0898607) q[3];
sx q[3];
rz(-1.4822019) q[3];
sx q[3];
rz(-0.54852329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.406245) q[2];
sx q[2];
rz(-0.9333868) q[2];
sx q[2];
rz(-2.1437342) q[2];
rz(-2.2570299) q[3];
sx q[3];
rz(-1.9622784) q[3];
sx q[3];
rz(-0.71200371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74580055) q[0];
sx q[0];
rz(-2.6048248) q[0];
sx q[0];
rz(-2.0641548) q[0];
rz(2.5224345) q[1];
sx q[1];
rz(-1.2666603) q[1];
sx q[1];
rz(-1.223863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1503939) q[0];
sx q[0];
rz(-1.7454892) q[0];
sx q[0];
rz(-1.2758486) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5406088) q[2];
sx q[2];
rz(-1.5689625) q[2];
sx q[2];
rz(0.48170127) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89586414) q[1];
sx q[1];
rz(-2.1460377) q[1];
sx q[1];
rz(1.0771345) q[1];
rz(0.72325752) q[3];
sx q[3];
rz(-0.92953909) q[3];
sx q[3];
rz(-2.0632405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9549442) q[2];
sx q[2];
rz(-2.2961605) q[2];
sx q[2];
rz(2.1050982) q[2];
rz(-1.9519818) q[3];
sx q[3];
rz(-1.2396953) q[3];
sx q[3];
rz(-0.069124393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.67967296) q[0];
sx q[0];
rz(-0.24149495) q[0];
sx q[0];
rz(1.6756206) q[0];
rz(-1.2566603) q[1];
sx q[1];
rz(-0.64777056) q[1];
sx q[1];
rz(0.57674903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38274469) q[0];
sx q[0];
rz(-0.85475105) q[0];
sx q[0];
rz(2.3780253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5336214) q[2];
sx q[2];
rz(-1.1122744) q[2];
sx q[2];
rz(2.9165845) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8654792) q[1];
sx q[1];
rz(-0.41655585) q[1];
sx q[1];
rz(1.4426484) q[1];
rz(-pi) q[2];
rz(1.0774319) q[3];
sx q[3];
rz(-0.18571346) q[3];
sx q[3];
rz(-0.41994219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3999148) q[2];
sx q[2];
rz(-0.33881131) q[2];
sx q[2];
rz(-0.76876387) q[2];
rz(-2.9703043) q[3];
sx q[3];
rz(-1.7007217) q[3];
sx q[3];
rz(0.35458529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6866368) q[0];
sx q[0];
rz(-2.2577715) q[0];
sx q[0];
rz(2.6052642) q[0];
rz(-0.80279154) q[1];
sx q[1];
rz(-1.2840459) q[1];
sx q[1];
rz(-3.0435496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2778846) q[0];
sx q[0];
rz(-1.4914762) q[0];
sx q[0];
rz(-1.2427928) q[0];
rz(-0.98020245) q[2];
sx q[2];
rz(-0.88104311) q[2];
sx q[2];
rz(-0.99926567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6373972) q[1];
sx q[1];
rz(-1.3679844) q[1];
sx q[1];
rz(-1.0493953) q[1];
rz(-0.81252247) q[3];
sx q[3];
rz(-2.5794263) q[3];
sx q[3];
rz(-2.6443329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1596277) q[2];
sx q[2];
rz(-1.030913) q[2];
sx q[2];
rz(2.8221596) q[2];
rz(-0.56644136) q[3];
sx q[3];
rz(-2.3947377) q[3];
sx q[3];
rz(1.9802861) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40120688) q[0];
sx q[0];
rz(-1.3196608) q[0];
sx q[0];
rz(0.55289406) q[0];
rz(0.7792019) q[1];
sx q[1];
rz(-0.82256493) q[1];
sx q[1];
rz(-2.353277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32892431) q[0];
sx q[0];
rz(-1.6397489) q[0];
sx q[0];
rz(-2.8806995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5469903) q[2];
sx q[2];
rz(-1.362065) q[2];
sx q[2];
rz(-3.0301651) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7692881) q[1];
sx q[1];
rz(-1.2944229) q[1];
sx q[1];
rz(-0.45120542) q[1];
x q[2];
rz(0.68627091) q[3];
sx q[3];
rz(-2.7618558) q[3];
sx q[3];
rz(0.59829955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.39764443) q[2];
sx q[2];
rz(-2.0136191) q[2];
sx q[2];
rz(-2.9359342) q[2];
rz(-1.3501984) q[3];
sx q[3];
rz(-2.4009027) q[3];
sx q[3];
rz(0.010644309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81269294) q[0];
sx q[0];
rz(-0.44726547) q[0];
sx q[0];
rz(-1.6269667) q[0];
rz(-2.336592) q[1];
sx q[1];
rz(-1.4470419) q[1];
sx q[1];
rz(-0.31594333) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63271967) q[0];
sx q[0];
rz(-2.609786) q[0];
sx q[0];
rz(0.72152941) q[0];
x q[1];
rz(2.5300643) q[2];
sx q[2];
rz(-1.2622764) q[2];
sx q[2];
rz(-2.6055544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.915357) q[1];
sx q[1];
rz(-1.8457883) q[1];
sx q[1];
rz(0.26696856) q[1];
rz(-2.831275) q[3];
sx q[3];
rz(-1.4253221) q[3];
sx q[3];
rz(-3.1273354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.463795) q[2];
sx q[2];
rz(-0.69910502) q[2];
sx q[2];
rz(2.9839436) q[2];
rz(-2.6080103) q[3];
sx q[3];
rz(-1.6505417) q[3];
sx q[3];
rz(1.5177479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59876281) q[0];
sx q[0];
rz(-0.49648008) q[0];
sx q[0];
rz(-0.98094034) q[0];
rz(-2.6988103) q[1];
sx q[1];
rz(-1.9552224) q[1];
sx q[1];
rz(-2.669899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2344246) q[0];
sx q[0];
rz(-1.8629854) q[0];
sx q[0];
rz(2.0227595) q[0];
rz(-2.8624257) q[2];
sx q[2];
rz(-1.4067603) q[2];
sx q[2];
rz(-3.0085473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16817465) q[1];
sx q[1];
rz(-1.7061632) q[1];
sx q[1];
rz(-1.5568887) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3863435) q[3];
sx q[3];
rz(-1.235699) q[3];
sx q[3];
rz(-1.9012251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9337351) q[2];
sx q[2];
rz(-0.8693049) q[2];
sx q[2];
rz(-3.0233439) q[2];
rz(0.56784981) q[3];
sx q[3];
rz(-1.7852781) q[3];
sx q[3];
rz(0.80287272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0826913) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(-0.68728224) q[0];
rz(1.2673238) q[1];
sx q[1];
rz(-1.9272389) q[1];
sx q[1];
rz(0.6257239) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95803496) q[0];
sx q[0];
rz(-1.5661217) q[0];
sx q[0];
rz(1.5918533) q[0];
x q[1];
rz(0.56628801) q[2];
sx q[2];
rz(-1.0121965) q[2];
sx q[2];
rz(1.2142912) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/4) q[1];
sx q[1];
rz(-0.57481261) q[1];
sx q[1];
rz(-1.5479799) q[1];
rz(-2.8318303) q[3];
sx q[3];
rz(-2.1283669) q[3];
sx q[3];
rz(-0.46771177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9132793) q[2];
sx q[2];
rz(-1.2484756) q[2];
sx q[2];
rz(1.0949562) q[2];
rz(-0.45371184) q[3];
sx q[3];
rz(-2.3504421) q[3];
sx q[3];
rz(2.9554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5472645) q[0];
sx q[0];
rz(-2.6481977) q[0];
sx q[0];
rz(0.067597978) q[0];
rz(2.1122081) q[1];
sx q[1];
rz(-1.8600347) q[1];
sx q[1];
rz(0.13551113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96582039) q[0];
sx q[0];
rz(-1.0403353) q[0];
sx q[0];
rz(0.6483174) q[0];
x q[1];
rz(1.4082528) q[2];
sx q[2];
rz(-0.83627273) q[2];
sx q[2];
rz(0.61406174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1780336) q[1];
sx q[1];
rz(-1.8521223) q[1];
sx q[1];
rz(2.540349) q[1];
x q[2];
rz(-1.3247847) q[3];
sx q[3];
rz(-1.8455077) q[3];
sx q[3];
rz(-1.5266071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4542666) q[2];
sx q[2];
rz(-2.8886075) q[2];
sx q[2];
rz(-2.6542286) q[2];
rz(-0.38960114) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(3.0756557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366632) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(1.7568463) q[0];
rz(1.0549649) q[1];
sx q[1];
rz(-0.87433785) q[1];
sx q[1];
rz(2.8816282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1459634) q[0];
sx q[0];
rz(-0.71486799) q[0];
sx q[0];
rz(-1.9674106) q[0];
rz(-2.2222338) q[2];
sx q[2];
rz(-0.16460379) q[2];
sx q[2];
rz(0.92606629) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1321226) q[1];
sx q[1];
rz(-2.1984221) q[1];
sx q[1];
rz(-3.0982137) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22180827) q[3];
sx q[3];
rz(-1.7647499) q[3];
sx q[3];
rz(1.2916331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7197623) q[2];
sx q[2];
rz(-2.7037342) q[2];
sx q[2];
rz(1.3502632) q[2];
rz(-2.0566025) q[3];
sx q[3];
rz(-1.9271873) q[3];
sx q[3];
rz(-1.1291645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0709025) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(1.6593973) q[1];
sx q[1];
rz(-1.8435602) q[1];
sx q[1];
rz(-1.6866121) q[1];
rz(-1.0537335) q[2];
sx q[2];
rz(-1.0521493) q[2];
sx q[2];
rz(-0.041368749) q[2];
rz(2.0856711) q[3];
sx q[3];
rz(-0.56031651) q[3];
sx q[3];
rz(0.65800695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
