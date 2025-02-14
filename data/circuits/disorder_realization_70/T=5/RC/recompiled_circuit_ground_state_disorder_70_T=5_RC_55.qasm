OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.31324759) q[0];
sx q[0];
rz(-2.3656443) q[0];
sx q[0];
rz(1.1021855) q[0];
rz(-1.518353) q[1];
sx q[1];
rz(-2.0145388) q[1];
sx q[1];
rz(2.9161646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0505053) q[0];
sx q[0];
rz(-2.4417672) q[0];
sx q[0];
rz(-0.096629337) q[0];
x q[1];
rz(-1.7182452) q[2];
sx q[2];
rz(-0.96012291) q[2];
sx q[2];
rz(1.5924647) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3833405) q[1];
sx q[1];
rz(-1.302725) q[1];
sx q[1];
rz(-1.9490543) q[1];
x q[2];
rz(0.59360154) q[3];
sx q[3];
rz(-0.57596016) q[3];
sx q[3];
rz(-2.0570448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6585091) q[2];
sx q[2];
rz(-0.67407125) q[2];
sx q[2];
rz(-0.020641208) q[2];
rz(2.9523197) q[3];
sx q[3];
rz(-0.34957379) q[3];
sx q[3];
rz(-1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610158) q[0];
sx q[0];
rz(-0.78106946) q[0];
sx q[0];
rz(-2.1625157) q[0];
rz(0.11115221) q[1];
sx q[1];
rz(-2.4039098) q[1];
sx q[1];
rz(-2.0806064) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7951668) q[0];
sx q[0];
rz(-1.1099263) q[0];
sx q[0];
rz(-0.55924846) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1355746) q[2];
sx q[2];
rz(-0.13745795) q[2];
sx q[2];
rz(2.6743367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6305914) q[1];
sx q[1];
rz(-1.7770635) q[1];
sx q[1];
rz(0.49623112) q[1];
x q[2];
rz(-1.4270328) q[3];
sx q[3];
rz(-1.9562436) q[3];
sx q[3];
rz(-0.072162554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7079118) q[2];
sx q[2];
rz(-1.8209063) q[2];
sx q[2];
rz(2.743538) q[2];
rz(-1.4237283) q[3];
sx q[3];
rz(-0.28702304) q[3];
sx q[3];
rz(0.48582336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5240391) q[0];
sx q[0];
rz(-0.51152873) q[0];
sx q[0];
rz(2.0892573) q[0];
rz(2.6984093) q[1];
sx q[1];
rz(-2.1801703) q[1];
sx q[1];
rz(2.4876432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1285514) q[0];
sx q[0];
rz(-1.6975901) q[0];
sx q[0];
rz(2.3953505) q[0];
rz(0.072418173) q[2];
sx q[2];
rz(-0.94885599) q[2];
sx q[2];
rz(-0.79889132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2005077) q[1];
sx q[1];
rz(-1.1533179) q[1];
sx q[1];
rz(2.3509376) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(2.7707278) q[2];
sx q[2];
rz(-2.3720522) q[2];
sx q[2];
rz(2.9097843) q[2];
rz(1.2109463) q[3];
sx q[3];
rz(-1.959266) q[3];
sx q[3];
rz(0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66832191) q[0];
sx q[0];
rz(-2.7372161) q[0];
sx q[0];
rz(-0.37687287) q[0];
rz(2.6301774) q[1];
sx q[1];
rz(-1.8835953) q[1];
sx q[1];
rz(0.84397856) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653042) q[0];
sx q[0];
rz(-1.6805959) q[0];
sx q[0];
rz(3.1240433) q[0];
rz(-pi) q[1];
rz(-0.99177166) q[2];
sx q[2];
rz(-1.3879773) q[2];
sx q[2];
rz(0.43695517) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5153) q[1];
sx q[1];
rz(-0.76071435) q[1];
sx q[1];
rz(-0.73360755) q[1];
rz(-pi) q[2];
rz(-2.6199201) q[3];
sx q[3];
rz(-1.5893916) q[3];
sx q[3];
rz(-2.2357495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1692928) q[2];
sx q[2];
rz(-1.0147213) q[2];
sx q[2];
rz(2.443552) q[2];
rz(1.9418779) q[3];
sx q[3];
rz(-0.90568393) q[3];
sx q[3];
rz(-2.5419295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19071628) q[0];
sx q[0];
rz(-2.8671725) q[0];
sx q[0];
rz(3.0158667) q[0];
rz(-2.114864) q[1];
sx q[1];
rz(-1.3871565) q[1];
sx q[1];
rz(0.52621192) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65186319) q[0];
sx q[0];
rz(-1.644977) q[0];
sx q[0];
rz(-0.47898888) q[0];
rz(-1.7650928) q[2];
sx q[2];
rz(-2.435413) q[2];
sx q[2];
rz(0.45748177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3354644) q[1];
sx q[1];
rz(-1.4683691) q[1];
sx q[1];
rz(-2.7598937) q[1];
rz(-0.79546572) q[3];
sx q[3];
rz(-1.0364) q[3];
sx q[3];
rz(-0.55783844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0937664) q[2];
sx q[2];
rz(-2.3036849) q[2];
sx q[2];
rz(2.8313336) q[2];
rz(-3.1066762) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(0.78970277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652311) q[0];
sx q[0];
rz(-1.4820453) q[0];
sx q[0];
rz(-0.37102997) q[0];
rz(-3.0767483) q[1];
sx q[1];
rz(-2.3193181) q[1];
sx q[1];
rz(1.8745905) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081535) q[0];
sx q[0];
rz(-1.2282208) q[0];
sx q[0];
rz(-2.2517363) q[0];
x q[1];
rz(-2.0376181) q[2];
sx q[2];
rz(-2.4577123) q[2];
sx q[2];
rz(3.0710222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53643287) q[1];
sx q[1];
rz(-2.7419232) q[1];
sx q[1];
rz(3.0529778) q[1];
rz(-pi) q[2];
x q[2];
rz(2.606845) q[3];
sx q[3];
rz(-0.52899299) q[3];
sx q[3];
rz(-2.3803866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12396699) q[2];
sx q[2];
rz(-1.664868) q[2];
sx q[2];
rz(1.6032479) q[2];
rz(2.7247143) q[3];
sx q[3];
rz(-2.6346801) q[3];
sx q[3];
rz(-2.9712334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040319547) q[0];
sx q[0];
rz(-2.9077001) q[0];
sx q[0];
rz(-2.7129569) q[0];
rz(2.0173232) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(0.7695778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807629) q[0];
sx q[0];
rz(-1.5676982) q[0];
sx q[0];
rz(-0.0018381434) q[0];
x q[1];
rz(-1.0753845) q[2];
sx q[2];
rz(-1.1009463) q[2];
sx q[2];
rz(0.097493492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0391072) q[1];
sx q[1];
rz(-1.8643446) q[1];
sx q[1];
rz(0.62112804) q[1];
rz(-0.34129391) q[3];
sx q[3];
rz(-1.0847632) q[3];
sx q[3];
rz(0.38859545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4838685) q[2];
sx q[2];
rz(-1.2843853) q[2];
sx q[2];
rz(-3.0460975) q[2];
rz(2.5996082) q[3];
sx q[3];
rz(-2.675455) q[3];
sx q[3];
rz(0.12921648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1456881) q[0];
sx q[0];
rz(-2.2112084) q[0];
sx q[0];
rz(-3.0297739) q[0];
rz(-1.3189141) q[1];
sx q[1];
rz(-1.8967862) q[1];
sx q[1];
rz(1.3056614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3506136) q[0];
sx q[0];
rz(-1.0270938) q[0];
sx q[0];
rz(-0.14474317) q[0];
rz(2.6839031) q[2];
sx q[2];
rz(-2.8297242) q[2];
sx q[2];
rz(1.0491187) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9247465) q[1];
sx q[1];
rz(-1.9439183) q[1];
sx q[1];
rz(-1.8736963) q[1];
rz(-pi) q[2];
rz(-1.1505125) q[3];
sx q[3];
rz(-1.0627565) q[3];
sx q[3];
rz(-0.86364323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8127415) q[2];
sx q[2];
rz(-1.988827) q[2];
sx q[2];
rz(-2.0218772) q[2];
rz(2.8730734) q[3];
sx q[3];
rz(-1.4041785) q[3];
sx q[3];
rz(1.9560811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28065228) q[0];
sx q[0];
rz(-0.76186162) q[0];
sx q[0];
rz(2.5326488) q[0];
rz(-2.2641585) q[1];
sx q[1];
rz(-2.139822) q[1];
sx q[1];
rz(0.61224365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5800047) q[0];
sx q[0];
rz(-1.1682434) q[0];
sx q[0];
rz(-0.7640362) q[0];
rz(-2.5721278) q[2];
sx q[2];
rz(-1.8761022) q[2];
sx q[2];
rz(-1.832163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7944736) q[1];
sx q[1];
rz(-1.11382) q[1];
sx q[1];
rz(-2.4040078) q[1];
rz(-3.0123467) q[3];
sx q[3];
rz(-2.5967483) q[3];
sx q[3];
rz(-2.2428577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.073254243) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(0.18745984) q[2];
rz(2.0295664) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(-2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336695) q[0];
sx q[0];
rz(-2.3275571) q[0];
sx q[0];
rz(2.8975876) q[0];
rz(-1.6581274) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(0.77267486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940154) q[0];
sx q[0];
rz(-2.4227067) q[0];
sx q[0];
rz(0.80475828) q[0];
x q[1];
rz(-1.7079321) q[2];
sx q[2];
rz(-1.8923645) q[2];
sx q[2];
rz(2.0069881) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15068842) q[1];
sx q[1];
rz(-0.3462308) q[1];
sx q[1];
rz(-0.33106403) q[1];
rz(-pi) q[2];
rz(0.78759191) q[3];
sx q[3];
rz(-2.1887472) q[3];
sx q[3];
rz(1.4026812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8030168) q[2];
sx q[2];
rz(-0.96119857) q[2];
sx q[2];
rz(2.4147066) q[2];
rz(-1.1072985) q[3];
sx q[3];
rz(-0.50873435) q[3];
sx q[3];
rz(-1.0072964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016639391) q[0];
sx q[0];
rz(-1.4062842) q[0];
sx q[0];
rz(1.9840065) q[0];
rz(-2.6955556) q[1];
sx q[1];
rz(-2.0560494) q[1];
sx q[1];
rz(2.5037419) q[1];
rz(-2.73861) q[2];
sx q[2];
rz(-1.4558515) q[2];
sx q[2];
rz(0.2414138) q[2];
rz(0.42556006) q[3];
sx q[3];
rz(-1.7055409) q[3];
sx q[3];
rz(2.9421229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
