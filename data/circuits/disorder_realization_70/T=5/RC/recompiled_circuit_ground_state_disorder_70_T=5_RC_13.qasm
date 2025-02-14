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
rz(3.917541) q[0];
sx q[0];
rz(10.526963) q[0];
rz(-1.518353) q[1];
sx q[1];
rz(-2.0145388) q[1];
sx q[1];
rz(2.9161646) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537215) q[0];
sx q[0];
rz(-1.632977) q[0];
sx q[0];
rz(2.4440701) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7182452) q[2];
sx q[2];
rz(-2.1814697) q[2];
sx q[2];
rz(-1.5491279) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3833405) q[1];
sx q[1];
rz(-1.302725) q[1];
sx q[1];
rz(1.1925384) q[1];
rz(-pi) q[2];
rz(-2.5479911) q[3];
sx q[3];
rz(-0.57596016) q[3];
sx q[3];
rz(-2.0570448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6585091) q[2];
sx q[2];
rz(-0.67407125) q[2];
sx q[2];
rz(0.020641208) q[2];
rz(-0.189273) q[3];
sx q[3];
rz(-0.34957379) q[3];
sx q[3];
rz(1.6674204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78057688) q[0];
sx q[0];
rz(-2.3605232) q[0];
sx q[0];
rz(0.97907698) q[0];
rz(-3.0304404) q[1];
sx q[1];
rz(-2.4039098) q[1];
sx q[1];
rz(1.0609863) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34642588) q[0];
sx q[0];
rz(-2.0316664) q[0];
sx q[0];
rz(0.55924846) q[0];
rz(-3.1355746) q[2];
sx q[2];
rz(-3.0041347) q[2];
sx q[2];
rz(0.46725592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0909522) q[1];
sx q[1];
rz(-1.0860069) q[1];
sx q[1];
rz(-1.3371972) q[1];
x q[2];
rz(0.33943601) q[3];
sx q[3];
rz(-0.41012479) q[3];
sx q[3];
rz(-2.8462178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4336808) q[2];
sx q[2];
rz(-1.3206864) q[2];
sx q[2];
rz(2.743538) q[2];
rz(1.7178644) q[3];
sx q[3];
rz(-2.8545696) q[3];
sx q[3];
rz(-0.48582336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(0.61755359) q[0];
sx q[0];
rz(-0.51152873) q[0];
sx q[0];
rz(-1.0523354) q[0];
rz(2.6984093) q[1];
sx q[1];
rz(-0.96142238) q[1];
sx q[1];
rz(0.65394941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4674462) q[0];
sx q[0];
rz(-2.3096414) q[0];
sx q[0];
rz(-1.7427) q[0];
x q[1];
rz(-0.94761272) q[2];
sx q[2];
rz(-1.511956) q[2];
sx q[2];
rz(-0.72966444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7591651) q[1];
sx q[1];
rz(-2.2780721) q[1];
sx q[1];
rz(-1.00818) q[1];
rz(1.3929358) q[3];
sx q[3];
rz(-1.7962958) q[3];
sx q[3];
rz(0.64060831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7707278) q[2];
sx q[2];
rz(-2.3720522) q[2];
sx q[2];
rz(2.9097843) q[2];
rz(-1.9306463) q[3];
sx q[3];
rz(-1.959266) q[3];
sx q[3];
rz(-2.8162091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66832191) q[0];
sx q[0];
rz(-2.7372161) q[0];
sx q[0];
rz(0.37687287) q[0];
rz(0.51141524) q[1];
sx q[1];
rz(-1.2579974) q[1];
sx q[1];
rz(0.84397856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653042) q[0];
sx q[0];
rz(-1.6805959) q[0];
sx q[0];
rz(0.017549347) q[0];
rz(-pi) q[1];
rz(-0.21739809) q[2];
sx q[2];
rz(-2.1389642) q[2];
sx q[2];
rz(1.252144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36605814) q[1];
sx q[1];
rz(-2.0506128) q[1];
sx q[1];
rz(2.5262031) q[1];
rz(-pi) q[2];
rz(1.5922437) q[3];
sx q[3];
rz(-1.0492232) q[3];
sx q[3];
rz(-0.65426588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1692928) q[2];
sx q[2];
rz(-2.1268714) q[2];
sx q[2];
rz(-0.69804066) q[2];
rz(-1.1997148) q[3];
sx q[3];
rz(-0.90568393) q[3];
sx q[3];
rz(0.5996632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19071628) q[0];
sx q[0];
rz(-0.27442014) q[0];
sx q[0];
rz(0.12572591) q[0];
rz(2.114864) q[1];
sx q[1];
rz(-1.7544361) q[1];
sx q[1];
rz(0.52621192) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4897295) q[0];
sx q[0];
rz(-1.4966156) q[0];
sx q[0];
rz(-0.47898888) q[0];
rz(-pi) q[1];
rz(2.9783812) q[2];
sx q[2];
rz(-2.2610352) q[2];
sx q[2];
rz(2.9371967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3359068) q[1];
sx q[1];
rz(-1.1912002) q[1];
sx q[1];
rz(1.6811045) q[1];
rz(-pi) q[2];
rz(-0.69198841) q[3];
sx q[3];
rz(-0.92433911) q[3];
sx q[3];
rz(-2.5916168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0478263) q[2];
sx q[2];
rz(-0.83790773) q[2];
sx q[2];
rz(2.8313336) q[2];
rz(0.034916498) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(-2.3518899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0763615) q[0];
sx q[0];
rz(-1.6595474) q[0];
sx q[0];
rz(-0.37102997) q[0];
rz(3.0767483) q[1];
sx q[1];
rz(-0.82227451) q[1];
sx q[1];
rz(1.8745905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081535) q[0];
sx q[0];
rz(-1.9133718) q[0];
sx q[0];
rz(-0.88985635) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94159884) q[2];
sx q[2];
rz(-1.8591188) q[2];
sx q[2];
rz(-1.1278111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6051598) q[1];
sx q[1];
rz(-0.39966941) q[1];
sx q[1];
rz(3.0529778) q[1];
rz(-pi) q[2];
rz(-0.53474769) q[3];
sx q[3];
rz(-0.52899299) q[3];
sx q[3];
rz(-2.3803866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0176257) q[2];
sx q[2];
rz(-1.4767246) q[2];
sx q[2];
rz(-1.6032479) q[2];
rz(-0.41687837) q[3];
sx q[3];
rz(-0.50691253) q[3];
sx q[3];
rz(-0.1703593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012731) q[0];
sx q[0];
rz(-0.23389255) q[0];
sx q[0];
rz(-0.42863578) q[0];
rz(-2.0173232) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(2.3720149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3316318) q[0];
sx q[0];
rz(-1.5689582) q[0];
sx q[0];
rz(1.5738945) q[0];
rz(1.0753845) q[2];
sx q[2];
rz(-1.1009463) q[2];
sx q[2];
rz(-0.097493492) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10248549) q[1];
sx q[1];
rz(-1.277248) q[1];
sx q[1];
rz(0.62112804) q[1];
rz(-0.34129391) q[3];
sx q[3];
rz(-2.0568295) q[3];
sx q[3];
rz(2.7529972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4838685) q[2];
sx q[2];
rz(-1.2843853) q[2];
sx q[2];
rz(-3.0460975) q[2];
rz(-2.5996082) q[3];
sx q[3];
rz(-2.675455) q[3];
sx q[3];
rz(3.0123762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.99590456) q[0];
sx q[0];
rz(-0.93038428) q[0];
sx q[0];
rz(3.0297739) q[0];
rz(1.3189141) q[1];
sx q[1];
rz(-1.8967862) q[1];
sx q[1];
rz(1.8359312) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3506136) q[0];
sx q[0];
rz(-1.0270938) q[0];
sx q[0];
rz(2.9968495) q[0];
rz(-0.4576896) q[2];
sx q[2];
rz(-2.8297242) q[2];
sx q[2];
rz(1.0491187) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9247465) q[1];
sx q[1];
rz(-1.1976744) q[1];
sx q[1];
rz(-1.8736963) q[1];
rz(-pi) q[2];
rz(-0.6324082) q[3];
sx q[3];
rz(-0.64738368) q[3];
sx q[3];
rz(1.5349015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32885113) q[2];
sx q[2];
rz(-1.1527656) q[2];
sx q[2];
rz(2.0218772) q[2];
rz(2.8730734) q[3];
sx q[3];
rz(-1.4041785) q[3];
sx q[3];
rz(-1.1855116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28065228) q[0];
sx q[0];
rz(-2.379731) q[0];
sx q[0];
rz(-2.5326488) q[0];
rz(-0.87743419) q[1];
sx q[1];
rz(-1.0017706) q[1];
sx q[1];
rz(0.61224365) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7626678) q[0];
sx q[0];
rz(-0.84419709) q[0];
sx q[0];
rz(0.55171497) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2127146) q[2];
sx q[2];
rz(-1.0306669) q[2];
sx q[2];
rz(-3.0703406) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.73758482) q[1];
rz(-pi) q[2];
rz(-1.4928451) q[3];
sx q[3];
rz(-2.1105937) q[3];
sx q[3];
rz(-1.0495561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0683384) q[2];
sx q[2];
rz(-1.6011668) q[2];
sx q[2];
rz(2.9541328) q[2];
rz(-2.0295664) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336695) q[0];
sx q[0];
rz(-0.81403553) q[0];
sx q[0];
rz(-2.8975876) q[0];
rz(-1.6581274) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(-2.3689178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84757728) q[0];
sx q[0];
rz(-2.4227067) q[0];
sx q[0];
rz(-0.80475828) q[0];
rz(-2.7521801) q[2];
sx q[2];
rz(-0.34865278) q[2];
sx q[2];
rz(1.5462923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15068842) q[1];
sx q[1];
rz(-2.7953618) q[1];
sx q[1];
rz(0.33106403) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3546702) q[3];
sx q[3];
rz(-0.95810181) q[3];
sx q[3];
rz(0.69132346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8030168) q[2];
sx q[2];
rz(-0.96119857) q[2];
sx q[2];
rz(0.72688603) q[2];
rz(-2.0342942) q[3];
sx q[3];
rz(-2.6328583) q[3];
sx q[3];
rz(2.1342962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1399287) q[0];
sx q[0];
rz(-1.4062842) q[0];
sx q[0];
rz(1.9840065) q[0];
rz(2.6955556) q[1];
sx q[1];
rz(-1.0855433) q[1];
sx q[1];
rz(-0.6378508) q[1];
rz(-0.40298265) q[2];
sx q[2];
rz(-1.6857412) q[2];
sx q[2];
rz(-2.9001789) q[2];
rz(-1.7185531) q[3];
sx q[3];
rz(-1.9922517) q[3];
sx q[3];
rz(-1.7094517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
