OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1096119) q[0];
sx q[0];
rz(-1.0053827) q[0];
sx q[0];
rz(0.96860743) q[0];
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(3.1117575) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.580509) q[0];
sx q[0];
rz(-2.5445815) q[0];
sx q[0];
rz(-2.719814) q[0];
rz(-1.6462684) q[2];
sx q[2];
rz(-2.3879037) q[2];
sx q[2];
rz(1.2686632) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3181495) q[1];
sx q[1];
rz(-1.9572721) q[1];
sx q[1];
rz(2.3837368) q[1];
rz(-pi) q[2];
rz(-1.6601059) q[3];
sx q[3];
rz(-1.4353283) q[3];
sx q[3];
rz(0.40218807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5030824) q[2];
sx q[2];
rz(-2.5677887) q[2];
sx q[2];
rz(0.32162515) q[2];
rz(-0.62421787) q[3];
sx q[3];
rz(-1.7523242) q[3];
sx q[3];
rz(-1.483884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44840789) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-0.93628991) q[0];
rz(-1.869465) q[1];
sx q[1];
rz(-1.5969758) q[1];
sx q[1];
rz(1.7535271) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427338) q[0];
sx q[0];
rz(-1.6279477) q[0];
sx q[0];
rz(-1.7903891) q[0];
x q[1];
rz(1.040561) q[2];
sx q[2];
rz(-2.2746646) q[2];
sx q[2];
rz(-2.3153697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85390988) q[1];
sx q[1];
rz(-2.1016995) q[1];
sx q[1];
rz(0.7002875) q[1];
rz(-2.2514492) q[3];
sx q[3];
rz(-0.33484866) q[3];
sx q[3];
rz(-1.0687506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0442514) q[2];
sx q[2];
rz(-2.1772431) q[2];
sx q[2];
rz(-0.67330366) q[2];
rz(-1.8246957) q[3];
sx q[3];
rz(-1.293106) q[3];
sx q[3];
rz(-0.84002686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3299385) q[0];
sx q[0];
rz(-1.5978156) q[0];
sx q[0];
rz(-0.04976186) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.754909) q[1];
sx q[1];
rz(1.6516364) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4987652) q[0];
sx q[0];
rz(-2.1826361) q[0];
sx q[0];
rz(-1.4239514) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0115205) q[2];
sx q[2];
rz(-2.080938) q[2];
sx q[2];
rz(2.5464692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7327642) q[1];
sx q[1];
rz(-1.8547426) q[1];
sx q[1];
rz(-2.8676621) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8036267) q[3];
sx q[3];
rz(-2.0062048) q[3];
sx q[3];
rz(2.0056279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9516307) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(0.0541617) q[2];
rz(-1.1545898) q[3];
sx q[3];
rz(-1.8535987) q[3];
sx q[3];
rz(-2.0258928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40256777) q[0];
sx q[0];
rz(-1.4250647) q[0];
sx q[0];
rz(-0.85790747) q[0];
rz(-0.67445451) q[1];
sx q[1];
rz(-2.1916788) q[1];
sx q[1];
rz(1.8106921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895578) q[0];
sx q[0];
rz(-1.4343795) q[0];
sx q[0];
rz(-2.9035764) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2670008) q[2];
sx q[2];
rz(-1.7335606) q[2];
sx q[2];
rz(-0.13653423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1787599) q[1];
sx q[1];
rz(-2.5542027) q[1];
sx q[1];
rz(1.2580832) q[1];
rz(-0.0033802948) q[3];
sx q[3];
rz(-1.1515318) q[3];
sx q[3];
rz(-1.5685035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7669507) q[2];
sx q[2];
rz(-0.16051126) q[2];
sx q[2];
rz(0.10425076) q[2];
rz(0.57779622) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(-0.92208636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.6660026) q[0];
sx q[0];
rz(-0.9444899) q[0];
rz(-1.7780444) q[1];
sx q[1];
rz(-1.8133769) q[1];
sx q[1];
rz(1.69467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3730225) q[0];
sx q[0];
rz(-1.3813586) q[0];
sx q[0];
rz(-0.85590881) q[0];
rz(3.1051272) q[2];
sx q[2];
rz(-1.8210095) q[2];
sx q[2];
rz(1.6584876) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80396485) q[1];
sx q[1];
rz(-1.9586438) q[1];
sx q[1];
rz(-2.4799073) q[1];
x q[2];
rz(0.57525697) q[3];
sx q[3];
rz(-1.0759531) q[3];
sx q[3];
rz(0.32512966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.038534433) q[2];
sx q[2];
rz(-0.47163042) q[2];
sx q[2];
rz(0.72059694) q[2];
rz(3.0125812) q[3];
sx q[3];
rz(-1.3611662) q[3];
sx q[3];
rz(-1.8425997) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593889) q[0];
sx q[0];
rz(-2.144618) q[0];
sx q[0];
rz(0.61656117) q[0];
rz(0.95245862) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(1.2430826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37694398) q[0];
sx q[0];
rz(-0.29452919) q[0];
sx q[0];
rz(1.0824979) q[0];
rz(-pi) q[1];
rz(1.0988192) q[2];
sx q[2];
rz(-2.3772559) q[2];
sx q[2];
rz(-0.42407045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7156429) q[1];
sx q[1];
rz(-2.2897807) q[1];
sx q[1];
rz(-1.4904316) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0674446) q[3];
sx q[3];
rz(-2.1909602) q[3];
sx q[3];
rz(0.95297232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9055966) q[2];
sx q[2];
rz(-0.2874898) q[2];
sx q[2];
rz(0.074404152) q[2];
rz(2.060804) q[3];
sx q[3];
rz(-1.4932884) q[3];
sx q[3];
rz(-2.1006445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3022795) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(-3.0795414) q[0];
rz(-1.7265559) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(0.089990377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9190529) q[0];
sx q[0];
rz(-2.3947859) q[0];
sx q[0];
rz(1.1928305) q[0];
rz(-pi) q[1];
rz(0.16570602) q[2];
sx q[2];
rz(-0.88962727) q[2];
sx q[2];
rz(-2.990304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5127232) q[1];
sx q[1];
rz(-2.1640393) q[1];
sx q[1];
rz(-0.42494659) q[1];
rz(2.2136647) q[3];
sx q[3];
rz(-2.0629598) q[3];
sx q[3];
rz(1.6097691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0061079582) q[2];
sx q[2];
rz(-1.8718448) q[2];
sx q[2];
rz(0.6832068) q[2];
rz(1.3225383) q[3];
sx q[3];
rz(-0.95219487) q[3];
sx q[3];
rz(1.7542084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4676062) q[0];
sx q[0];
rz(-1.4961996) q[0];
sx q[0];
rz(2.6737387) q[0];
rz(2.4521008) q[1];
sx q[1];
rz(-2.4263224) q[1];
sx q[1];
rz(-0.47324866) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9625585) q[0];
sx q[0];
rz(-2.9518572) q[0];
sx q[0];
rz(-1.9535854) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6976349) q[2];
sx q[2];
rz(-1.4251764) q[2];
sx q[2];
rz(-0.99010003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2040904) q[1];
sx q[1];
rz(-0.023462208) q[1];
sx q[1];
rz(-0.84659441) q[1];
rz(0.28041081) q[3];
sx q[3];
rz(-1.3254998) q[3];
sx q[3];
rz(0.51285686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0143934) q[2];
sx q[2];
rz(-1.8202123) q[2];
sx q[2];
rz(-2.7180706) q[2];
rz(-1.315377) q[3];
sx q[3];
rz(-2.6500621) q[3];
sx q[3];
rz(-1.8486842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2143283) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.3342791) q[0];
rz(-1.3182053) q[1];
sx q[1];
rz(-0.91000906) q[1];
sx q[1];
rz(3.1277025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0158169) q[0];
sx q[0];
rz(-2.7770574) q[0];
sx q[0];
rz(-1.4992072) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35160323) q[2];
sx q[2];
rz(-1.4228063) q[2];
sx q[2];
rz(2.6949712) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.909643) q[1];
sx q[1];
rz(-1.0000537) q[1];
sx q[1];
rz(-2.8464223) q[1];
x q[2];
rz(-3.0180279) q[3];
sx q[3];
rz(-0.84513226) q[3];
sx q[3];
rz(2.7853876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5615329) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(2.0884936) q[2];
rz(2.2140391) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(0.48228669) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841566) q[0];
sx q[0];
rz(-0.3967537) q[0];
sx q[0];
rz(2.5411153) q[0];
rz(-2.8720169) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(-1.762766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.775009) q[0];
sx q[0];
rz(-1.8407158) q[0];
sx q[0];
rz(1.2054514) q[0];
rz(-2.7711033) q[2];
sx q[2];
rz(-1.8571808) q[2];
sx q[2];
rz(0.011025393) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32374183) q[1];
sx q[1];
rz(-1.3881589) q[1];
sx q[1];
rz(-0.43437965) q[1];
rz(-pi) q[2];
rz(-2.5254978) q[3];
sx q[3];
rz(-1.2047676) q[3];
sx q[3];
rz(-0.6835365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1158585) q[2];
sx q[2];
rz(-0.33418843) q[2];
sx q[2];
rz(-2.1791229) q[2];
rz(2.1879503) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(-1.3941221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3772603) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(-1.3537021) q[1];
sx q[1];
rz(-1.1693015) q[1];
sx q[1];
rz(2.404626) q[1];
rz(-2.4719292) q[2];
sx q[2];
rz(-1.0367994) q[2];
sx q[2];
rz(-1.8112469) q[2];
rz(2.6346281) q[3];
sx q[3];
rz(-1.4486113) q[3];
sx q[3];
rz(2.2511528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
