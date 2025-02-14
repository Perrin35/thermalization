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
rz(2.0394072) q[0];
rz(1.6232396) q[1];
sx q[1];
rz(2.0145388) q[1];
sx q[1];
rz(9.1993499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537215) q[0];
sx q[0];
rz(-1.5086156) q[0];
sx q[0];
rz(-0.6975226) q[0];
x q[1];
rz(-2.5257843) q[2];
sx q[2];
rz(-1.4501418) q[2];
sx q[2];
rz(-3.0349611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7582522) q[1];
sx q[1];
rz(-1.8388677) q[1];
sx q[1];
rz(-1.1925384) q[1];
rz(-1.9192237) q[3];
sx q[3];
rz(-1.1023848) q[3];
sx q[3];
rz(-2.7345776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48308358) q[2];
sx q[2];
rz(-0.67407125) q[2];
sx q[2];
rz(-3.1209514) q[2];
rz(-0.189273) q[3];
sx q[3];
rz(-0.34957379) q[3];
sx q[3];
rz(-1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610158) q[0];
sx q[0];
rz(-0.78106946) q[0];
sx q[0];
rz(-2.1625157) q[0];
rz(-0.11115221) q[1];
sx q[1];
rz(-0.73768288) q[1];
sx q[1];
rz(-2.0806064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1887061) q[0];
sx q[0];
rz(-2.0659819) q[0];
sx q[0];
rz(1.0409036) q[0];
rz(-0.13745549) q[2];
sx q[2];
rz(-1.5699717) q[2];
sx q[2];
rz(2.0440136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6305914) q[1];
sx q[1];
rz(-1.7770635) q[1];
sx q[1];
rz(-2.6453615) q[1];
rz(1.7145598) q[3];
sx q[3];
rz(-1.9562436) q[3];
sx q[3];
rz(3.0694301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7079118) q[2];
sx q[2];
rz(-1.8209063) q[2];
sx q[2];
rz(-0.39805463) q[2];
rz(1.7178644) q[3];
sx q[3];
rz(-0.28702304) q[3];
sx q[3];
rz(0.48582336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240391) q[0];
sx q[0];
rz(-0.51152873) q[0];
sx q[0];
rz(1.0523354) q[0];
rz(2.6984093) q[1];
sx q[1];
rz(-2.1801703) q[1];
sx q[1];
rz(-0.65394941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67414647) q[0];
sx q[0];
rz(-2.3096414) q[0];
sx q[0];
rz(-1.7427) q[0];
x q[1];
rz(-0.072418173) q[2];
sx q[2];
rz(-2.1927367) q[2];
sx q[2];
rz(2.3427013) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3824275) q[1];
sx q[1];
rz(-2.2780721) q[1];
sx q[1];
rz(2.1334126) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3929358) q[3];
sx q[3];
rz(-1.3452969) q[3];
sx q[3];
rz(-0.64060831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7707278) q[2];
sx q[2];
rz(-0.76954049) q[2];
sx q[2];
rz(-0.23180836) q[2];
rz(1.9306463) q[3];
sx q[3];
rz(-1.1823267) q[3];
sx q[3];
rz(0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66832191) q[0];
sx q[0];
rz(-2.7372161) q[0];
sx q[0];
rz(0.37687287) q[0];
rz(-0.51141524) q[1];
sx q[1];
rz(-1.8835953) q[1];
sx q[1];
rz(-2.2976141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6451615) q[0];
sx q[0];
rz(-1.58824) q[0];
sx q[0];
rz(-1.6806127) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9241946) q[2];
sx q[2];
rz(-1.0026284) q[2];
sx q[2];
rz(1.8894486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36605814) q[1];
sx q[1];
rz(-2.0506128) q[1];
sx q[1];
rz(-0.61538954) q[1];
x q[2];
rz(-2.6199201) q[3];
sx q[3];
rz(-1.552201) q[3];
sx q[3];
rz(2.2357495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9722998) q[2];
sx q[2];
rz(-2.1268714) q[2];
sx q[2];
rz(-2.443552) q[2];
rz(1.1997148) q[3];
sx q[3];
rz(-0.90568393) q[3];
sx q[3];
rz(2.5419295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19071628) q[0];
sx q[0];
rz(-0.27442014) q[0];
sx q[0];
rz(3.0158667) q[0];
rz(1.0267286) q[1];
sx q[1];
rz(-1.3871565) q[1];
sx q[1];
rz(-2.6153807) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4897295) q[0];
sx q[0];
rz(-1.644977) q[0];
sx q[0];
rz(-0.47898888) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87398087) q[2];
sx q[2];
rz(-1.6964198) q[2];
sx q[2];
rz(-1.2619293) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51523415) q[1];
sx q[1];
rz(-0.39455104) q[1];
sx q[1];
rz(2.872353) q[1];
rz(-2.3461269) q[3];
sx q[3];
rz(-1.0364) q[3];
sx q[3];
rz(-2.5837542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0937664) q[2];
sx q[2];
rz(-2.3036849) q[2];
sx q[2];
rz(0.31025904) q[2];
rz(-0.034916498) q[3];
sx q[3];
rz(-1.755545) q[3];
sx q[3];
rz(0.78970277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0652311) q[0];
sx q[0];
rz(-1.6595474) q[0];
sx q[0];
rz(2.7705627) q[0];
rz(-0.064844355) q[1];
sx q[1];
rz(-0.82227451) q[1];
sx q[1];
rz(1.8745905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112087) q[0];
sx q[0];
rz(-0.74978854) q[0];
sx q[0];
rz(-2.0862338) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0376181) q[2];
sx q[2];
rz(-2.4577123) q[2];
sx q[2];
rz(0.070570408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53643287) q[1];
sx q[1];
rz(-0.39966941) q[1];
sx q[1];
rz(3.0529778) q[1];
rz(-pi) q[2];
rz(-2.606845) q[3];
sx q[3];
rz(-0.52899299) q[3];
sx q[3];
rz(2.3803866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0176257) q[2];
sx q[2];
rz(-1.664868) q[2];
sx q[2];
rz(-1.6032479) q[2];
rz(0.41687837) q[3];
sx q[3];
rz(-0.50691253) q[3];
sx q[3];
rz(0.1703593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040319547) q[0];
sx q[0];
rz(-2.9077001) q[0];
sx q[0];
rz(2.7129569) q[0];
rz(-2.0173232) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(2.3720149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76082976) q[0];
sx q[0];
rz(-1.5676982) q[0];
sx q[0];
rz(-0.0018381434) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0662082) q[2];
sx q[2];
rz(-2.0406463) q[2];
sx q[2];
rz(3.0440992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0391072) q[1];
sx q[1];
rz(-1.8643446) q[1];
sx q[1];
rz(0.62112804) q[1];
rz(-pi) q[2];
rz(-2.0817716) q[3];
sx q[3];
rz(-1.8712412) q[3];
sx q[3];
rz(-1.0177801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.6577242) q[2];
sx q[2];
rz(-1.2843853) q[2];
sx q[2];
rz(-0.095495187) q[2];
rz(2.5996082) q[3];
sx q[3];
rz(-2.675455) q[3];
sx q[3];
rz(-3.0123762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99590456) q[0];
sx q[0];
rz(-2.2112084) q[0];
sx q[0];
rz(-0.11181871) q[0];
rz(1.3189141) q[1];
sx q[1];
rz(-1.2448064) q[1];
sx q[1];
rz(1.3056614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163286) q[0];
sx q[0];
rz(-2.5808307) q[0];
sx q[0];
rz(1.8050342) q[0];
rz(2.6839031) q[2];
sx q[2];
rz(-0.31186843) q[2];
sx q[2];
rz(2.092474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2160141) q[1];
sx q[1];
rz(-0.47607254) q[1];
sx q[1];
rz(-2.4904597) q[1];
rz(2.5939529) q[3];
sx q[3];
rz(-1.2062916) q[3];
sx q[3];
rz(-2.648516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8127415) q[2];
sx q[2];
rz(-1.988827) q[2];
sx q[2];
rz(2.0218772) q[2];
rz(2.8730734) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(1.1855116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609404) q[0];
sx q[0];
rz(-2.379731) q[0];
sx q[0];
rz(0.60894388) q[0];
rz(2.2641585) q[1];
sx q[1];
rz(-1.0017706) q[1];
sx q[1];
rz(-2.529349) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5615879) q[0];
sx q[0];
rz(-1.9733493) q[0];
sx q[0];
rz(2.3775565) q[0];
x q[1];
rz(-2.5721278) q[2];
sx q[2];
rz(-1.8761022) q[2];
sx q[2];
rz(1.3094297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5366) q[1];
sx q[1];
rz(-0.92291622) q[1];
sx q[1];
rz(-0.98439321) q[1];
rz(-0.12924592) q[3];
sx q[3];
rz(-0.5448444) q[3];
sx q[3];
rz(0.89873492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.073254243) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(-2.9541328) q[2];
rz(2.0295664) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(-2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40792313) q[0];
sx q[0];
rz(-2.3275571) q[0];
sx q[0];
rz(-2.8975876) q[0];
rz(1.4834652) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(0.77267486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2940154) q[0];
sx q[0];
rz(-2.4227067) q[0];
sx q[0];
rz(-0.80475828) q[0];
rz(2.7521801) q[2];
sx q[2];
rz(-0.34865278) q[2];
sx q[2];
rz(1.5953003) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9909042) q[1];
sx q[1];
rz(-2.7953618) q[1];
sx q[1];
rz(2.8105286) q[1];
x q[2];
rz(-2.3540007) q[3];
sx q[3];
rz(-0.95284546) q[3];
sx q[3];
rz(1.7389115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3385758) q[2];
sx q[2];
rz(-2.1803941) q[2];
sx q[2];
rz(-2.4147066) q[2];
rz(1.1072985) q[3];
sx q[3];
rz(-0.50873435) q[3];
sx q[3];
rz(-2.1342962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-3.1399287) q[0];
sx q[0];
rz(-1.4062842) q[0];
sx q[0];
rz(1.9840065) q[0];
rz(-0.4460371) q[1];
sx q[1];
rz(-1.0855433) q[1];
sx q[1];
rz(-0.6378508) q[1];
rz(-1.4459417) q[2];
sx q[2];
rz(-1.9709675) q[2];
sx q[2];
rz(1.7633543) q[2];
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
