OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7115294) q[0];
sx q[0];
rz(-1.9646514) q[0];
sx q[0];
rz(0.060056134) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(1.905029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3136183) q[1];
sx q[1];
rz(-0.59578124) q[1];
sx q[1];
rz(1.6525364) q[1];
x q[2];
rz(-2.0446159) q[3];
sx q[3];
rz(-1.5690656) q[3];
sx q[3];
rz(-2.2497821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(-1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.21355024) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.7785633) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33107685) q[0];
sx q[0];
rz(-1.1805981) q[0];
sx q[0];
rz(-0.66731989) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.734191) q[2];
sx q[2];
rz(-0.50440895) q[2];
sx q[2];
rz(0.87770578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0568309) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(-0.057573307) q[1];
x q[2];
rz(0.95008534) q[3];
sx q[3];
rz(-0.32334298) q[3];
sx q[3];
rz(0.33085535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.8691241) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0740046) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(2.8787676) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(-1.2353108) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10997406) q[0];
sx q[0];
rz(-1.585916) q[0];
sx q[0];
rz(-1.581024) q[0];
rz(-pi) q[1];
rz(-3.1149408) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(-1.9453366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23535145) q[1];
sx q[1];
rz(-1.6760577) q[1];
sx q[1];
rz(0.32354849) q[1];
x q[2];
rz(-0.87576207) q[3];
sx q[3];
rz(-1.2998298) q[3];
sx q[3];
rz(-0.7489487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(0.18033218) q[2];
rz(0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.3746049) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(1.0346574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2186787) q[0];
sx q[0];
rz(-1.5905587) q[0];
sx q[0];
rz(3.0826871) q[0];
rz(-1.4356394) q[2];
sx q[2];
rz(-1.5237336) q[2];
sx q[2];
rz(-1.0263718) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6369578) q[1];
sx q[1];
rz(-0.89244288) q[1];
sx q[1];
rz(-1.9370609) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6684768) q[3];
sx q[3];
rz(-1.4149727) q[3];
sx q[3];
rz(-1.5161878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(-2.9053524) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(0.62228084) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21475269) q[0];
sx q[0];
rz(-1.8536957) q[0];
sx q[0];
rz(-0.19006417) q[0];
rz(-pi) q[1];
rz(3.0201966) q[2];
sx q[2];
rz(-2.4897794) q[2];
sx q[2];
rz(0.35169841) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9363074) q[1];
sx q[1];
rz(-0.19166066) q[1];
sx q[1];
rz(-1.081341) q[1];
x q[2];
rz(-0.47982521) q[3];
sx q[3];
rz(-1.2180426) q[3];
sx q[3];
rz(1.477581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92419147) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(0.021281555) q[2];
rz(1.6993258) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570046) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-0.044629991) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32193004) q[0];
sx q[0];
rz(-0.18580431) q[0];
sx q[0];
rz(1.6447322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0738515) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(1.5397132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1457445) q[1];
sx q[1];
rz(-1.2461582) q[1];
sx q[1];
rz(1.0041298) q[1];
rz(-pi) q[2];
rz(-1.7354292) q[3];
sx q[3];
rz(-0.87222404) q[3];
sx q[3];
rz(0.64100953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(2.1748523) q[2];
rz(-3.1294075) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(-2.1475041) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(2.4023043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4261739) q[0];
sx q[0];
rz(-2.9746375) q[0];
sx q[0];
rz(-0.963361) q[0];
x q[1];
rz(-0.57619608) q[2];
sx q[2];
rz(-1.9857166) q[2];
sx q[2];
rz(2.8232099) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.90034396) q[1];
sx q[1];
rz(-2.5785891) q[1];
sx q[1];
rz(-0.36069936) q[1];
x q[2];
rz(-2.7474653) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69650841) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-0.27841321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42785545) q[0];
sx q[0];
rz(-0.58361485) q[0];
sx q[0];
rz(-1.7996656) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7114696) q[2];
sx q[2];
rz(-0.42818907) q[2];
sx q[2];
rz(-0.45454121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6555357) q[1];
sx q[1];
rz(-1.9209407) q[1];
sx q[1];
rz(-2.9085367) q[1];
x q[2];
rz(-1.9592459) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(-2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2723096) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(-0.38796866) q[2];
rz(1.3502454) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.2824654) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-0.073721185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026422231) q[0];
sx q[0];
rz(-2.0048884) q[0];
sx q[0];
rz(0.11154453) q[0];
x q[1];
rz(2.164052) q[2];
sx q[2];
rz(-2.9974555) q[2];
sx q[2];
rz(1.6050715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8158055) q[1];
sx q[1];
rz(-2.1539638) q[1];
sx q[1];
rz(1.6156242) q[1];
rz(0.1707465) q[3];
sx q[3];
rz(-2.5617544) q[3];
sx q[3];
rz(2.8454526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59447294) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(-0.42868844) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73228943) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(-1.5111142) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139075) q[0];
sx q[0];
rz(-1.5965009) q[0];
sx q[0];
rz(-1.684042) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0610602) q[2];
sx q[2];
rz(-2.3057612) q[2];
sx q[2];
rz(0.9860952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1765602) q[1];
sx q[1];
rz(-1.5000048) q[1];
sx q[1];
rz(-1.895158) q[1];
rz(-pi) q[2];
rz(0.24095778) q[3];
sx q[3];
rz(-1.3853419) q[3];
sx q[3];
rz(-2.0897739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80031359) q[2];
sx q[2];
rz(-2.6276402) q[2];
sx q[2];
rz(-0.094816118) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6288347) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(-1.9376471) q[2];
sx q[2];
rz(-2.9137712) q[2];
sx q[2];
rz(-3.055345) q[2];
rz(2.6125828) q[3];
sx q[3];
rz(-2.257824) q[3];
sx q[3];
rz(-2.7469225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
