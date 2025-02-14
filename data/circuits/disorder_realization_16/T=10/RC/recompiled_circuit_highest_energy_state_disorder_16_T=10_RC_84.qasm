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
rz(1.9749405) q[0];
sx q[0];
rz(-2.6986172) q[0];
sx q[0];
rz(-2.7583581) q[0];
rz(1.8546328) q[1];
sx q[1];
rz(-1.9998963) q[1];
sx q[1];
rz(0.72007522) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135277) q[0];
sx q[0];
rz(-1.1733455) q[0];
sx q[0];
rz(1.9883519) q[0];
rz(-pi) q[1];
rz(-1.7699446) q[2];
sx q[2];
rz(-1.7590283) q[2];
sx q[2];
rz(-3.0668099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9125058) q[1];
sx q[1];
rz(-2.6709963) q[1];
sx q[1];
rz(-2.5060593) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3377417) q[3];
sx q[3];
rz(-0.62352288) q[3];
sx q[3];
rz(-0.064176809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8522475) q[2];
sx q[2];
rz(-2.2699558) q[2];
sx q[2];
rz(2.3343425) q[2];
rz(-1.4489669) q[3];
sx q[3];
rz(-0.92985409) q[3];
sx q[3];
rz(0.71678954) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221231) q[0];
sx q[0];
rz(-1.8417646) q[0];
sx q[0];
rz(1.8362554) q[0];
rz(2.3985825) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(1.8441127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8332307) q[0];
sx q[0];
rz(-1.9675281) q[0];
sx q[0];
rz(0.19773592) q[0];
x q[1];
rz(-2.1427311) q[2];
sx q[2];
rz(-1.9543658) q[2];
sx q[2];
rz(0.072173031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0875133) q[1];
sx q[1];
rz(-2.1554408) q[1];
sx q[1];
rz(-2.784009) q[1];
x q[2];
rz(-0.36139523) q[3];
sx q[3];
rz(-1.6797429) q[3];
sx q[3];
rz(1.4055085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3798736) q[2];
sx q[2];
rz(-2.1642809) q[2];
sx q[2];
rz(-2.1011815) q[2];
rz(-2.8271293) q[3];
sx q[3];
rz(-1.209126) q[3];
sx q[3];
rz(1.8837121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8085025) q[0];
sx q[0];
rz(-2.1603778) q[0];
sx q[0];
rz(-1.8044949) q[0];
rz(0.62726504) q[1];
sx q[1];
rz(-1.7341055) q[1];
sx q[1];
rz(-1.6076535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21625265) q[0];
sx q[0];
rz(-2.0006764) q[0];
sx q[0];
rz(2.2735808) q[0];
x q[1];
rz(2.5467403) q[2];
sx q[2];
rz(-2.640758) q[2];
sx q[2];
rz(-0.55933096) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.99974224) q[1];
sx q[1];
rz(-2.6313884) q[1];
sx q[1];
rz(2.7562618) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19948761) q[3];
sx q[3];
rz(-1.9428245) q[3];
sx q[3];
rz(-0.92687273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4072121) q[2];
sx q[2];
rz(-0.89443365) q[2];
sx q[2];
rz(-1.9341932) q[2];
rz(2.1842128) q[3];
sx q[3];
rz(-1.9060241) q[3];
sx q[3];
rz(0.68937075) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44416881) q[0];
sx q[0];
rz(-0.95132315) q[0];
sx q[0];
rz(-0.091766894) q[0];
rz(0.16608876) q[1];
sx q[1];
rz(-1.2867915) q[1];
sx q[1];
rz(-1.0413569) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7163775) q[0];
sx q[0];
rz(-1.5272523) q[0];
sx q[0];
rz(-2.3460044) q[0];
x q[1];
rz(1.9729593) q[2];
sx q[2];
rz(-2.3170174) q[2];
sx q[2];
rz(-0.74092016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3709621) q[1];
sx q[1];
rz(-0.47415552) q[1];
sx q[1];
rz(-1.6459792) q[1];
rz(3.0001767) q[3];
sx q[3];
rz(-1.4849471) q[3];
sx q[3];
rz(-0.72535634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.912821) q[2];
sx q[2];
rz(-0.81835881) q[2];
sx q[2];
rz(-1.5405601) q[2];
rz(1.0176963) q[3];
sx q[3];
rz(-1.8120268) q[3];
sx q[3];
rz(-1.2525919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3200662) q[0];
sx q[0];
rz(-1.3952661) q[0];
sx q[0];
rz(-2.9030002) q[0];
rz(0.78856167) q[1];
sx q[1];
rz(-2.8193654) q[1];
sx q[1];
rz(-2.2408748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64553211) q[0];
sx q[0];
rz(-2.0010316) q[0];
sx q[0];
rz(1.2853464) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9600257) q[2];
sx q[2];
rz(-1.4837564) q[2];
sx q[2];
rz(0.6861251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8713995) q[1];
sx q[1];
rz(-2.7084623) q[1];
sx q[1];
rz(-0.30850839) q[1];
x q[2];
rz(-0.78910335) q[3];
sx q[3];
rz(-2.2450617) q[3];
sx q[3];
rz(-1.2885119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0476524) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(-1.368604) q[2];
rz(-2.7641344) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(-1.9774168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222328) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(-2.9121616) q[0];
rz(0.58652985) q[1];
sx q[1];
rz(-0.49085453) q[1];
sx q[1];
rz(-1.4803001) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23585715) q[0];
sx q[0];
rz(-0.93556306) q[0];
sx q[0];
rz(-2.6978289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56917135) q[2];
sx q[2];
rz(-1.7915951) q[2];
sx q[2];
rz(-2.6463531) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52119285) q[1];
sx q[1];
rz(-1.344465) q[1];
sx q[1];
rz(1.5160376) q[1];
x q[2];
rz(1.3020206) q[3];
sx q[3];
rz(-0.69221262) q[3];
sx q[3];
rz(-0.48066329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71969879) q[2];
sx q[2];
rz(-2.5514558) q[2];
sx q[2];
rz(-1.7936919) q[2];
rz(2.9501996) q[3];
sx q[3];
rz(-0.32262155) q[3];
sx q[3];
rz(-1.2389244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557068) q[0];
sx q[0];
rz(-1.5140336) q[0];
sx q[0];
rz(-0.53644449) q[0];
rz(-0.084511936) q[1];
sx q[1];
rz(-2.5630496) q[1];
sx q[1];
rz(1.9523778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41302785) q[0];
sx q[0];
rz(-2.4997093) q[0];
sx q[0];
rz(-1.8299787) q[0];
rz(-pi) q[1];
rz(0.68185735) q[2];
sx q[2];
rz(-1.8769662) q[2];
sx q[2];
rz(0.078753565) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.099965) q[1];
sx q[1];
rz(-2.065127) q[1];
sx q[1];
rz(-2.7848141) q[1];
x q[2];
rz(-2.9042894) q[3];
sx q[3];
rz(-2.0198698) q[3];
sx q[3];
rz(-2.9782691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1875923) q[2];
sx q[2];
rz(-1.2728929) q[2];
sx q[2];
rz(-0.083258955) q[2];
rz(1.8795053) q[3];
sx q[3];
rz(-2.3048461) q[3];
sx q[3];
rz(0.73778233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5543723) q[0];
sx q[0];
rz(-1.502259) q[0];
sx q[0];
rz(-2.3578405) q[0];
rz(-1.4844249) q[1];
sx q[1];
rz(-0.52394167) q[1];
sx q[1];
rz(2.8481683) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2121289) q[0];
sx q[0];
rz(-2.0269454) q[0];
sx q[0];
rz(-0.41123234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90250315) q[2];
sx q[2];
rz(-1.74729) q[2];
sx q[2];
rz(0.89418156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3345016) q[1];
sx q[1];
rz(-1.1576022) q[1];
sx q[1];
rz(2.9838802) q[1];
rz(-pi) q[2];
rz(1.8752304) q[3];
sx q[3];
rz(-2.4633411) q[3];
sx q[3];
rz(1.1224318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7744814) q[2];
sx q[2];
rz(-1.7256871) q[2];
sx q[2];
rz(1.6488546) q[2];
rz(-1.1501009) q[3];
sx q[3];
rz(-0.29905683) q[3];
sx q[3];
rz(-2.2404631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8915326) q[0];
sx q[0];
rz(-1.0279259) q[0];
sx q[0];
rz(0.0082536396) q[0];
rz(3.1390269) q[1];
sx q[1];
rz(-2.6305514) q[1];
sx q[1];
rz(2.0770238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7637007) q[0];
sx q[0];
rz(-2.3449909) q[0];
sx q[0];
rz(-2.5515351) q[0];
rz(-pi) q[1];
rz(1.6029036) q[2];
sx q[2];
rz(-2.8920482) q[2];
sx q[2];
rz(1.7966934) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4618116) q[1];
sx q[1];
rz(-1.0625328) q[1];
sx q[1];
rz(-2.630788) q[1];
rz(-0.58925924) q[3];
sx q[3];
rz(-1.6943101) q[3];
sx q[3];
rz(2.6753787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0801487) q[2];
sx q[2];
rz(-1.4465569) q[2];
sx q[2];
rz(-2.9385369) q[2];
rz(1.2114581) q[3];
sx q[3];
rz(-1.0695499) q[3];
sx q[3];
rz(-3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9515297) q[0];
sx q[0];
rz(-2.0385346) q[0];
sx q[0];
rz(-0.67361012) q[0];
rz(-0.31189648) q[1];
sx q[1];
rz(-1.2047647) q[1];
sx q[1];
rz(-1.9685251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6098227) q[0];
sx q[0];
rz(-1.8635475) q[0];
sx q[0];
rz(1.0171153) q[0];
x q[1];
rz(-2.3854957) q[2];
sx q[2];
rz(-1.399516) q[2];
sx q[2];
rz(1.1403699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58880776) q[1];
sx q[1];
rz(-1.9845595) q[1];
sx q[1];
rz(0.082774712) q[1];
rz(2.0579751) q[3];
sx q[3];
rz(-0.77278944) q[3];
sx q[3];
rz(-0.98675283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6761026) q[2];
sx q[2];
rz(-2.5238621) q[2];
sx q[2];
rz(2.9205196) q[2];
rz(0.57340932) q[3];
sx q[3];
rz(-2.2677877) q[3];
sx q[3];
rz(-1.1055841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40972805) q[0];
sx q[0];
rz(-1.8479713) q[0];
sx q[0];
rz(0.027298409) q[0];
rz(1.191054) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(1.2192192) q[2];
sx q[2];
rz(-0.58813358) q[2];
sx q[2];
rz(1.0434601) q[2];
rz(2.7577362) q[3];
sx q[3];
rz(-1.0105269) q[3];
sx q[3];
rz(-0.15906048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
