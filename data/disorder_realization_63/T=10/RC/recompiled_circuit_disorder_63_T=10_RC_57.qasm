OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(-0.35342616) q[0];
sx q[0];
rz(1.0647635) q[0];
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6807841) q[0];
sx q[0];
rz(-0.55235282) q[0];
sx q[0];
rz(1.1171364) q[0];
x q[1];
rz(0.067692368) q[2];
sx q[2];
rz(-0.86172416) q[2];
sx q[2];
rz(-0.46082218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55878996) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(2.5488528) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73379559) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(-2.2759329) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.3902364) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-2.4332411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10193292) q[0];
sx q[0];
rz(-1.205737) q[0];
sx q[0];
rz(-1.870116) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3588261) q[2];
sx q[2];
rz(-2.5296629) q[2];
sx q[2];
rz(1.3509392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0185512) q[1];
sx q[1];
rz(-2.0611144) q[1];
sx q[1];
rz(2.3398188) q[1];
rz(-pi) q[2];
rz(-2.0833011) q[3];
sx q[3];
rz(-1.1304434) q[3];
sx q[3];
rz(3.1194221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(-1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-0.57166878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0535307) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(-0.17266973) q[0];
rz(-pi) q[1];
rz(2.7823206) q[2];
sx q[2];
rz(-1.1604939) q[2];
sx q[2];
rz(-1.8112195) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0464222) q[1];
sx q[1];
rz(-1.7032402) q[1];
sx q[1];
rz(-0.73352791) q[1];
rz(-pi) q[2];
rz(-2.6895371) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(0.053754036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-0.70880115) q[2];
rz(-0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(0.79950142) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-2.9611011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9468964) q[0];
sx q[0];
rz(-1.6616271) q[0];
sx q[0];
rz(-2.7541861) q[0];
x q[1];
rz(0.36914354) q[2];
sx q[2];
rz(-1.7338848) q[2];
sx q[2];
rz(-0.8085608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0536249) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(1.4474523) q[1];
rz(-pi) q[2];
rz(3.1293731) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(-1.5578516) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809526) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(2.8900237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0753206) q[0];
sx q[0];
rz(-2.2499654) q[0];
sx q[0];
rz(2.9156296) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45995633) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(0.83173448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1269762) q[1];
sx q[1];
rz(-1.3754305) q[1];
sx q[1];
rz(1.5429392) q[1];
rz(-2.6336446) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.75491607) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(2.7485671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0133936) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(1.8830995) q[0];
rz(-pi) q[1];
rz(0.82627798) q[2];
sx q[2];
rz(-1.3958566) q[2];
sx q[2];
rz(0.28184055) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0361573) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(2.2219707) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.844448) q[3];
sx q[3];
rz(-1.2840464) q[3];
sx q[3];
rz(-1.8991889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-2.1940103) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-0.0059676776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80124827) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(-1.9968541) q[0];
rz(2.6628394) q[2];
sx q[2];
rz(-0.81214777) q[2];
sx q[2];
rz(0.85613814) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2318864) q[1];
sx q[1];
rz(-1.5071553) q[1];
sx q[1];
rz(-2.1257036) q[1];
x q[2];
rz(2.7695914) q[3];
sx q[3];
rz(-1.5270546) q[3];
sx q[3];
rz(1.7895167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(-2.2082632) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(2.1933864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78806879) q[0];
sx q[0];
rz(-2.3234962) q[0];
sx q[0];
rz(-1.3226932) q[0];
rz(-pi) q[1];
rz(0.87585978) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(-1.1832331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58029786) q[1];
sx q[1];
rz(-1.5082238) q[1];
sx q[1];
rz(-0.22217521) q[1];
rz(-0.35685278) q[3];
sx q[3];
rz(-1.6423422) q[3];
sx q[3];
rz(2.1944291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.3482288) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76294476) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(0.25094029) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-2.7889263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6610722) q[0];
sx q[0];
rz(-1.4561704) q[0];
sx q[0];
rz(1.4474439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7151095) q[2];
sx q[2];
rz(-1.2735954) q[2];
sx q[2];
rz(0.65629634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7949617) q[1];
sx q[1];
rz(-1.7966086) q[1];
sx q[1];
rz(2.5717616) q[1];
rz(2.114186) q[3];
sx q[3];
rz(-2.3282442) q[3];
sx q[3];
rz(-3.120595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(0.68332589) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(-1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-2.2050819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6652128) q[0];
sx q[0];
rz(-1.8139651) q[0];
sx q[0];
rz(2.5961844) q[0];
x q[1];
rz(-1.0087183) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(2.9269232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22589707) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(0.26030635) q[1];
x q[2];
rz(0.31302932) q[3];
sx q[3];
rz(-2.1306681) q[3];
sx q[3];
rz(2.5525639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-0.014952095) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(2.0832534) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(-2.6111504) q[3];
sx q[3];
rz(-1.3417288) q[3];
sx q[3];
rz(-0.19977278) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
