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
rz(-2.0768291) q[0];
rz(3.937768) q[1];
sx q[1];
rz(1.9328971) q[1];
sx q[1];
rz(9.9608496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(-1.0648849) q[0];
rz(1.6494765) q[2];
sx q[2];
rz(-0.71173758) q[2];
sx q[2];
rz(0.56456883) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5828027) q[1];
sx q[1];
rz(-1.7585808) q[1];
sx q[1];
rz(0.59273984) q[1];
rz(-pi) q[2];
rz(0.53341289) q[3];
sx q[3];
rz(-2.3331113) q[3];
sx q[3];
rz(-1.7318788) q[3];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.3902364) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-2.4332411) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10193292) q[0];
sx q[0];
rz(-1.9358557) q[0];
sx q[0];
rz(-1.870116) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96946851) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(-3.0960992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(2.3398188) q[1];
rz(-pi) q[2];
rz(-0.80531081) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(-2.2413072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-2.34989) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(2.5699239) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5105195) q[0];
sx q[0];
rz(-1.4820931) q[0];
sx q[0];
rz(2.6064794) q[0];
rz(-pi) q[1];
rz(-2.0057726) q[2];
sx q[2];
rz(-1.8990714) q[2];
sx q[2];
rz(3.0498743) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5204822) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(-0.1964257) q[1];
rz(-0.27922697) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(-0.79950142) q[0];
rz(3.0917621) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66089345) q[0];
sx q[0];
rz(-1.1850712) q[0];
sx q[0];
rz(1.6688523) q[0];
x q[1];
rz(-2.7137202) q[2];
sx q[2];
rz(-0.40204907) q[2];
sx q[2];
rz(2.7767162) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6199477) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(2.820693) q[1];
x q[2];
rz(-3.1293731) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(0.11158768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-0.25156897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35206301) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(-2.2625838) q[0];
x q[1];
rz(-2.6816363) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(0.83173448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.561589) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(-0.1954397) q[1];
rz(-pi) q[2];
rz(2.2734363) q[3];
sx q[3];
rz(-1.1690825) q[3];
sx q[3];
rz(2.274184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5315006) q[0];
sx q[0];
rz(-1.8787662) q[0];
sx q[0];
rz(-2.969335) q[0];
x q[1];
rz(0.23586258) q[2];
sx q[2];
rz(-0.84025192) q[2];
sx q[2];
rz(-1.1299396) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0361573) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(-2.2219707) q[1];
x q[2];
rz(-1.844448) q[3];
sx q[3];
rz(-1.2840464) q[3];
sx q[3];
rz(-1.8991889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3605911) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(-2.1684516) q[2];
rz(0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25516605) q[0];
sx q[0];
rz(-2.2212941) q[0];
sx q[0];
rz(0.35264539) q[0];
rz(-2.0231831) q[2];
sx q[2];
rz(-0.87086073) q[2];
sx q[2];
rz(1.6391022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90970627) q[1];
sx q[1];
rz(-1.6344374) q[1];
sx q[1];
rz(2.1257036) q[1];
rz(-pi) q[2];
rz(-3.0217516) q[3];
sx q[3];
rz(-0.37444515) q[3];
sx q[3];
rz(-0.10711174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(-0.73474187) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(-0.93332943) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(2.1933864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535239) q[0];
sx q[0];
rz(-2.3234962) q[0];
sx q[0];
rz(1.8188994) q[0];
rz(-0.44342946) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(-0.37920096) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5612948) q[1];
sx q[1];
rz(-1.5082238) q[1];
sx q[1];
rz(-0.22217521) q[1];
rz(-1.6471345) q[3];
sx q[3];
rz(-1.2148972) q[3];
sx q[3];
rz(0.59698856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(2.7889263) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48052045) q[0];
sx q[0];
rz(-1.6854223) q[0];
sx q[0];
rz(1.4474439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7151095) q[2];
sx q[2];
rz(-1.2735954) q[2];
sx q[2];
rz(0.65629634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11201227) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(-2.7390202) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50036939) q[3];
sx q[3];
rz(-2.2420068) q[3];
sx q[3];
rz(-2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3744366) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(1.3002243) q[0];
rz(0.70872778) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-0.93651071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6652128) q[0];
sx q[0];
rz(-1.8139651) q[0];
sx q[0];
rz(-0.54540821) q[0];
rz(-pi) q[1];
rz(0.18239393) q[2];
sx q[2];
rz(-1.2904022) q[2];
sx q[2];
rz(0.80255752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5156538) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(-0.84968062) q[1];
x q[2];
rz(1.1141206) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(-0.041681899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54291723) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(-1.5785718) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(-2.9984152) q[2];
sx q[2];
rz(-1.3224052) q[2];
sx q[2];
rz(0.094766141) q[2];
rz(2.6111504) q[3];
sx q[3];
rz(-1.7998639) q[3];
sx q[3];
rz(2.9418199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
