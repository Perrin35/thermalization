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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(1.0648849) q[0];
x q[1];
rz(-1.4921161) q[2];
sx q[2];
rz(-0.71173758) q[2];
sx q[2];
rz(0.56456883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8592035) q[1];
sx q[1];
rz(-0.61835641) q[1];
sx q[1];
rz(-2.8137141) q[1];
rz(-0.73379559) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(2.9154582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(2.1146963) q[2];
rz(-0.25201592) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(1.3902364) q[0];
rz(-0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-0.70835152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81472155) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(0.65713709) q[0];
rz(0.14658908) q[2];
sx q[2];
rz(-2.1671038) q[2];
sx q[2];
rz(1.5335611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0201455) q[1];
sx q[1];
rz(-0.91031204) q[1];
sx q[1];
rz(2.5026908) q[1];
rz(-pi) q[2];
rz(0.80531081) q[3];
sx q[3];
rz(-0.66262965) q[3];
sx q[3];
rz(-2.2413072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(0.29176816) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.5244938) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3290688) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(0.57166878) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088061995) q[0];
sx q[0];
rz(-0.54170875) q[0];
sx q[0];
rz(0.17266973) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2505635) q[2];
sx q[2];
rz(-0.53855145) q[2];
sx q[2];
rz(1.0559168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7844312) q[1];
sx q[1];
rz(-0.845134) q[1];
sx q[1];
rz(-1.3933338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4324576) q[3];
sx q[3];
rz(-1.1225015) q[3];
sx q[3];
rz(1.5642779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(-0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4335094) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(0.18049151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9468964) q[0];
sx q[0];
rz(-1.6616271) q[0];
sx q[0];
rz(2.7541861) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7137202) q[2];
sx q[2];
rz(-2.7395436) q[2];
sx q[2];
rz(2.7767162) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0879678) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(-1.4474523) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5482076) q[3];
sx q[3];
rz(-0.49593192) q[3];
sx q[3];
rz(0.13726928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(-1.5578516) q[2];
rz(1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(1.1622693) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(0.25156897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7154327) q[0];
sx q[0];
rz(-0.71007219) q[0];
sx q[0];
rz(-1.8415113) q[0];
rz(2.6816363) q[2];
sx q[2];
rz(-1.6486247) q[2];
sx q[2];
rz(0.83173448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.561589) q[1];
sx q[1];
rz(-1.5434693) q[1];
sx q[1];
rz(0.1954397) q[1];
x q[2];
rz(-0.86815636) q[3];
sx q[3];
rz(-1.9725102) q[3];
sx q[3];
rz(0.86740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(2.2914698) q[0];
rz(-1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(2.7485671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0133936) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(-1.2584932) q[0];
rz(-pi) q[1];
rz(1.8259465) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(-1.4756502) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6720851) q[1];
sx q[1];
rz(-0.91963327) q[1];
sx q[1];
rz(-3.1347515) q[1];
x q[2];
rz(1.2971446) q[3];
sx q[3];
rz(-1.8575462) q[3];
sx q[3];
rz(-1.2424038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(2.1684516) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(3.135625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0963421) q[0];
sx q[0];
rz(-1.8492286) q[0];
sx q[0];
rz(-2.2521426) q[0];
rz(-pi) q[1];
rz(-0.75255021) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(-0.37170751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5827427) q[1];
sx q[1];
rz(-0.55816459) q[1];
sx q[1];
rz(1.4504257) q[1];
x q[2];
rz(-3.0217516) q[3];
sx q[3];
rz(-2.7671475) q[3];
sx q[3];
rz(-3.0344809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(-0.93332943) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873916) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(-2.3733634) q[0];
x q[1];
rz(0.44342946) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(0.37920096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1369689) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(1.6349413) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20235297) q[3];
sx q[3];
rz(-2.7779397) q[3];
sx q[3];
rz(2.3285151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(0.25094029) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964325) q[0];
sx q[0];
rz(-2.9734018) q[0];
sx q[0];
rz(0.81859421) q[0];
rz(-2.7025931) q[2];
sx q[2];
rz(-0.3294496) q[2];
sx q[2];
rz(-2.0246558) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11201227) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(2.7390202) q[1];
rz(-0.8351164) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(1.3002243) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47637981) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(0.54540821) q[0];
rz(-pi) q[1];
rz(2.1328743) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(2.9269232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22589707) q[1];
sx q[1];
rz(-0.86663336) q[1];
sx q[1];
rz(2.8812863) q[1];
rz(0.31302932) q[3];
sx q[3];
rz(-2.1306681) q[3];
sx q[3];
rz(2.5525639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-3.1266406) q[2];
rz(0.22710083) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.2623513) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1508355) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(1.8216495) q[2];
sx q[2];
rz(-1.7095507) q[2];
sx q[2];
rz(1.6301353) q[2];
rz(-0.43185497) q[3];
sx q[3];
rz(-2.5681744) q[3];
sx q[3];
rz(-1.4011866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
