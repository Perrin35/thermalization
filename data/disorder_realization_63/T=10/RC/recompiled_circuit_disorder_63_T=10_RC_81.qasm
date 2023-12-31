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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4608085) q[0];
sx q[0];
rz(-0.55235282) q[0];
sx q[0];
rz(-2.0244563) q[0];
rz(-pi) q[1];
rz(3.0739003) q[2];
sx q[2];
rz(-2.2798685) q[2];
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
rz(1.2823892) q[1];
sx q[1];
rz(-0.61835641) q[1];
sx q[1];
rz(-2.8137141) q[1];
x q[2];
rz(-2.4077971) q[3];
sx q[3];
rz(-1.194209) q[3];
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
rz(-1.0268964) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.3902364) q[0];
rz(-0.83579666) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-2.4332411) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7824552) q[0];
sx q[0];
rz(-1.8498427) q[0];
sx q[0];
rz(2.7611087) q[0];
x q[1];
rz(1.7827665) q[2];
sx q[2];
rz(-0.6119298) q[2];
sx q[2];
rz(1.7906534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0201455) q[1];
sx q[1];
rz(-2.2312806) q[1];
sx q[1];
rz(-2.5026908) q[1];
x q[2];
rz(2.3362818) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.3290688) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(0.30971757) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(-0.57166878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6310731) q[0];
sx q[0];
rz(-1.6594995) q[0];
sx q[0];
rz(-0.53511329) q[0];
rz(0.35927202) q[2];
sx q[2];
rz(-1.1604939) q[2];
sx q[2];
rz(1.8112195) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.5773147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-0.70880115) q[2];
rz(-0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66089345) q[0];
sx q[0];
rz(-1.1850712) q[0];
sx q[0];
rz(1.6688523) q[0];
rz(2.7724491) q[2];
sx q[2];
rz(-1.7338848) q[2];
sx q[2];
rz(0.8085608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4304639) q[1];
sx q[1];
rz(-2.8006878) q[1];
sx q[1];
rz(-0.35699637) q[1];
x q[2];
rz(1.5482076) q[3];
sx q[3];
rz(-0.49593192) q[3];
sx q[3];
rz(-0.13726928) q[3];
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
rz(1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(0.25156897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35206301) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(2.2625838) q[0];
rz(-pi) q[1];
rz(-1.6576084) q[2];
sx q[2];
rz(-1.112339) q[2];
sx q[2];
rz(2.364033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1269762) q[1];
sx q[1];
rz(-1.3754305) q[1];
sx q[1];
rz(-1.5986534) q[1];
rz(-pi) q[2];
rz(-2.6336446) q[3];
sx q[3];
rz(-0.93379279) q[3];
sx q[3];
rz(-0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.183737) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(-2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0133936) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(-1.2584932) q[0];
rz(-pi) q[1];
rz(2.3153147) q[2];
sx q[2];
rz(-1.3958566) q[2];
sx q[2];
rz(-0.28184055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45822083) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(-1.5618192) q[1];
rz(-pi) q[2];
rz(2.8443908) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3605911) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(0.0059676776) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25516605) q[0];
sx q[0];
rz(-2.2212941) q[0];
sx q[0];
rz(2.7889473) q[0];
rz(2.0231831) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(-1.5024904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90970627) q[1];
sx q[1];
rz(-1.5071553) q[1];
sx q[1];
rz(1.015889) q[1];
rz(-pi) q[2];
rz(-1.5238477) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(-0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(2.9515284) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(-2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-2.1933864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873916) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(2.3733634) q[0];
x q[1];
rz(-2.6981632) q[2];
sx q[2];
rz(-1.0955053) q[2];
sx q[2];
rz(2.7623917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0046237) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(1.6349413) q[1];
rz(-pi) q[2];
rz(2.9392397) q[3];
sx q[3];
rz(-0.36365299) q[3];
sx q[3];
rz(2.3285151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-0.25094029) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7964325) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(0.81859421) q[0];
rz(-pi) q[1];
rz(1.7151095) q[2];
sx q[2];
rz(-1.2735954) q[2];
sx q[2];
rz(-0.65629634) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7749412) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(1.8371546) q[1];
x q[2];
rz(2.6412233) q[3];
sx q[3];
rz(-0.89958588) q[3];
sx q[3];
rz(2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-0.93651071) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6694326) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(2.6955312) q[0];
x q[1];
rz(2.9591987) q[2];
sx q[2];
rz(-1.2904022) q[2];
sx q[2];
rz(-0.80255752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22589707) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(2.8812863) q[1];
x q[2];
rz(-2.1533222) q[3];
sx q[3];
rz(-1.3068145) q[3];
sx q[3];
rz(-1.1519983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.3199432) q[2];
sx q[2];
rz(-1.432042) q[2];
sx q[2];
rz(-1.5114573) q[2];
rz(-1.3068009) q[3];
sx q[3];
rz(-2.0859857) q[3];
sx q[3];
rz(1.2386238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
