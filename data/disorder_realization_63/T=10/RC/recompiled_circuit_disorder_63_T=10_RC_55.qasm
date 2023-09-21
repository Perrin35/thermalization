OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9665943) q[0];
sx q[0];
rz(-2.7881665) q[0];
sx q[0];
rz(2.0768291) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4608085) q[0];
sx q[0];
rz(-0.55235282) q[0];
sx q[0];
rz(-2.0244563) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0739003) q[2];
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
rz(0.88692611) q[1];
sx q[1];
rz(-0.98985043) q[1];
sx q[1];
rz(1.3455774) q[1];
x q[2];
rz(-2.4077971) q[3];
sx q[3];
rz(-1.194209) q[3];
sx q[3];
rz(-0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(-0.25201592) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(0.70835152) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7824552) q[0];
sx q[0];
rz(-1.8498427) q[0];
sx q[0];
rz(2.7611087) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96946851) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(-0.045493424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9006151) q[1];
sx q[1];
rz(-0.88417378) q[1];
sx q[1];
rz(-2.2254506) q[1];
rz(-2.0833011) q[3];
sx q[3];
rz(-2.0111492) q[3];
sx q[3];
rz(-3.1194221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(0.29176816) q[2];
rz(-0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3290688) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(-0.57166878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088061995) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(0.17266973) q[0];
x q[1];
rz(2.2505635) q[2];
sx q[2];
rz(-2.6030412) q[2];
sx q[2];
rz(-2.0856759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7844312) q[1];
sx q[1];
rz(-0.845134) q[1];
sx q[1];
rz(-1.7482589) q[1];
rz(-pi) q[2];
rz(2.8623657) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(2.4327915) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.18049151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7364396) q[0];
sx q[0];
rz(-2.7442051) q[0];
sx q[0];
rz(-2.9050164) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36914354) q[2];
sx q[2];
rz(-1.7338848) q[2];
sx q[2];
rz(0.8085608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0879678) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(-1.4474523) q[1];
x q[2];
rz(-2.0666215) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(-1.4533952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.8245274) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809526) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-2.8900237) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0753206) q[0];
sx q[0];
rz(-0.89162725) q[0];
sx q[0];
rz(0.22596304) q[0];
x q[1];
rz(-1.6576084) q[2];
sx q[2];
rz(-2.0292536) q[2];
sx q[2];
rz(-2.364033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.561589) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(2.9461529) q[1];
rz(-pi) q[2];
rz(2.6336446) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(2.7579443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-1.0305369) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(-0.39302557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61009208) q[0];
sx q[0];
rz(-1.8787662) q[0];
sx q[0];
rz(-2.969335) q[0];
x q[1];
rz(1.3156462) q[2];
sx q[2];
rz(-0.7609376) q[2];
sx q[2];
rz(1.6659425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46950754) q[1];
sx q[1];
rz(-0.91963327) q[1];
sx q[1];
rz(3.1347515) q[1];
x q[2];
rz(-2.399746) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(2.0239995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(-0.97314107) q[2];
rz(0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-0.0059676776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(1.1447385) q[0];
rz(-pi) q[1];
rz(-2.6628394) q[2];
sx q[2];
rz(-2.3294449) q[2];
sx q[2];
rz(-2.2854545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70049268) q[1];
sx q[1];
rz(-1.0171434) q[1];
sx q[1];
rz(0.074837491) q[1];
rz(-pi) q[2];
rz(1.5238477) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(-2.9058128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(-2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95420102) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(2.3733634) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6981632) q[2];
sx q[2];
rz(-1.0955053) q[2];
sx q[2];
rz(-0.37920096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72045418) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(-2.8645664) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35685278) q[3];
sx q[3];
rz(-1.4992504) q[3];
sx q[3];
rz(-2.1944291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-0.25094029) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1044554) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(-3.0260968) q[0];
x q[1];
rz(1.4264832) q[2];
sx q[2];
rz(-1.8679973) q[2];
sx q[2];
rz(-0.65629634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.346631) q[1];
sx q[1];
rz(-1.3449841) q[1];
sx q[1];
rz(-2.5717616) q[1];
rz(-pi) q[2];
rz(0.8351164) q[3];
sx q[3];
rz(-1.1856688) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6694326) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(-0.44606146) q[0];
rz(-pi) q[1];
rz(-1.2859225) q[2];
sx q[2];
rz(-1.3956009) q[2];
sx q[2];
rz(0.81923649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22589707) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(2.8812863) q[1];
rz(-pi) q[2];
rz(-2.0274721) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(-3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(-3.1266406) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(2.0832534) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(1.3068009) q[3];
sx q[3];
rz(-1.055607) q[3];
sx q[3];
rz(-1.9029688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];