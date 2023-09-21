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
rz(3.937768) q[1];
sx q[1];
rz(1.9328971) q[1];
sx q[1];
rz(9.9608496) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(-2.0767077) q[0];
rz(2.2810018) q[2];
sx q[2];
rz(-1.5194367) q[2];
sx q[2];
rz(-1.0658588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5828027) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(2.5488528) q[1];
rz(0.73379559) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(-2.9154582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(1.3902364) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(0.70835152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81472155) q[0];
sx q[0];
rz(-2.673809) q[0];
sx q[0];
rz(0.65713709) q[0];
rz(2.1721241) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(0.045493424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1214472) q[1];
sx q[1];
rz(-0.91031204) q[1];
sx q[1];
rz(0.63890181) q[1];
rz(2.6459341) q[3];
sx q[3];
rz(-2.0303876) q[3];
sx q[3];
rz(-1.3575777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-2.8498245) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(-1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(0.30971757) q[0];
rz(1.4801056) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(-2.5699239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6310731) q[0];
sx q[0];
rz(-1.4820931) q[0];
sx q[0];
rz(-0.53511329) q[0];
rz(-pi) q[1];
rz(2.7823206) q[2];
sx q[2];
rz(-1.9810988) q[2];
sx q[2];
rz(1.8112195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.09517041) q[1];
sx q[1];
rz(-1.4383525) q[1];
sx q[1];
rz(-2.4080647) q[1];
rz(-0.27922697) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(3.0917621) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-2.9611011) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364396) q[0];
sx q[0];
rz(-2.7442051) q[0];
sx q[0];
rz(-0.23657628) q[0];
rz(2.7724491) q[2];
sx q[2];
rz(-1.4077079) q[2];
sx q[2];
rz(-0.8085608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4304639) q[1];
sx q[1];
rz(-2.8006878) q[1];
sx q[1];
rz(2.7845963) q[1];
rz(-3.1293731) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(0.11158768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(1.7662988) q[3];
sx q[3];
rz(-1.8245274) q[3];
sx q[3];
rz(0.93262514) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809526) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-0.25156897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4261599) q[0];
sx q[0];
rz(-2.4315205) q[0];
sx q[0];
rz(-1.8415113) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45995633) q[2];
sx q[2];
rz(-1.6486247) q[2];
sx q[2];
rz(2.3098582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5800036) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(-2.9461529) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1523347) q[3];
sx q[3];
rz(-0.79205081) q[3];
sx q[3];
rz(2.0056412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.183737) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(-1.7306227) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(-1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(2.7485671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.052554) q[0];
sx q[0];
rz(-0.35152838) q[0];
sx q[0];
rz(-1.0765443) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3156462) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(1.6659425) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45822083) q[1];
sx q[1];
rz(-0.65119377) q[1];
sx q[1];
rz(1.5797735) q[1];
x q[2];
rz(2.8443908) q[3];
sx q[3];
rz(-1.8330049) q[3];
sx q[3];
rz(-0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3605911) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(-2.1684516) q[2];
rz(-2.1940103) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(-2.7667926) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(0.0059676776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(-1.9968541) q[0];
x q[1];
rz(1.1184095) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(1.5024904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4411) q[1];
sx q[1];
rz(-2.1244493) q[1];
sx q[1];
rz(3.0667552) q[1];
rz(-pi) q[2];
rz(-0.11984101) q[3];
sx q[3];
rz(-2.7671475) q[3];
sx q[3];
rz(3.0344809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.5678762) q[1];
sx q[1];
rz(0.94820625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084003) q[0];
sx q[0];
rz(-0.78484479) q[0];
sx q[0];
rz(2.8851896) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0886705) q[2];
sx q[2];
rz(-1.9621984) q[2];
sx q[2];
rz(2.1640167) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1369689) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(1.5066513) q[1];
rz(2.9392397) q[3];
sx q[3];
rz(-2.7779397) q[3];
sx q[3];
rz(0.81307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.3482288) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-0.33319107) q[0];
sx q[0];
rz(-0.25094029) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1044554) q[0];
sx q[0];
rz(-1.6933352) q[0];
sx q[0];
rz(-0.11549581) q[0];
rz(-0.43899957) q[2];
sx q[2];
rz(-0.3294496) q[2];
sx q[2];
rz(-1.1169369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7949617) q[1];
sx q[1];
rz(-1.7966086) q[1];
sx q[1];
rz(-2.5717616) q[1];
rz(-pi) q[2];
rz(-0.8351164) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(0.68332589) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47637981) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(2.5961844) q[0];
x q[1];
rz(-1.0087183) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(2.9269232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(0.84968062) q[1];
rz(-pi) q[2];
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
rz(0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-0.014952095) q[2];
rz(0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.8792413) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1508355) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(-1.5630209) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(0.14317748) q[2];
sx q[2];
rz(-1.3224052) q[2];
sx q[2];
rz(0.094766141) q[2];
rz(-0.53044226) q[3];
sx q[3];
rz(-1.7998639) q[3];
sx q[3];
rz(2.9418199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
