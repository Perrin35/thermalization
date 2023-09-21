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
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(-0.53607166) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85815) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(-1.0648849) q[0];
rz(-0.067692368) q[2];
sx q[2];
rz(-2.2798685) q[2];
sx q[2];
rz(-0.46082218) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2823892) q[1];
sx q[1];
rz(-2.5232362) q[1];
sx q[1];
rz(-2.8137141) q[1];
rz(2.6081798) q[3];
sx q[3];
rz(-2.3331113) q[3];
sx q[3];
rz(-1.4097139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7321695) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.3902364) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-2.4332411) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81472155) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(-2.4844556) q[0];
rz(-2.9950036) q[2];
sx q[2];
rz(-0.97448889) q[2];
sx q[2];
rz(-1.5335611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1214472) q[1];
sx q[1];
rz(-0.91031204) q[1];
sx q[1];
rz(2.5026908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0833011) q[3];
sx q[3];
rz(-1.1304434) q[3];
sx q[3];
rz(-0.022170598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(-3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(-1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(2.5699239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6310731) q[0];
sx q[0];
rz(-1.4820931) q[0];
sx q[0];
rz(-0.53511329) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.13582) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(3.0498743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5204822) q[1];
sx q[1];
rz(-0.74319327) q[1];
sx q[1];
rz(2.945167) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45205558) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(-0.053754036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16033515) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(-3.0917621) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-2.9611011) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364396) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(0.23657628) q[0];
rz(2.7724491) q[2];
sx q[2];
rz(-1.4077079) q[2];
sx q[2];
rz(-0.8085608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71112878) q[1];
sx q[1];
rz(-0.34090484) q[1];
sx q[1];
rz(-2.7845963) q[1];
rz(-3.1293731) q[3];
sx q[3];
rz(-1.0750024) q[3];
sx q[3];
rz(3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(-0.98651648) q[0];
rz(1.1622693) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-2.8900237) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0662721) q[0];
sx q[0];
rz(-0.89162725) q[0];
sx q[0];
rz(2.9156296) q[0];
rz(-pi) q[1];
rz(-0.17390522) q[2];
sx q[2];
rz(-2.675563) q[2];
sx q[2];
rz(-0.58338651) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5800036) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(2.9461529) q[1];
rz(2.2734363) q[3];
sx q[3];
rz(-1.1690825) q[3];
sx q[3];
rz(-0.86740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(-0.39302557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5315006) q[0];
sx q[0];
rz(-1.8787662) q[0];
sx q[0];
rz(2.969335) q[0];
rz(-2.3153147) q[2];
sx q[2];
rz(-1.7457361) q[2];
sx q[2];
rz(2.8597521) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46950754) q[1];
sx q[1];
rz(-2.2219594) q[1];
sx q[1];
rz(0.0068411946) q[1];
x q[2];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.6181035) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-0.0059676776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0452506) q[0];
sx q[0];
rz(-1.8492286) q[0];
sx q[0];
rz(-2.2521426) q[0];
x q[1];
rz(0.75255021) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(-2.7698851) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4411) q[1];
sx q[1];
rz(-1.0171434) q[1];
sx q[1];
rz(-3.0667552) q[1];
rz(1.617745) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(-0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(2.9515284) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0311325) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95420102) q[0];
sx q[0];
rz(-1.750995) q[0];
sx q[0];
rz(2.3733634) q[0];
rz(-0.87585978) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(1.1832331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0046237) q[1];
sx q[1];
rz(-1.7925295) q[1];
sx q[1];
rz(-1.5066513) q[1];
rz(-1.6471345) q[3];
sx q[3];
rz(-1.9266955) q[3];
sx q[3];
rz(-0.59698856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.640921) q[2];
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
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-0.35266638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0371373) q[0];
sx q[0];
rz(-1.6933352) q[0];
sx q[0];
rz(3.0260968) q[0];
rz(-pi) q[1];
rz(0.30013957) q[2];
sx q[2];
rz(-1.4328513) q[2];
sx q[2];
rz(-0.87196841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.346631) q[1];
sx q[1];
rz(-1.3449841) q[1];
sx q[1];
rz(-2.5717616) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50036939) q[3];
sx q[3];
rz(-0.89958588) q[3];
sx q[3];
rz(0.74218111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(-2.4328649) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47637981) q[0];
sx q[0];
rz(-1.8139651) q[0];
sx q[0];
rz(2.5961844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18239393) q[2];
sx q[2];
rz(-1.8511904) q[2];
sx q[2];
rz(-0.80255752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(-2.291912) q[1];
x q[2];
rz(2.0274721) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
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
rz(0.53044226) q[3];
sx q[3];
rz(-1.3417288) q[3];
sx q[3];
rz(-0.19977278) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];