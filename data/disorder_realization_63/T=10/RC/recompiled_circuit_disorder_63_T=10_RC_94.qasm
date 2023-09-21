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
rz(2.605521) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9809496) q[0];
sx q[0];
rz(-2.0619218) q[0];
sx q[0];
rz(0.26382291) q[0];
rz(-pi) q[1];
rz(2.2810018) q[2];
sx q[2];
rz(-1.5194367) q[2];
sx q[2];
rz(2.0757338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8592035) q[1];
sx q[1];
rz(-2.5232362) q[1];
sx q[1];
rz(-0.32787852) q[1];
rz(-pi) q[2];
rz(-1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(1.0244474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7321695) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.3902364) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3591374) q[0];
sx q[0];
rz(-1.29175) q[0];
sx q[0];
rz(2.7611087) q[0];
x q[1];
rz(0.14658908) q[2];
sx q[2];
rz(-2.1671038) q[2];
sx q[2];
rz(1.5335611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(0.80177387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3362818) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(-2.2413072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-0.57166878) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088061995) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(0.17266973) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.13582) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(-0.091718397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5204822) q[1];
sx q[1];
rz(-0.74319327) q[1];
sx q[1];
rz(2.945167) q[1];
x q[2];
rz(2.8623657) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(-1.8750909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-2.3420912) q[0];
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
rz(2.7364396) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(-2.9050164) q[0];
x q[1];
rz(-2.7724491) q[2];
sx q[2];
rz(-1.4077079) q[2];
sx q[2];
rz(-2.3330319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6199477) q[1];
sx q[1];
rz(-1.6879028) q[1];
sx q[1];
rz(0.32089969) q[1];
rz(-pi) q[2];
rz(-2.0666215) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(1.6881975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0522456) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(-1.583741) q[2];
rz(1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(2.2089675) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(1.1622693) q[1];
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
rz(-1.7154327) q[0];
sx q[0];
rz(-0.71007219) q[0];
sx q[0];
rz(-1.3000814) q[0];
x q[1];
rz(0.17390522) q[2];
sx q[2];
rz(-0.46602962) q[2];
sx q[2];
rz(-0.58338651) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0136452) q[1];
sx q[1];
rz(-2.9442759) q[1];
sx q[1];
rz(0.13983388) q[1];
x q[2];
rz(2.1523347) q[3];
sx q[3];
rz(-2.3495418) q[3];
sx q[3];
rz(2.0056412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(2.1110558) q[2];
rz(-1.7306227) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(2.2914698) q[0];
rz(-1.9592346) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(-2.7485671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5315006) q[0];
sx q[0];
rz(-1.8787662) q[0];
sx q[0];
rz(-2.969335) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82627798) q[2];
sx q[2];
rz(-1.7457361) q[2];
sx q[2];
rz(0.28184055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45822083) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(1.5797735) q[1];
rz(2.399746) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(-2.0239995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3605911) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(2.1684516) q[2];
rz(2.1940103) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(3.135625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(1.9968541) q[0];
x q[1];
rz(2.3890424) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(2.7698851) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5827427) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(-1.6911669) q[1];
x q[2];
rz(-2.7695914) q[3];
sx q[3];
rz(-1.614538) q[3];
sx q[3];
rz(-1.352076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5696047) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(-2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(0.93332943) q[1];
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
rz(-1.750995) q[0];
sx q[0];
rz(-2.3733634) q[0];
x q[1];
rz(-0.44342946) q[2];
sx q[2];
rz(-1.0955053) q[2];
sx q[2];
rz(-2.7623917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72045418) q[1];
sx q[1];
rz(-2.9109143) q[1];
sx q[1];
rz(2.8645664) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4944581) q[3];
sx q[3];
rz(-1.9266955) q[3];
sx q[3];
rz(0.59698856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(-2.9013157) q[2];
rz(0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371373) q[0];
sx q[0];
rz(-1.6933352) q[0];
sx q[0];
rz(-3.0260968) q[0];
x q[1];
rz(1.4264832) q[2];
sx q[2];
rz(-1.2735954) q[2];
sx q[2];
rz(0.65629634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7949617) q[1];
sx q[1];
rz(-1.7966086) q[1];
sx q[1];
rz(0.56983106) q[1];
rz(-pi) q[2];
rz(-1.0274067) q[3];
sx q[3];
rz(-0.81334844) q[3];
sx q[3];
rz(3.120595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(0.68332589) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(-1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(2.4328649) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6652128) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(-0.54540821) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2859225) q[2];
sx q[2];
rz(-1.7459918) q[2];
sx q[2];
rz(-2.3223562) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5156538) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(2.291912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1141206) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(1.5785718) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(2.0832534) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
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
