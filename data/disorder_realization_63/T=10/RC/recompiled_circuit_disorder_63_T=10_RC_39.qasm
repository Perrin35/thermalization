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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6807841) q[0];
sx q[0];
rz(-2.5892398) q[0];
sx q[0];
rz(2.0244563) q[0];
rz(-pi) q[1];
rz(-0.86059086) q[2];
sx q[2];
rz(-1.6221559) q[2];
sx q[2];
rz(1.0658588) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8592035) q[1];
sx q[1];
rz(-2.5232362) q[1];
sx q[1];
rz(-2.8137141) q[1];
rz(-pi) q[2];
rz(1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(-1.0244474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(0.25201592) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.7513562) q[0];
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
rz(-3.0396597) q[0];
sx q[0];
rz(-1.9358557) q[0];
sx q[0];
rz(-1.2714766) q[0];
x q[1];
rz(-2.9950036) q[2];
sx q[2];
rz(-2.1671038) q[2];
sx q[2];
rz(-1.6080315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0201455) q[1];
sx q[1];
rz(-2.2312806) q[1];
sx q[1];
rz(-2.5026908) q[1];
rz(-0.80531081) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-0.57166878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0535307) q[0];
sx q[0];
rz(-0.54170875) q[0];
sx q[0];
rz(-0.17266973) q[0];
rz(-pi) q[1];
rz(1.13582) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(0.091718397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7844312) q[1];
sx q[1];
rz(-2.2964587) q[1];
sx q[1];
rz(-1.7482589) q[1];
rz(-1.709135) q[3];
sx q[3];
rz(-1.1225015) q[3];
sx q[3];
rz(-1.5773147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(0.18049151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4806992) q[0];
sx q[0];
rz(-1.1850712) q[0];
sx q[0];
rz(-1.6688523) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3961584) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(0.82496914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52164493) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(2.820693) q[1];
x q[2];
rz(-0.012219592) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.583741) q[2];
rz(1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3809526) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(1.1622693) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(2.8900237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4261599) q[0];
sx q[0];
rz(-0.71007219) q[0];
sx q[0];
rz(-1.3000814) q[0];
rz(0.45995633) q[2];
sx q[2];
rz(-1.6486247) q[2];
sx q[2];
rz(-0.83173448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12794749) q[1];
sx q[1];
rz(-0.19731678) q[1];
sx q[1];
rz(0.13983388) q[1];
x q[2];
rz(-0.86815636) q[3];
sx q[3];
rz(-1.1690825) q[3];
sx q[3];
rz(-0.86740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(2.1110558) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(-2.7485671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0133936) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(-1.8830995) q[0];
rz(-pi) q[1];
rz(0.23586258) q[2];
sx q[2];
rz(-0.84025192) q[2];
sx q[2];
rz(2.0116531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6720851) q[1];
sx q[1];
rz(-0.91963327) q[1];
sx q[1];
rz(-0.0068411946) q[1];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.8330049) q[3];
sx q[3];
rz(-0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(-0.97314107) q[2];
rz(0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.1943641) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457552) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-0.0059676776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0963421) q[0];
sx q[0];
rz(-1.8492286) q[0];
sx q[0];
rz(-0.88945008) q[0];
rz(-pi) q[1];
rz(-2.0231831) q[2];
sx q[2];
rz(-0.87086073) q[2];
sx q[2];
rz(-1.5024904) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90970627) q[1];
sx q[1];
rz(-1.6344374) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(2.9515284) q[2];
rz(-0.41641411) q[3];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(-2.2082632) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78806879) q[0];
sx q[0];
rz(-0.81809645) q[0];
sx q[0];
rz(1.8188994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0886705) q[2];
sx q[2];
rz(-1.1793943) q[2];
sx q[2];
rz(-2.1640167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1369689) q[1];
sx q[1];
rz(-1.7925295) q[1];
sx q[1];
rz(-1.5066513) q[1];
x q[2];
rz(2.9392397) q[3];
sx q[3];
rz(-2.7779397) q[3];
sx q[3];
rz(0.81307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(-2.9013157) q[2];
rz(-2.7622973) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.7933638) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(0.25094029) q[0];
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
rz(-0.34516016) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(2.3229984) q[0];
rz(0.30013957) q[2];
sx q[2];
rz(-1.7087414) q[2];
sx q[2];
rz(-2.2696242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0295804) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(-2.7390202) q[1];
rz(-pi) q[2];
rz(2.3064763) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(-1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(0.70872778) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4721601) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(-0.44606146) q[0];
rz(-pi) q[1];
rz(-2.1328743) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(0.21466941) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6259389) q[1];
sx q[1];
rz(-1.3733555) q[1];
sx q[1];
rz(-2.291912) q[1];
x q[2];
rz(-1.1141206) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(0.041681899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99075714) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(1.5785718) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(1.0583393) q[2];
sx q[2];
rz(-0.28596157) q[2];
sx q[2];
rz(2.7059976) q[2];
rz(-1.8347918) q[3];
sx q[3];
rz(-1.055607) q[3];
sx q[3];
rz(-1.9029688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
