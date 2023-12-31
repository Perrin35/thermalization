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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85815) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(2.0767077) q[0];
x q[1];
rz(-0.067692368) q[2];
sx q[2];
rz(-2.2798685) q[2];
sx q[2];
rz(-0.46082218) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5828027) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(-0.59273984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6081798) q[3];
sx q[3];
rz(-0.80848137) q[3];
sx q[3];
rz(1.7318788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(2.4332411) q[1];
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
rz(-pi) q[1];
x q[1];
rz(2.1721241) q[2];
sx q[2];
rz(-1.6919486) q[2];
sx q[2];
rz(3.0960992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2409776) q[1];
sx q[1];
rz(-0.88417378) q[1];
sx q[1];
rz(2.2254506) q[1];
x q[2];
rz(1.0582916) q[3];
sx q[3];
rz(-1.1304434) q[3];
sx q[3];
rz(3.1194221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(-1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-2.5699239) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5105195) q[0];
sx q[0];
rz(-1.4820931) q[0];
sx q[0];
rz(2.6064794) q[0];
x q[1];
rz(0.89102913) q[2];
sx q[2];
rz(-0.53855145) q[2];
sx q[2];
rz(-2.0856759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7844312) q[1];
sx q[1];
rz(-0.845134) q[1];
sx q[1];
rz(1.7482589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6895371) q[3];
sx q[3];
rz(-1.6953903) q[3];
sx q[3];
rz(0.053754036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.16033515) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(-0.30250868) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.4335094) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7364396) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(2.9050164) q[0];
rz(0.36914354) q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0879678) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(-1.4474523) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5933851) q[3];
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
rz(-1.9767438) q[2];
sx q[2];
rz(1.5578516) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(2.8900237) q[1];
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
rz(0.45995633) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(0.83173448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12794749) q[1];
sx q[1];
rz(-0.19731678) q[1];
sx q[1];
rz(3.0017588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98925791) q[3];
sx q[3];
rz(-0.79205081) q[3];
sx q[3];
rz(-1.1359515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(-1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(2.7485671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128199) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(1.2584932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9057301) q[2];
sx q[2];
rz(-0.84025192) q[2];
sx q[2];
rz(-2.0116531) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6833718) q[1];
sx q[1];
rz(-0.65119377) q[1];
sx q[1];
rz(-1.5797735) q[1];
rz(-pi) q[2];
rz(2.8443908) q[3];
sx q[3];
rz(-1.8330049) q[3];
sx q[3];
rz(-0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7810016) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(-0.97314107) q[2];
rz(-2.1940103) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-0.0059676776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0963421) q[0];
sx q[0];
rz(-1.2923641) q[0];
sx q[0];
rz(-2.2521426) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1184095) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(1.6391022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4411) q[1];
sx q[1];
rz(-1.0171434) q[1];
sx q[1];
rz(-3.0667552) q[1];
rz(-pi) q[2];
rz(-0.11984101) q[3];
sx q[3];
rz(-2.7671475) q[3];
sx q[3];
rz(-0.10711174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(2.1933864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535239) q[0];
sx q[0];
rz(-2.3234962) q[0];
sx q[0];
rz(-1.3226932) q[0];
rz(2.2657329) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(-1.9583595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1369689) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(1.5066513) q[1];
rz(-pi) q[2];
rz(1.4944581) q[3];
sx q[3];
rz(-1.9266955) q[3];
sx q[3];
rz(2.5446041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371373) q[0];
sx q[0];
rz(-1.6933352) q[0];
sx q[0];
rz(-3.0260968) q[0];
rz(-pi) q[1];
rz(-0.30013957) q[2];
sx q[2];
rz(-1.4328513) q[2];
sx q[2];
rz(0.87196841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36665146) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(-1.8371546) q[1];
rz(-pi) q[2];
rz(-2.6412233) q[3];
sx q[3];
rz(-2.2420068) q[3];
sx q[3];
rz(2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76715604) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(-2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-0.93651071) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721601) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(0.44606146) q[0];
rz(-2.1328743) q[2];
sx q[2];
rz(-2.808411) q[2];
sx q[2];
rz(2.9269232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22589707) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(-2.8812863) q[1];
x q[2];
rz(2.0274721) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(-3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(1.5785718) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(2.9984152) q[2];
sx q[2];
rz(-1.8191874) q[2];
sx q[2];
rz(-3.0468265) q[2];
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
