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
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.160643) q[0];
sx q[0];
rz(-2.0619218) q[0];
sx q[0];
rz(0.26382291) q[0];
x q[1];
rz(-0.067692368) q[2];
sx q[2];
rz(-0.86172416) q[2];
sx q[2];
rz(0.46082218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2546665) q[1];
sx q[1];
rz(-2.1517422) q[1];
sx q[1];
rz(-1.7960153) q[1];
x q[2];
rz(2.4077971) q[3];
sx q[3];
rz(-1.194209) q[3];
sx q[3];
rz(0.22613444) q[3];
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
rz(-2.8895767) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-0.70835152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3268711) q[0];
sx q[0];
rz(-2.673809) q[0];
sx q[0];
rz(-0.65713709) q[0];
x q[1];
rz(1.3588261) q[2];
sx q[2];
rz(-0.6119298) q[2];
sx q[2];
rz(1.3509392) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(-0.80177387) q[1];
rz(-pi) q[2];
rz(0.49565855) q[3];
sx q[3];
rz(-1.111205) q[3];
sx q[3];
rz(1.7840149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0248802) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(-2.8498245) q[2];
rz(-3.0388888) q[3];
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
rz(pi/2) q[3];
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
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(-1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(2.5699239) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6310731) q[0];
sx q[0];
rz(-1.6594995) q[0];
sx q[0];
rz(-2.6064794) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35927202) q[2];
sx q[2];
rz(-1.9810988) q[2];
sx q[2];
rz(1.8112195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6211105) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(-0.1964257) q[1];
x q[2];
rz(-0.27922697) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(-1.8750909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40515306) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(0.23657628) q[0];
rz(-pi) q[1];
rz(1.7454342) q[2];
sx q[2];
rz(-1.2067814) q[2];
sx q[2];
rz(0.82496914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71112878) q[1];
sx q[1];
rz(-0.34090484) q[1];
sx q[1];
rz(2.7845963) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0666215) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(1.6881975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809526) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(2.8900237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35206301) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(2.2625838) q[0];
x q[1];
rz(-2.6816363) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(0.83173448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12794749) q[1];
sx q[1];
rz(-2.9442759) q[1];
sx q[1];
rz(-0.13983388) q[1];
rz(2.1523347) q[3];
sx q[3];
rz(-0.79205081) q[3];
sx q[3];
rz(-2.0056412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.183737) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(2.1110558) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(-0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(-1.1823581) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0133936) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(-1.8830995) q[0];
x q[1];
rz(-1.3156462) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(1.6659425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0361573) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(-0.91962199) q[1];
x q[2];
rz(-0.7418467) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(1.1175931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3605911) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(-2.1684516) q[2];
rz(2.1940103) q[3];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(-2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(3.135625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3403444) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(-1.1447385) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47875328) q[2];
sx q[2];
rz(-0.81214777) q[2];
sx q[2];
rz(-0.85613814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55885) q[1];
sx q[1];
rz(-0.55816459) q[1];
sx q[1];
rz(-1.6911669) q[1];
x q[2];
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
rz(pi/2) q[1];
rz(-0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0311325) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(-2.2082632) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
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
rz(-2.6981632) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(-2.7623917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72045418) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(-2.8645664) q[1];
rz(2.9392397) q[3];
sx q[3];
rz(-0.36365299) q[3];
sx q[3];
rz(-0.81307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964325) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(2.3229984) q[0];
x q[1];
rz(1.4264832) q[2];
sx q[2];
rz(-1.8679973) q[2];
sx q[2];
rz(2.4852963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7949617) q[1];
sx q[1];
rz(-1.3449841) q[1];
sx q[1];
rz(-0.56983106) q[1];
x q[2];
rz(1.0274067) q[3];
sx q[3];
rz(-2.3282442) q[3];
sx q[3];
rz(-0.020997626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(0.68332589) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4721601) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(2.6955312) q[0];
x q[1];
rz(-1.8556701) q[2];
sx q[2];
rz(-1.3956009) q[2];
sx q[2];
rz(-0.81923649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22589707) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(-0.26030635) q[1];
rz(-0.98827045) q[3];
sx q[3];
rz(-1.3068145) q[3];
sx q[3];
rz(1.1519983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(0.14317748) q[2];
sx q[2];
rz(-1.3224052) q[2];
sx q[2];
rz(0.094766141) q[2];
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
