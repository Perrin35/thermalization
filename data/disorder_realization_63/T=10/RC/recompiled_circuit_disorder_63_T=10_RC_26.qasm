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
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(-1.0648849) q[0];
rz(3.0739003) q[2];
sx q[2];
rz(-0.86172416) q[2];
sx q[2];
rz(0.46082218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5828027) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(2.5488528) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4077971) q[3];
sx q[3];
rz(-1.194209) q[3];
sx q[3];
rz(-0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7321695) q[2];
sx q[2];
rz(-2.7004341) q[2];
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
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3591374) q[0];
sx q[0];
rz(-1.29175) q[0];
sx q[0];
rz(2.7611087) q[0];
rz(-pi) q[1];
rz(-0.96946851) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(-3.0960992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9006151) q[1];
sx q[1];
rz(-2.2574189) q[1];
sx q[1];
rz(-2.2254506) q[1];
rz(-2.3362818) q[3];
sx q[3];
rz(-0.66262965) q[3];
sx q[3];
rz(-2.2413072) q[3];
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
rz(-3.0388888) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(2.8318751) q[0];
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
rz(3.0535307) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(0.17266973) q[0];
rz(-2.7823206) q[2];
sx q[2];
rz(-1.1604939) q[2];
sx q[2];
rz(1.8112195) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.09517041) q[1];
sx q[1];
rz(-1.4383525) q[1];
sx q[1];
rz(-0.73352791) q[1];
rz(-pi) q[2];
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
rz(-0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.9992874) q[3];
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
sx q[0];
rz(-pi) q[0];
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
rz(-0.79950142) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1946963) q[0];
sx q[0];
rz(-1.4799656) q[0];
sx q[0];
rz(-0.38740654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7454342) q[2];
sx q[2];
rz(-1.2067814) q[2];
sx q[2];
rz(-2.3166235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0536249) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(-1.4474523) q[1];
rz(-pi) q[2];
x q[2];
rz(0.012219592) q[3];
sx q[3];
rz(-1.0750024) q[3];
sx q[3];
rz(-0.11158768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(0.93262514) q[3];
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
rz(pi/2) q[0];
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
rz(-0.33712688) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(0.25156897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7895296) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(2.2625838) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6816363) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(-2.3098582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12794749) q[1];
sx q[1];
rz(-2.9442759) q[1];
sx q[1];
rz(-3.0017588) q[1];
x q[2];
rz(2.6336446) q[3];
sx q[3];
rz(-0.93379279) q[3];
sx q[3];
rz(0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
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
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3866766) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(2.2914698) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5315006) q[0];
sx q[0];
rz(-1.2628265) q[0];
sx q[0];
rz(2.969335) q[0];
x q[1];
rz(-2.3153147) q[2];
sx q[2];
rz(-1.7457361) q[2];
sx q[2];
rz(-0.28184055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45822083) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(-1.5797735) q[1];
x q[2];
rz(2.8443908) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(-2.8924243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
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
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958375) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(-2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(3.135625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0963421) q[0];
sx q[0];
rz(-1.2923641) q[0];
sx q[0];
rz(-0.88945008) q[0];
rz(-0.75255021) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(-0.37170751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5827427) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(1.6911669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.617745) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(2.9058128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-2.5964952) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-0.53623143) q[0];
rz(2.2082632) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(0.94820625) q[1];
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
rz(1.8188994) q[0];
rz(0.87585978) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(-1.1832331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72045418) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(-2.8645664) q[1];
rz(-0.20235297) q[3];
sx q[3];
rz(-2.7779397) q[3];
sx q[3];
rz(0.81307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-2.7622973) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76294476) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964325) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(2.3229984) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4264832) q[2];
sx q[2];
rz(-1.8679973) q[2];
sx q[2];
rz(2.4852963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0295804) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(-2.7390202) q[1];
x q[2];
rz(-2.3064763) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(2.4582668) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(-1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721601) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(-2.6955312) q[0];
x q[1];
rz(-1.2859225) q[2];
sx q[2];
rz(-1.7459918) q[2];
sx q[2];
rz(-0.81923649) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9771802) q[1];
sx q[1];
rz(-2.3986448) q[1];
sx q[1];
rz(-1.8650024) q[1];
rz(2.0274721) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(0.041681899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-3.1266406) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1508355) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(1.0583393) q[2];
sx q[2];
rz(-0.28596157) q[2];
sx q[2];
rz(2.7059976) q[2];
rz(-2.6111504) q[3];
sx q[3];
rz(-1.3417288) q[3];
sx q[3];
rz(-0.19977278) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
