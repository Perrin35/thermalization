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
rz(2.7881665) q[0];
sx q[0];
rz(8.3600144) q[0];
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.160643) q[0];
sx q[0];
rz(-1.0796709) q[0];
sx q[0];
rz(-0.26382291) q[0];
rz(-pi) q[1];
rz(1.4921161) q[2];
sx q[2];
rz(-2.4298551) q[2];
sx q[2];
rz(0.56456883) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5828027) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(-0.59273984) q[1];
x q[2];
rz(2.6081798) q[3];
sx q[3];
rz(-0.80848137) q[3];
sx q[3];
rz(-1.7318788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7321695) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(0.25201592) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(2.2759329) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(2.4332411) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3268711) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(0.65713709) q[0];
rz(-pi) q[1];
rz(2.1721241) q[2];
sx q[2];
rz(-1.6919486) q[2];
sx q[2];
rz(-0.045493424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12304141) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(-2.3398188) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49565855) q[3];
sx q[3];
rz(-2.0303876) q[3];
sx q[3];
rz(1.7840149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(0.29176816) q[2];
rz(-3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(2.5699239) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0288491) q[0];
sx q[0];
rz(-1.0380121) q[0];
sx q[0];
rz(-1.4677731) q[0];
rz(-2.2505635) q[2];
sx q[2];
rz(-2.6030412) q[2];
sx q[2];
rz(2.0856759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6211105) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(2.945167) q[1];
rz(-0.27922697) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9812575) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4806992) q[0];
sx q[0];
rz(-1.9565214) q[0];
sx q[0];
rz(1.4727403) q[0];
rz(-0.36914354) q[2];
sx q[2];
rz(-1.7338848) q[2];
sx q[2];
rz(-2.3330319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52164493) q[1];
sx q[1];
rz(-1.6879028) q[1];
sx q[1];
rz(2.820693) q[1];
rz(3.1293731) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(-0.11158768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.583741) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809526) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(0.25156897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0753206) q[0];
sx q[0];
rz(-0.89162725) q[0];
sx q[0];
rz(-0.22596304) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9676874) q[2];
sx q[2];
rz(-0.46602962) q[2];
sx q[2];
rz(0.58338651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.014616414) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(-1.5429392) q[1];
rz(-2.6336446) q[3];
sx q[3];
rz(-0.93379279) q[3];
sx q[3];
rz(-0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-1.0305369) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-2.2914698) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(0.39302557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128199) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(-1.2584932) q[0];
x q[1];
rz(2.3153147) q[2];
sx q[2];
rz(-1.3958566) q[2];
sx q[2];
rz(-0.28184055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45822083) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(1.5618192) q[1];
rz(-pi) q[2];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(-2.8924243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457552) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-3.135625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80124827) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(1.1447385) q[0];
rz(-0.47875328) q[2];
sx q[2];
rz(-2.3294449) q[2];
sx q[2];
rz(-0.85613814) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55885) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(1.6911669) q[1];
x q[2];
rz(-1.5238477) q[3];
sx q[3];
rz(-1.1991683) q[3];
sx q[3];
rz(0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(-2.7251785) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(2.1933864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535239) q[0];
sx q[0];
rz(-0.81809645) q[0];
sx q[0];
rz(1.3226932) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87585978) q[2];
sx q[2];
rz(-2.5033853) q[2];
sx q[2];
rz(-1.1832331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5612948) q[1];
sx q[1];
rz(-1.6333688) q[1];
sx q[1];
rz(2.9194174) q[1];
rz(-pi) q[2];
rz(-0.35685278) q[3];
sx q[3];
rz(-1.6423422) q[3];
sx q[3];
rz(-0.94716351) q[3];
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
rz(0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-0.35266638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34516016) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(-0.81859421) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30013957) q[2];
sx q[2];
rz(-1.4328513) q[2];
sx q[2];
rz(2.2696242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0295804) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(0.40257247) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8351164) q[3];
sx q[3];
rz(-1.1856688) q[3];
sx q[3];
rz(-1.9851828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.654489) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(-2.4328649) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1922558) q[0];
sx q[0];
rz(-2.098447) q[0];
sx q[0];
rz(-1.2883745) q[0];
x q[1];
rz(2.9591987) q[2];
sx q[2];
rz(-1.2904022) q[2];
sx q[2];
rz(2.3390351) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22589707) q[1];
sx q[1];
rz(-2.2749593) q[1];
sx q[1];
rz(-2.8812863) q[1];
x q[2];
rz(1.1141206) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(-3.0999108) q[3];
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
rz(-0.014952095) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(1.2623513) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
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
