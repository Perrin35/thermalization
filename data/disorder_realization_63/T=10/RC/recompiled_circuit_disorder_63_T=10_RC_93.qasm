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
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(1.0648849) q[0];
x q[1];
rz(-1.6494765) q[2];
sx q[2];
rz(-2.4298551) q[2];
sx q[2];
rz(0.56456883) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55878996) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(2.5488528) q[1];
rz(-pi) q[2];
rz(-1.0814813) q[3];
sx q[3];
rz(-0.89852528) q[3];
sx q[3];
rz(2.1171452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(0.25201592) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
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
rz(-0.73671571) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81472155) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(-0.65713709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1721241) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(3.0960992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(-0.80177387) q[1];
rz(-pi) q[2];
rz(-0.80531081) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(-2.2413072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-0.29176816) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(0.30971757) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(0.57166878) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5105195) q[0];
sx q[0];
rz(-1.4820931) q[0];
sx q[0];
rz(-0.53511329) q[0];
x q[1];
rz(0.89102913) q[2];
sx q[2];
rz(-0.53855145) q[2];
sx q[2];
rz(1.0559168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3571615) q[1];
sx q[1];
rz(-0.845134) q[1];
sx q[1];
rz(-1.3933338) q[1];
rz(-pi) q[2];
rz(1.709135) q[3];
sx q[3];
rz(-1.1225015) q[3];
sx q[3];
rz(1.5773147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16033515) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-0.79950142) q[0];
rz(-3.0917621) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40515306) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(2.9050164) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3961584) q[2];
sx q[2];
rz(-1.2067814) q[2];
sx q[2];
rz(0.82496914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6199477) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(0.32089969) q[1];
rz(-pi) q[2];
rz(-0.012219592) q[3];
sx q[3];
rz(-1.0750024) q[3];
sx q[3];
rz(-3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.583741) q[2];
rz(1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-0.93262514) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(2.8900237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0662721) q[0];
sx q[0];
rz(-2.2499654) q[0];
sx q[0];
rz(0.22596304) q[0];
x q[1];
rz(2.6816363) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(2.3098582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12794749) q[1];
sx q[1];
rz(-2.9442759) q[1];
sx q[1];
rz(-3.0017588) q[1];
rz(0.86815636) q[3];
sx q[3];
rz(-1.9725102) q[3];
sx q[3];
rz(-0.86740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(-1.1823581) q[1];
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
rz(0.61009208) q[0];
sx q[0];
rz(-1.8787662) q[0];
sx q[0];
rz(-0.17225762) q[0];
rz(1.3156462) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(1.4756502) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46950754) q[1];
sx q[1];
rz(-2.2219594) q[1];
sx q[1];
rz(3.1347515) q[1];
x q[2];
rz(2.399746) q[3];
sx q[3];
rz(-2.7478455) q[3];
sx q[3];
rz(2.0239995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(-2.1684516) q[2];
rz(-2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958375) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(-1.6181035) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25516605) q[0];
sx q[0];
rz(-2.2212941) q[0];
sx q[0];
rz(2.7889473) q[0];
x q[1];
rz(-2.6628394) q[2];
sx q[2];
rz(-2.3294449) q[2];
sx q[2];
rz(0.85613814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70049268) q[1];
sx q[1];
rz(-2.1244493) q[1];
sx q[1];
rz(-0.074837491) q[1];
rz(0.11984101) q[3];
sx q[3];
rz(-2.7671475) q[3];
sx q[3];
rz(0.10711174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(-2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(-0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-2.1933864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95420102) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(-0.76822922) q[0];
rz(-0.87585978) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(1.1832331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72045418) q[1];
sx q[1];
rz(-2.9109143) q[1];
sx q[1];
rz(-2.8645664) q[1];
rz(0.35685278) q[3];
sx q[3];
rz(-1.6423422) q[3];
sx q[3];
rz(0.94716351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3809985) q[1];
sx q[1];
rz(0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0371373) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(-3.0260968) q[0];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-1.7949617) q[1];
sx q[1];
rz(-1.7966086) q[1];
sx q[1];
rz(-0.56983106) q[1];
rz(-pi) q[2];
rz(1.0274067) q[3];
sx q[3];
rz(-0.81334844) q[3];
sx q[3];
rz(-3.120595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(0.68332589) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6694326) q[0];
sx q[0];
rz(-2.5494808) q[0];
sx q[0];
rz(0.44606146) q[0];
rz(-pi) q[1];
rz(-1.2859225) q[2];
sx q[2];
rz(-1.3956009) q[2];
sx q[2];
rz(0.81923649) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9771802) q[1];
sx q[1];
rz(-2.3986448) q[1];
sx q[1];
rz(1.2765902) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1533222) q[3];
sx q[3];
rz(-1.8347782) q[3];
sx q[3];
rz(1.9895944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(-0.014952095) q[2];
rz(-0.22710083) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.0583393) q[2];
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
