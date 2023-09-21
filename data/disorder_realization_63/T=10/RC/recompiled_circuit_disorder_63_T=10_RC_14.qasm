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
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(-2.0767077) q[0];
rz(-pi) q[1];
rz(1.6494765) q[2];
sx q[2];
rz(-2.4298551) q[2];
sx q[2];
rz(2.5770238) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55878996) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(0.59273984) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0814813) q[3];
sx q[3];
rz(-0.89852528) q[3];
sx q[3];
rz(-2.1171452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7321695) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(-0.25201592) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(2.2759329) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-0.70835152) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81472155) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(-2.4844556) q[0];
rz(-pi) q[1];
rz(1.7827665) q[2];
sx q[2];
rz(-0.6119298) q[2];
sx q[2];
rz(1.7906534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0201455) q[1];
sx q[1];
rz(-0.91031204) q[1];
sx q[1];
rz(2.5026908) q[1];
rz(-1.0582916) q[3];
sx q[3];
rz(-2.0111492) q[3];
sx q[3];
rz(-0.022170598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
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
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-2.5699239) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0288491) q[0];
sx q[0];
rz(-1.0380121) q[0];
sx q[0];
rz(1.4677731) q[0];
rz(-pi) q[1];
rz(-1.13582) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(3.0498743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5204822) q[1];
sx q[1];
rz(-0.74319327) q[1];
sx q[1];
rz(-0.1964257) q[1];
rz(-pi) q[2];
rz(1.4324576) q[3];
sx q[3];
rz(-2.0190911) q[3];
sx q[3];
rz(-1.5642779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16033515) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-0.70880115) q[2];
rz(0.30250868) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-0.049830534) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-0.18049151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9468964) q[0];
sx q[0];
rz(-1.4799656) q[0];
sx q[0];
rz(-2.7541861) q[0];
x q[1];
rz(1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(2.3166235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0536249) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(-1.6941403) q[1];
x q[2];
rz(3.1293731) q[3];
sx q[3];
rz(-1.0750024) q[3];
sx q[3];
rz(-3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(-1.583741) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809526) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(-0.98651648) q[0];
rz(1.1622693) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-0.25156897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0753206) q[0];
sx q[0];
rz(-2.2499654) q[0];
sx q[0];
rz(2.9156296) q[0];
rz(-1.4839843) q[2];
sx q[2];
rz(-1.112339) q[2];
sx q[2];
rz(0.77755962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.014616414) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(-1.5429392) q[1];
rz(0.86815636) q[3];
sx q[3];
rz(-1.9725102) q[3];
sx q[3];
rz(2.274184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(-2.7485671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052554) q[0];
sx q[0];
rz(-2.7900643) q[0];
sx q[0];
rz(-1.0765443) q[0];
rz(-2.9057301) q[2];
sx q[2];
rz(-0.84025192) q[2];
sx q[2];
rz(2.0116531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6720851) q[1];
sx q[1];
rz(-2.2219594) q[1];
sx q[1];
rz(0.0068411946) q[1];
x q[2];
rz(-2.399746) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(-1.1175931) q[3];
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
rz(0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(3.135625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80124827) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(-1.1447385) q[0];
rz(-pi) q[1];
rz(1.1184095) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(1.5024904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5827427) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(-1.4504257) q[1];
rz(-pi) q[2];
rz(0.11984101) q[3];
sx q[3];
rz(-0.37444515) q[3];
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
rz(-2.5964952) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(2.1933864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084003) q[0];
sx q[0];
rz(-0.78484479) q[0];
sx q[0];
rz(-0.25640304) q[0];
rz(-2.6981632) q[2];
sx q[2];
rz(-1.0955053) q[2];
sx q[2];
rz(2.7623917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.72045418) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(2.8645664) q[1];
rz(-pi) q[2];
rz(-0.35685278) q[3];
sx q[3];
rz(-1.6423422) q[3];
sx q[3];
rz(-0.94716351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(0.24027696) q[2];
rz(-2.7622973) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(0.25094029) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(0.35266638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0371373) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(0.11549581) q[0];
x q[1];
rz(1.7151095) q[2];
sx q[2];
rz(-1.8679973) q[2];
sx q[2];
rz(0.65629634) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7749412) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(-1.8371546) q[1];
rz(-pi) q[2];
rz(2.6412233) q[3];
sx q[3];
rz(-2.2420068) q[3];
sx q[3];
rz(-2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(-1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.93651071) q[1];
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
rz(-1.0087183) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(2.9269232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5156538) q[1];
sx q[1];
rz(-1.3733555) q[1];
sx q[1];
rz(2.291912) q[1];
x q[2];
rz(-2.8285633) q[3];
sx q[3];
rz(-1.0109245) q[3];
sx q[3];
rz(0.58902878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(0.014952095) q[2];
rz(2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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