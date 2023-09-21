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
rz(-0.53607166) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4608085) q[0];
sx q[0];
rz(-0.55235282) q[0];
sx q[0];
rz(-2.0244563) q[0];
rz(-pi) q[1];
rz(-3.0739003) q[2];
sx q[2];
rz(-0.86172416) q[2];
sx q[2];
rz(-0.46082218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55878996) q[1];
sx q[1];
rz(-1.7585808) q[1];
sx q[1];
rz(0.59273984) q[1];
rz(2.4077971) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(2.9154582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7321695) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(2.1146963) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(0.70835152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10193292) q[0];
sx q[0];
rz(-1.205737) q[0];
sx q[0];
rz(-1.2714766) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9950036) q[2];
sx q[2];
rz(-0.97448889) q[2];
sx q[2];
rz(-1.6080315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2409776) q[1];
sx q[1];
rz(-0.88417378) q[1];
sx q[1];
rz(-0.91614206) q[1];
rz(-pi) q[2];
rz(2.6459341) q[3];
sx q[3];
rz(-2.0303876) q[3];
sx q[3];
rz(1.7840149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-2.8498245) q[2];
rz(-3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(0.57166878) q[1];
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
rz(-0.89102913) q[2];
sx q[2];
rz(-2.6030412) q[2];
sx q[2];
rz(1.0559168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0464222) q[1];
sx q[1];
rz(-1.7032402) q[1];
sx q[1];
rz(2.4080647) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27922697) q[3];
sx q[3];
rz(-2.673827) q[3];
sx q[3];
rz(1.8750909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.1423053) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-2.9611011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40515306) q[0];
sx q[0];
rz(-2.7442051) q[0];
sx q[0];
rz(-0.23657628) q[0];
rz(1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(2.3166235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0536249) q[1];
sx q[1];
rz(-1.2521724) q[1];
sx q[1];
rz(1.6941403) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1293731) q[3];
sx q[3];
rz(-1.0750024) q[3];
sx q[3];
rz(-0.11158768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.583741) q[2];
rz(1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(0.25156897) q[1];
rz(-pi) q[2];
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
rz(-0.45995633) q[2];
sx q[2];
rz(-1.6486247) q[2];
sx q[2];
rz(0.83173448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.561589) q[1];
sx q[1];
rz(-1.5434693) q[1];
sx q[1];
rz(0.1954397) q[1];
x q[2];
rz(-2.6336446) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(-2.7579443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.183737) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(2.7485671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052554) q[0];
sx q[0];
rz(-2.7900643) q[0];
sx q[0];
rz(-1.0765443) q[0];
rz(-pi) q[1];
rz(-2.9057301) q[2];
sx q[2];
rz(-2.3013407) q[2];
sx q[2];
rz(1.1299396) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6833718) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(1.5618192) q[1];
x q[2];
rz(-1.844448) q[3];
sx q[3];
rz(-1.2840464) q[3];
sx q[3];
rz(1.2424038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3605911) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(2.1684516) q[2];
rz(-2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.0059676776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(1.9968541) q[0];
x q[1];
rz(-2.3890424) q[2];
sx q[2];
rz(-1.9117022) q[2];
sx q[2];
rz(-0.37170751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70049268) q[1];
sx q[1];
rz(-1.0171434) q[1];
sx q[1];
rz(3.0667552) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5238477) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5696047) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.0311325) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(-2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(0.94820625) q[1];
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
rz(-1.3226932) q[0];
rz(-pi) q[1];
rz(-0.87585978) q[2];
sx q[2];
rz(-2.5033853) q[2];
sx q[2];
rz(-1.1832331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58029786) q[1];
sx q[1];
rz(-1.5082238) q[1];
sx q[1];
rz(0.22217521) q[1];
rz(-0.35685278) q[3];
sx q[3];
rz(-1.4992504) q[3];
sx q[3];
rz(-2.1944291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-0.25094029) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-2.7889263) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371373) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(-0.11549581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7025931) q[2];
sx q[2];
rz(-0.3294496) q[2];
sx q[2];
rz(-2.0246558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7749412) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(-1.304438) q[1];
x q[2];
rz(2.3064763) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94933687) q[0];
sx q[0];
rz(-2.098447) q[0];
sx q[0];
rz(1.2883745) q[0];
rz(-pi) q[1];
rz(-1.0087183) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(2.9269232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22589707) q[1];
sx q[1];
rz(-0.86663336) q[1];
sx q[1];
rz(-2.8812863) q[1];
x q[2];
rz(2.0274721) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(-0.041681899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(0.014952095) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(2.7097377) q[3];
sx q[3];
rz(-2.5681744) q[3];
sx q[3];
rz(-1.4011866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];