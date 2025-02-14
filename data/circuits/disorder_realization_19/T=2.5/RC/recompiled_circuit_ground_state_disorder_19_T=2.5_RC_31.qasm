OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1524393) q[0];
sx q[0];
rz(-1.5751155) q[0];
sx q[0];
rz(-2.0164665) q[0];
rz(1.0220802) q[1];
sx q[1];
rz(-0.66749579) q[1];
sx q[1];
rz(0.8134841) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9620243) q[0];
sx q[0];
rz(-1.429707) q[0];
sx q[0];
rz(-1.2015094) q[0];
x q[1];
rz(1.2101233) q[2];
sx q[2];
rz(-1.9934335) q[2];
sx q[2];
rz(-0.3061184) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5578545) q[1];
sx q[1];
rz(-2.0101133) q[1];
sx q[1];
rz(1.735461) q[1];
x q[2];
rz(2.0735246) q[3];
sx q[3];
rz(-0.18618551) q[3];
sx q[3];
rz(-1.8753827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6525314) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(0.92864621) q[2];
rz(-1.5422025) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20392513) q[0];
sx q[0];
rz(-1.3805905) q[0];
sx q[0];
rz(-3.0145338) q[0];
rz(-0.98310414) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(-2.3703221) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5707023) q[0];
sx q[0];
rz(-1.2178019) q[0];
sx q[0];
rz(2.3285026) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0808582) q[2];
sx q[2];
rz(-1.2337451) q[2];
sx q[2];
rz(3.0063546) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3202618) q[1];
sx q[1];
rz(-1.2148569) q[1];
sx q[1];
rz(1.9411646) q[1];
rz(-0.57621376) q[3];
sx q[3];
rz(-2.0348843) q[3];
sx q[3];
rz(2.4903542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9942921) q[2];
sx q[2];
rz(-1.9376126) q[2];
sx q[2];
rz(3.1331983) q[2];
rz(-0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(-2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10107772) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(0.44152942) q[0];
rz(2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-3.006014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39865935) q[0];
sx q[0];
rz(-0.017983111) q[0];
sx q[0];
rz(-2.6989486) q[0];
rz(-pi) q[1];
rz(-0.71696059) q[2];
sx q[2];
rz(-1.2768942) q[2];
sx q[2];
rz(-0.64168054) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6171268) q[1];
sx q[1];
rz(-1.9248665) q[1];
sx q[1];
rz(-1.0431784) q[1];
x q[2];
rz(1.7884004) q[3];
sx q[3];
rz(-0.74624589) q[3];
sx q[3];
rz(-1.7640132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2077937) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(-0.03820339) q[2];
rz(2.6162052) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(-0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4811089) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(-2.5307181) q[0];
rz(-1.3350217) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(2.9023721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611899) q[0];
sx q[0];
rz(-2.1227269) q[0];
sx q[0];
rz(0.30776382) q[0];
rz(-pi) q[1];
rz(0.47232136) q[2];
sx q[2];
rz(-1.9237674) q[2];
sx q[2];
rz(-1.2833088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94915056) q[1];
sx q[1];
rz(-0.56862967) q[1];
sx q[1];
rz(2.5123764) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(0.1399006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4227582) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(0.0017496721) q[2];
rz(-2.8478029) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(0.26688117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2413498) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(-0.84306651) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(1.5379803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3030515) q[0];
sx q[0];
rz(-1.3796796) q[0];
sx q[0];
rz(2.3421846) q[0];
rz(-2.9672253) q[2];
sx q[2];
rz(-1.6233168) q[2];
sx q[2];
rz(-0.32584056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2258721) q[1];
sx q[1];
rz(-0.32494007) q[1];
sx q[1];
rz(-2.48808) q[1];
rz(-pi) q[2];
rz(2.9318453) q[3];
sx q[3];
rz(-1.4919859) q[3];
sx q[3];
rz(2.1926853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43188492) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(-1.0464) q[2];
rz(-0.31878582) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980314) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(1.3290149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15547046) q[0];
sx q[0];
rz(-2.8948445) q[0];
sx q[0];
rz(-2.0339436) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29000303) q[2];
sx q[2];
rz(-1.3652186) q[2];
sx q[2];
rz(0.21943352) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15550286) q[1];
sx q[1];
rz(-0.8900607) q[1];
sx q[1];
rz(-0.40386856) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51402199) q[3];
sx q[3];
rz(-1.2482621) q[3];
sx q[3];
rz(-3.1399825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7602188) q[2];
sx q[2];
rz(-1.2083283) q[2];
sx q[2];
rz(1.5488497) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(-2.2729592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(0.94183952) q[0];
rz(-2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036445905) q[0];
sx q[0];
rz(-2.9147286) q[0];
sx q[0];
rz(-1.731621) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0925557) q[2];
sx q[2];
rz(-1.5360439) q[2];
sx q[2];
rz(0.032973789) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3698764) q[1];
sx q[1];
rz(-1.2315005) q[1];
sx q[1];
rz(2.3849065) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0683769) q[3];
sx q[3];
rz(-1.5514152) q[3];
sx q[3];
rz(1.3598702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3840702) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(2.7080217) q[2];
rz(1.5571669) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(0.43237329) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(-1.3442511) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(1.7787836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8936854) q[0];
sx q[0];
rz(-3.1046668) q[0];
sx q[0];
rz(-0.93494995) q[0];
rz(-pi) q[1];
rz(1.3819456) q[2];
sx q[2];
rz(-2.3879693) q[2];
sx q[2];
rz(-3.0830887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40920112) q[1];
sx q[1];
rz(-0.26136569) q[1];
sx q[1];
rz(-2.638326) q[1];
rz(-1.9523432) q[3];
sx q[3];
rz(-2.4182712) q[3];
sx q[3];
rz(1.2398694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5726996) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(-2.1577238) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(-2.4510395) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47138658) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(0.27467003) q[0];
rz(1.7999016) q[1];
sx q[1];
rz(-2.5119753) q[1];
sx q[1];
rz(-2.0955657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90993308) q[0];
sx q[0];
rz(-2.446736) q[0];
sx q[0];
rz(-2.6207663) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9555904) q[2];
sx q[2];
rz(-0.21285393) q[2];
sx q[2];
rz(-1.4215046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6089692) q[1];
sx q[1];
rz(-1.5057505) q[1];
sx q[1];
rz(0.58397997) q[1];
rz(-pi) q[2];
rz(0.041954354) q[3];
sx q[3];
rz(-2.0734412) q[3];
sx q[3];
rz(1.6142538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2785953) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(0.39946237) q[2];
rz(-0.56600371) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(-2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28867662) q[0];
sx q[0];
rz(-1.4194019) q[0];
sx q[0];
rz(0.5823108) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(0.80642548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18412735) q[0];
sx q[0];
rz(-0.13617198) q[0];
sx q[0];
rz(2.2720112) q[0];
rz(2.8060032) q[2];
sx q[2];
rz(-1.3177455) q[2];
sx q[2];
rz(-2.0584681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2576221) q[1];
sx q[1];
rz(-0.7695573) q[1];
sx q[1];
rz(-0.29410024) q[1];
x q[2];
rz(2.8424758) q[3];
sx q[3];
rz(-2.7795305) q[3];
sx q[3];
rz(-0.55940926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9091984) q[2];
sx q[2];
rz(-2.4269203) q[2];
sx q[2];
rz(2.3363028) q[2];
rz(-0.14690873) q[3];
sx q[3];
rz(-0.78607905) q[3];
sx q[3];
rz(-0.73827353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88937) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(1.5994785) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-3.0270544) q[2];
sx q[2];
rz(-2.2907612) q[2];
sx q[2];
rz(-1.9775122) q[2];
rz(-0.64503786) q[3];
sx q[3];
rz(-1.8492263) q[3];
sx q[3];
rz(3.1254569) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
