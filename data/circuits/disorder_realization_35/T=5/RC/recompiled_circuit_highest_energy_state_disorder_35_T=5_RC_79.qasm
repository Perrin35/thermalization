OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.063989446) q[0];
sx q[0];
rz(4.0539157) q[0];
sx q[0];
rz(11.231448) q[0];
rz(1.9995243) q[1];
sx q[1];
rz(-2.6231782) q[1];
sx q[1];
rz(-1.6496744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27150422) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(-1.2369878) q[0];
rz(-0.69500287) q[2];
sx q[2];
rz(-2.0796513) q[2];
sx q[2];
rz(-2.5093504) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8644997) q[1];
sx q[1];
rz(-2.1756209) q[1];
sx q[1];
rz(2.7983309) q[1];
x q[2];
rz(1.434695) q[3];
sx q[3];
rz(-0.24472642) q[3];
sx q[3];
rz(-0.32258209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3492744) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(-0.054923687) q[2];
rz(1.6867636) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(-1.5168064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216264) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(-2.8935905) q[0];
rz(-1.2451046) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(-1.9128333) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7354483) q[0];
sx q[0];
rz(-1.7559253) q[0];
sx q[0];
rz(-1.6044751) q[0];
rz(-pi) q[1];
rz(0.51313674) q[2];
sx q[2];
rz(-2.0433132) q[2];
sx q[2];
rz(3.0233011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23763021) q[1];
sx q[1];
rz(-1.3591237) q[1];
sx q[1];
rz(0.71118506) q[1];
rz(-pi) q[2];
x q[2];
rz(3.008083) q[3];
sx q[3];
rz(-0.4236246) q[3];
sx q[3];
rz(-0.51160073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9241141) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(-0.58383101) q[2];
rz(-2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9871224) q[0];
sx q[0];
rz(-0.38726375) q[0];
sx q[0];
rz(1.467147) q[0];
rz(-1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(-2.0416226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8407134) q[0];
sx q[0];
rz(-1.7689118) q[0];
sx q[0];
rz(1.4639356) q[0];
x q[1];
rz(-1.2312382) q[2];
sx q[2];
rz(-1.2969839) q[2];
sx q[2];
rz(1.6664315) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4277122) q[1];
sx q[1];
rz(-1.1888479) q[1];
sx q[1];
rz(1.2101342) q[1];
rz(-0.99697014) q[3];
sx q[3];
rz(-1.1862635) q[3];
sx q[3];
rz(-1.1763193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.786342) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(-1.5477017) q[2];
rz(-1.330438) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(2.0681341) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610157) q[0];
sx q[0];
rz(-1.3207734) q[0];
sx q[0];
rz(0.15783489) q[0];
rz(0.24179587) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(0.48113021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8932228) q[0];
sx q[0];
rz(-0.71740323) q[0];
sx q[0];
rz(0.1347085) q[0];
rz(2.782576) q[2];
sx q[2];
rz(-1.331777) q[2];
sx q[2];
rz(-2.0109107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2896636) q[1];
sx q[1];
rz(-0.73411513) q[1];
sx q[1];
rz(2.0925702) q[1];
x q[2];
rz(-1.8577544) q[3];
sx q[3];
rz(-1.7718414) q[3];
sx q[3];
rz(-1.6359117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87159291) q[2];
sx q[2];
rz(-0.81120482) q[2];
sx q[2];
rz(-2.8672186) q[2];
rz(1.4659878) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(-0.93332851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45193732) q[0];
sx q[0];
rz(-0.87757293) q[0];
sx q[0];
rz(2.080132) q[0];
rz(-0.44867107) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(1.4400858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3502304) q[0];
sx q[0];
rz(-1.4141948) q[0];
sx q[0];
rz(-2.3788405) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8605804) q[2];
sx q[2];
rz(-2.7216879) q[2];
sx q[2];
rz(-0.29422255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0764112) q[1];
sx q[1];
rz(-1.5440327) q[1];
sx q[1];
rz(1.3127432) q[1];
x q[2];
rz(-0.96025709) q[3];
sx q[3];
rz(-1.4818076) q[3];
sx q[3];
rz(1.8833835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3638641) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(1.1931922) q[2];
rz(0.39919546) q[3];
sx q[3];
rz(-1.4123071) q[3];
sx q[3];
rz(-3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4611918) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(-2.4611018) q[0];
rz(-0.83241278) q[1];
sx q[1];
rz(-1.8871555) q[1];
sx q[1];
rz(-2.9852273) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6844821) q[0];
sx q[0];
rz(-1.6231939) q[0];
sx q[0];
rz(-2.8657317) q[0];
rz(-0.84951289) q[2];
sx q[2];
rz(-2.2597183) q[2];
sx q[2];
rz(2.920814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5030843) q[1];
sx q[1];
rz(-1.2933532) q[1];
sx q[1];
rz(1.9126519) q[1];
rz(1.353305) q[3];
sx q[3];
rz(-1.6902349) q[3];
sx q[3];
rz(-0.34298204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51123315) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(0.88999256) q[2];
rz(0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(-2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.8222734) q[0];
sx q[0];
rz(-0.61221468) q[0];
rz(-1.7828434) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(0.30119687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4801203) q[0];
sx q[0];
rz(-1.299322) q[0];
sx q[0];
rz(0.090711509) q[0];
x q[1];
rz(-1.1015755) q[2];
sx q[2];
rz(-0.71327268) q[2];
sx q[2];
rz(0.81339806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14872486) q[1];
sx q[1];
rz(-0.5827671) q[1];
sx q[1];
rz(2.4087743) q[1];
rz(-pi) q[2];
rz(2.4047818) q[3];
sx q[3];
rz(-2.7932146) q[3];
sx q[3];
rz(-0.37955561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23413868) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(-2.0937199) q[2];
rz(-0.35510865) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(0.32383305) q[0];
rz(0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-2.8388265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9186818) q[0];
sx q[0];
rz(-2.7245173) q[0];
sx q[0];
rz(-0.66584754) q[0];
rz(0.19925929) q[2];
sx q[2];
rz(-2.8997921) q[2];
sx q[2];
rz(-1.3346145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1582022) q[1];
sx q[1];
rz(-1.9333955) q[1];
sx q[1];
rz(2.1561978) q[1];
rz(-pi) q[2];
rz(-2.1015424) q[3];
sx q[3];
rz(-0.47785366) q[3];
sx q[3];
rz(-2.8213906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4964464) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(2.5808064) q[2];
rz(1.0049817) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(-2.0463478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305369) q[0];
sx q[0];
rz(-2.472214) q[0];
sx q[0];
rz(0.74991599) q[0];
rz(-0.38052446) q[1];
sx q[1];
rz(-0.98697248) q[1];
sx q[1];
rz(2.1662625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.014054) q[0];
sx q[0];
rz(-0.97146266) q[0];
sx q[0];
rz(-2.4571193) q[0];
rz(3.0815711) q[2];
sx q[2];
rz(-1.2916471) q[2];
sx q[2];
rz(2.699083) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1854461) q[1];
sx q[1];
rz(-2.5393894) q[1];
sx q[1];
rz(-2.5700997) q[1];
rz(0.94176835) q[3];
sx q[3];
rz(-0.36702752) q[3];
sx q[3];
rz(2.0218771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9425977) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(-3.0723451) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-2.8157225) q[0];
sx q[0];
rz(-0.3057873) q[0];
rz(0.30793134) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(-1.1526795) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15241218) q[0];
sx q[0];
rz(-2.3539641) q[0];
sx q[0];
rz(-0.37352011) q[0];
rz(-1.8434373) q[2];
sx q[2];
rz(-2.4783588) q[2];
sx q[2];
rz(-0.014201268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.252723) q[1];
sx q[1];
rz(-2.6111341) q[1];
sx q[1];
rz(2.2525851) q[1];
rz(-pi) q[2];
rz(-3.1265101) q[3];
sx q[3];
rz(-1.507457) q[3];
sx q[3];
rz(-0.36288211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5311188) q[2];
sx q[2];
rz(-2.0788772) q[2];
sx q[2];
rz(-1.6884241) q[2];
rz(0.48062634) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(1.0427955) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(-1.3026793) q[2];
sx q[2];
rz(-2.0901879) q[2];
sx q[2];
rz(-1.2218634) q[2];
rz(-0.35814169) q[3];
sx q[3];
rz(-1.609057) q[3];
sx q[3];
rz(-1.9032794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
