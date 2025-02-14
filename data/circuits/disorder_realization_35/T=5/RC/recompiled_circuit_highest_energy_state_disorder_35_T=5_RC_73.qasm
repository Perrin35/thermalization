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
rz(-2.2292697) q[0];
sx q[0];
rz(1.8066701) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700884) q[0];
sx q[0];
rz(-3.0850924) q[0];
sx q[0];
rz(1.9046049) q[0];
x q[1];
rz(-2.4249627) q[2];
sx q[2];
rz(-0.83558768) q[2];
sx q[2];
rz(2.731833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0483889) q[1];
sx q[1];
rz(-1.2902765) q[1];
sx q[1];
rz(0.93757052) q[1];
x q[2];
rz(3.1077216) q[3];
sx q[3];
rz(-1.3283786) q[3];
sx q[3];
rz(-2.9592379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7923183) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(0.054923687) q[2];
rz(1.4548291) q[3];
sx q[3];
rz(-2.7456386) q[3];
sx q[3];
rz(-1.5168064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216264) q[0];
sx q[0];
rz(-1.8333789) q[0];
sx q[0];
rz(0.24800214) q[0];
rz(-1.2451046) q[1];
sx q[1];
rz(-0.73687941) q[1];
sx q[1];
rz(-1.2287593) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7354483) q[0];
sx q[0];
rz(-1.3856674) q[0];
sx q[0];
rz(1.5371176) q[0];
rz(2.1013772) q[2];
sx q[2];
rz(-2.0231915) q[2];
sx q[2];
rz(-1.7035521) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23763021) q[1];
sx q[1];
rz(-1.3591237) q[1];
sx q[1];
rz(-2.4304076) q[1];
x q[2];
rz(1.6307488) q[3];
sx q[3];
rz(-1.9904117) q[3];
sx q[3];
rz(-2.7762716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2174786) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(2.5577616) q[2];
rz(2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(1.0113641) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-0.38726375) q[0];
sx q[0];
rz(1.6744457) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.4091622) q[1];
sx q[1];
rz(-2.0416226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29102688) q[0];
sx q[0];
rz(-1.675559) q[0];
sx q[0];
rz(2.9423703) q[0];
rz(-0.87022484) q[2];
sx q[2];
rz(-2.7087475) q[2];
sx q[2];
rz(-0.74898042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4277122) q[1];
sx q[1];
rz(-1.9527447) q[1];
sx q[1];
rz(-1.2101342) q[1];
x q[2];
rz(-2.1446225) q[3];
sx q[3];
rz(-1.9553292) q[3];
sx q[3];
rz(-1.1763193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.786342) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(-1.5477017) q[2];
rz(1.8111546) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(-1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.380577) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(-0.15783489) q[0];
rz(-0.24179587) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(-2.6604624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07040992) q[0];
sx q[0];
rz(-2.2803223) q[0];
sx q[0];
rz(-1.6874403) q[0];
x q[1];
rz(0.60637668) q[2];
sx q[2];
rz(-0.42841347) q[2];
sx q[2];
rz(-2.1386752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0193023) q[1];
sx q[1];
rz(-1.230352) q[1];
sx q[1];
rz(-0.90695088) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8577544) q[3];
sx q[3];
rz(-1.7718414) q[3];
sx q[3];
rz(1.6359117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87159291) q[2];
sx q[2];
rz(-0.81120482) q[2];
sx q[2];
rz(-2.8672186) q[2];
rz(-1.4659878) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(0.93332851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6896553) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(-1.0614606) q[0];
rz(-0.44867107) q[1];
sx q[1];
rz(-2.6866388) q[1];
sx q[1];
rz(1.7015069) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7730624) q[0];
sx q[0];
rz(-0.8196747) q[0];
sx q[0];
rz(1.7858265) q[0];
rz(-pi) q[1];
rz(1.8605804) q[2];
sx q[2];
rz(-0.41990478) q[2];
sx q[2];
rz(-2.8473701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0764112) q[1];
sx q[1];
rz(-1.5440327) q[1];
sx q[1];
rz(-1.3127432) q[1];
x q[2];
rz(3.0331221) q[3];
sx q[3];
rz(-2.1785695) q[3];
sx q[3];
rz(2.8911107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.8455467) q[2];
sx q[2];
rz(1.9484005) q[2];
rz(2.7423972) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(-3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6804009) q[0];
sx q[0];
rz(-1.8541279) q[0];
sx q[0];
rz(-2.4611018) q[0];
rz(-2.3091799) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(0.15636538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45711058) q[0];
sx q[0];
rz(-1.6231939) q[0];
sx q[0];
rz(-2.8657317) q[0];
rz(-0.6758718) q[2];
sx q[2];
rz(-0.9524494) q[2];
sx q[2];
rz(2.4175274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.63850832) q[1];
sx q[1];
rz(-1.8482395) q[1];
sx q[1];
rz(1.9126519) q[1];
x q[2];
rz(-1.353305) q[3];
sx q[3];
rz(-1.6902349) q[3];
sx q[3];
rz(0.34298204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51123315) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(-2.2516001) q[2];
rz(2.6651799) q[3];
sx q[3];
rz(-0.92372957) q[3];
sx q[3];
rz(0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(0.61221468) q[0];
rz(-1.7828434) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(-2.8403958) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66147236) q[0];
sx q[0];
rz(-1.299322) q[0];
sx q[0];
rz(3.0508811) q[0];
rz(-pi) q[1];
rz(-0.91353784) q[2];
sx q[2];
rz(-1.2704256) q[2];
sx q[2];
rz(-1.1235088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4674356) q[1];
sx q[1];
rz(-1.9922207) q[1];
sx q[1];
rz(-1.155505) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8101569) q[3];
sx q[3];
rz(-1.3151957) q[3];
sx q[3];
rz(1.1472697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.907454) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(0.35510865) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(-2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2316786) q[0];
sx q[0];
rz(-2.7993918) q[0];
sx q[0];
rz(0.32383305) q[0];
rz(-0.58770761) q[1];
sx q[1];
rz(-2.4617742) q[1];
sx q[1];
rz(-2.8388265) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1708831) q[0];
sx q[0];
rz(-1.3178749) q[0];
sx q[0];
rz(-2.80632) q[0];
rz(-pi) q[1];
rz(-1.6195756) q[2];
sx q[2];
rz(-1.3338727) q[2];
sx q[2];
rz(1.1295484) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3564612) q[1];
sx q[1];
rz(-2.11368) q[1];
sx q[1];
rz(-2.7144542) q[1];
rz(-pi) q[2];
rz(-0.25637324) q[3];
sx q[3];
rz(-1.9786165) q[3];
sx q[3];
rz(0.90424171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6451463) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(-2.5808064) q[2];
rz(1.0049817) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(2.0463478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305369) q[0];
sx q[0];
rz(-2.472214) q[0];
sx q[0];
rz(0.74991599) q[0];
rz(2.7610682) q[1];
sx q[1];
rz(-0.98697248) q[1];
sx q[1];
rz(2.1662625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.153424) q[0];
sx q[0];
rz(-2.1200709) q[0];
sx q[0];
rz(-0.84814056) q[0];
rz(-pi) q[1];
rz(1.2911694) q[2];
sx q[2];
rz(-1.5131009) q[2];
sx q[2];
rz(-1.1448432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6187894) q[1];
sx q[1];
rz(-2.0674043) q[1];
sx q[1];
rz(-1.2148395) q[1];
rz(1.2694025) q[3];
sx q[3];
rz(-1.7835254) q[3];
sx q[3];
rz(-1.0476867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9425977) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(2.060176) q[2];
rz(3.0723451) q[3];
sx q[3];
rz(-1.4360177) q[3];
sx q[3];
rz(-0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0104495) q[0];
sx q[0];
rz(-2.8157225) q[0];
sx q[0];
rz(0.3057873) q[0];
rz(-0.30793134) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(-1.9889132) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9891805) q[0];
sx q[0];
rz(-2.3539641) q[0];
sx q[0];
rz(2.7680725) q[0];
x q[1];
rz(-2.2158578) q[2];
sx q[2];
rz(-1.7373499) q[2];
sx q[2];
rz(-1.8018617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0704502) q[1];
sx q[1];
rz(-1.8952888) q[1];
sx q[1];
rz(-1.1434511) q[1];
rz(-1.3373371) q[3];
sx q[3];
rz(-0.065107927) q[3];
sx q[3];
rz(-2.5447735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61047381) q[2];
sx q[2];
rz(-2.0788772) q[2];
sx q[2];
rz(1.6884241) q[2];
rz(-0.48062634) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(-1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65344812) q[0];
sx q[0];
rz(-1.6076037) q[0];
sx q[0];
rz(-1.6758767) q[0];
rz(1.0427955) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(-1.8389134) q[2];
sx q[2];
rz(-1.0514048) q[2];
sx q[2];
rz(1.9197293) q[2];
rz(-1.5299464) q[3];
sx q[3];
rz(-1.9286641) q[3];
sx q[3];
rz(2.7947938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
