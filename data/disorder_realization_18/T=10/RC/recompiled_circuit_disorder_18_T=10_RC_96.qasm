OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(1.7083038) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(1.1448316) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0776318) q[0];
sx q[0];
rz(-1.504717) q[0];
sx q[0];
rz(1.328701) q[0];
x q[1];
rz(0.28540622) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(-1.0813431) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.634234) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(-1.5276315) q[1];
rz(-0.47031109) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(-2.2944962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23628274) q[0];
sx q[0];
rz(-2.5295527) q[0];
sx q[0];
rz(0.73487868) q[0];
rz(-pi) q[1];
rz(-0.77913021) q[2];
sx q[2];
rz(-1.3574294) q[2];
sx q[2];
rz(1.9386292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3165247) q[1];
sx q[1];
rz(-1.0267338) q[1];
sx q[1];
rz(1.6506509) q[1];
rz(-2.2291525) q[3];
sx q[3];
rz(-0.53386253) q[3];
sx q[3];
rz(-2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(-2.4222597) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(-2.5701994) q[0];
rz(0.96673036) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8674832) q[0];
sx q[0];
rz(-1.880078) q[0];
sx q[0];
rz(-1.5779737) q[0];
rz(0.40076077) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(1.00373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1151162) q[1];
sx q[1];
rz(-0.83362245) q[1];
sx q[1];
rz(1.2181746) q[1];
rz(-pi) q[2];
rz(1.2445883) q[3];
sx q[3];
rz(-2.2766114) q[3];
sx q[3];
rz(2.0446442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6986615) q[0];
sx q[0];
rz(-2.1974265) q[0];
sx q[0];
rz(-2.4311964) q[0];
rz(-pi) q[1];
rz(-0.063886558) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(-2.137616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5570045) q[1];
sx q[1];
rz(-1.6465721) q[1];
sx q[1];
rz(1.9083379) q[1];
rz(-pi) q[2];
rz(-2.5097333) q[3];
sx q[3];
rz(-1.251289) q[3];
sx q[3];
rz(-0.39998049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(-3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3591946) q[0];
sx q[0];
rz(-0.054464666) q[0];
sx q[0];
rz(2.1468303) q[0];
x q[1];
rz(1.1281625) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(0.15304676) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6981814) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(-2.0425668) q[1];
rz(3.0693552) q[3];
sx q[3];
rz(-0.17599711) q[3];
sx q[3];
rz(0.76398677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(0.0897952) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-0.18009137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371532) q[0];
sx q[0];
rz(-2.3817872) q[0];
sx q[0];
rz(-0.8888437) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96713709) q[2];
sx q[2];
rz(-1.7269616) q[2];
sx q[2];
rz(2.9784163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.084701531) q[1];
sx q[1];
rz(-1.4086205) q[1];
sx q[1];
rz(-0.2918891) q[1];
x q[2];
rz(1.6052386) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(-2.8480414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.95773762) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9366074) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(1.3315014) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1646541) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(-0.050780642) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0446339) q[1];
sx q[1];
rz(-0.63475906) q[1];
sx q[1];
rz(1.3332913) q[1];
rz(-pi) q[2];
rz(0.76544806) q[3];
sx q[3];
rz(-2.7141889) q[3];
sx q[3];
rz(0.039507341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-0.17383943) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-2.2824536) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.404495) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.6361902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508142) q[0];
sx q[0];
rz(-0.53508004) q[0];
sx q[0];
rz(-0.65061609) q[0];
rz(-0.6263528) q[2];
sx q[2];
rz(-2.2144631) q[2];
sx q[2];
rz(-2.6477637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7561188) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(0.69617747) q[1];
x q[2];
rz(1.1106311) q[3];
sx q[3];
rz(-2.4880829) q[3];
sx q[3];
rz(1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(2.495893) q[0];
rz(0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.8766778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7495959) q[0];
sx q[0];
rz(-0.8526593) q[0];
sx q[0];
rz(1.0438265) q[0];
rz(2.4931156) q[2];
sx q[2];
rz(-2.4816374) q[2];
sx q[2];
rz(-2.6476268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7281108) q[1];
sx q[1];
rz(-2.9915447) q[1];
sx q[1];
rz(-2.8996182) q[1];
rz(-pi) q[2];
rz(-0.42322741) q[3];
sx q[3];
rz(-1.3864577) q[3];
sx q[3];
rz(-1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(-2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(-0.24542228) q[0];
rz(2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.6773178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368162) q[0];
sx q[0];
rz(-1.4602772) q[0];
sx q[0];
rz(0.21893455) q[0];
rz(-2.9173031) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(1.6155417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-0.34983695) q[1];
rz(-pi) q[2];
rz(-1.7467935) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(-1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650919) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-0.71184288) q[2];
sx q[2];
rz(-1.1887475) q[2];
sx q[2];
rz(2.9954994) q[2];
rz(-2.5916354) q[3];
sx q[3];
rz(-2.4933542) q[3];
sx q[3];
rz(-1.8174432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
