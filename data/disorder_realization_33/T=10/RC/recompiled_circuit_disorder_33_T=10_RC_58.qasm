OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(-1.5712665) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9309064) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(-2.7772285) q[0];
rz(-pi) q[1];
rz(-2.3762796) q[2];
sx q[2];
rz(-1.0368477) q[2];
sx q[2];
rz(-0.050616654) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.03304122) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(0.34717314) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0899815) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(-1.4663565) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(2.9247608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8685535) q[0];
sx q[0];
rz(-0.93398636) q[0];
sx q[0];
rz(0.36006948) q[0];
rz(2.1255323) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(-2.1330657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5237907) q[1];
sx q[1];
rz(-2.4317867) q[1];
sx q[1];
rz(1.8989423) q[1];
x q[2];
rz(-0.64036815) q[3];
sx q[3];
rz(-2.159517) q[3];
sx q[3];
rz(1.1191739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(-2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771616) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(0.60423869) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55721012) q[0];
sx q[0];
rz(-0.32214468) q[0];
sx q[0];
rz(1.4003217) q[0];
x q[1];
rz(0.18205299) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(1.455866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8285117) q[1];
sx q[1];
rz(-0.57144895) q[1];
sx q[1];
rz(-2.9213195) q[1];
x q[2];
rz(0.91498418) q[3];
sx q[3];
rz(-2.4768156) q[3];
sx q[3];
rz(-1.2802326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-2.1742163) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.6436228) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.362975) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(1.0233364) q[0];
x q[1];
rz(1.285032) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(-0.49516585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3264309) q[1];
sx q[1];
rz(-1.8306499) q[1];
sx q[1];
rz(1.8111147) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1180531) q[3];
sx q[3];
rz(-1.1467883) q[3];
sx q[3];
rz(-0.31183576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(-2.4397819) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410626) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5532903) q[0];
sx q[0];
rz(-0.36670812) q[0];
sx q[0];
rz(-2.328863) q[0];
rz(1.5830718) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(2.2001681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.293562) q[1];
sx q[1];
rz(-1.8082431) q[1];
sx q[1];
rz(2.8192987) q[1];
x q[2];
rz(2.1225554) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(-2.2321731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(2.2560789) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-0.12983233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3467305) q[0];
sx q[0];
rz(-0.71160337) q[0];
sx q[0];
rz(2.5559588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0340704) q[2];
sx q[2];
rz(-2.495129) q[2];
sx q[2];
rz(1.549364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0778724) q[1];
sx q[1];
rz(-1.93396) q[1];
sx q[1];
rz(-0.6086463) q[1];
rz(-pi) q[2];
rz(2.8270709) q[3];
sx q[3];
rz(-2.5701227) q[3];
sx q[3];
rz(-0.86626569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3123902) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(-0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-0.68626219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2080363) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(0.74525381) q[0];
rz(-pi) q[1];
rz(-0.83001901) q[2];
sx q[2];
rz(-2.0247211) q[2];
sx q[2];
rz(-1.1656851) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9069179) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(-3.0659552) q[1];
x q[2];
rz(2.5043082) q[3];
sx q[3];
rz(-0.50989671) q[3];
sx q[3];
rz(-2.2989458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(-2.5799675) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.6954533) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.5015645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088889) q[0];
sx q[0];
rz(-2.6484657) q[0];
sx q[0];
rz(-1.5296442) q[0];
rz(-1.6325475) q[2];
sx q[2];
rz(-1.5442863) q[2];
sx q[2];
rz(2.3960631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9403119) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(-2.8708354) q[1];
rz(-pi) q[2];
rz(2.1498508) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(0.21608298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(0.35167545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6707014) q[0];
sx q[0];
rz(-1.1180709) q[0];
sx q[0];
rz(-0.41505138) q[0];
rz(-pi) q[1];
rz(-1.8750538) q[2];
sx q[2];
rz(-0.90196246) q[2];
sx q[2];
rz(-0.63389102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55191509) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(-3.0685436) q[1];
x q[2];
rz(-2.6978108) q[3];
sx q[3];
rz(-3.0203331) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7982771) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-0.19432755) q[0];
rz(1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(2.1077572) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68034222) q[0];
sx q[0];
rz(-2.470812) q[0];
sx q[0];
rz(1.781342) q[0];
rz(-pi) q[1];
rz(2.1543703) q[2];
sx q[2];
rz(-0.7910896) q[2];
sx q[2];
rz(0.9466048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49627134) q[1];
sx q[1];
rz(-1.6722582) q[1];
sx q[1];
rz(-2.1437777) q[1];
rz(-pi) q[2];
rz(-0.71318993) q[3];
sx q[3];
rz(-1.7367559) q[3];
sx q[3];
rz(0.99456577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(-1.4355961) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(0.031899115) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(1.3676436) q[3];
sx q[3];
rz(-2.805134) q[3];
sx q[3];
rz(0.72611879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
