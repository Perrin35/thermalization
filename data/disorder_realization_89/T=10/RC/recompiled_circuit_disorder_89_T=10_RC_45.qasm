OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5428107) q[0];
sx q[0];
rz(-0.069617696) q[0];
sx q[0];
rz(1.9570062) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4550081) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(0.1498915) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1301073) q[1];
sx q[1];
rz(-1.752578) q[1];
sx q[1];
rz(2.6891522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5476417) q[3];
sx q[3];
rz(-1.6868601) q[3];
sx q[3];
rz(-0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(8*pi/11) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0533957) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-0.85900599) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.3536317) q[0];
sx q[0];
rz(0.3152245) q[0];
rz(-1.2334521) q[2];
sx q[2];
rz(-1.0609027) q[2];
sx q[2];
rz(0.22491977) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6444781) q[1];
sx q[1];
rz(-1.2116355) q[1];
sx q[1];
rz(-2.4317846) q[1];
x q[2];
rz(-2.6724085) q[3];
sx q[3];
rz(-0.57933925) q[3];
sx q[3];
rz(0.34512305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(2.583288) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(-1.7680426) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2620619) q[0];
sx q[0];
rz(-0.52723215) q[0];
sx q[0];
rz(-3.0920045) q[0];
rz(-pi) q[1];
rz(-1.9203556) q[2];
sx q[2];
rz(-2.3974843) q[2];
sx q[2];
rz(2.0958054) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9436685) q[1];
sx q[1];
rz(-2.1197882) q[1];
sx q[1];
rz(-0.82358574) q[1];
x q[2];
rz(2.3397066) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(1.2847881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(2.1219357) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85128879) q[0];
sx q[0];
rz(-1.5595058) q[0];
sx q[0];
rz(-1.6800866) q[0];
rz(-pi) q[1];
rz(2.3158014) q[2];
sx q[2];
rz(-0.50160393) q[2];
sx q[2];
rz(2.290291) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8205386) q[1];
sx q[1];
rz(-0.74756261) q[1];
sx q[1];
rz(-1.9588406) q[1];
x q[2];
rz(-2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(2.9042517) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(-1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.6089815) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.3132494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4603235) q[0];
sx q[0];
rz(-1.2460421) q[0];
sx q[0];
rz(-2.2446509) q[0];
rz(-pi) q[1];
rz(3.1057642) q[2];
sx q[2];
rz(-1.5140859) q[2];
sx q[2];
rz(-2.5113475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6254127) q[1];
sx q[1];
rz(-1.9429632) q[1];
sx q[1];
rz(-0.401293) q[1];
rz(-0.89574121) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(-2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1398754) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076684549) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(0.48318133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1937712) q[0];
sx q[0];
rz(-1.7521439) q[0];
sx q[0];
rz(0.062747196) q[0];
rz(-pi) q[1];
rz(-3.1088164) q[2];
sx q[2];
rz(-1.1083229) q[2];
sx q[2];
rz(1.0810766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15258458) q[1];
sx q[1];
rz(-2.2629988) q[1];
sx q[1];
rz(-0.89905507) q[1];
x q[2];
rz(-0.11404927) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(1.9412458) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(2.8894997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68483401) q[0];
sx q[0];
rz(-1.4258458) q[0];
sx q[0];
rz(0.42379871) q[0];
x q[1];
rz(-0.78105314) q[2];
sx q[2];
rz(-2.1456246) q[2];
sx q[2];
rz(0.99036723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.568087) q[1];
sx q[1];
rz(-1.9316257) q[1];
sx q[1];
rz(2.2682857) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0310568) q[3];
sx q[3];
rz(-1.2292976) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.4138387) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.3716912) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1641952) q[0];
sx q[0];
rz(-1.482555) q[0];
sx q[0];
rz(2.1910358) q[0];
rz(-pi) q[1];
rz(0.76099446) q[2];
sx q[2];
rz(-1.4174263) q[2];
sx q[2];
rz(-0.93190565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5270556) q[1];
sx q[1];
rz(-2.136288) q[1];
sx q[1];
rz(1.2177699) q[1];
rz(-pi) q[2];
rz(-0.38512226) q[3];
sx q[3];
rz(-2.754854) q[3];
sx q[3];
rz(0.56798565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82683212) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-0.27059069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9070248) q[0];
sx q[0];
rz(-1.3493291) q[0];
sx q[0];
rz(-1.7937167) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0653535) q[2];
sx q[2];
rz(-2.7521172) q[2];
sx q[2];
rz(-1.7324093) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2968263) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(-1.5674577) q[1];
rz(-pi) q[2];
rz(-3.0398265) q[3];
sx q[3];
rz(-1.093285) q[3];
sx q[3];
rz(2.4448642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.6516997) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001227) q[0];
sx q[0];
rz(-2.1614657) q[0];
sx q[0];
rz(2.3175879) q[0];
rz(-pi) q[1];
rz(0.97147271) q[2];
sx q[2];
rz(-2.4571672) q[2];
sx q[2];
rz(-1.3542092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91126373) q[1];
sx q[1];
rz(-1.2250369) q[1];
sx q[1];
rz(1.5589578) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9418482) q[3];
sx q[3];
rz(-2.9962857) q[3];
sx q[3];
rz(-2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(2.3201597) q[2];
sx q[2];
rz(-1.3315622) q[2];
sx q[2];
rz(-0.17377725) q[2];
rz(3.001694) q[3];
sx q[3];
rz(-1.7053053) q[3];
sx q[3];
rz(-2.7444621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
