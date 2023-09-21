OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(0.4217622) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5428107) q[0];
sx q[0];
rz(-0.069617696) q[0];
sx q[0];
rz(1.9570062) q[0];
rz(-pi) q[1];
rz(2.4550081) q[2];
sx q[2];
rz(-1.540544) q[2];
sx q[2];
rz(2.9917012) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64695839) q[1];
sx q[1];
rz(-2.0152433) q[1];
sx q[1];
rz(-1.7723945) q[1];
x q[2];
rz(-0.20538879) q[3];
sx q[3];
rz(-2.5377512) q[3];
sx q[3];
rz(2.3103034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(3*pi/11) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.4650311) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37061781) q[0];
sx q[0];
rz(-0.3807225) q[0];
sx q[0];
rz(0.61852635) q[0];
x q[1];
rz(-2.6066166) q[2];
sx q[2];
rz(-1.8638532) q[2];
sx q[2];
rz(-1.5154293) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4620004) q[1];
sx q[1];
rz(-2.3604217) q[1];
sx q[1];
rz(2.6189234) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46918418) q[3];
sx q[3];
rz(-2.5622534) q[3];
sx q[3];
rz(-2.7964696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(-1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(-0.71584654) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(-1.7680426) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64840245) q[0];
sx q[0];
rz(-1.5957386) q[0];
sx q[0];
rz(0.52669749) q[0];
rz(-2.8361514) q[2];
sx q[2];
rz(-0.88103308) q[2];
sx q[2];
rz(-1.6357712) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92257567) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(0.69505691) q[1];
rz(-pi) q[2];
rz(1.8009637) q[3];
sx q[3];
rz(-1.3535415) q[3];
sx q[3];
rz(2.1118856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(-2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90137988) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(2.531321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82204098) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(1.6739474) q[0];
x q[1];
rz(2.3158014) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(0.85130168) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1000881) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(-0.86160223) q[1];
rz(-pi) q[2];
rz(-2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(-1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(1.8283432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816313) q[0];
sx q[0];
rz(-0.93802035) q[0];
sx q[0];
rz(0.40681337) q[0];
x q[1];
rz(3.1057642) q[2];
sx q[2];
rz(-1.6275068) q[2];
sx q[2];
rz(2.5113475) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65391958) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(0.78507702) q[1];
x q[2];
rz(-2.5562708) q[3];
sx q[3];
rz(-0.96671852) q[3];
sx q[3];
rz(-2.4620591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076684549) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(2.6584113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85855243) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(-1.900308) q[0];
x q[1];
rz(-1.6364355) q[2];
sx q[2];
rz(-2.6780431) q[2];
sx q[2];
rz(2.133873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9890081) q[1];
sx q[1];
rz(-0.87859381) q[1];
sx q[1];
rz(2.2425376) q[1];
rz(-pi) q[2];
rz(0.11404927) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(1.9432817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(0.085263578) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(-2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(2.8894997) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68483401) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(2.7177939) q[0];
rz(0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(0.99036723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8512307) q[1];
sx q[1];
rz(-0.92612672) q[1];
sx q[1];
rz(2.6840997) q[1];
x q[2];
rz(-1.1105359) q[3];
sx q[3];
rz(-1.2292976) q[3];
sx q[3];
rz(2.5550492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(-2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.7699014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773974) q[0];
sx q[0];
rz(-1.6590377) q[0];
sx q[0];
rz(-0.95055687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76099446) q[2];
sx q[2];
rz(-1.7241663) q[2];
sx q[2];
rz(-2.209687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5270556) q[1];
sx q[1];
rz(-1.0053047) q[1];
sx q[1];
rz(1.9238227) q[1];
x q[2];
rz(-2.7564704) q[3];
sx q[3];
rz(-0.38673863) q[3];
sx q[3];
rz(-2.573607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(-1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(0.23432215) q[0];
rz(-1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-0.27059069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.035707) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(-2.3653415) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9491896) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(2.2601688) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27442676) q[1];
sx q[1];
rz(-1.5741036) q[1];
sx q[1];
rz(3.0043688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0911646) q[3];
sx q[3];
rz(-1.661146) q[3];
sx q[3];
rz(2.3144212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(-1.489893) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437317) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(0.74238366) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61481793) q[0];
sx q[0];
rz(-0.97133884) q[0];
sx q[0];
rz(-0.74032289) q[0];
x q[1];
rz(-0.43138357) q[2];
sx q[2];
rz(-1.0215534) q[2];
sx q[2];
rz(1.0647578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66354499) q[1];
sx q[1];
rz(-1.5819342) q[1];
sx q[1];
rz(0.34578172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5998245) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(-0.27424973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29006526) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-0.82143299) q[2];
sx q[2];
rz(-1.3315622) q[2];
sx q[2];
rz(-0.17377725) q[2];
rz(-0.7704173) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
