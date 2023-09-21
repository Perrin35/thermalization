OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(-1.4810286) q[0];
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45863736) q[0];
sx q[0];
rz(-1.4541172) q[0];
sx q[0];
rz(1.7200243) q[0];
x q[1];
rz(0.68732287) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-0.99422115) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8926156) q[1];
sx q[1];
rz(-1.2274582) q[1];
sx q[1];
rz(2.9825319) q[1];
x q[2];
rz(0.25306563) q[3];
sx q[3];
rz(-2.1742637) q[3];
sx q[3];
rz(2.8367708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(2.375405) q[2];
rz(-0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8121174) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-3.0564953) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21762411) q[0];
sx q[0];
rz(-0.94857615) q[0];
sx q[0];
rz(1.5909965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9688103) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(1.3041376) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47463402) q[1];
sx q[1];
rz(-2.9999795) q[1];
sx q[1];
rz(2.2255564) q[1];
rz(-0.3195023) q[3];
sx q[3];
rz(-1.1211683) q[3];
sx q[3];
rz(3.086123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(-1.7571626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2565838) q[0];
sx q[0];
rz(-2.6589767) q[0];
sx q[0];
rz(-2.9628739) q[0];
x q[1];
rz(-1.0648107) q[2];
sx q[2];
rz(-1.2227321) q[2];
sx q[2];
rz(-0.51007523) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3123734) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(-1.6014598) q[1];
x q[2];
rz(0.82377394) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(0.87574524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(1.7975851) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-0.20733325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79065381) q[0];
sx q[0];
rz(-1.5293984) q[0];
sx q[0];
rz(-0.008965094) q[0];
rz(-pi) q[1];
rz(-1.7447128) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(-2.7001691) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.748121) q[1];
sx q[1];
rz(-1.8078184) q[1];
sx q[1];
rz(2.9841828) q[1];
rz(-pi) q[2];
rz(-1.5563577) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(-0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(2.6079544) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1874439) q[0];
sx q[0];
rz(-1.2866409) q[0];
sx q[0];
rz(2.8548293) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9075378) q[2];
sx q[2];
rz(-2.3599527) q[2];
sx q[2];
rz(-2.3603338) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4644949) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(2.5527843) q[1];
x q[2];
rz(-2.2002033) q[3];
sx q[3];
rz(-1.7377059) q[3];
sx q[3];
rz(-3.133579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(2.1369381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83197901) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(-1.2020822) q[0];
rz(-pi) q[1];
x q[1];
rz(1.651628) q[2];
sx q[2];
rz(-1.9664552) q[2];
sx q[2];
rz(-2.253502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.172563) q[1];
sx q[1];
rz(-1.8168212) q[1];
sx q[1];
rz(0.56117705) q[1];
x q[2];
rz(-0.011990487) q[3];
sx q[3];
rz(-1.9846091) q[3];
sx q[3];
rz(3.0179265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7375609) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7710966) q[0];
sx q[0];
rz(-1.578195) q[0];
sx q[0];
rz(-1.7367944) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80588801) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(-2.9758331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8183648) q[1];
sx q[1];
rz(-1.0867026) q[1];
sx q[1];
rz(-1.1546043) q[1];
x q[2];
rz(-0.70127212) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(0.14227223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(2.9833802) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(0.80668443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5935532) q[0];
sx q[0];
rz(-1.1386477) q[0];
sx q[0];
rz(-0.12950626) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47862349) q[2];
sx q[2];
rz(-2.9106986) q[2];
sx q[2];
rz(-0.17639562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.403277) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(-2.7110389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8716378) q[3];
sx q[3];
rz(-1.6013718) q[3];
sx q[3];
rz(-1.2310864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(2.1419443) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(-1.1243533) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(0.65752423) q[0];
rz(-0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(0.34067571) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1239615) q[0];
sx q[0];
rz(-0.88273662) q[0];
sx q[0];
rz(0.4171564) q[0];
rz(-pi) q[1];
rz(1.0057783) q[2];
sx q[2];
rz(-1.7592332) q[2];
sx q[2];
rz(-2.4609158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.162519) q[1];
sx q[1];
rz(-0.58774978) q[1];
sx q[1];
rz(1.9025004) q[1];
x q[2];
rz(0.19823234) q[3];
sx q[3];
rz(-1.7361589) q[3];
sx q[3];
rz(0.73393047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(-0.07671193) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4519276) q[0];
sx q[0];
rz(-1.2566914) q[0];
sx q[0];
rz(3.0118224) q[0];
rz(-pi) q[1];
rz(1.3086583) q[2];
sx q[2];
rz(-2.5524676) q[2];
sx q[2];
rz(0.72074705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0127276) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(0.058458316) q[1];
rz(-pi) q[2];
rz(-1.4823227) q[3];
sx q[3];
rz(-1.3946597) q[3];
sx q[3];
rz(0.74151553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(2.0417487) q[2];
sx q[2];
rz(-1.8269314) q[2];
sx q[2];
rz(-1.6614428) q[2];
rz(-1.205411) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
