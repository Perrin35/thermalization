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
rz(0.29851222) q[0];
sx q[0];
rz(-1.8520344) q[0];
sx q[0];
rz(0.02331743) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(-2.4040931) q[1];
sx q[1];
rz(1.6657383) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6029089) q[0];
sx q[0];
rz(-1.5011461) q[0];
sx q[0];
rz(-3.0979741) q[0];
x q[1];
rz(-1.453773) q[2];
sx q[2];
rz(-0.3323148) q[2];
sx q[2];
rz(1.3261283) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1920584) q[1];
sx q[1];
rz(-0.84334521) q[1];
sx q[1];
rz(-2.6019168) q[1];
x q[2];
rz(2.2470583) q[3];
sx q[3];
rz(-1.9893551) q[3];
sx q[3];
rz(-0.422131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6797356) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(-2.8409345) q[2];
rz(0.65827185) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(-0.58832735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268426) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(-2.9657189) q[0];
rz(-2.580592) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(0.19634136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6602621) q[0];
sx q[0];
rz(-2.5721533) q[0];
sx q[0];
rz(-0.91668949) q[0];
x q[1];
rz(0.11949338) q[2];
sx q[2];
rz(-1.6867078) q[2];
sx q[2];
rz(-1.3136065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2113116) q[1];
sx q[1];
rz(-1.0553816) q[1];
sx q[1];
rz(3.1392908) q[1];
rz(-pi) q[2];
rz(-2.7982233) q[3];
sx q[3];
rz(-1.4447304) q[3];
sx q[3];
rz(2.8186563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0520797) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(1.6717795) q[2];
rz(-2.6164264) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(2.4770881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1216275) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(-2.5110733) q[0];
rz(1.9969253) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(0.05680457) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81679427) q[0];
sx q[0];
rz(-1.9787119) q[0];
sx q[0];
rz(-2.7403031) q[0];
rz(-2.2571083) q[2];
sx q[2];
rz(-1.9504042) q[2];
sx q[2];
rz(1.1432709) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89598846) q[1];
sx q[1];
rz(-1.1088246) q[1];
sx q[1];
rz(1.4126591) q[1];
x q[2];
rz(-2.9288559) q[3];
sx q[3];
rz(-1.9895344) q[3];
sx q[3];
rz(1.1513082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28222617) q[2];
sx q[2];
rz(-1.6702009) q[2];
sx q[2];
rz(-2.7024506) q[2];
rz(-1.3965083) q[3];
sx q[3];
rz(-1.2625932) q[3];
sx q[3];
rz(-0.5784353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0948931) q[0];
sx q[0];
rz(-0.72143227) q[0];
sx q[0];
rz(-2.1266345) q[0];
rz(3.0789442) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(0.28422022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009722) q[0];
sx q[0];
rz(-0.56864244) q[0];
sx q[0];
rz(2.4533662) q[0];
x q[1];
rz(-0.49782984) q[2];
sx q[2];
rz(-1.4477028) q[2];
sx q[2];
rz(0.13409889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38317162) q[1];
sx q[1];
rz(-1.9291506) q[1];
sx q[1];
rz(1.5580295) q[1];
x q[2];
rz(1.3174414) q[3];
sx q[3];
rz(-1.5572963) q[3];
sx q[3];
rz(2.6576192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3011498) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(-0.76048771) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-1.1287929) q[3];
sx q[3];
rz(1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29192057) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(-2.736295) q[0];
rz(-2.6536476) q[1];
sx q[1];
rz(-2.0351724) q[1];
sx q[1];
rz(2.778756) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0340509) q[0];
sx q[0];
rz(-2.0011609) q[0];
sx q[0];
rz(1.2314046) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53952257) q[2];
sx q[2];
rz(-2.4837257) q[2];
sx q[2];
rz(-2.4488329) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2606973) q[1];
sx q[1];
rz(-2.915307) q[1];
sx q[1];
rz(-0.21842167) q[1];
x q[2];
rz(2.7801301) q[3];
sx q[3];
rz(-2.4130954) q[3];
sx q[3];
rz(2.0036774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1399347) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(0.72251594) q[2];
rz(2.6288988) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743415) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(2.784506) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(1.9836609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0561075) q[0];
sx q[0];
rz(-0.72626136) q[0];
sx q[0];
rz(1.9421436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59717565) q[2];
sx q[2];
rz(-0.9998601) q[2];
sx q[2];
rz(-1.2025637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14640644) q[1];
sx q[1];
rz(-0.9764834) q[1];
sx q[1];
rz(-1.3239149) q[1];
x q[2];
rz(-0.66863184) q[3];
sx q[3];
rz(-0.9210081) q[3];
sx q[3];
rz(-1.970552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57924119) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-2.5518899) q[2];
rz(-1.3127182) q[3];
sx q[3];
rz(-1.6629985) q[3];
sx q[3];
rz(-2.5780799) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967459) q[0];
sx q[0];
rz(-0.99269301) q[0];
sx q[0];
rz(0.27031159) q[0];
rz(-2.7145794) q[1];
sx q[1];
rz(-1.7662363) q[1];
sx q[1];
rz(2.7332773) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1751375) q[0];
sx q[0];
rz(-1.9388559) q[0];
sx q[0];
rz(-1.4428663) q[0];
rz(-pi) q[1];
rz(0.87323453) q[2];
sx q[2];
rz(-1.4931275) q[2];
sx q[2];
rz(-0.49688646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9198614) q[1];
sx q[1];
rz(-1.1581005) q[1];
sx q[1];
rz(-2.2151715) q[1];
x q[2];
rz(2.4189878) q[3];
sx q[3];
rz(-1.4935023) q[3];
sx q[3];
rz(2.3376436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2030877) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(-2.2173524) q[2];
rz(-2.3067394) q[3];
sx q[3];
rz(-2.3122841) q[3];
sx q[3];
rz(2.1448962) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860109) q[0];
sx q[0];
rz(-0.38014933) q[0];
sx q[0];
rz(2.7741449) q[0];
rz(2.5732749) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(2.7265991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4545373) q[0];
sx q[0];
rz(-2.6090985) q[0];
sx q[0];
rz(2.103602) q[0];
rz(0.52936036) q[2];
sx q[2];
rz(-1.0870013) q[2];
sx q[2];
rz(-2.8279378) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69232443) q[1];
sx q[1];
rz(-1.0639166) q[1];
sx q[1];
rz(0.70517069) q[1];
x q[2];
rz(1.2250118) q[3];
sx q[3];
rz(-1.5756851) q[3];
sx q[3];
rz(1.1556311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7830398) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(-2.3084194) q[2];
rz(0.35946515) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.6370714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9678765) q[0];
sx q[0];
rz(-0.40206566) q[0];
sx q[0];
rz(-3.1267082) q[0];
rz(1.124565) q[1];
sx q[1];
rz(-0.24528565) q[1];
sx q[1];
rz(0.10717779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10725645) q[0];
sx q[0];
rz(-1.7198623) q[0];
sx q[0];
rz(2.3019522) q[0];
x q[1];
rz(-0.82294269) q[2];
sx q[2];
rz(-1.2185893) q[2];
sx q[2];
rz(0.65300452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31871492) q[1];
sx q[1];
rz(-1.4221622) q[1];
sx q[1];
rz(-0.48245247) q[1];
rz(-pi) q[2];
rz(2.8802683) q[3];
sx q[3];
rz(-1.3991941) q[3];
sx q[3];
rz(-1.2444519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(2.2819819) q[2];
rz(-0.25935069) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(-0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6421826) q[0];
sx q[0];
rz(-0.12633093) q[0];
sx q[0];
rz(-1.9835749) q[0];
rz(0.49531373) q[1];
sx q[1];
rz(-1.9751578) q[1];
sx q[1];
rz(-0.29702979) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659781) q[0];
sx q[0];
rz(-0.31904623) q[0];
sx q[0];
rz(-1.5837529) q[0];
rz(0.30029617) q[2];
sx q[2];
rz(-1.5321333) q[2];
sx q[2];
rz(3.0186332) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1446643) q[1];
sx q[1];
rz(-1.7283617) q[1];
sx q[1];
rz(-2.8466606) q[1];
rz(3.1252485) q[3];
sx q[3];
rz(-2.2333003) q[3];
sx q[3];
rz(-1.7733044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26934066) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(2.2807109) q[2];
rz(-2.6779029) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(1.1369107) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9853482) q[0];
sx q[0];
rz(-1.5736268) q[0];
sx q[0];
rz(-1.9851984) q[0];
rz(-1.3337749) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(2.479066) q[2];
sx q[2];
rz(-1.1077729) q[2];
sx q[2];
rz(-0.33453579) q[2];
rz(-1.9577338) q[3];
sx q[3];
rz(-1.1087742) q[3];
sx q[3];
rz(-2.5840989) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
