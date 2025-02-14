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
rz(-1.3402101) q[0];
sx q[0];
rz(-2.8647381) q[0];
sx q[0];
rz(1.0387596) q[0];
rz(2.7998595) q[1];
sx q[1];
rz(-0.83655292) q[1];
sx q[1];
rz(2.7247782) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57274103) q[0];
sx q[0];
rz(-2.4089455) q[0];
sx q[0];
rz(-2.9824663) q[0];
rz(-1.8174581) q[2];
sx q[2];
rz(-2.1132872) q[2];
sx q[2];
rz(0.88298029) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8326679) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(0.69242386) q[1];
rz(-pi) q[2];
rz(-0.19617041) q[3];
sx q[3];
rz(-2.9769398) q[3];
sx q[3];
rz(-2.1712554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(1.8876342) q[2];
rz(-1.5818671) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9480243) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(0.76876202) q[0];
rz(1.3668775) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(-2.1034525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6383253) q[0];
sx q[0];
rz(-1.5621788) q[0];
sx q[0];
rz(-1.5805954) q[0];
rz(-2.4203796) q[2];
sx q[2];
rz(-1.8497457) q[2];
sx q[2];
rz(-0.21706377) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86377108) q[1];
sx q[1];
rz(-1.9367332) q[1];
sx q[1];
rz(2.4863431) q[1];
rz(2.059405) q[3];
sx q[3];
rz(-1.587579) q[3];
sx q[3];
rz(0.29573655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.495503) q[2];
sx q[2];
rz(-1.2868519) q[2];
sx q[2];
rz(0.82026473) q[2];
rz(-0.25343728) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(2.7804815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609817) q[0];
sx q[0];
rz(-2.6237223) q[0];
sx q[0];
rz(-2.2542727) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(1.1393772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761918) q[0];
sx q[0];
rz(-1.7125107) q[0];
sx q[0];
rz(-2.7928674) q[0];
x q[1];
rz(1.0372889) q[2];
sx q[2];
rz(-2.5441558) q[2];
sx q[2];
rz(3.1410599) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92147747) q[1];
sx q[1];
rz(-2.8954828) q[1];
sx q[1];
rz(-3.0819478) q[1];
rz(1.5185131) q[3];
sx q[3];
rz(-1.0526592) q[3];
sx q[3];
rz(-0.62925807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1246216) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-0.43928453) q[2];
rz(0.47075054) q[3];
sx q[3];
rz(-0.079340383) q[3];
sx q[3];
rz(-1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(-3.0261107) q[0];
rz(3.0237517) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(1.3062564) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4233526) q[0];
sx q[0];
rz(-2.6246492) q[0];
sx q[0];
rz(2.1812308) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9491862) q[2];
sx q[2];
rz(-0.94589627) q[2];
sx q[2];
rz(-1.6922127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0769121) q[1];
sx q[1];
rz(-2.2222328) q[1];
sx q[1];
rz(-2.3802735) q[1];
rz(0.99808295) q[3];
sx q[3];
rz(-1.6340165) q[3];
sx q[3];
rz(1.8259189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8370886) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(2.9869249) q[2];
rz(-1.6020487) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(2.535517) q[3];
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
rz(-2.5681169) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(0.83531761) q[0];
rz(-0.71594816) q[1];
sx q[1];
rz(-2.7138111) q[1];
sx q[1];
rz(0.11944019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855921) q[0];
sx q[0];
rz(-1.5876666) q[0];
sx q[0];
rz(-1.5149679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7238462) q[2];
sx q[2];
rz(-2.0156392) q[2];
sx q[2];
rz(-2.9314624) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24158289) q[1];
sx q[1];
rz(-1.1092343) q[1];
sx q[1];
rz(-0.66630967) q[1];
x q[2];
rz(-2.6089657) q[3];
sx q[3];
rz(-1.2029193) q[3];
sx q[3];
rz(2.0834578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2401838) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(0.048246233) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(0.72859305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.401684) q[0];
sx q[0];
rz(-0.70867276) q[0];
sx q[0];
rz(2.0106864) q[0];
rz(-2.6629958) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(1.7344249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3597108) q[0];
sx q[0];
rz(-3.0880083) q[0];
sx q[0];
rz(-3.1126541) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2532151) q[2];
sx q[2];
rz(-0.7478315) q[2];
sx q[2];
rz(-1.795447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6428522) q[1];
sx q[1];
rz(-1.7987218) q[1];
sx q[1];
rz(-3.0251316) q[1];
rz(-pi) q[2];
rz(-2.9286372) q[3];
sx q[3];
rz(-1.7476488) q[3];
sx q[3];
rz(2.9573754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5799134) q[2];
sx q[2];
rz(-2.734197) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(-2.0922349) q[3];
sx q[3];
rz(-0.5368084) q[3];
sx q[3];
rz(1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9310164) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(1.806102) q[0];
rz(-1.3212851) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(-2.8335422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0548693) q[0];
sx q[0];
rz(-2.6941774) q[0];
sx q[0];
rz(-2.6347876) q[0];
rz(-0.51527889) q[2];
sx q[2];
rz(-1.2211498) q[2];
sx q[2];
rz(-2.6315029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8561473) q[1];
sx q[1];
rz(-2.3458907) q[1];
sx q[1];
rz(-0.85983069) q[1];
x q[2];
rz(1.9568029) q[3];
sx q[3];
rz(-1.6148477) q[3];
sx q[3];
rz(-3.0795627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41219741) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(-0.13128734) q[2];
rz(2.4042551) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(-0.15784119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8312296) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(0.53037733) q[0];
rz(1.4250379) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(-0.72149611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3193766) q[0];
sx q[0];
rz(-1.1255029) q[0];
sx q[0];
rz(0.67586835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7185831) q[2];
sx q[2];
rz(-2.035454) q[2];
sx q[2];
rz(0.3275268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65039448) q[1];
sx q[1];
rz(-1.7868687) q[1];
sx q[1];
rz(-2.9465527) q[1];
rz(0.11105342) q[3];
sx q[3];
rz(-0.78563443) q[3];
sx q[3];
rz(1.812404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43712744) q[2];
sx q[2];
rz(-0.4054873) q[2];
sx q[2];
rz(-1.6638157) q[2];
rz(1.9780698) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(-1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32996938) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(2.5323618) q[0];
rz(-1.5902663) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.685166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1066172) q[0];
sx q[0];
rz(-1.9272695) q[0];
sx q[0];
rz(-2.8750505) q[0];
x q[1];
rz(1.6316192) q[2];
sx q[2];
rz(-2.0552539) q[2];
sx q[2];
rz(-1.4549507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85553414) q[1];
sx q[1];
rz(-2.0288336) q[1];
sx q[1];
rz(0.18064255) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2942469) q[3];
sx q[3];
rz(-2.06978) q[3];
sx q[3];
rz(2.8049198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41813254) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(0.88031236) q[2];
rz(-0.20600016) q[3];
sx q[3];
rz(-2.7166631) q[3];
sx q[3];
rz(-0.4900842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708165) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(3.1296375) q[0];
rz(-1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(2.5451122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947114) q[0];
sx q[0];
rz(-0.46976006) q[0];
sx q[0];
rz(-0.042094783) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61052236) q[2];
sx q[2];
rz(-0.43402616) q[2];
sx q[2];
rz(-0.69930017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47190753) q[1];
sx q[1];
rz(-1.0958671) q[1];
sx q[1];
rz(-0.088598786) q[1];
x q[2];
rz(-1.7285198) q[3];
sx q[3];
rz(-1.8650569) q[3];
sx q[3];
rz(2.074948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.97800469) q[2];
sx q[2];
rz(-0.96077335) q[2];
sx q[2];
rz(2.9445924) q[2];
rz(1.9893076) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(-2.6939189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86354179) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(0.20881431) q[1];
sx q[1];
rz(-1.6046235) q[1];
sx q[1];
rz(-1.0135289) q[1];
rz(2.2603358) q[2];
sx q[2];
rz(-1.4980346) q[2];
sx q[2];
rz(-3.0537506) q[2];
rz(2.1590334) q[3];
sx q[3];
rz(-2.0053902) q[3];
sx q[3];
rz(-2.7872661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
