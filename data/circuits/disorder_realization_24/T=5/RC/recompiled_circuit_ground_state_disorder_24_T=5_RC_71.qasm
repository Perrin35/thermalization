OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.71819031) q[0];
sx q[0];
rz(3.1443449) q[0];
sx q[0];
rz(10.436463) q[0];
rz(0.45971316) q[1];
sx q[1];
rz(-1.8399532) q[1];
sx q[1];
rz(-0.17170061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7778439) q[0];
sx q[0];
rz(-2.1743589) q[0];
sx q[0];
rz(-1.4715172) q[0];
rz(-pi) q[1];
rz(0.19440513) q[2];
sx q[2];
rz(-1.462359) q[2];
sx q[2];
rz(-0.51314236) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57331177) q[1];
sx q[1];
rz(-2.3900095) q[1];
sx q[1];
rz(2.4922396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0954451) q[3];
sx q[3];
rz(-0.10766115) q[3];
sx q[3];
rz(0.3357418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1823938) q[2];
sx q[2];
rz(-2.6540519) q[2];
sx q[2];
rz(1.7215151) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(1.0678631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(3.0576393) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(-0.49348304) q[0];
rz(-2.3656942) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(0.50481558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579266) q[0];
sx q[0];
rz(-2.2424881) q[0];
sx q[0];
rz(0.85606411) q[0];
x q[1];
rz(-2.6146982) q[2];
sx q[2];
rz(-1.3119446) q[2];
sx q[2];
rz(-0.35824725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2097834) q[1];
sx q[1];
rz(-1.9701013) q[1];
sx q[1];
rz(-2.4834391) q[1];
x q[2];
rz(0.22498954) q[3];
sx q[3];
rz(-2.8691022) q[3];
sx q[3];
rz(-0.30411965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8043171) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(-2.6709225) q[2];
rz(-0.63052952) q[3];
sx q[3];
rz(-2.1157406) q[3];
sx q[3];
rz(1.7414909) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0640963) q[0];
sx q[0];
rz(-0.12501669) q[0];
sx q[0];
rz(0.65573829) q[0];
rz(2.8302622) q[1];
sx q[1];
rz(-2.2475524) q[1];
sx q[1];
rz(0.071455926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27256672) q[0];
sx q[0];
rz(-1.9211968) q[0];
sx q[0];
rz(-0.053215543) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33768968) q[2];
sx q[2];
rz(-1.2215541) q[2];
sx q[2];
rz(1.2472635) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.25117043) q[1];
sx q[1];
rz(-0.4931207) q[1];
sx q[1];
rz(-1.7095196) q[1];
rz(-0.35122613) q[3];
sx q[3];
rz(-0.5595419) q[3];
sx q[3];
rz(-0.21332394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91325703) q[2];
sx q[2];
rz(-1.5568638) q[2];
sx q[2];
rz(0.060001686) q[2];
rz(1.9289198) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(2.5643964) q[0];
rz(0.82950854) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(-2.3037691) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55623193) q[0];
sx q[0];
rz(-2.5149226) q[0];
sx q[0];
rz(0.13154948) q[0];
rz(-pi) q[1];
rz(0.060537593) q[2];
sx q[2];
rz(-2.1601387) q[2];
sx q[2];
rz(-2.4174487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29463331) q[1];
sx q[1];
rz(-0.75645743) q[1];
sx q[1];
rz(1.0119757) q[1];
rz(-2.6959723) q[3];
sx q[3];
rz(-2.1975448) q[3];
sx q[3];
rz(2.8901951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1159346) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(-1.9177829) q[2];
rz(-2.6754248) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(2.536072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25662988) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(-0.10109854) q[0];
rz(-2.7283607) q[1];
sx q[1];
rz(-0.44094545) q[1];
sx q[1];
rz(2.0535927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013070949) q[0];
sx q[0];
rz(-2.5422342) q[0];
sx q[0];
rz(0.96512633) q[0];
x q[1];
rz(2.3067683) q[2];
sx q[2];
rz(-1.2918279) q[2];
sx q[2];
rz(-2.5503412) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1428296) q[1];
sx q[1];
rz(-1.3266449) q[1];
sx q[1];
rz(-1.2296089) q[1];
rz(-pi) q[2];
rz(-0.66941485) q[3];
sx q[3];
rz(-1.6261887) q[3];
sx q[3];
rz(0.15633571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92855144) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(1.5555596) q[2];
rz(-1.0726311) q[3];
sx q[3];
rz(-1.6173247) q[3];
sx q[3];
rz(-0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17305408) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(-1.435085) q[0];
rz(-0.75812078) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(3.039956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1124737) q[0];
sx q[0];
rz(-2.2879507) q[0];
sx q[0];
rz(2.4957335) q[0];
x q[1];
rz(-0.36113895) q[2];
sx q[2];
rz(-1.1419347) q[2];
sx q[2];
rz(-2.7782235) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1730301) q[1];
sx q[1];
rz(-0.30245879) q[1];
sx q[1];
rz(0.36424251) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6061344) q[3];
sx q[3];
rz(-2.1369918) q[3];
sx q[3];
rz(-2.9009852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8018084) q[2];
sx q[2];
rz(-1.1016176) q[2];
sx q[2];
rz(-1.3327117) q[2];
rz(1.0821651) q[3];
sx q[3];
rz(-0.75282955) q[3];
sx q[3];
rz(1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69671714) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(-2.3984997) q[0];
rz(0.7630868) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(-0.9185763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6800623) q[0];
sx q[0];
rz(-1.5627075) q[0];
sx q[0];
rz(-1.5957521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0985561) q[2];
sx q[2];
rz(-1.8794606) q[2];
sx q[2];
rz(-0.93738467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1700701) q[1];
sx q[1];
rz(-2.5405209) q[1];
sx q[1];
rz(-2.5843589) q[1];
x q[2];
rz(-0.72414805) q[3];
sx q[3];
rz(-1.482943) q[3];
sx q[3];
rz(-1.7334565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7243728) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(0.43583885) q[2];
rz(-3.0271652) q[3];
sx q[3];
rz(-0.2468214) q[3];
sx q[3];
rz(0.59823263) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33335394) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(2.5575141) q[0];
rz(-0.43560478) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(1.5199419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1849862) q[0];
sx q[0];
rz(-1.2756366) q[0];
sx q[0];
rz(1.0888238) q[0];
x q[1];
rz(0.21824117) q[2];
sx q[2];
rz(-1.9710961) q[2];
sx q[2];
rz(-1.2117653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4055109) q[1];
sx q[1];
rz(-2.9353432) q[1];
sx q[1];
rz(-3.095605) q[1];
rz(0.55840839) q[3];
sx q[3];
rz(-0.30277751) q[3];
sx q[3];
rz(-1.6807792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0674151) q[2];
sx q[2];
rz(-1.3520853) q[2];
sx q[2];
rz(1.1161067) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-2.0383056) q[3];
sx q[3];
rz(-0.11837676) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(2.351601) q[0];
rz(3.0780011) q[1];
sx q[1];
rz(-2.5345232) q[1];
sx q[1];
rz(-2.7008609) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9960956) q[0];
sx q[0];
rz(-2.7005368) q[0];
sx q[0];
rz(-0.8121374) q[0];
rz(-1.9750544) q[2];
sx q[2];
rz(-1.1102997) q[2];
sx q[2];
rz(1.782589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1524618) q[1];
sx q[1];
rz(-0.50639443) q[1];
sx q[1];
rz(-1.9146862) q[1];
x q[2];
rz(-3.1001631) q[3];
sx q[3];
rz(-1.3815615) q[3];
sx q[3];
rz(1.3845058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8792087) q[2];
sx q[2];
rz(-2.3449506) q[2];
sx q[2];
rz(-0.66514307) q[2];
rz(-2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(-1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.7803698) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(2.1369456) q[0];
rz(1.4077582) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(-1.2900603) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697851) q[0];
sx q[0];
rz(-1.3243073) q[0];
sx q[0];
rz(2.3355961) q[0];
rz(-pi) q[1];
rz(2.9384841) q[2];
sx q[2];
rz(-0.98047335) q[2];
sx q[2];
rz(3.1121569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5651144) q[1];
sx q[1];
rz(-1.8449515) q[1];
sx q[1];
rz(1.3512011) q[1];
rz(-2.9674888) q[3];
sx q[3];
rz(-0.977036) q[3];
sx q[3];
rz(1.6792149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1349692) q[2];
sx q[2];
rz(-1.1639736) q[2];
sx q[2];
rz(0.25035614) q[2];
rz(2.1641459) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(1.6710963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(-2.3969338) q[1];
sx q[1];
rz(-0.79569334) q[1];
sx q[1];
rz(-2.4462499) q[1];
rz(2.9982243) q[2];
sx q[2];
rz(-1.5833387) q[2];
sx q[2];
rz(0.34064731) q[2];
rz(0.86359371) q[3];
sx q[3];
rz(-0.61932388) q[3];
sx q[3];
rz(-1.2854734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
