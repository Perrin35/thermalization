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
rz(1.5155091) q[0];
sx q[0];
rz(3.4184366) q[0];
sx q[0];
rz(9.3079216) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(-0.35419551) q[1];
sx q[1];
rz(0.49122223) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15165937) q[0];
sx q[0];
rz(-1.0826246) q[0];
sx q[0];
rz(-3.0753646) q[0];
rz(-pi) q[1];
rz(-2.6543219) q[2];
sx q[2];
rz(-1.7585764) q[2];
sx q[2];
rz(-1.0010546) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3957805) q[1];
sx q[1];
rz(-2.2015101) q[1];
sx q[1];
rz(1.4480026) q[1];
rz(3.1163636) q[3];
sx q[3];
rz(-1.2275891) q[3];
sx q[3];
rz(1.1017088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7626071) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(-2.911705) q[2];
rz(-3.0843132) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(0.91276401) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463992) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-0.65509534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0948308) q[0];
sx q[0];
rz(-0.74709409) q[0];
sx q[0];
rz(-1.2741035) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35327618) q[2];
sx q[2];
rz(-0.52918079) q[2];
sx q[2];
rz(-2.8904817) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8959143) q[1];
sx q[1];
rz(-1.0686456) q[1];
sx q[1];
rz(-2.4163626) q[1];
x q[2];
rz(2.9981587) q[3];
sx q[3];
rz(-2.3000882) q[3];
sx q[3];
rz(-2.7021331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(1.8435271) q[2];
rz(1.2201747) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(-0.65742457) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87797457) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(-2.5275912) q[0];
rz(1.7738495) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(0.49555379) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14392631) q[0];
sx q[0];
rz(-0.94519061) q[0];
sx q[0];
rz(0.99220522) q[0];
rz(-pi) q[1];
x q[1];
rz(2.173893) q[2];
sx q[2];
rz(-0.7843547) q[2];
sx q[2];
rz(1.7627782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6537955) q[1];
sx q[1];
rz(-1.3843517) q[1];
sx q[1];
rz(-1.3273147) q[1];
x q[2];
rz(-1.9248149) q[3];
sx q[3];
rz(-1.1292158) q[3];
sx q[3];
rz(-1.3819249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76501781) q[2];
sx q[2];
rz(-1.0708662) q[2];
sx q[2];
rz(0.59256727) q[2];
rz(-2.9060034) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4950824) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(2.5616772) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(0.52111202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082588293) q[0];
sx q[0];
rz(-0.0085235611) q[0];
sx q[0];
rz(-0.99446358) q[0];
rz(-pi) q[1];
rz(0.02946202) q[2];
sx q[2];
rz(-0.97004393) q[2];
sx q[2];
rz(3.1371869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5856984) q[1];
sx q[1];
rz(-1.1194849) q[1];
sx q[1];
rz(-2.0183619) q[1];
x q[2];
rz(-0.24346015) q[3];
sx q[3];
rz(-2.820558) q[3];
sx q[3];
rz(2.3784172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(-2.8233675) q[2];
rz(0.71277726) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37079674) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(-3.1053012) q[0];
rz(-2.6151784) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(-2.859419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68809915) q[0];
sx q[0];
rz(-2.4465313) q[0];
sx q[0];
rz(-1.1316677) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87558117) q[2];
sx q[2];
rz(-1.1291847) q[2];
sx q[2];
rz(1.3087866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6101004) q[1];
sx q[1];
rz(-2.2740915) q[1];
sx q[1];
rz(-1.6628357) q[1];
rz(-2.3822278) q[3];
sx q[3];
rz(-2.1753484) q[3];
sx q[3];
rz(-2.5370363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37461773) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(3.0275893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58979708) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.9675323) q[0];
rz(-0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(3.1016268) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55520457) q[0];
sx q[0];
rz(-0.38643906) q[0];
sx q[0];
rz(-2.1512335) q[0];
rz(-pi) q[1];
rz(2.8103175) q[2];
sx q[2];
rz(-0.9563947) q[2];
sx q[2];
rz(-1.5671052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3921622) q[1];
sx q[1];
rz(-2.431173) q[1];
sx q[1];
rz(2.602052) q[1];
rz(-2.9257366) q[3];
sx q[3];
rz(-0.83789589) q[3];
sx q[3];
rz(-2.6449316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93969718) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(-0.8197909) q[2];
rz(2.36006) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(-1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49067295) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(-1.8248722) q[0];
rz(-1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(2.7046611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1986982) q[0];
sx q[0];
rz(-1.6075396) q[0];
sx q[0];
rz(-1.7612639) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4605791) q[2];
sx q[2];
rz(-1.4535365) q[2];
sx q[2];
rz(2.3667468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93561427) q[1];
sx q[1];
rz(-2.2956478) q[1];
sx q[1];
rz(1.0826675) q[1];
rz(-pi) q[2];
rz(2.1070621) q[3];
sx q[3];
rz(-1.9428184) q[3];
sx q[3];
rz(-2.0780502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74066073) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(1.9756165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400035) q[0];
sx q[0];
rz(-2.8497301) q[0];
sx q[0];
rz(0.18184161) q[0];
rz(2.8002411) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(-0.15886074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10578622) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(0.44802702) q[0];
x q[1];
rz(-2.7188825) q[2];
sx q[2];
rz(-1.8313932) q[2];
sx q[2];
rz(-0.24525951) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5103858) q[1];
sx q[1];
rz(-2.800731) q[1];
sx q[1];
rz(-2.4426961) q[1];
rz(2.587593) q[3];
sx q[3];
rz(-1.4408198) q[3];
sx q[3];
rz(-2.8838571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86026704) q[2];
sx q[2];
rz(-0.93026668) q[2];
sx q[2];
rz(2.2567828) q[2];
rz(-1.7224711) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(-0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849843) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(0.9935317) q[0];
rz(-1.7811071) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(2.5885168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4489514) q[0];
sx q[0];
rz(-1.7990489) q[0];
sx q[0];
rz(-2.7518163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.07358609) q[2];
sx q[2];
rz(-2.4017757) q[2];
sx q[2];
rz(0.70895665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7874541) q[1];
sx q[1];
rz(-2.1864357) q[1];
sx q[1];
rz(1.8880185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2129799) q[3];
sx q[3];
rz(-0.67104895) q[3];
sx q[3];
rz(-2.8007361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-0.23703144) q[2];
rz(1.6009181) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700664) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-0.1189098) q[0];
rz(0.43926829) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-0.063442245) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3812696) q[0];
sx q[0];
rz(-2.1043692) q[0];
sx q[0];
rz(-2.3593581) q[0];
rz(-pi) q[1];
x q[1];
rz(1.756606) q[2];
sx q[2];
rz(-2.1868621) q[2];
sx q[2];
rz(2.1587929) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0902151) q[1];
sx q[1];
rz(-1.8974638) q[1];
sx q[1];
rz(-0.0089545687) q[1];
rz(-pi) q[2];
rz(-0.88313734) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(0.36206743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0110317) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(-1.6434742) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.3907923) q[0];
sx q[0];
rz(-1.1815429) q[0];
sx q[0];
rz(-1.0829157) q[0];
rz(-0.87286585) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(-0.68442576) q[2];
sx q[2];
rz(-1.050907) q[2];
sx q[2];
rz(2.4466865) q[2];
rz(-0.11029966) q[3];
sx q[3];
rz(-2.5613789) q[3];
sx q[3];
rz(-1.2753244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
