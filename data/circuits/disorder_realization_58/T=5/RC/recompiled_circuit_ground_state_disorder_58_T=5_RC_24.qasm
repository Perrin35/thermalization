OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4492884) q[0];
sx q[0];
rz(3.6078499) q[0];
sx q[0];
rz(10.719263) q[0];
rz(2.8757088) q[1];
sx q[1];
rz(-0.20800132) q[1];
sx q[1];
rz(-1.2535569) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3304015) q[0];
sx q[0];
rz(-1.6881144) q[0];
sx q[0];
rz(-2.6091924) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1325344) q[2];
sx q[2];
rz(-1.9423368) q[2];
sx q[2];
rz(-1.8105992) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8575864) q[1];
sx q[1];
rz(-2.2145971) q[1];
sx q[1];
rz(0.30217742) q[1];
rz(-pi) q[2];
rz(-1.6276609) q[3];
sx q[3];
rz(-2.8167776) q[3];
sx q[3];
rz(1.1991833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0365389) q[2];
sx q[2];
rz(-1.7221071) q[2];
sx q[2];
rz(-1.4744021) q[2];
rz(-2.8640532) q[3];
sx q[3];
rz(-1.2801291) q[3];
sx q[3];
rz(0.92777073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9934746) q[0];
sx q[0];
rz(-2.9463705) q[0];
sx q[0];
rz(-1.3268205) q[0];
rz(0.38816342) q[1];
sx q[1];
rz(-1.5048985) q[1];
sx q[1];
rz(-0.51663748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.680707) q[0];
sx q[0];
rz(-0.26038489) q[0];
sx q[0];
rz(-2.5234114) q[0];
x q[1];
rz(0.20262589) q[2];
sx q[2];
rz(-1.1595402) q[2];
sx q[2];
rz(-0.52792462) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0320466) q[1];
sx q[1];
rz(-1.418772) q[1];
sx q[1];
rz(-0.8117453) q[1];
rz(-pi) q[2];
rz(-1.850892) q[3];
sx q[3];
rz(-1.3317654) q[3];
sx q[3];
rz(1.355442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8119729) q[2];
sx q[2];
rz(-2.6361578) q[2];
sx q[2];
rz(2.2774515) q[2];
rz(-2.4498074) q[3];
sx q[3];
rz(-1.1312048) q[3];
sx q[3];
rz(-2.2085021) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8253887) q[0];
sx q[0];
rz(-2.3360257) q[0];
sx q[0];
rz(1.7387996) q[0];
rz(0.56651506) q[1];
sx q[1];
rz(-1.6367876) q[1];
sx q[1];
rz(-1.8623955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9815807) q[0];
sx q[0];
rz(-1.564508) q[0];
sx q[0];
rz(-1.8304197) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98033078) q[2];
sx q[2];
rz(-1.1644018) q[2];
sx q[2];
rz(-0.097336285) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0356939) q[1];
sx q[1];
rz(-2.0146684) q[1];
sx q[1];
rz(-3.1015293) q[1];
rz(-pi) q[2];
rz(0.62554977) q[3];
sx q[3];
rz(-1.2569142) q[3];
sx q[3];
rz(-0.63077273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0684356) q[2];
sx q[2];
rz(-0.34978875) q[2];
sx q[2];
rz(-1.5901828) q[2];
rz(2.6326211) q[3];
sx q[3];
rz(-0.83566982) q[3];
sx q[3];
rz(-1.6905748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029025404) q[0];
sx q[0];
rz(-1.4968766) q[0];
sx q[0];
rz(-2.577884) q[0];
rz(1.6911223) q[1];
sx q[1];
rz(-1.1553973) q[1];
sx q[1];
rz(-2.3203826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9552069) q[0];
sx q[0];
rz(-1.5816551) q[0];
sx q[0];
rz(-0.71594724) q[0];
x q[1];
rz(2.8217373) q[2];
sx q[2];
rz(-1.2884239) q[2];
sx q[2];
rz(-0.47114633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6159284) q[1];
sx q[1];
rz(-2.2771336) q[1];
sx q[1];
rz(0.86752992) q[1];
x q[2];
rz(2.1139007) q[3];
sx q[3];
rz(-1.538435) q[3];
sx q[3];
rz(-1.1278271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44068367) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(-1.7636501) q[2];
rz(-0.16436973) q[3];
sx q[3];
rz(-1.5956968) q[3];
sx q[3];
rz(-0.36851287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.565777) q[0];
sx q[0];
rz(-1.4482647) q[0];
sx q[0];
rz(2.6337295) q[0];
rz(-2.2841618) q[1];
sx q[1];
rz(-1.6897759) q[1];
sx q[1];
rz(2.6618777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21151152) q[0];
sx q[0];
rz(-0.9418315) q[0];
sx q[0];
rz(-1.1674561) q[0];
x q[1];
rz(1.6417333) q[2];
sx q[2];
rz(-1.1878769) q[2];
sx q[2];
rz(-1.4986582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3277726) q[1];
sx q[1];
rz(-1.7952807) q[1];
sx q[1];
rz(-0.95275615) q[1];
rz(1.2913088) q[3];
sx q[3];
rz(-1.5429792) q[3];
sx q[3];
rz(1.3007377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63626426) q[2];
sx q[2];
rz(-2.0365066) q[2];
sx q[2];
rz(0.21978933) q[2];
rz(-0.2291186) q[3];
sx q[3];
rz(-2.2380405) q[3];
sx q[3];
rz(-0.98178274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.583928) q[0];
sx q[0];
rz(-2.5983577) q[0];
sx q[0];
rz(1.1274717) q[0];
rz(2.8358031) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(0.19010273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0738433) q[0];
sx q[0];
rz(-2.0408568) q[0];
sx q[0];
rz(-2.7795634) q[0];
x q[1];
rz(0.79391329) q[2];
sx q[2];
rz(-1.8474692) q[2];
sx q[2];
rz(0.45396462) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2711004) q[1];
sx q[1];
rz(-1.2855347) q[1];
sx q[1];
rz(2.3149957) q[1];
x q[2];
rz(-0.25214002) q[3];
sx q[3];
rz(-1.6448955) q[3];
sx q[3];
rz(-2.4197297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.019913435) q[2];
sx q[2];
rz(-1.4608258) q[2];
sx q[2];
rz(-1.100568) q[2];
rz(1.3044283) q[3];
sx q[3];
rz(-1.5521939) q[3];
sx q[3];
rz(-1.3531551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562427) q[0];
sx q[0];
rz(-0.42735639) q[0];
sx q[0];
rz(-0.57619488) q[0];
rz(2.088749) q[1];
sx q[1];
rz(-2.5710227) q[1];
sx q[1];
rz(1.0594692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5813025) q[0];
sx q[0];
rz(-1.6972739) q[0];
sx q[0];
rz(2.4471388) q[0];
rz(-pi) q[1];
rz(1.1722819) q[2];
sx q[2];
rz(-1.5969807) q[2];
sx q[2];
rz(-2.6790533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6061343) q[1];
sx q[1];
rz(-0.48608695) q[1];
sx q[1];
rz(-0.41907678) q[1];
rz(0.93848159) q[3];
sx q[3];
rz(-0.96358591) q[3];
sx q[3];
rz(-0.069165088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8656371) q[2];
sx q[2];
rz(-0.58791462) q[2];
sx q[2];
rz(2.0666583) q[2];
rz(-2.2771207) q[3];
sx q[3];
rz(-1.8755707) q[3];
sx q[3];
rz(-1.5629432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7445755) q[0];
sx q[0];
rz(-2.1730142) q[0];
sx q[0];
rz(0.50719914) q[0];
rz(1.1811258) q[1];
sx q[1];
rz(-1.6543417) q[1];
sx q[1];
rz(2.2307253) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9932803) q[0];
sx q[0];
rz(-2.1851843) q[0];
sx q[0];
rz(-2.5117158) q[0];
rz(-pi) q[1];
rz(-1.2664001) q[2];
sx q[2];
rz(-0.60303973) q[2];
sx q[2];
rz(1.3040486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0720328) q[1];
sx q[1];
rz(-1.7246453) q[1];
sx q[1];
rz(0.01158684) q[1];
rz(2.0924545) q[3];
sx q[3];
rz(-2.6781265) q[3];
sx q[3];
rz(-1.3090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7562423) q[2];
sx q[2];
rz(-1.664398) q[2];
sx q[2];
rz(-0.64794668) q[2];
rz(0.09662763) q[3];
sx q[3];
rz(-2.1164618) q[3];
sx q[3];
rz(-2.1881762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18428093) q[0];
sx q[0];
rz(-2.1537557) q[0];
sx q[0];
rz(-1.8010358) q[0];
rz(2.6181009) q[1];
sx q[1];
rz(-1.5308056) q[1];
sx q[1];
rz(1.9370105) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8891994) q[0];
sx q[0];
rz(-2.3601818) q[0];
sx q[0];
rz(-1.8535421) q[0];
rz(-2.9926489) q[2];
sx q[2];
rz(-1.7686378) q[2];
sx q[2];
rz(2.1959675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5960418) q[1];
sx q[1];
rz(-1.7298455) q[1];
sx q[1];
rz(0.56175128) q[1];
rz(-2.7369237) q[3];
sx q[3];
rz(-2.0467374) q[3];
sx q[3];
rz(2.1309545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1022819) q[2];
sx q[2];
rz(-2.1151586) q[2];
sx q[2];
rz(0.2956051) q[2];
rz(-0.11232703) q[3];
sx q[3];
rz(-1.5887518) q[3];
sx q[3];
rz(0.86265341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.2304147) q[0];
sx q[0];
rz(-0.4998315) q[0];
sx q[0];
rz(-2.3590132) q[0];
rz(-2.8453907) q[1];
sx q[1];
rz(-1.9007416) q[1];
sx q[1];
rz(0.20015073) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3302932) q[0];
sx q[0];
rz(-1.7331373) q[0];
sx q[0];
rz(-0.77438023) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1364983) q[2];
sx q[2];
rz(-1.3315233) q[2];
sx q[2];
rz(2.4636961) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1910341) q[1];
sx q[1];
rz(-2.28015) q[1];
sx q[1];
rz(2.2851431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93882046) q[3];
sx q[3];
rz(-2.0858313) q[3];
sx q[3];
rz(2.7523628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.079411) q[2];
sx q[2];
rz(-0.14180413) q[2];
sx q[2];
rz(1.0395435) q[2];
rz(0.027035106) q[3];
sx q[3];
rz(-2.1578433) q[3];
sx q[3];
rz(2.1496617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1886002) q[0];
sx q[0];
rz(-2.47692) q[0];
sx q[0];
rz(1.1823786) q[0];
rz(-1.7091119) q[1];
sx q[1];
rz(-1.2624337) q[1];
sx q[1];
rz(1.591325) q[1];
rz(2.9903632) q[2];
sx q[2];
rz(-1.4561903) q[2];
sx q[2];
rz(-1.0753808) q[2];
rz(0.19743528) q[3];
sx q[3];
rz(-2.6876269) q[3];
sx q[3];
rz(-0.85314565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
