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
rz(1.9361629) q[0];
sx q[0];
rz(-0.054447629) q[0];
sx q[0];
rz(0.86583889) q[0];
rz(-1.539433) q[1];
sx q[1];
rz(-1.2442234) q[1];
sx q[1];
rz(-3.08334) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7476539) q[0];
sx q[0];
rz(-2.9705715) q[0];
sx q[0];
rz(-1.291541) q[0];
rz(-pi) q[1];
rz(-1.7069346) q[2];
sx q[2];
rz(-2.6420569) q[2];
sx q[2];
rz(-1.4972715) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.733898) q[1];
sx q[1];
rz(-1.0419894) q[1];
sx q[1];
rz(-0.66441345) q[1];
rz(-2.8179161) q[3];
sx q[3];
rz(-1.9012682) q[3];
sx q[3];
rz(0.068745384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3129348) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(-0.9642967) q[2];
rz(-0.079484552) q[3];
sx q[3];
rz(-1.3026404) q[3];
sx q[3];
rz(2.3704119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013414772) q[0];
sx q[0];
rz(-2.7981813) q[0];
sx q[0];
rz(1.6012023) q[0];
rz(-2.3004498) q[1];
sx q[1];
rz(-1.4079739) q[1];
sx q[1];
rz(-2.2841563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682789) q[0];
sx q[0];
rz(-0.43860093) q[0];
sx q[0];
rz(-1.5579786) q[0];
rz(-pi) q[1];
rz(-0.964999) q[2];
sx q[2];
rz(-0.73179663) q[2];
sx q[2];
rz(0.082187637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7341566) q[1];
sx q[1];
rz(-2.4333471) q[1];
sx q[1];
rz(-0.70869653) q[1];
rz(1.4748092) q[3];
sx q[3];
rz(-1.848683) q[3];
sx q[3];
rz(0.052770719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2555344) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(-1.9319755) q[2];
rz(1.3813193) q[3];
sx q[3];
rz(-1.8807024) q[3];
sx q[3];
rz(-2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.34476122) q[0];
sx q[0];
rz(-1.3153356) q[0];
sx q[0];
rz(-1.3321846) q[0];
rz(-1.4765129) q[1];
sx q[1];
rz(-1.8107199) q[1];
sx q[1];
rz(2.5217893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3416416) q[0];
sx q[0];
rz(-1.5551621) q[0];
sx q[0];
rz(-0.077341103) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99360433) q[2];
sx q[2];
rz(-0.30221488) q[2];
sx q[2];
rz(-0.72590088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42592749) q[1];
sx q[1];
rz(-1.8879379) q[1];
sx q[1];
rz(-1.4955382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2215106) q[3];
sx q[3];
rz(-0.6699282) q[3];
sx q[3];
rz(1.0494029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1406113) q[2];
sx q[2];
rz(-2.8162214) q[2];
sx q[2];
rz(0.78835431) q[2];
rz(-1.8654478) q[3];
sx q[3];
rz(-1.9143462) q[3];
sx q[3];
rz(-2.2787794) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9901504) q[0];
sx q[0];
rz(-1.7553512) q[0];
sx q[0];
rz(1.066712) q[0];
rz(1.3534631) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(1.1563168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3301834) q[0];
sx q[0];
rz(-1.5519987) q[0];
sx q[0];
rz(1.4292595) q[0];
x q[1];
rz(-1.5300445) q[2];
sx q[2];
rz(-1.5446071) q[2];
sx q[2];
rz(-1.9350236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9383685) q[1];
sx q[1];
rz(-2.6611106) q[1];
sx q[1];
rz(2.234747) q[1];
rz(-pi) q[2];
rz(1.636999) q[3];
sx q[3];
rz(-2.2341223) q[3];
sx q[3];
rz(-2.5732793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(3.0827674) q[2];
rz(-2.5022653) q[3];
sx q[3];
rz(-1.2213734) q[3];
sx q[3];
rz(0.85734573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095379742) q[0];
sx q[0];
rz(-1.7261427) q[0];
sx q[0];
rz(3.131102) q[0];
rz(1.563021) q[1];
sx q[1];
rz(-0.69067162) q[1];
sx q[1];
rz(2.6522327) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90411579) q[0];
sx q[0];
rz(-1.4129708) q[0];
sx q[0];
rz(-0.62741168) q[0];
rz(2.7362664) q[2];
sx q[2];
rz(-2.2203682) q[2];
sx q[2];
rz(-2.6774466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8243049) q[1];
sx q[1];
rz(-2.2154059) q[1];
sx q[1];
rz(2.0177474) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5602091) q[3];
sx q[3];
rz(-1.6015617) q[3];
sx q[3];
rz(-0.72584541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8420777) q[2];
sx q[2];
rz(-2.044951) q[2];
sx q[2];
rz(-0.00046029885) q[2];
rz(-0.67356235) q[3];
sx q[3];
rz(-0.898415) q[3];
sx q[3];
rz(2.4323997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8507268) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(-1.3036183) q[0];
rz(0.29779008) q[1];
sx q[1];
rz(-0.89556634) q[1];
sx q[1];
rz(-0.29388014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0654945) q[0];
sx q[0];
rz(-0.84653234) q[0];
sx q[0];
rz(-0.80418555) q[0];
rz(-pi) q[1];
rz(2.0877327) q[2];
sx q[2];
rz(-0.74956761) q[2];
sx q[2];
rz(-0.38656879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9781295) q[1];
sx q[1];
rz(-2.4678964) q[1];
sx q[1];
rz(-0.81836318) q[1];
rz(-pi) q[2];
rz(-2.5081283) q[3];
sx q[3];
rz(-0.38486275) q[3];
sx q[3];
rz(0.72609392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46875724) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(2.920816) q[2];
rz(-1.8686434) q[3];
sx q[3];
rz(-1.3947398) q[3];
sx q[3];
rz(-1.9934191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4918168) q[0];
sx q[0];
rz(-2.5204372) q[0];
sx q[0];
rz(-2.5671) q[0];
rz(-2.5993787) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(0.39074674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6764561) q[0];
sx q[0];
rz(-2.1628597) q[0];
sx q[0];
rz(-1.3964064) q[0];
rz(1.0263881) q[2];
sx q[2];
rz(-1.6154535) q[2];
sx q[2];
rz(0.80151973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.407651) q[1];
sx q[1];
rz(-1.9875257) q[1];
sx q[1];
rz(-1.5768361) q[1];
x q[2];
rz(3.0454759) q[3];
sx q[3];
rz(-0.36727723) q[3];
sx q[3];
rz(2.8678081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6812402) q[2];
sx q[2];
rz(-1.1972903) q[2];
sx q[2];
rz(2.528842) q[2];
rz(-0.98226205) q[3];
sx q[3];
rz(-1.8270315) q[3];
sx q[3];
rz(2.6764892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9923582) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(-3.0862578) q[0];
rz(1.5845567) q[1];
sx q[1];
rz(-1.0880071) q[1];
sx q[1];
rz(0.43992821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2429139) q[0];
sx q[0];
rz(-2.0283594) q[0];
sx q[0];
rz(-1.8458864) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7715334) q[2];
sx q[2];
rz(-2.6208813) q[2];
sx q[2];
rz(-1.0974018) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59938731) q[1];
sx q[1];
rz(-2.4894307) q[1];
sx q[1];
rz(-2.60531) q[1];
rz(-2.8596419) q[3];
sx q[3];
rz(-1.022517) q[3];
sx q[3];
rz(1.5632526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3736734) q[2];
sx q[2];
rz(-1.3481216) q[2];
sx q[2];
rz(-2.8295753) q[2];
rz(2.2219374) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31174082) q[0];
sx q[0];
rz(-2.1599202) q[0];
sx q[0];
rz(-3.1134636) q[0];
rz(-1.6955388) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(-0.45949724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7336373) q[0];
sx q[0];
rz(-1.5351632) q[0];
sx q[0];
rz(0.028546988) q[0];
x q[1];
rz(1.2120281) q[2];
sx q[2];
rz(-0.78132403) q[2];
sx q[2];
rz(1.5272905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7920835) q[1];
sx q[1];
rz(-2.4448423) q[1];
sx q[1];
rz(1.3443483) q[1];
rz(-0.13925456) q[3];
sx q[3];
rz(-1.9495268) q[3];
sx q[3];
rz(-1.1564573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0405937) q[2];
sx q[2];
rz(-1.7969635) q[2];
sx q[2];
rz(0.25516587) q[2];
rz(2.1549639) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(-1.423098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3808909) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(-1.3421407) q[0];
rz(0.0019207151) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(-1.049918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59219071) q[0];
sx q[0];
rz(-1.4225385) q[0];
sx q[0];
rz(1.3561983) q[0];
x q[1];
rz(-3.1403818) q[2];
sx q[2];
rz(-1.7467611) q[2];
sx q[2];
rz(-0.75226417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44651664) q[1];
sx q[1];
rz(-2.7088266) q[1];
sx q[1];
rz(2.5997735) q[1];
rz(-pi) q[2];
rz(-2.9258419) q[3];
sx q[3];
rz(-2.2892244) q[3];
sx q[3];
rz(1.6529447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92901015) q[2];
sx q[2];
rz(-1.2139823) q[2];
sx q[2];
rz(-2.6455961) q[2];
rz(-1.6811194) q[3];
sx q[3];
rz(-1.3417599) q[3];
sx q[3];
rz(-1.9546485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37954189) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(2.5024391) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(2.6668702) q[2];
sx q[2];
rz(-1.3480777) q[2];
sx q[2];
rz(2.011275) q[2];
rz(-2.4927327) q[3];
sx q[3];
rz(-0.93001233) q[3];
sx q[3];
rz(-0.42540023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
