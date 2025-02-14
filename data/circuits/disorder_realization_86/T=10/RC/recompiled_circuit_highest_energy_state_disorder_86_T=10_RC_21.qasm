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
rz(1.9011693) q[0];
sx q[0];
rz(-1.9440396) q[0];
sx q[0];
rz(-2.9252606) q[0];
rz(-2.1842015) q[1];
sx q[1];
rz(-0.67617813) q[1];
sx q[1];
rz(1.6300936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0974329) q[0];
sx q[0];
rz(-0.55632797) q[0];
sx q[0];
rz(3.1233643) q[0];
rz(1.2443107) q[2];
sx q[2];
rz(-1.1190345) q[2];
sx q[2];
rz(-1.2690074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2547467) q[1];
sx q[1];
rz(-2.0664082) q[1];
sx q[1];
rz(-1.9504471) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39325289) q[3];
sx q[3];
rz(-2.234708) q[3];
sx q[3];
rz(-2.8877986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.034885255) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(1.0484288) q[2];
rz(-0.18167051) q[3];
sx q[3];
rz(-2.1766267) q[3];
sx q[3];
rz(0.65339965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6492017) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(-0.44678584) q[0];
rz(-2.2611639) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(-0.78278881) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7715817) q[0];
sx q[0];
rz(-0.25213045) q[0];
sx q[0];
rz(-1.3229516) q[0];
x q[1];
rz(-1.9982463) q[2];
sx q[2];
rz(-0.99787092) q[2];
sx q[2];
rz(-1.3818416) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0228717) q[1];
sx q[1];
rz(-2.0346271) q[1];
sx q[1];
rz(1.2341586) q[1];
x q[2];
rz(1.0083593) q[3];
sx q[3];
rz(-2.376412) q[3];
sx q[3];
rz(-0.4972813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.030152628) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(0.97935575) q[2];
rz(-2.7566946) q[3];
sx q[3];
rz(-1.220547) q[3];
sx q[3];
rz(-0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.6556743) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(-0.78980494) q[0];
rz(0.18665953) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(2.1479215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6826166) q[0];
sx q[0];
rz(-2.1135215) q[0];
sx q[0];
rz(1.2786464) q[0];
x q[1];
rz(2.5285401) q[2];
sx q[2];
rz(-1.384735) q[2];
sx q[2];
rz(-1.2836054) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6216065) q[1];
sx q[1];
rz(-2.4944759) q[1];
sx q[1];
rz(-0.81955975) q[1];
rz(-1.2993811) q[3];
sx q[3];
rz(-0.86835734) q[3];
sx q[3];
rz(-2.0347119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32187244) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(-0.37128386) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(-2.546052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6868941) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(1.7373079) q[0];
rz(-0.67717254) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(1.6489395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4801475) q[0];
sx q[0];
rz(-2.8146525) q[0];
sx q[0];
rz(-0.72823712) q[0];
x q[1];
rz(2.0022961) q[2];
sx q[2];
rz(-1.8427263) q[2];
sx q[2];
rz(1.6587342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3918152) q[1];
sx q[1];
rz(-1.3484553) q[1];
sx q[1];
rz(2.8676377) q[1];
rz(-pi) q[2];
rz(-1.221232) q[3];
sx q[3];
rz(-2.3858983) q[3];
sx q[3];
rz(0.44103482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45381418) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(2.7933534) q[2];
rz(-1.6866775) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763879) q[0];
sx q[0];
rz(-0.56469733) q[0];
sx q[0];
rz(-2.2494466) q[0];
rz(0.46420321) q[1];
sx q[1];
rz(-1.8966388) q[1];
sx q[1];
rz(1.2947327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97335739) q[0];
sx q[0];
rz(-2.0414314) q[0];
sx q[0];
rz(-0.55460978) q[0];
rz(-pi) q[1];
rz(0.21315766) q[2];
sx q[2];
rz(-2.5937754) q[2];
sx q[2];
rz(2.167706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.059493493) q[1];
sx q[1];
rz(-2.0096742) q[1];
sx q[1];
rz(-2.782269) q[1];
rz(2.5312349) q[3];
sx q[3];
rz(-2.1564283) q[3];
sx q[3];
rz(-0.40025362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(-0.56582212) q[2];
rz(-3.0602509) q[3];
sx q[3];
rz(-0.9809202) q[3];
sx q[3];
rz(0.62059039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.464798) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(1.5555405) q[0];
rz(2.0809035) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(-2.5097844) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2151129) q[0];
sx q[0];
rz(-0.2479015) q[0];
sx q[0];
rz(0.32899022) q[0];
rz(-2.6423694) q[2];
sx q[2];
rz(-2.2167335) q[2];
sx q[2];
rz(1.8725841) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28367701) q[1];
sx q[1];
rz(-1.2322958) q[1];
sx q[1];
rz(-1.0033016) q[1];
x q[2];
rz(1.7177714) q[3];
sx q[3];
rz(-0.89255652) q[3];
sx q[3];
rz(0.036546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98171988) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(2.6427606) q[2];
rz(1.8286797) q[3];
sx q[3];
rz(-0.73176089) q[3];
sx q[3];
rz(-1.4341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53443921) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(2.4420807) q[0];
rz(-0.46547678) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(-2.2043998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995183) q[0];
sx q[0];
rz(-2.8694911) q[0];
sx q[0];
rz(2.4689552) q[0];
rz(-pi) q[1];
rz(1.6752536) q[2];
sx q[2];
rz(-1.3323931) q[2];
sx q[2];
rz(1.4739715) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4027462) q[1];
sx q[1];
rz(-0.35772309) q[1];
sx q[1];
rz(-2.2132232) q[1];
rz(-1.2770416) q[3];
sx q[3];
rz(-1.86092) q[3];
sx q[3];
rz(-1.4690831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84404868) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(-2.76827) q[2];
rz(1.1897872) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3968286) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(-2.3936791) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(-3.1386197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82389861) q[0];
sx q[0];
rz(-0.71461535) q[0];
sx q[0];
rz(1.8665642) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6417129) q[2];
sx q[2];
rz(-1.5902963) q[2];
sx q[2];
rz(1.8766581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21196762) q[1];
sx q[1];
rz(-0.16008437) q[1];
sx q[1];
rz(-1.2950241) q[1];
rz(-pi) q[2];
rz(0.19488867) q[3];
sx q[3];
rz(-1.7656544) q[3];
sx q[3];
rz(-3.0611567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.030674) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(-2.0745011) q[2];
rz(-0.083960697) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(-2.3626204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.344051) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(3.0352266) q[0];
rz(2.1616409) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(-0.76593691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.879012) q[0];
sx q[0];
rz(-1.1404788) q[0];
sx q[0];
rz(2.5613214) q[0];
rz(1.6607051) q[2];
sx q[2];
rz(-1.769763) q[2];
sx q[2];
rz(-1.0294017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3801743) q[1];
sx q[1];
rz(-1.6055709) q[1];
sx q[1];
rz(-2.4654287) q[1];
x q[2];
rz(3.1274904) q[3];
sx q[3];
rz(-1.5670793) q[3];
sx q[3];
rz(2.2110018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47436675) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-2.4412947) q[2];
rz(-2.4750366) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(-1.0880281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661082) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(-1.3960557) q[0];
rz(-1.667977) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(-2.4748763) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7233492) q[0];
sx q[0];
rz(-2.1921792) q[0];
sx q[0];
rz(-0.80051144) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2856977) q[2];
sx q[2];
rz(-0.93109967) q[2];
sx q[2];
rz(0.49298795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4426431) q[1];
sx q[1];
rz(-1.8869699) q[1];
sx q[1];
rz(1.433028) q[1];
x q[2];
rz(-0.57228831) q[3];
sx q[3];
rz(-2.0438571) q[3];
sx q[3];
rz(-2.8561887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0746158) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(0.89078772) q[2];
rz(-2.7095419) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4074832) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(1.2474077) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(2.9802889) q[2];
sx q[2];
rz(-1.2304753) q[2];
sx q[2];
rz(-0.58438042) q[2];
rz(3.0790764) q[3];
sx q[3];
rz(-2.6517031) q[3];
sx q[3];
rz(-0.22417886) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
