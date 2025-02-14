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
rz(-1.7492548) q[0];
sx q[0];
rz(-0.37819401) q[0];
sx q[0];
rz(0.35051546) q[0];
rz(-2.2502083) q[1];
sx q[1];
rz(-0.79217029) q[1];
sx q[1];
rz(1.9686735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25305155) q[0];
sx q[0];
rz(-1.3810147) q[0];
sx q[0];
rz(-1.3773328) q[0];
rz(-pi) q[1];
rz(-2.1831398) q[2];
sx q[2];
rz(-0.26587379) q[2];
sx q[2];
rz(2.079351) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4306372) q[1];
sx q[1];
rz(-1.4078377) q[1];
sx q[1];
rz(2.8204172) q[1];
x q[2];
rz(0.32781847) q[3];
sx q[3];
rz(-0.49428764) q[3];
sx q[3];
rz(-1.0456628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9709836) q[2];
sx q[2];
rz(-1.8656518) q[2];
sx q[2];
rz(-0.21273908) q[2];
rz(3.0563266) q[3];
sx q[3];
rz(-1.8906967) q[3];
sx q[3];
rz(-1.2712449) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3474715) q[0];
sx q[0];
rz(-0.81231064) q[0];
sx q[0];
rz(1.043327) q[0];
rz(3.0087545) q[1];
sx q[1];
rz(-0.21505198) q[1];
sx q[1];
rz(-1.9827693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67135982) q[0];
sx q[0];
rz(-1.7152107) q[0];
sx q[0];
rz(0.40035046) q[0];
x q[1];
rz(0.01387502) q[2];
sx q[2];
rz(-1.4094947) q[2];
sx q[2];
rz(-0.32729766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4606884) q[1];
sx q[1];
rz(-1.0863606) q[1];
sx q[1];
rz(1.2699732) q[1];
rz(-pi) q[2];
rz(-1.1400677) q[3];
sx q[3];
rz(-1.8716806) q[3];
sx q[3];
rz(1.4804076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1966689) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(-2.4169253) q[2];
rz(0.66257462) q[3];
sx q[3];
rz(-0.96810883) q[3];
sx q[3];
rz(2.841943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1835943) q[0];
sx q[0];
rz(-1.2677001) q[0];
sx q[0];
rz(-0.9147574) q[0];
rz(2.5547408) q[1];
sx q[1];
rz(-1.7210759) q[1];
sx q[1];
rz(0.14952001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5356333) q[0];
sx q[0];
rz(-2.0894755) q[0];
sx q[0];
rz(-2.9796322) q[0];
rz(-pi) q[1];
rz(1.9376041) q[2];
sx q[2];
rz(-1.2292394) q[2];
sx q[2];
rz(1.2836518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5897474) q[1];
sx q[1];
rz(-2.0974297) q[1];
sx q[1];
rz(2.2174775) q[1];
x q[2];
rz(-0.87823509) q[3];
sx q[3];
rz(-2.1562169) q[3];
sx q[3];
rz(-2.6458322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2933423) q[2];
sx q[2];
rz(-1.9363656) q[2];
sx q[2];
rz(-2.3411574) q[2];
rz(-0.46659255) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(0.65565562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7993497) q[0];
sx q[0];
rz(-1.8514587) q[0];
sx q[0];
rz(-2.2829862) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-2.938439) q[1];
sx q[1];
rz(2.1536749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31269095) q[0];
sx q[0];
rz(-0.36417555) q[0];
sx q[0];
rz(2.0065424) q[0];
rz(-pi) q[1];
rz(1.3484389) q[2];
sx q[2];
rz(-0.50560617) q[2];
sx q[2];
rz(0.33481516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6776442) q[1];
sx q[1];
rz(-1.8419187) q[1];
sx q[1];
rz(1.8810924) q[1];
rz(2.6060054) q[3];
sx q[3];
rz(-0.37112826) q[3];
sx q[3];
rz(2.8674026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19114384) q[2];
sx q[2];
rz(-1.3942275) q[2];
sx q[2];
rz(-2.6960755) q[2];
rz(2.4281003) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(-2.2215686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54818654) q[0];
sx q[0];
rz(-2.0596518) q[0];
sx q[0];
rz(-1.5997546) q[0];
rz(-1.2708739) q[1];
sx q[1];
rz(-1.7393232) q[1];
sx q[1];
rz(0.52938968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75225509) q[0];
sx q[0];
rz(-1.2239972) q[0];
sx q[0];
rz(2.8777468) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7725792) q[2];
sx q[2];
rz(-2.1786961) q[2];
sx q[2];
rz(2.4948134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7449977) q[1];
sx q[1];
rz(-2.0740119) q[1];
sx q[1];
rz(-1.2835591) q[1];
x q[2];
rz(1.9233723) q[3];
sx q[3];
rz(-1.7099172) q[3];
sx q[3];
rz(-1.2253882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9993837) q[2];
sx q[2];
rz(-1.9500407) q[2];
sx q[2];
rz(1.8184398) q[2];
rz(-0.38149825) q[3];
sx q[3];
rz(-1.1161209) q[3];
sx q[3];
rz(-0.39314666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8489654) q[0];
sx q[0];
rz(-0.75858527) q[0];
sx q[0];
rz(1.0211771) q[0];
rz(-1.5317597) q[1];
sx q[1];
rz(-2.1949218) q[1];
sx q[1];
rz(-2.3381332) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96413104) q[0];
sx q[0];
rz(-1.2636856) q[0];
sx q[0];
rz(-1.403128) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7293936) q[2];
sx q[2];
rz(-1.4934818) q[2];
sx q[2];
rz(-0.22082034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4804621) q[1];
sx q[1];
rz(-1.4509038) q[1];
sx q[1];
rz(0.57173034) q[1];
rz(-1.4705974) q[3];
sx q[3];
rz(-0.9462983) q[3];
sx q[3];
rz(-1.1794832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44567406) q[2];
sx q[2];
rz(-0.58738223) q[2];
sx q[2];
rz(-2.7109801) q[2];
rz(-2.5066091) q[3];
sx q[3];
rz(-3.1228784) q[3];
sx q[3];
rz(1.2682605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9950614) q[0];
sx q[0];
rz(-2.596031) q[0];
sx q[0];
rz(-1.5978285) q[0];
rz(-2.457288) q[1];
sx q[1];
rz(-1.6938208) q[1];
sx q[1];
rz(-2.3421471) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7921126) q[0];
sx q[0];
rz(-1.6474012) q[0];
sx q[0];
rz(1.5653866) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2861023) q[2];
sx q[2];
rz(-1.5134619) q[2];
sx q[2];
rz(2.3279026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23161665) q[1];
sx q[1];
rz(-1.5808531) q[1];
sx q[1];
rz(-2.6811428) q[1];
x q[2];
rz(0.83752172) q[3];
sx q[3];
rz(-1.9068516) q[3];
sx q[3];
rz(-0.94810644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6058309) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(-2.3568995) q[2];
rz(-0.11624087) q[3];
sx q[3];
rz(-1.5250456) q[3];
sx q[3];
rz(-1.8893265) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6004953) q[0];
sx q[0];
rz(-1.2299812) q[0];
sx q[0];
rz(-2.1121209) q[0];
rz(-3.0211499) q[1];
sx q[1];
rz(-1.30013) q[1];
sx q[1];
rz(-2.5708503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.127205) q[0];
sx q[0];
rz(-0.61349166) q[0];
sx q[0];
rz(-2.2741689) q[0];
x q[1];
rz(-3.013754) q[2];
sx q[2];
rz(-0.72261506) q[2];
sx q[2];
rz(1.3764718) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2201669) q[1];
sx q[1];
rz(-1.1296185) q[1];
sx q[1];
rz(-2.7814072) q[1];
rz(-pi) q[2];
rz(1.086497) q[3];
sx q[3];
rz(-0.37621337) q[3];
sx q[3];
rz(-2.0813326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8735147) q[2];
sx q[2];
rz(-0.87468481) q[2];
sx q[2];
rz(2.6178005) q[2];
rz(-1.9889471) q[3];
sx q[3];
rz(-0.58061424) q[3];
sx q[3];
rz(1.4279648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607587) q[0];
sx q[0];
rz(-1.9402215) q[0];
sx q[0];
rz(-2.7177366) q[0];
rz(2.2747874) q[1];
sx q[1];
rz(-1.024217) q[1];
sx q[1];
rz(-2.9383235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05860672) q[0];
sx q[0];
rz(-1.5241677) q[0];
sx q[0];
rz(-2.4871422) q[0];
rz(-pi) q[1];
rz(1.5166984) q[2];
sx q[2];
rz(-1.2612169) q[2];
sx q[2];
rz(0.54540173) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79938526) q[1];
sx q[1];
rz(-0.3179271) q[1];
sx q[1];
rz(2.6683776) q[1];
x q[2];
rz(0.8747845) q[3];
sx q[3];
rz(-2.751707) q[3];
sx q[3];
rz(2.3830151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.086143494) q[2];
sx q[2];
rz(-1.92675) q[2];
sx q[2];
rz(0.27935585) q[2];
rz(1.2934359) q[3];
sx q[3];
rz(-1.3118298) q[3];
sx q[3];
rz(-1.2156585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348895) q[0];
sx q[0];
rz(-0.80253974) q[0];
sx q[0];
rz(1.7247024) q[0];
rz(-0.92719999) q[1];
sx q[1];
rz(-1.5798774) q[1];
sx q[1];
rz(1.4097811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7977236) q[0];
sx q[0];
rz(-2.2316405) q[0];
sx q[0];
rz(3.0190574) q[0];
x q[1];
rz(2.1883606) q[2];
sx q[2];
rz(-0.57170638) q[2];
sx q[2];
rz(0.99105159) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7076787) q[1];
sx q[1];
rz(-2.0053889) q[1];
sx q[1];
rz(0.7264002) q[1];
x q[2];
rz(0.042155592) q[3];
sx q[3];
rz(-1.32844) q[3];
sx q[3];
rz(-0.66756304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7659144) q[2];
sx q[2];
rz(-1.8199074) q[2];
sx q[2];
rz(-0.970617) q[2];
rz(0.64540234) q[3];
sx q[3];
rz(-2.161945) q[3];
sx q[3];
rz(-2.1360883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29393016) q[0];
sx q[0];
rz(-1.6456589) q[0];
sx q[0];
rz(1.8346067) q[0];
rz(-0.38181276) q[1];
sx q[1];
rz(-1.3419071) q[1];
sx q[1];
rz(0.76795427) q[1];
rz(1.0316331) q[2];
sx q[2];
rz(-0.8349541) q[2];
sx q[2];
rz(-0.63944774) q[2];
rz(2.0313203) q[3];
sx q[3];
rz(-1.1694179) q[3];
sx q[3];
rz(-2.0760134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
