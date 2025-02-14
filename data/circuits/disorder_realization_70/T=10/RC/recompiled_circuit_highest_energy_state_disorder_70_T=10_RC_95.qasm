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
rz(1.7477859) q[0];
sx q[0];
rz(-2.2414247) q[0];
sx q[0];
rz(1.409344) q[0];
rz(-2.4318168) q[1];
sx q[1];
rz(-0.33325279) q[1];
sx q[1];
rz(-0.38212734) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9242121) q[0];
sx q[0];
rz(-1.4367625) q[0];
sx q[0];
rz(-2.4045375) q[0];
rz(-2.0517778) q[2];
sx q[2];
rz(-0.28533563) q[2];
sx q[2];
rz(1.3621881) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.00013) q[1];
sx q[1];
rz(-1.2034109) q[1];
sx q[1];
rz(0.39333435) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7149147) q[3];
sx q[3];
rz(-0.69732795) q[3];
sx q[3];
rz(2.5643045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0742775) q[2];
sx q[2];
rz(-1.5019608) q[2];
sx q[2];
rz(0.074946694) q[2];
rz(1.8986757) q[3];
sx q[3];
rz(-0.26151812) q[3];
sx q[3];
rz(-2.7459131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13260929) q[0];
sx q[0];
rz(-2.7217396) q[0];
sx q[0];
rz(0.62927759) q[0];
rz(-0.83726007) q[1];
sx q[1];
rz(-0.75051296) q[1];
sx q[1];
rz(-2.5619521) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73635676) q[0];
sx q[0];
rz(-1.7742298) q[0];
sx q[0];
rz(-0.74958165) q[0];
x q[1];
rz(-2.7461065) q[2];
sx q[2];
rz(-1.3278393) q[2];
sx q[2];
rz(2.0045351) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77387735) q[1];
sx q[1];
rz(-1.2396149) q[1];
sx q[1];
rz(0.41172387) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0218815) q[3];
sx q[3];
rz(-0.16489794) q[3];
sx q[3];
rz(2.4462857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2418182) q[2];
sx q[2];
rz(-2.9782229) q[2];
sx q[2];
rz(0.74431288) q[2];
rz(2.7742079) q[3];
sx q[3];
rz(-1.8495879) q[3];
sx q[3];
rz(-1.415409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71017569) q[0];
sx q[0];
rz(-0.95040584) q[0];
sx q[0];
rz(2.1344192) q[0];
rz(-0.011064359) q[1];
sx q[1];
rz(-0.31115752) q[1];
sx q[1];
rz(0.99753582) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5704114) q[0];
sx q[0];
rz(-1.7083252) q[0];
sx q[0];
rz(3.1155706) q[0];
x q[1];
rz(1.6901971) q[2];
sx q[2];
rz(-0.7660256) q[2];
sx q[2];
rz(-1.3090145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.026873206) q[1];
sx q[1];
rz(-2.623154) q[1];
sx q[1];
rz(0.015109574) q[1];
rz(-pi) q[2];
rz(0.51181958) q[3];
sx q[3];
rz(-1.5835754) q[3];
sx q[3];
rz(-0.46841533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1204388) q[2];
sx q[2];
rz(-2.8595371) q[2];
sx q[2];
rz(2.2998478) q[2];
rz(0.65819955) q[3];
sx q[3];
rz(-2.2704312) q[3];
sx q[3];
rz(-0.021520821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-2.4206674) q[0];
sx q[0];
rz(-0.1249211) q[0];
sx q[0];
rz(2.4125873) q[0];
rz(-2.3782102) q[1];
sx q[1];
rz(-0.63012505) q[1];
sx q[1];
rz(-2.7785832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7621146) q[0];
sx q[0];
rz(-1.2804556) q[0];
sx q[0];
rz(-1.3069673) q[0];
rz(-pi) q[1];
rz(-2.5914089) q[2];
sx q[2];
rz(-0.90298684) q[2];
sx q[2];
rz(-3.0444403) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91528401) q[1];
sx q[1];
rz(-1.7123509) q[1];
sx q[1];
rz(-2.841921) q[1];
rz(-pi) q[2];
rz(0.12673817) q[3];
sx q[3];
rz(-1.6235545) q[3];
sx q[3];
rz(-2.3873752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9558941) q[2];
sx q[2];
rz(-2.318435) q[2];
sx q[2];
rz(1.4996747) q[2];
rz(-0.58756346) q[3];
sx q[3];
rz(-0.9891808) q[3];
sx q[3];
rz(0.51830083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8374306) q[0];
sx q[0];
rz(-0.70697933) q[0];
sx q[0];
rz(-0.28513232) q[0];
rz(2.888491) q[1];
sx q[1];
rz(-1.0292091) q[1];
sx q[1];
rz(2.0957799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8469197) q[0];
sx q[0];
rz(-0.6631279) q[0];
sx q[0];
rz(2.4756433) q[0];
rz(-pi) q[1];
x q[1];
rz(1.990647) q[2];
sx q[2];
rz(-1.1641181) q[2];
sx q[2];
rz(-2.5620154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6599413) q[1];
sx q[1];
rz(-1.627198) q[1];
sx q[1];
rz(-0.39404558) q[1];
rz(-pi) q[2];
rz(2.4638484) q[3];
sx q[3];
rz(-1.4137005) q[3];
sx q[3];
rz(1.7435297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4046341) q[2];
sx q[2];
rz(-1.0598695) q[2];
sx q[2];
rz(-0.57445478) q[2];
rz(-2.5308841) q[3];
sx q[3];
rz(-2.6210531) q[3];
sx q[3];
rz(-1.0569388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33559281) q[0];
sx q[0];
rz(-1.3845504) q[0];
sx q[0];
rz(-0.35032508) q[0];
rz(-0.2925182) q[1];
sx q[1];
rz(-3.0151093) q[1];
sx q[1];
rz(0.77936053) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2823863) q[0];
sx q[0];
rz(-0.13257688) q[0];
sx q[0];
rz(-1.1135121) q[0];
rz(-pi) q[1];
rz(0.47996491) q[2];
sx q[2];
rz(-1.4128608) q[2];
sx q[2];
rz(-1.2564893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23664973) q[1];
sx q[1];
rz(-2.7588435) q[1];
sx q[1];
rz(-0.56884258) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52042689) q[3];
sx q[3];
rz(-1.4311858) q[3];
sx q[3];
rz(2.0922144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.902035) q[2];
sx q[2];
rz(-1.6263447) q[2];
sx q[2];
rz(2.1774192) q[2];
rz(0.17298175) q[3];
sx q[3];
rz(-1.0123342) q[3];
sx q[3];
rz(0.13121901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5445589) q[0];
sx q[0];
rz(-1.1976765) q[0];
sx q[0];
rz(-0.069393754) q[0];
rz(1.3736877) q[1];
sx q[1];
rz(-1.2488139) q[1];
sx q[1];
rz(-0.51838851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.501058) q[0];
sx q[0];
rz(-0.69307166) q[0];
sx q[0];
rz(1.3035167) q[0];
x q[1];
rz(2.3396427) q[2];
sx q[2];
rz(-2.5031002) q[2];
sx q[2];
rz(-2.0921767) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1063819) q[1];
sx q[1];
rz(-1.7089881) q[1];
sx q[1];
rz(1.5968678) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33643548) q[3];
sx q[3];
rz(-2.9671768) q[3];
sx q[3];
rz(-0.96272308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9531276) q[2];
sx q[2];
rz(-1.943482) q[2];
sx q[2];
rz(0.24448621) q[2];
rz(-1.2178347) q[3];
sx q[3];
rz(-0.094450258) q[3];
sx q[3];
rz(0.34734669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090488) q[0];
sx q[0];
rz(-0.29260391) q[0];
sx q[0];
rz(2.4421413) q[0];
rz(-2.4199016) q[1];
sx q[1];
rz(-1.7803918) q[1];
sx q[1];
rz(1.9500505) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22868294) q[0];
sx q[0];
rz(-1.7828724) q[0];
sx q[0];
rz(1.6985083) q[0];
rz(-pi) q[1];
rz(1.7315242) q[2];
sx q[2];
rz(-1.3597915) q[2];
sx q[2];
rz(2.4461022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2649768) q[1];
sx q[1];
rz(-0.54538762) q[1];
sx q[1];
rz(-1.7343069) q[1];
x q[2];
rz(0.64719836) q[3];
sx q[3];
rz(-2.5997926) q[3];
sx q[3];
rz(1.7447646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2193853) q[2];
sx q[2];
rz(-0.81874138) q[2];
sx q[2];
rz(-1.4156263) q[2];
rz(2.7042232) q[3];
sx q[3];
rz(-2.8919817) q[3];
sx q[3];
rz(-0.50284809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1943844) q[0];
sx q[0];
rz(-1.9358862) q[0];
sx q[0];
rz(2.6956287) q[0];
rz(-1.3392316) q[1];
sx q[1];
rz(-0.65473348) q[1];
sx q[1];
rz(1.1599249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8561607) q[0];
sx q[0];
rz(-0.92575476) q[0];
sx q[0];
rz(-2.0061532) q[0];
rz(0.075322876) q[2];
sx q[2];
rz(-1.8315856) q[2];
sx q[2];
rz(1.0950077) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7801108) q[1];
sx q[1];
rz(-1.0358397) q[1];
sx q[1];
rz(1.8423716) q[1];
rz(0.91083093) q[3];
sx q[3];
rz(-2.0487222) q[3];
sx q[3];
rz(-0.35199374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40598536) q[2];
sx q[2];
rz(-1.0286101) q[2];
sx q[2];
rz(-2.410991) q[2];
rz(-0.90100151) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(0.034339529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5071097) q[0];
sx q[0];
rz(-2.6431838) q[0];
sx q[0];
rz(-0.46257567) q[0];
rz(-1.0628465) q[1];
sx q[1];
rz(-1.7741508) q[1];
sx q[1];
rz(-0.06632334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6505867) q[0];
sx q[0];
rz(-2.1438476) q[0];
sx q[0];
rz(1.9377524) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0991597) q[2];
sx q[2];
rz(-0.37274578) q[2];
sx q[2];
rz(-2.3605704) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8005261) q[1];
sx q[1];
rz(-2.1459346) q[1];
sx q[1];
rz(2.2209206) q[1];
rz(-pi) q[2];
rz(2.0420426) q[3];
sx q[3];
rz(-0.46898919) q[3];
sx q[3];
rz(2.2820306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.037420951) q[2];
sx q[2];
rz(-2.6080242) q[2];
sx q[2];
rz(1.2332234) q[2];
rz(-0.45796606) q[3];
sx q[3];
rz(-2.870324) q[3];
sx q[3];
rz(-0.39877322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0290699) q[0];
sx q[0];
rz(-1.4334913) q[0];
sx q[0];
rz(1.7423472) q[0];
rz(2.9872672) q[1];
sx q[1];
rz(-1.4230774) q[1];
sx q[1];
rz(1.9565061) q[1];
rz(2.828601) q[2];
sx q[2];
rz(-1.9421158) q[2];
sx q[2];
rz(-2.4696642) q[2];
rz(-2.2521491) q[3];
sx q[3];
rz(-2.4460197) q[3];
sx q[3];
rz(-1.0424436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
