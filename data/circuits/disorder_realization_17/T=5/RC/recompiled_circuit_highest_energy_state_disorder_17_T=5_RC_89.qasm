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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(0.71001473) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(1.8097872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79373097) q[0];
sx q[0];
rz(-1.2958741) q[0];
sx q[0];
rz(1.8040276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3178145) q[2];
sx q[2];
rz(-1.903878) q[2];
sx q[2];
rz(0.87424226) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0441452) q[1];
sx q[1];
rz(-2.4821979) q[1];
sx q[1];
rz(-1.161518) q[1];
x q[2];
rz(2.2043742) q[3];
sx q[3];
rz(-1.211418) q[3];
sx q[3];
rz(-1.0868974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4189202) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(-2.326272) q[2];
rz(-0.80960649) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(-2.1965006) q[1];
sx q[1];
rz(-0.4393622) q[1];
sx q[1];
rz(1.5623215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5961869) q[0];
sx q[0];
rz(-0.89841398) q[0];
sx q[0];
rz(-2.5491712) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1223824) q[2];
sx q[2];
rz(-1.1875936) q[2];
sx q[2];
rz(-0.84898432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72772273) q[1];
sx q[1];
rz(-1.4663823) q[1];
sx q[1];
rz(-0.48522471) q[1];
rz(0.82206313) q[3];
sx q[3];
rz(-1.9371607) q[3];
sx q[3];
rz(1.0613393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14757806) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(0.4501403) q[2];
rz(-2.0077997) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-2.1639737) q[3];
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
rz(2.5728773) q[0];
sx q[0];
rz(-2.1935538) q[0];
sx q[0];
rz(-0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(2.5414355) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900213) q[0];
sx q[0];
rz(-0.29415392) q[0];
sx q[0];
rz(0.16071975) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96397206) q[2];
sx q[2];
rz(-1.6168878) q[2];
sx q[2];
rz(2.5986236) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.197613) q[1];
sx q[1];
rz(-2.463732) q[1];
sx q[1];
rz(-0.53450905) q[1];
rz(-pi) q[2];
rz(-1.3998109) q[3];
sx q[3];
rz(-1.0043084) q[3];
sx q[3];
rz(1.2681707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0448138) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(-2.0951994) q[2];
rz(-0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9699049) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(2.0890253) q[0];
rz(1.9453847) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(-2.722091) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013166817) q[0];
sx q[0];
rz(-1.8343529) q[0];
sx q[0];
rz(-0.58851425) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9151462) q[2];
sx q[2];
rz(-0.94880494) q[2];
sx q[2];
rz(-2.7055307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0177815) q[1];
sx q[1];
rz(-0.7582802) q[1];
sx q[1];
rz(2.1348544) q[1];
x q[2];
rz(1.8617083) q[3];
sx q[3];
rz(-1.0802764) q[3];
sx q[3];
rz(-1.377493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18998751) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(-2.0513963) q[2];
rz(-2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200318) q[0];
sx q[0];
rz(-1.5721385) q[0];
sx q[0];
rz(1.74362) q[0];
rz(-0.13294237) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(-0.001210777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053731669) q[0];
sx q[0];
rz(-2.7606761) q[0];
sx q[0];
rz(0.043405966) q[0];
rz(-pi) q[1];
rz(-1.8563849) q[2];
sx q[2];
rz(-1.4175709) q[2];
sx q[2];
rz(-2.0375348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.397314) q[1];
sx q[1];
rz(-0.71956735) q[1];
sx q[1];
rz(-2.6707763) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4210971) q[3];
sx q[3];
rz(-1.544501) q[3];
sx q[3];
rz(-0.61248518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2289537) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(0.21052989) q[2];
rz(2.1052965) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(-1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(0.92232651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3137445) q[0];
sx q[0];
rz(-1.2822312) q[0];
sx q[0];
rz(-3.0305578) q[0];
rz(-3.0386557) q[2];
sx q[2];
rz(-2.2349226) q[2];
sx q[2];
rz(0.31598202) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90627024) q[1];
sx q[1];
rz(-0.962469) q[1];
sx q[1];
rz(-0.13923213) q[1];
rz(-pi) q[2];
rz(0.57862307) q[3];
sx q[3];
rz(-2.4994183) q[3];
sx q[3];
rz(0.35143055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4313844) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(-2.7396438) q[2];
rz(1.8534144) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(1.2791546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61889082) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(0.064362137) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.3305957) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272328) q[0];
sx q[0];
rz(-1.4454675) q[0];
sx q[0];
rz(-2.5663239) q[0];
rz(2.0688017) q[2];
sx q[2];
rz(-2.118022) q[2];
sx q[2];
rz(0.71820532) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45695282) q[1];
sx q[1];
rz(-1.350953) q[1];
sx q[1];
rz(-0.5348285) q[1];
rz(-pi) q[2];
rz(-2.9693611) q[3];
sx q[3];
rz(-0.96559286) q[3];
sx q[3];
rz(-1.8021405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7360721) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(-0.68354496) q[2];
rz(-0.73379597) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7965294) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(-1.9684568) q[0];
rz(2.7660811) q[1];
sx q[1];
rz(-0.48986062) q[1];
sx q[1];
rz(-2.1048022) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4597367) q[0];
sx q[0];
rz(-1.5516743) q[0];
sx q[0];
rz(-0.023650344) q[0];
rz(-pi) q[1];
rz(1.2944968) q[2];
sx q[2];
rz(-1.7299998) q[2];
sx q[2];
rz(2.4197848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45932508) q[1];
sx q[1];
rz(-0.81574355) q[1];
sx q[1];
rz(-1.7641279) q[1];
x q[2];
rz(0.49077197) q[3];
sx q[3];
rz(-1.0939071) q[3];
sx q[3];
rz(-0.27974883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5658687) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(-0.73053989) q[2];
rz(2.9324487) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(-1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.4850979) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(-3.094161) q[0];
rz(2.9196396) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(-0.91112959) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.52235) q[0];
sx q[0];
rz(-1.6043092) q[0];
sx q[0];
rz(1.5977809) q[0];
x q[1];
rz(-3.1268551) q[2];
sx q[2];
rz(-2.212425) q[2];
sx q[2];
rz(-0.82889885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67942372) q[1];
sx q[1];
rz(-2.7738681) q[1];
sx q[1];
rz(-3.1064139) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0044596) q[3];
sx q[3];
rz(-1.3977504) q[3];
sx q[3];
rz(-0.53849788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86965108) q[2];
sx q[2];
rz(-2.7615774) q[2];
sx q[2];
rz(-0.16858777) q[2];
rz(-0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.8144089) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1472226) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(-0.45743531) q[0];
rz(-1.2262729) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(1.3630684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5209893) q[0];
sx q[0];
rz(-2.531157) q[0];
sx q[0];
rz(-0.10535289) q[0];
x q[1];
rz(2.8122452) q[2];
sx q[2];
rz(-2.4852356) q[2];
sx q[2];
rz(-0.48329566) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7292028) q[1];
sx q[1];
rz(-2.7250184) q[1];
sx q[1];
rz(2.7135486) q[1];
rz(-pi) q[2];
rz(3.0567597) q[3];
sx q[3];
rz(-0.89003497) q[3];
sx q[3];
rz(-2.5670402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3739796) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(2.6895831) q[2];
rz(2.7376145) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(1.1261136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6954738) q[0];
sx q[0];
rz(-1.6071381) q[0];
sx q[0];
rz(0.18679609) q[0];
rz(2.6192464) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(-1.3138554) q[2];
sx q[2];
rz(-1.8443454) q[2];
sx q[2];
rz(0.92583427) q[2];
rz(1.4215076) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
