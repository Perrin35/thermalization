OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96149093) q[0];
sx q[0];
rz(3.7082727) q[0];
sx q[0];
rz(10.215496) q[0];
rz(-0.92791954) q[1];
sx q[1];
rz(-0.035049573) q[1];
sx q[1];
rz(-2.8383281) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7105616) q[0];
sx q[0];
rz(-0.91791422) q[0];
sx q[0];
rz(1.8588572) q[0];
rz(-0.35776414) q[2];
sx q[2];
rz(-1.9310037) q[2];
sx q[2];
rz(1.2586762) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7590847) q[1];
sx q[1];
rz(-1.0945935) q[1];
sx q[1];
rz(-0.97530968) q[1];
x q[2];
rz(-3.0500134) q[3];
sx q[3];
rz(-0.73741961) q[3];
sx q[3];
rz(2.9980286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7683893) q[2];
sx q[2];
rz(-0.36036569) q[2];
sx q[2];
rz(-2.4599794) q[2];
rz(2.5810589) q[3];
sx q[3];
rz(-0.78442854) q[3];
sx q[3];
rz(-2.4973629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.949837) q[0];
sx q[0];
rz(-0.90832174) q[0];
sx q[0];
rz(-2.7857696) q[0];
rz(0.25335723) q[1];
sx q[1];
rz(-0.8496049) q[1];
sx q[1];
rz(1.2454978) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863128) q[0];
sx q[0];
rz(-1.7944424) q[0];
sx q[0];
rz(-2.6564084) q[0];
rz(-pi) q[1];
rz(-2.5077057) q[2];
sx q[2];
rz(-1.1340464) q[2];
sx q[2];
rz(0.89981198) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.624234) q[1];
sx q[1];
rz(-1.9962203) q[1];
sx q[1];
rz(2.5555771) q[1];
x q[2];
rz(-0.83578556) q[3];
sx q[3];
rz(-2.4525149) q[3];
sx q[3];
rz(-3.0398577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7443098) q[2];
sx q[2];
rz(-2.3645568) q[2];
sx q[2];
rz(-3.0139319) q[2];
rz(2.4658261) q[3];
sx q[3];
rz(-0.45857576) q[3];
sx q[3];
rz(-0.42135409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0965939) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(1.1340207) q[0];
rz(-2.2484312) q[1];
sx q[1];
rz(-0.90407073) q[1];
sx q[1];
rz(-1.2718511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24090996) q[0];
sx q[0];
rz(-2.5782452) q[0];
sx q[0];
rz(-1.5213837) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8804466) q[2];
sx q[2];
rz(-1.4952588) q[2];
sx q[2];
rz(2.6713918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0772897) q[1];
sx q[1];
rz(-1.3659119) q[1];
sx q[1];
rz(0.85475342) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92471931) q[3];
sx q[3];
rz(-2.8689403) q[3];
sx q[3];
rz(0.70875185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21358718) q[2];
sx q[2];
rz(-2.2760133) q[2];
sx q[2];
rz(-1.340284) q[2];
rz(-1.853893) q[3];
sx q[3];
rz(-1.4475334) q[3];
sx q[3];
rz(0.48937669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61557788) q[0];
sx q[0];
rz(-2.5647793) q[0];
sx q[0];
rz(2.4082129) q[0];
rz(0.79542696) q[1];
sx q[1];
rz(-2.9454102) q[1];
sx q[1];
rz(1.9733852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8775032) q[0];
sx q[0];
rz(-1.6764219) q[0];
sx q[0];
rz(1.9838404) q[0];
x q[1];
rz(-2.204037) q[2];
sx q[2];
rz(-2.0481234) q[2];
sx q[2];
rz(3.1138495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2760022) q[1];
sx q[1];
rz(-1.6011213) q[1];
sx q[1];
rz(0.042976168) q[1];
rz(-0.9549035) q[3];
sx q[3];
rz(-2.1400937) q[3];
sx q[3];
rz(-2.3950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1802133) q[2];
sx q[2];
rz(-0.81575477) q[2];
sx q[2];
rz(-0.95480314) q[2];
rz(1.0500326) q[3];
sx q[3];
rz(-0.78767109) q[3];
sx q[3];
rz(-0.55154705) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0410205) q[0];
sx q[0];
rz(-2.6028778) q[0];
sx q[0];
rz(2.4464497) q[0];
rz(-1.2160542) q[1];
sx q[1];
rz(-2.1247517) q[1];
sx q[1];
rz(3.1207808) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7835158) q[0];
sx q[0];
rz(-2.4672271) q[0];
sx q[0];
rz(-1.3410617) q[0];
rz(-pi) q[1];
rz(1.2131507) q[2];
sx q[2];
rz(-2.3721173) q[2];
sx q[2];
rz(-2.3703252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4592308) q[1];
sx q[1];
rz(-1.3317079) q[1];
sx q[1];
rz(3.1296631) q[1];
rz(-pi) q[2];
rz(2.5567804) q[3];
sx q[3];
rz(-3.0303114) q[3];
sx q[3];
rz(2.5113784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9713328) q[2];
sx q[2];
rz(-0.33997619) q[2];
sx q[2];
rz(0.60410947) q[2];
rz(2.4938834) q[3];
sx q[3];
rz(-0.52102399) q[3];
sx q[3];
rz(-2.6807396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1243509) q[0];
sx q[0];
rz(-1.0931953) q[0];
sx q[0];
rz(0.35853115) q[0];
rz(-2.5784946) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(2.8360352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222851) q[0];
sx q[0];
rz(-0.74757517) q[0];
sx q[0];
rz(-2.5966552) q[0];
rz(-1.7420058) q[2];
sx q[2];
rz(-1.7504901) q[2];
sx q[2];
rz(3.0269738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87898743) q[1];
sx q[1];
rz(-1.078558) q[1];
sx q[1];
rz(-0.61779554) q[1];
rz(1.5005388) q[3];
sx q[3];
rz(-1.4936325) q[3];
sx q[3];
rz(2.7479395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.87865889) q[2];
sx q[2];
rz(-0.96076751) q[2];
sx q[2];
rz(3.089454) q[2];
rz(-2.7122998) q[3];
sx q[3];
rz(-1.7310127) q[3];
sx q[3];
rz(2.8894292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4572064) q[0];
sx q[0];
rz(-2.6173213) q[0];
sx q[0];
rz(-2.9065409) q[0];
rz(0.60335195) q[1];
sx q[1];
rz(-0.82913202) q[1];
sx q[1];
rz(2.5650909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7991195) q[0];
sx q[0];
rz(-2.2469889) q[0];
sx q[0];
rz(-2.4732117) q[0];
x q[1];
rz(0.7754972) q[2];
sx q[2];
rz(-1.8172872) q[2];
sx q[2];
rz(-1.6131439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53945069) q[1];
sx q[1];
rz(-1.2509402) q[1];
sx q[1];
rz(1.9906688) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1436809) q[3];
sx q[3];
rz(-2.0485123) q[3];
sx q[3];
rz(-2.9766072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0779886) q[2];
sx q[2];
rz(-0.84090191) q[2];
sx q[2];
rz(-1.3464751) q[2];
rz(-2.6320362) q[3];
sx q[3];
rz(-1.2600803) q[3];
sx q[3];
rz(-0.30663651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6154196) q[0];
sx q[0];
rz(-1.0831447) q[0];
sx q[0];
rz(2.7857067) q[0];
rz(2.4399759) q[1];
sx q[1];
rz(-2.9569148) q[1];
sx q[1];
rz(2.1772749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50011694) q[0];
sx q[0];
rz(-1.4371514) q[0];
sx q[0];
rz(2.0290613) q[0];
rz(-pi) q[1];
rz(-3.1213254) q[2];
sx q[2];
rz(-1.8635443) q[2];
sx q[2];
rz(2.0890638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0135994) q[1];
sx q[1];
rz(-2.0656054) q[1];
sx q[1];
rz(1.5484936) q[1];
x q[2];
rz(1.4889433) q[3];
sx q[3];
rz(-1.8813895) q[3];
sx q[3];
rz(1.6037343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4926766) q[2];
sx q[2];
rz(-2.4339088) q[2];
sx q[2];
rz(0.81563449) q[2];
rz(2.2260769) q[3];
sx q[3];
rz(-2.2736277) q[3];
sx q[3];
rz(-0.22667949) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30064073) q[0];
sx q[0];
rz(-0.17423593) q[0];
sx q[0];
rz(-1.9006282) q[0];
rz(-3.0366483) q[1];
sx q[1];
rz(-0.34067708) q[1];
sx q[1];
rz(0.45936432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3973839) q[0];
sx q[0];
rz(-2.7046015) q[0];
sx q[0];
rz(2.1676012) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2392339) q[2];
sx q[2];
rz(-1.3156097) q[2];
sx q[2];
rz(-2.55748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33973673) q[1];
sx q[1];
rz(-1.1366399) q[1];
sx q[1];
rz(1.3375907) q[1];
rz(-pi) q[2];
rz(2.5673994) q[3];
sx q[3];
rz(-1.4530655) q[3];
sx q[3];
rz(1.2494465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0910054) q[2];
sx q[2];
rz(-0.688474) q[2];
sx q[2];
rz(-0.090593226) q[2];
rz(2.374384) q[3];
sx q[3];
rz(-1.6588666) q[3];
sx q[3];
rz(-1.6902794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42152688) q[0];
sx q[0];
rz(-2.2179715) q[0];
sx q[0];
rz(-1.5371171) q[0];
rz(0.9556669) q[1];
sx q[1];
rz(-1.4612528) q[1];
sx q[1];
rz(2.59424) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794649) q[0];
sx q[0];
rz(-1.5766931) q[0];
sx q[0];
rz(1.4910758) q[0];
x q[1];
rz(-0.53031594) q[2];
sx q[2];
rz(-1.2387453) q[2];
sx q[2];
rz(0.9584934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9319084) q[1];
sx q[1];
rz(-2.3386777) q[1];
sx q[1];
rz(1.241339) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56741107) q[3];
sx q[3];
rz(-2.3443522) q[3];
sx q[3];
rz(2.286943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1110111) q[2];
sx q[2];
rz(-2.9647398) q[2];
sx q[2];
rz(-2.492823) q[2];
rz(0.23160058) q[3];
sx q[3];
rz(-2.2571371) q[3];
sx q[3];
rz(-3.054936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(-1.9301013) q[0];
sx q[0];
rz(1.7420266) q[0];
rz(2.4823785) q[1];
sx q[1];
rz(-1.8623687) q[1];
sx q[1];
rz(1.6946793) q[1];
rz(-1.7805889) q[2];
sx q[2];
rz(-0.2332415) q[2];
sx q[2];
rz(2.94549) q[2];
rz(-1.4249887) q[3];
sx q[3];
rz(-1.5167699) q[3];
sx q[3];
rz(-1.9595893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
