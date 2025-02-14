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
rz(-2.4315779) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(11.234565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0667257) q[0];
sx q[0];
rz(-2.7829889) q[0];
sx q[0];
rz(-2.4551366) q[0];
rz(-2.7983886) q[2];
sx q[2];
rz(-1.3319974) q[2];
sx q[2];
rz(-2.3607139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59939042) q[1];
sx q[1];
rz(-0.97394651) q[1];
sx q[1];
rz(-2.8423897) q[1];
rz(-pi) q[2];
rz(0.93721849) q[3];
sx q[3];
rz(-1.211418) q[3];
sx q[3];
rz(-2.0546953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72267246) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(0.81532064) q[2];
rz(0.80960649) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(-2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.078449) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(3.0310042) q[0];
rz(0.94509205) q[1];
sx q[1];
rz(-0.4393622) q[1];
sx q[1];
rz(1.5623215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62842436) q[0];
sx q[0];
rz(-1.1187176) q[0];
sx q[0];
rz(-0.80597178) q[0];
rz(2.0192102) q[2];
sx q[2];
rz(-1.1875936) q[2];
sx q[2];
rz(2.2926083) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4936714) q[1];
sx q[1];
rz(-2.646138) q[1];
sx q[1];
rz(-0.22101553) q[1];
rz(-2.0840629) q[3];
sx q[3];
rz(-2.3239414) q[3];
sx q[3];
rz(3.0000837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9940146) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(-2.6914524) q[2];
rz(1.133793) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728773) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(0.6001572) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068014) q[0];
sx q[0];
rz(-1.8610483) q[0];
sx q[0];
rz(-1.5223548) q[0];
x q[1];
rz(-0.96397206) q[2];
sx q[2];
rz(-1.6168878) q[2];
sx q[2];
rz(-2.5986236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0827209) q[1];
sx q[1];
rz(-1.8959672) q[1];
sx q[1];
rz(-0.60589686) q[1];
rz(-pi) q[2];
rz(2.8801877) q[3];
sx q[3];
rz(-0.58901826) q[3];
sx q[3];
rz(0.95688577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0448138) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(2.0951994) q[2];
rz(-0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-3.0090295) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1716877) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(-2.0890253) q[0];
rz(-1.9453847) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(-2.722091) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3854806) q[0];
sx q[0];
rz(-1.0051553) q[0];
sx q[0];
rz(-1.2570981) q[0];
x q[1];
rz(1.8742895) q[2];
sx q[2];
rz(-0.65676531) q[2];
sx q[2];
rz(3.0820897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12246015) q[1];
sx q[1];
rz(-1.9472709) q[1];
sx q[1];
rz(2.2458162) q[1];
rz(-pi) q[2];
rz(2.64873) q[3];
sx q[3];
rz(-0.56418428) q[3];
sx q[3];
rz(1.197937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9516051) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(-1.0901964) q[2];
rz(-2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215609) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.74362) q[0];
rz(-0.13294237) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(3.1403819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346136) q[0];
sx q[0];
rz(-1.190257) q[0];
sx q[0];
rz(-1.5881722) q[0];
rz(-1.8563849) q[2];
sx q[2];
rz(-1.7240217) q[2];
sx q[2];
rz(-1.1040579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9494755) q[1];
sx q[1];
rz(-1.8743974) q[1];
sx q[1];
rz(0.66302256) q[1];
rz(-pi) q[2];
rz(-1.7204956) q[3];
sx q[3];
rz(-1.544501) q[3];
sx q[3];
rz(0.61248518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(-0.21052989) q[2];
rz(1.0362961) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(-1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(-2.7358828) q[0];
rz(-1.3544719) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(0.92232651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2887696) q[0];
sx q[0];
rz(-1.6772224) q[0];
sx q[0];
rz(-1.861051) q[0];
rz(-pi) q[1];
rz(-1.4402661) q[2];
sx q[2];
rz(-2.4707327) q[2];
sx q[2];
rz(2.6595569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90627024) q[1];
sx q[1];
rz(-2.1791237) q[1];
sx q[1];
rz(3.0023605) q[1];
rz(-pi) q[2];
rz(-0.57862307) q[3];
sx q[3];
rz(-2.4994183) q[3];
sx q[3];
rz(-0.35143055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(0.40194884) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(-1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(0.064362137) q[0];
rz(-0.83802682) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.8109969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9884686) q[0];
sx q[0];
rz(-0.58725427) q[0];
sx q[0];
rz(2.9140317) q[0];
x q[1];
rz(2.5352201) q[2];
sx q[2];
rz(-1.9909711) q[2];
sx q[2];
rz(2.5647031) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8992675) q[1];
sx q[1];
rz(-1.0501672) q[1];
sx q[1];
rz(1.8249056) q[1];
x q[2];
rz(-1.3279541) q[3];
sx q[3];
rz(-2.5153219) q[3];
sx q[3];
rz(1.0427208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4055206) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(-2.4580477) q[2];
rz(-0.73379597) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(0.84764135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(-1.1731359) q[0];
rz(-2.7660811) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(-2.1048022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6818559) q[0];
sx q[0];
rz(-1.5899183) q[0];
sx q[0];
rz(3.1179423) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16536819) q[2];
sx q[2];
rz(-1.2980808) q[2];
sx q[2];
rz(0.80406666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97809659) q[1];
sx q[1];
rz(-1.4304203) q[1];
sx q[1];
rz(2.3771493) q[1];
rz(-pi) q[2];
rz(1.0408786) q[3];
sx q[3];
rz(-1.138759) q[3];
sx q[3];
rz(-2.0910859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.57572395) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(0.73053989) q[2];
rz(2.9324487) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(-1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65649477) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(-3.094161) q[0];
rz(2.9196396) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(2.2304631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1909433) q[0];
sx q[0];
rz(-1.5438269) q[0];
sx q[0];
rz(-3.1080676) q[0];
rz(-pi) q[1];
rz(-1.59052) q[2];
sx q[2];
rz(-2.4998186) q[2];
sx q[2];
rz(2.3373147) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67942372) q[1];
sx q[1];
rz(-0.36772455) q[1];
sx q[1];
rz(-3.1064139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8857314) q[3];
sx q[3];
rz(-0.58940998) q[3];
sx q[3];
rz(-2.3737597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2719416) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(-2.9730049) q[2];
rz(2.9224959) q[3];
sx q[3];
rz(-2.2897661) q[3];
sx q[3];
rz(1.8144089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.99437) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(2.6841573) q[0];
rz(-1.2262729) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(1.7785243) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4922614) q[0];
sx q[0];
rz(-2.1773585) q[0];
sx q[0];
rz(-1.6442293) q[0];
x q[1];
rz(-0.3293475) q[2];
sx q[2];
rz(-0.65635704) q[2];
sx q[2];
rz(0.48329566) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7292028) q[1];
sx q[1];
rz(-2.7250184) q[1];
sx q[1];
rz(-0.42804407) q[1];
rz(-pi) q[2];
rz(-1.4665589) q[3];
sx q[3];
rz(-0.6851894) q[3];
sx q[3];
rz(0.70885056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7676131) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(0.45200959) q[2];
rz(-0.40397817) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(-1.1261136) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(0.52234621) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(1.8277373) q[2];
sx q[2];
rz(-1.8443454) q[2];
sx q[2];
rz(0.92583427) q[2];
rz(-2.8780739) q[3];
sx q[3];
rz(-1.4265887) q[3];
sx q[3];
rz(-1.9909457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
