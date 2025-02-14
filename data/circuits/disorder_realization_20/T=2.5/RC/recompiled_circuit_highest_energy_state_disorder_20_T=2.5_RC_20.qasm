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
rz(-0.72467726) q[0];
sx q[0];
rz(4.4530498) q[0];
sx q[0];
rz(11.761576) q[0];
rz(0.59504396) q[1];
sx q[1];
rz(-2.9437328) q[1];
sx q[1];
rz(-1.9207538) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0752995) q[0];
sx q[0];
rz(-1.4508712) q[0];
sx q[0];
rz(-1.6060702) q[0];
rz(-0.27313613) q[2];
sx q[2];
rz(-2.663586) q[2];
sx q[2];
rz(0.82551685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2455622) q[1];
sx q[1];
rz(-0.82791019) q[1];
sx q[1];
rz(1.5106535) q[1];
rz(-pi) q[2];
rz(1.4335459) q[3];
sx q[3];
rz(-2.1992634) q[3];
sx q[3];
rz(-0.31627218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4294942) q[2];
sx q[2];
rz(-0.50769371) q[2];
sx q[2];
rz(-1.3362159) q[2];
rz(1.7441689) q[3];
sx q[3];
rz(-1.5066527) q[3];
sx q[3];
rz(-1.5857182) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3812688) q[0];
sx q[0];
rz(-1.5070494) q[0];
sx q[0];
rz(0.67833483) q[0];
rz(-0.61596576) q[1];
sx q[1];
rz(-1.0266285) q[1];
sx q[1];
rz(2.9982627) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6852138) q[0];
sx q[0];
rz(-2.7806578) q[0];
sx q[0];
rz(-1.1120615) q[0];
rz(-pi) q[1];
rz(3.0264261) q[2];
sx q[2];
rz(-2.2858742) q[2];
sx q[2];
rz(0.63181704) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7166087) q[1];
sx q[1];
rz(-1.0657446) q[1];
sx q[1];
rz(-1.5165971) q[1];
rz(-pi) q[2];
rz(-0.42558221) q[3];
sx q[3];
rz(-1.2493361) q[3];
sx q[3];
rz(1.6265353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6501288) q[2];
sx q[2];
rz(-0.77892196) q[2];
sx q[2];
rz(1.9219363) q[2];
rz(-0.31798142) q[3];
sx q[3];
rz(-0.82428437) q[3];
sx q[3];
rz(2.4021751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3984482) q[0];
sx q[0];
rz(-0.92107451) q[0];
sx q[0];
rz(2.5019124) q[0];
rz(-1.3321053) q[1];
sx q[1];
rz(-0.94146856) q[1];
sx q[1];
rz(2.926362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23519653) q[0];
sx q[0];
rz(-1.0825038) q[0];
sx q[0];
rz(1.8129691) q[0];
rz(-1.2864743) q[2];
sx q[2];
rz(-1.5657565) q[2];
sx q[2];
rz(1.7667631) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58956026) q[1];
sx q[1];
rz(-1.560324) q[1];
sx q[1];
rz(1.9661436) q[1];
rz(3.0212175) q[3];
sx q[3];
rz(-2.0248499) q[3];
sx q[3];
rz(1.9714485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3053863) q[2];
sx q[2];
rz(-0.47577327) q[2];
sx q[2];
rz(-2.9212908) q[2];
rz(2.3588755) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(0.14557423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43133217) q[0];
sx q[0];
rz(-2.1556957) q[0];
sx q[0];
rz(1.8248935) q[0];
rz(-0.45007625) q[1];
sx q[1];
rz(-1.7066259) q[1];
sx q[1];
rz(3.0302474) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25791439) q[0];
sx q[0];
rz(-0.92843572) q[0];
sx q[0];
rz(1.8877554) q[0];
rz(-1.8315122) q[2];
sx q[2];
rz(-0.33277724) q[2];
sx q[2];
rz(2.8944601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3831454) q[1];
sx q[1];
rz(-2.1229593) q[1];
sx q[1];
rz(-0.32917413) q[1];
rz(-2.4347859) q[3];
sx q[3];
rz(-0.68663678) q[3];
sx q[3];
rz(1.5070908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7340118) q[2];
sx q[2];
rz(-1.13052) q[2];
sx q[2];
rz(0.13047516) q[2];
rz(-2.981251) q[3];
sx q[3];
rz(-0.66157833) q[3];
sx q[3];
rz(-2.9714238) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57123667) q[0];
sx q[0];
rz(-0.2061051) q[0];
sx q[0];
rz(0.11827949) q[0];
rz(-0.89625278) q[1];
sx q[1];
rz(-1.6629442) q[1];
sx q[1];
rz(-0.55312696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5251585) q[0];
sx q[0];
rz(-1.2431974) q[0];
sx q[0];
rz(2.6382951) q[0];
x q[1];
rz(-2.011492) q[2];
sx q[2];
rz(-0.88437786) q[2];
sx q[2];
rz(-2.7946607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.432888) q[1];
sx q[1];
rz(-0.9917534) q[1];
sx q[1];
rz(-0.20198532) q[1];
rz(-pi) q[2];
rz(-0.40377577) q[3];
sx q[3];
rz(-2.194768) q[3];
sx q[3];
rz(-1.4194907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3502189) q[2];
sx q[2];
rz(-2.4198136) q[2];
sx q[2];
rz(2.1283894) q[2];
rz(-0.31747097) q[3];
sx q[3];
rz(-1.443913) q[3];
sx q[3];
rz(2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233699) q[0];
sx q[0];
rz(-0.88339266) q[0];
sx q[0];
rz(-0.50514847) q[0];
rz(2.743424) q[1];
sx q[1];
rz(-1.265637) q[1];
sx q[1];
rz(-1.8123951) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58513858) q[0];
sx q[0];
rz(-0.49665652) q[0];
sx q[0];
rz(-1.9209548) q[0];
rz(-2.6077725) q[2];
sx q[2];
rz(-1.8835734) q[2];
sx q[2];
rz(2.52616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8230629) q[1];
sx q[1];
rz(-1.4612056) q[1];
sx q[1];
rz(-2.1353882) q[1];
x q[2];
rz(-0.81327) q[3];
sx q[3];
rz(-0.36618847) q[3];
sx q[3];
rz(-1.0367108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.053981606) q[2];
sx q[2];
rz(-1.3302646) q[2];
sx q[2];
rz(2.3114253) q[2];
rz(-1.4812482) q[3];
sx q[3];
rz(-2.0989428) q[3];
sx q[3];
rz(-1.9771077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2723349) q[0];
sx q[0];
rz(-2.2255958) q[0];
sx q[0];
rz(-1.9328312) q[0];
rz(2.8973978) q[1];
sx q[1];
rz(-1.1714275) q[1];
sx q[1];
rz(-0.29921439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1580762) q[0];
sx q[0];
rz(-1.7868446) q[0];
sx q[0];
rz(-0.026834994) q[0];
x q[1];
rz(2.1019548) q[2];
sx q[2];
rz(-2.1130145) q[2];
sx q[2];
rz(-1.7249267) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0639887) q[1];
sx q[1];
rz(-0.65404679) q[1];
sx q[1];
rz(1.6507504) q[1];
rz(0.47513606) q[3];
sx q[3];
rz(-2.3804602) q[3];
sx q[3];
rz(2.6360398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19305912) q[2];
sx q[2];
rz(-1.7391354) q[2];
sx q[2];
rz(-0.94686568) q[2];
rz(1.3777422) q[3];
sx q[3];
rz(-0.54072127) q[3];
sx q[3];
rz(-0.64259678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4544025) q[0];
sx q[0];
rz(-2.7070422) q[0];
sx q[0];
rz(2.6026671) q[0];
rz(-2.9680805) q[1];
sx q[1];
rz(-0.75650802) q[1];
sx q[1];
rz(2.7573746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464041) q[0];
sx q[0];
rz(-2.6676819) q[0];
sx q[0];
rz(-2.7907811) q[0];
rz(-pi) q[1];
rz(0.60224763) q[2];
sx q[2];
rz(-1.7302697) q[2];
sx q[2];
rz(-0.98724706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.780336) q[1];
sx q[1];
rz(-1.5118264) q[1];
sx q[1];
rz(-0.026193729) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.027111091) q[3];
sx q[3];
rz(-1.016992) q[3];
sx q[3];
rz(1.1983271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3197799) q[2];
sx q[2];
rz(-0.76924789) q[2];
sx q[2];
rz(2.5452851) q[2];
rz(-1.0373621) q[3];
sx q[3];
rz(-0.77330247) q[3];
sx q[3];
rz(-1.1843225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7647112) q[0];
sx q[0];
rz(-2.3887971) q[0];
sx q[0];
rz(2.9994614) q[0];
rz(-2.0243952) q[1];
sx q[1];
rz(-1.486472) q[1];
sx q[1];
rz(-1.9140917) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3820597) q[0];
sx q[0];
rz(-2.5798424) q[0];
sx q[0];
rz(-2.6047833) q[0];
x q[1];
rz(1.4332716) q[2];
sx q[2];
rz(-1.5548717) q[2];
sx q[2];
rz(-0.28057306) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.092526768) q[1];
sx q[1];
rz(-0.62916189) q[1];
sx q[1];
rz(1.4607304) q[1];
rz(-pi) q[2];
rz(-2.7015721) q[3];
sx q[3];
rz(-0.87910637) q[3];
sx q[3];
rz(1.7753851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0907937) q[2];
sx q[2];
rz(-2.8992081) q[2];
sx q[2];
rz(2.3432815) q[2];
rz(-1.5618886) q[3];
sx q[3];
rz(-1.3825994) q[3];
sx q[3];
rz(-0.17670512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0853737) q[0];
sx q[0];
rz(-0.64301411) q[0];
sx q[0];
rz(3.1363078) q[0];
rz(-0.082848631) q[1];
sx q[1];
rz(-2.0382035) q[1];
sx q[1];
rz(-2.8740035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0326313) q[0];
sx q[0];
rz(-1.259045) q[0];
sx q[0];
rz(0.70272081) q[0];
x q[1];
rz(-1.1830953) q[2];
sx q[2];
rz(-1.299721) q[2];
sx q[2];
rz(-0.77927206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53658518) q[1];
sx q[1];
rz(-1.7844024) q[1];
sx q[1];
rz(1.2475292) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6676099) q[3];
sx q[3];
rz(-1.1055531) q[3];
sx q[3];
rz(-1.7758235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5679607) q[2];
sx q[2];
rz(-2.6601514) q[2];
sx q[2];
rz(1.2130223) q[2];
rz(-1.2447478) q[3];
sx q[3];
rz(-2.6601578) q[3];
sx q[3];
rz(1.1904233) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5382814) q[0];
sx q[0];
rz(-0.9700226) q[0];
sx q[0];
rz(0.097401311) q[0];
rz(3.1311323) q[1];
sx q[1];
rz(-0.86581007) q[1];
sx q[1];
rz(2.5444358) q[1];
rz(-0.53562989) q[2];
sx q[2];
rz(-1.6471072) q[2];
sx q[2];
rz(-2.9181388) q[2];
rz(2.1599471) q[3];
sx q[3];
rz(-1.4245778) q[3];
sx q[3];
rz(2.1733976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
