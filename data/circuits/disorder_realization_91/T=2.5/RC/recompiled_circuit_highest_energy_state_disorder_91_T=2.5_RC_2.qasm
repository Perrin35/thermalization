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
rz(-2.7741127) q[0];
sx q[0];
rz(-1.8125266) q[0];
sx q[0];
rz(-1.4867866) q[0];
rz(1.5005255) q[1];
sx q[1];
rz(-2.5379116) q[1];
sx q[1];
rz(-0.85890213) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4444111) q[0];
sx q[0];
rz(-1.0134122) q[0];
sx q[0];
rz(1.0943221) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98318451) q[2];
sx q[2];
rz(-1.1074083) q[2];
sx q[2];
rz(-2.7962229) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9382883) q[1];
sx q[1];
rz(-1.4611725) q[1];
sx q[1];
rz(2.2548864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8789419) q[3];
sx q[3];
rz(-1.7761782) q[3];
sx q[3];
rz(-2.8440203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9438802) q[2];
sx q[2];
rz(-2.3494425) q[2];
sx q[2];
rz(1.0666749) q[2];
rz(-3.0471622) q[3];
sx q[3];
rz(-2.5240199) q[3];
sx q[3];
rz(0.42816952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9601032) q[0];
sx q[0];
rz(-2.8747989) q[0];
sx q[0];
rz(1.2130523) q[0];
rz(-2.6890697) q[1];
sx q[1];
rz(-0.25233832) q[1];
sx q[1];
rz(-2.4620893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5951341) q[0];
sx q[0];
rz(-1.6041557) q[0];
sx q[0];
rz(-2.8921158) q[0];
rz(0.5084796) q[2];
sx q[2];
rz(-2.8114308) q[2];
sx q[2];
rz(-2.1521558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22917381) q[1];
sx q[1];
rz(-2.3245735) q[1];
sx q[1];
rz(-0.55025834) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3704471) q[3];
sx q[3];
rz(-0.89290038) q[3];
sx q[3];
rz(-1.6163449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6015893) q[2];
sx q[2];
rz(-2.3175779) q[2];
sx q[2];
rz(-0.69671112) q[2];
rz(1.8981029) q[3];
sx q[3];
rz(-1.9466629) q[3];
sx q[3];
rz(2.5534326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884647) q[0];
sx q[0];
rz(-2.1051814) q[0];
sx q[0];
rz(0.82365197) q[0];
rz(-0.15692391) q[1];
sx q[1];
rz(-1.4361959) q[1];
sx q[1];
rz(2.7395111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4068459) q[0];
sx q[0];
rz(-1.1752793) q[0];
sx q[0];
rz(2.9989373) q[0];
x q[1];
rz(0.04871647) q[2];
sx q[2];
rz(-1.8557544) q[2];
sx q[2];
rz(1.6866956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3797219) q[1];
sx q[1];
rz(-0.58062299) q[1];
sx q[1];
rz(2.1942744) q[1];
rz(0.22076729) q[3];
sx q[3];
rz(-1.1896092) q[3];
sx q[3];
rz(0.38635283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0154401) q[2];
sx q[2];
rz(-2.2045279) q[2];
sx q[2];
rz(-0.13119571) q[2];
rz(1.0570071) q[3];
sx q[3];
rz(-1.1934692) q[3];
sx q[3];
rz(1.9317651) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7159202) q[0];
sx q[0];
rz(-1.965006) q[0];
sx q[0];
rz(1.7536989) q[0];
rz(-1.7790986) q[1];
sx q[1];
rz(-1.250993) q[1];
sx q[1];
rz(-1.9349792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20509991) q[0];
sx q[0];
rz(-1.1171891) q[0];
sx q[0];
rz(-2.4587247) q[0];
rz(-0.6529867) q[2];
sx q[2];
rz(-2.3781502) q[2];
sx q[2];
rz(-0.12432822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8959103) q[1];
sx q[1];
rz(-2.8295662) q[1];
sx q[1];
rz(-0.25319497) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73829262) q[3];
sx q[3];
rz(-1.2311934) q[3];
sx q[3];
rz(-2.5552354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6917307) q[2];
sx q[2];
rz(-1.8064156) q[2];
sx q[2];
rz(-0.94592363) q[2];
rz(-2.6876884) q[3];
sx q[3];
rz(-2.4920521) q[3];
sx q[3];
rz(-0.18979931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7396963) q[0];
sx q[0];
rz(-1.8296158) q[0];
sx q[0];
rz(-2.358118) q[0];
rz(-2.2354194) q[1];
sx q[1];
rz(-2.2515191) q[1];
sx q[1];
rz(-2.5772212) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6221478) q[0];
sx q[0];
rz(-2.9542202) q[0];
sx q[0];
rz(-2.2617784) q[0];
rz(-pi) q[1];
x q[1];
rz(0.059877574) q[2];
sx q[2];
rz(-1.9765325) q[2];
sx q[2];
rz(-2.6889888) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5003779) q[1];
sx q[1];
rz(-2.1754773) q[1];
sx q[1];
rz(0.7577618) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4268812) q[3];
sx q[3];
rz(-1.6922975) q[3];
sx q[3];
rz(2.7845259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.924661) q[2];
sx q[2];
rz(-1.7380119) q[2];
sx q[2];
rz(-3.0972287) q[2];
rz(-1.7313322) q[3];
sx q[3];
rz(-2.9328465) q[3];
sx q[3];
rz(-1.9467719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57589543) q[0];
sx q[0];
rz(-0.79224753) q[0];
sx q[0];
rz(-2.3727697) q[0];
rz(-0.63111758) q[1];
sx q[1];
rz(-1.2080071) q[1];
sx q[1];
rz(2.9879976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0868743) q[0];
sx q[0];
rz(-1.8970034) q[0];
sx q[0];
rz(0.98814641) q[0];
rz(-2.4903712) q[2];
sx q[2];
rz(-1.3488646) q[2];
sx q[2];
rz(0.62276182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2163917) q[1];
sx q[1];
rz(-1.6782624) q[1];
sx q[1];
rz(2.8587384) q[1];
rz(2.3702341) q[3];
sx q[3];
rz(-1.5017209) q[3];
sx q[3];
rz(1.815064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.067387335) q[2];
sx q[2];
rz(-1.5978483) q[2];
sx q[2];
rz(-2.325861) q[2];
rz(-0.054281209) q[3];
sx q[3];
rz(-1.3868325) q[3];
sx q[3];
rz(-1.3666216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62189198) q[0];
sx q[0];
rz(-1.0655572) q[0];
sx q[0];
rz(0.52072293) q[0];
rz(-2.0479274) q[1];
sx q[1];
rz(-2.2126074) q[1];
sx q[1];
rz(1.3826133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1102229) q[0];
sx q[0];
rz(-0.58833226) q[0];
sx q[0];
rz(-1.3084685) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6135234) q[2];
sx q[2];
rz(-1.1748399) q[2];
sx q[2];
rz(-1.89324) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3471038) q[1];
sx q[1];
rz(-2.1784049) q[1];
sx q[1];
rz(-2.6332864) q[1];
x q[2];
rz(2.6723271) q[3];
sx q[3];
rz(-1.7466144) q[3];
sx q[3];
rz(-1.3202536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19892056) q[2];
sx q[2];
rz(-0.3519381) q[2];
sx q[2];
rz(1.4658296) q[2];
rz(0.28907019) q[3];
sx q[3];
rz(-1.7710268) q[3];
sx q[3];
rz(1.8606961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11956231) q[0];
sx q[0];
rz(-0.95198315) q[0];
sx q[0];
rz(2.6493678) q[0];
rz(2.4960663) q[1];
sx q[1];
rz(-1.5954834) q[1];
sx q[1];
rz(2.2135977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6715901) q[0];
sx q[0];
rz(-0.3151463) q[0];
sx q[0];
rz(1.1130087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7828752) q[2];
sx q[2];
rz(-1.2041098) q[2];
sx q[2];
rz(-2.3334954) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7592589) q[1];
sx q[1];
rz(-1.302726) q[1];
sx q[1];
rz(-2.2925582) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8230439) q[3];
sx q[3];
rz(-1.2179228) q[3];
sx q[3];
rz(-3.1016289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32138985) q[2];
sx q[2];
rz(-1.5265744) q[2];
sx q[2];
rz(-2.3825633) q[2];
rz(1.2718893) q[3];
sx q[3];
rz(-1.8153518) q[3];
sx q[3];
rz(2.429764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7976545) q[0];
sx q[0];
rz(-2.5929218) q[0];
sx q[0];
rz(2.3691673) q[0];
rz(1.3683569) q[1];
sx q[1];
rz(-1.1027579) q[1];
sx q[1];
rz(-0.30430421) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9979663) q[0];
sx q[0];
rz(-1.3990825) q[0];
sx q[0];
rz(2.986326) q[0];
x q[1];
rz(-1.7965806) q[2];
sx q[2];
rz(-1.5743557) q[2];
sx q[2];
rz(2.4709216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.180025) q[1];
sx q[1];
rz(-1.7503382) q[1];
sx q[1];
rz(-0.7431598) q[1];
rz(-pi) q[2];
rz(-1.6966882) q[3];
sx q[3];
rz(-1.1790118) q[3];
sx q[3];
rz(0.34989935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7336537) q[2];
sx q[2];
rz(-2.759178) q[2];
sx q[2];
rz(-1.1747053) q[2];
rz(2.5041653) q[3];
sx q[3];
rz(-1.8791684) q[3];
sx q[3];
rz(-1.2229961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680854) q[0];
sx q[0];
rz(-0.33184505) q[0];
sx q[0];
rz(-2.3688431) q[0];
rz(0.4575153) q[1];
sx q[1];
rz(-1.3950149) q[1];
sx q[1];
rz(-1.7083907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59638176) q[0];
sx q[0];
rz(-2.0186252) q[0];
sx q[0];
rz(-2.2660072) q[0];
rz(-pi) q[1];
rz(1.9338701) q[2];
sx q[2];
rz(-1.6350897) q[2];
sx q[2];
rz(-0.63426547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6683176) q[1];
sx q[1];
rz(-0.5328446) q[1];
sx q[1];
rz(-2.1168296) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3328934) q[3];
sx q[3];
rz(-2.3876841) q[3];
sx q[3];
rz(-1.8659908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2210803) q[2];
sx q[2];
rz(-2.1996193) q[2];
sx q[2];
rz(-2.4648049) q[2];
rz(1.5306728) q[3];
sx q[3];
rz(-0.30083209) q[3];
sx q[3];
rz(-2.0570741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21657011) q[0];
sx q[0];
rz(-1.8393479) q[0];
sx q[0];
rz(1.7787697) q[0];
rz(0.88381797) q[1];
sx q[1];
rz(-2.4088036) q[1];
sx q[1];
rz(2.1688681) q[1];
rz(-0.45223805) q[2];
sx q[2];
rz(-2.0932719) q[2];
sx q[2];
rz(-1.1050638) q[2];
rz(-2.6780405) q[3];
sx q[3];
rz(-1.7854073) q[3];
sx q[3];
rz(2.2525595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
