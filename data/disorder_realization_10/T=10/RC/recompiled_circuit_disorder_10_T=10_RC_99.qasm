OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(1.2115275) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(0.29830631) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878108) q[0];
sx q[0];
rz(-0.90667533) q[0];
sx q[0];
rz(1.5009297) q[0];
rz(1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(0.36117902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.28042291) q[1];
sx q[1];
rz(-1.0958332) q[1];
sx q[1];
rz(2.1850987) q[1];
rz(-pi) q[2];
rz(1.6739507) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(-0.32675693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6478708) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(-0.018218536) q[0];
rz(-2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(0.47168628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6149711) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(-1.4814266) q[0];
rz(-pi) q[1];
rz(-1.8505487) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(-0.56088698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7447409) q[1];
sx q[1];
rz(-2.5137915) q[1];
sx q[1];
rz(0.46698924) q[1];
x q[2];
rz(1.7945812) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(-2.1686045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(2.1726051) q[2];
rz(-2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(1.0415174) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-0.18181268) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(1.6859432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239379) q[0];
sx q[0];
rz(-2.0485989) q[0];
sx q[0];
rz(-0.91207232) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9934898) q[2];
sx q[2];
rz(-2.2659677) q[2];
sx q[2];
rz(0.19902755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1959343) q[1];
sx q[1];
rz(-3.0155026) q[1];
sx q[1];
rz(0.083421589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5163243) q[3];
sx q[3];
rz(-1.4843974) q[3];
sx q[3];
rz(0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(1.2329873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15105948) q[0];
sx q[0];
rz(-0.49976832) q[0];
sx q[0];
rz(-2.997056) q[0];
rz(-pi) q[1];
rz(2.1190676) q[2];
sx q[2];
rz(-1.0523445) q[2];
sx q[2];
rz(3.0419635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14703688) q[1];
sx q[1];
rz(-1.5347267) q[1];
sx q[1];
rz(-1.3922763) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2635157) q[3];
sx q[3];
rz(-1.3948166) q[3];
sx q[3];
rz(-0.39719492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91810742) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-0.89170757) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(2.9290501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5355797) q[0];
sx q[0];
rz(-1.5870464) q[0];
sx q[0];
rz(-1.5910801) q[0];
rz(-1.6264621) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(-0.84673131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37413874) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(-0.9064845) q[1];
rz(-2.0810633) q[3];
sx q[3];
rz(-0.79499309) q[3];
sx q[3];
rz(-1.0802964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(-0.5212211) q[2];
rz(1.3850348) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(-0.22512063) q[0];
rz(-1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-2.7640142) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7901944) q[0];
sx q[0];
rz(-1.5488008) q[0];
sx q[0];
rz(-1.3047332) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16212459) q[2];
sx q[2];
rz(-2.8505278) q[2];
sx q[2];
rz(-1.7578917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4695417) q[1];
sx q[1];
rz(-0.44476032) q[1];
sx q[1];
rz(-1.3252392) q[1];
rz(-0.029929786) q[3];
sx q[3];
rz(-2.3464591) q[3];
sx q[3];
rz(3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1569415) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(2.7379819) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(-0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(2.9220707) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(0.25442466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5146778) q[0];
sx q[0];
rz(-1.2509545) q[0];
sx q[0];
rz(-0.11492782) q[0];
rz(-2.1830325) q[2];
sx q[2];
rz(-1.6974291) q[2];
sx q[2];
rz(0.74778344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88528819) q[1];
sx q[1];
rz(-1.8418152) q[1];
sx q[1];
rz(-0.46896743) q[1];
x q[2];
rz(0.80612225) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1640132) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(-1.3158201) q[2];
rz(0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(-1.4550495) q[0];
rz(-0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5664068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4644509) q[0];
sx q[0];
rz(-1.8832708) q[0];
sx q[0];
rz(-2.421445) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0303866) q[2];
sx q[2];
rz(-1.724913) q[2];
sx q[2];
rz(0.87481462) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7944784) q[1];
sx q[1];
rz(-1.9273259) q[1];
sx q[1];
rz(-1.1785248) q[1];
rz(-pi) q[2];
rz(2.6353587) q[3];
sx q[3];
rz(-1.9444379) q[3];
sx q[3];
rz(-2.1831499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(2.8640462) q[2];
rz(-1.4510441) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(-1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9797416) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(1.6850527) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(-1.1368407) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6900078) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(0.33126979) q[0];
rz(-pi) q[1];
rz(1.1486263) q[2];
sx q[2];
rz(-1.2264226) q[2];
sx q[2];
rz(0.070377199) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1946823) q[1];
sx q[1];
rz(-1.2543251) q[1];
sx q[1];
rz(3.0419993) q[1];
x q[2];
rz(2.7342019) q[3];
sx q[3];
rz(-2.3630777) q[3];
sx q[3];
rz(0.22637573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4920766) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(-1.0409522) q[2];
rz(0.0020290931) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(0.94605207) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8004868) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(-0.97408803) q[0];
rz(-0.59801306) q[2];
sx q[2];
rz(-0.41053718) q[2];
sx q[2];
rz(1.41278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(-0.65995364) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53491433) q[3];
sx q[3];
rz(-1.2423008) q[3];
sx q[3];
rz(1.2319777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(2.7907794) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012797) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-2.4621261) q[2];
sx q[2];
rz(-1.1107399) q[2];
sx q[2];
rz(-2.6723292) q[2];
rz(1.894542) q[3];
sx q[3];
rz(-1.6508045) q[3];
sx q[3];
rz(1.9423021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
