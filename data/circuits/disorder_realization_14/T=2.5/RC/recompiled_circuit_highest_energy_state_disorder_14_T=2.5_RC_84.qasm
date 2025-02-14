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
rz(-0.7402339) q[0];
sx q[0];
rz(-1.482168) q[0];
sx q[0];
rz(2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(-0.60751539) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4273194) q[0];
sx q[0];
rz(-1.3392942) q[0];
sx q[0];
rz(2.0718859) q[0];
rz(1.7781177) q[2];
sx q[2];
rz(-2.6229515) q[2];
sx q[2];
rz(2.867709) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.059587) q[1];
sx q[1];
rz(-1.4114037) q[1];
sx q[1];
rz(0.48377796) q[1];
x q[2];
rz(2.0161765) q[3];
sx q[3];
rz(-1.9102671) q[3];
sx q[3];
rz(-0.31472963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6174378) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(0.080862008) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(-1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2422159) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(2.2826165) q[0];
rz(1.2902749) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(1.7346409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1340356) q[0];
sx q[0];
rz(-1.7985744) q[0];
sx q[0];
rz(2.0640949) q[0];
x q[1];
rz(-0.036769899) q[2];
sx q[2];
rz(-1.9027793) q[2];
sx q[2];
rz(-0.22184243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21995658) q[1];
sx q[1];
rz(-2.0825279) q[1];
sx q[1];
rz(0.10351609) q[1];
rz(0.59515335) q[3];
sx q[3];
rz(-0.66616026) q[3];
sx q[3];
rz(-0.43597886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0824288) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(-1.3272237) q[2];
rz(-1.0446769) q[3];
sx q[3];
rz(-0.71377126) q[3];
sx q[3];
rz(2.7248342) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5663719) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(2.7161993) q[0];
rz(-1.3771903) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(-1.4345217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4017556) q[0];
sx q[0];
rz(-0.65261894) q[0];
sx q[0];
rz(-0.1952066) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.038952053) q[2];
sx q[2];
rz(-2.4329081) q[2];
sx q[2];
rz(-1.5295636) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4340448) q[1];
sx q[1];
rz(-2.2128989) q[1];
sx q[1];
rz(-1.0471859) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1528963) q[3];
sx q[3];
rz(-2.4800081) q[3];
sx q[3];
rz(1.0583056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(1.2909935) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62035471) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(2.3714016) q[0];
rz(2.1417446) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(1.8005449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65085852) q[0];
sx q[0];
rz(-0.58877173) q[0];
sx q[0];
rz(-2.7522699) q[0];
x q[1];
rz(-2.2571428) q[2];
sx q[2];
rz(-1.9362861) q[2];
sx q[2];
rz(-0.7892424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.136678) q[1];
sx q[1];
rz(-1.305849) q[1];
sx q[1];
rz(2.8317004) q[1];
rz(-pi) q[2];
rz(-1.8561268) q[3];
sx q[3];
rz(-1.4754681) q[3];
sx q[3];
rz(0.25432977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3406713) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(-3.0411804) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95555821) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(1.7452128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3811548) q[0];
sx q[0];
rz(-1.210307) q[0];
sx q[0];
rz(1.1973778) q[0];
rz(1.7368083) q[2];
sx q[2];
rz(-1.1653656) q[2];
sx q[2];
rz(-2.977598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.760031) q[1];
sx q[1];
rz(-2.2027822) q[1];
sx q[1];
rz(1.6595592) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5377858) q[3];
sx q[3];
rz(-0.74652687) q[3];
sx q[3];
rz(-1.3804264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9153626) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(3.0653595) q[2];
rz(1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(1.7414198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48080322) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.3355108) q[0];
rz(2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(-2.9551771) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6680697) q[0];
sx q[0];
rz(-2.0041564) q[0];
sx q[0];
rz(2.1587203) q[0];
rz(-pi) q[1];
rz(-0.30419402) q[2];
sx q[2];
rz(-0.59202164) q[2];
sx q[2];
rz(1.6006058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0708145) q[1];
sx q[1];
rz(-2.1970587) q[1];
sx q[1];
rz(2.8423944) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0746434) q[3];
sx q[3];
rz(-2.0130139) q[3];
sx q[3];
rz(2.7520455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4621801) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(1.3235486) q[2];
rz(-1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(2.2511258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776176) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(0.63419813) q[0];
rz(2.0967261) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(-0.96010906) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2506367) q[0];
sx q[0];
rz(-1.5433558) q[0];
sx q[0];
rz(1.9695884) q[0];
rz(-1.5627925) q[2];
sx q[2];
rz(-1.2697313) q[2];
sx q[2];
rz(2.0640399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83852488) q[1];
sx q[1];
rz(-2.6160598) q[1];
sx q[1];
rz(1.5477033) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84585692) q[3];
sx q[3];
rz(-2.1921033) q[3];
sx q[3];
rz(1.814807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94379696) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(-1.724285) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(-1.5203169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.5328131) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(-2.1667495) q[0];
rz(0.10803647) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.1688165) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10442142) q[0];
sx q[0];
rz(-1.9871431) q[0];
sx q[0];
rz(-2.092931) q[0];
rz(-2.37974) q[2];
sx q[2];
rz(-1.5708062) q[2];
sx q[2];
rz(-0.61864432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7119693) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(-1.8046677) q[1];
rz(0.036691477) q[3];
sx q[3];
rz(-2.6495669) q[3];
sx q[3];
rz(1.8560611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.03269) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(-2.4857793) q[2];
rz(0.10410318) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(0.18812215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.6498097) q[0];
rz(1.0393633) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(2.6731491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85850785) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(-1.5751189) q[0];
rz(-pi) q[1];
rz(-0.55267398) q[2];
sx q[2];
rz(-2.6330593) q[2];
sx q[2];
rz(-2.7486211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5891582) q[1];
sx q[1];
rz(-0.8048519) q[1];
sx q[1];
rz(2.9167152) q[1];
rz(-pi) q[2];
rz(-0.78448589) q[3];
sx q[3];
rz(-2.753559) q[3];
sx q[3];
rz(-1.3875543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.046772) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(0.92791933) q[2];
rz(1.8038484) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(-1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054166404) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(2.0637276) q[0];
rz(-0.36733356) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7736762) q[0];
sx q[0];
rz(-1.319029) q[0];
sx q[0];
rz(-1.0836224) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4201035) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(-0.53021741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91032797) q[1];
sx q[1];
rz(-0.18016768) q[1];
sx q[1];
rz(0.50039165) q[1];
rz(-pi) q[2];
rz(-1.359109) q[3];
sx q[3];
rz(-0.84722391) q[3];
sx q[3];
rz(-2.7016957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1401356) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(0.4293116) q[2];
rz(2.9863827) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(-1.8163053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406381) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(0.86391972) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(-0.39045329) q[2];
sx q[2];
rz(-0.57885546) q[2];
sx q[2];
rz(-0.59138966) q[2];
rz(0.33959099) q[3];
sx q[3];
rz(-2.17008) q[3];
sx q[3];
rz(-1.9934987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
