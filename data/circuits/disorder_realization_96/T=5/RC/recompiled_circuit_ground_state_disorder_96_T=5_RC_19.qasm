OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95028967) q[0];
sx q[0];
rz(6.0130881) q[0];
sx q[0];
rz(10.313378) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(4.7409952) q[1];
sx q[1];
rz(10.805605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7021877) q[0];
sx q[0];
rz(-0.26964615) q[0];
sx q[0];
rz(0.7650956) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51989748) q[2];
sx q[2];
rz(-2.2714104) q[2];
sx q[2];
rz(-1.1288647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9819698) q[1];
sx q[1];
rz(-1.2352984) q[1];
sx q[1];
rz(-1.5606176) q[1];
rz(1.2065776) q[3];
sx q[3];
rz(-1.3602339) q[3];
sx q[3];
rz(2.5488146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8218653) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-2.3925609) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-2.1059683) q[3];
sx q[3];
rz(-3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6946436) q[0];
sx q[0];
rz(-1.6634989) q[0];
sx q[0];
rz(-1.9248167) q[0];
rz(-2.079839) q[1];
sx q[1];
rz(-2.3371425) q[1];
sx q[1];
rz(0.42713508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0674853) q[0];
sx q[0];
rz(-0.68746131) q[0];
sx q[0];
rz(-0.56337728) q[0];
x q[1];
rz(-2.2160276) q[2];
sx q[2];
rz(-1.3750018) q[2];
sx q[2];
rz(0.93709556) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6500351) q[1];
sx q[1];
rz(-1.8582398) q[1];
sx q[1];
rz(0.14658714) q[1];
rz(-1.7747545) q[3];
sx q[3];
rz(-1.5195623) q[3];
sx q[3];
rz(0.27328396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(-1.3986577) q[2];
rz(-0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(-1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1200714) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(1.1687763) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(-2.5659836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2164477) q[0];
sx q[0];
rz(-2.5024104) q[0];
sx q[0];
rz(3.1228288) q[0];
rz(-pi) q[1];
rz(2.2314638) q[2];
sx q[2];
rz(-1.5572539) q[2];
sx q[2];
rz(-2.8000268) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82390416) q[1];
sx q[1];
rz(-2.1145193) q[1];
sx q[1];
rz(0.25441092) q[1];
x q[2];
rz(-0.024954114) q[3];
sx q[3];
rz(-0.91812274) q[3];
sx q[3];
rz(-0.11634532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1316954) q[2];
sx q[2];
rz(-2.3991149) q[2];
sx q[2];
rz(2.1843145) q[2];
rz(0.57473985) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9855373) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(-1.2611457) q[0];
rz(-0.082322923) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(1.7877158) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.533723) q[0];
sx q[0];
rz(-1.3114616) q[0];
sx q[0];
rz(1.3277256) q[0];
rz(-pi) q[1];
rz(-0.33892314) q[2];
sx q[2];
rz(-0.51033516) q[2];
sx q[2];
rz(1.5938544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9685843) q[1];
sx q[1];
rz(-0.49793303) q[1];
sx q[1];
rz(3.0281316) q[1];
rz(-pi) q[2];
rz(2.1534377) q[3];
sx q[3];
rz(-2.372962) q[3];
sx q[3];
rz(-2.5386794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(-1.2565695) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.310815) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-0.34314439) q[0];
rz(0.061773069) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.7247346) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34300229) q[0];
sx q[0];
rz(-1.2932284) q[0];
sx q[0];
rz(2.9647602) q[0];
rz(2.4409962) q[2];
sx q[2];
rz(-1.8237517) q[2];
sx q[2];
rz(-2.7299936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38450634) q[1];
sx q[1];
rz(-2.050274) q[1];
sx q[1];
rz(-2.9887385) q[1];
x q[2];
rz(2.9067578) q[3];
sx q[3];
rz(-1.1560625) q[3];
sx q[3];
rz(-2.107634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(-1.0859547) q[2];
rz(-0.91935277) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(2.5277188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3980961) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(0.79291517) q[0];
rz(2.4010557) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(-0.83121306) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34666592) q[0];
sx q[0];
rz(-1.7788609) q[0];
sx q[0];
rz(-1.9545593) q[0];
x q[1];
rz(-1.5446051) q[2];
sx q[2];
rz(-2.8644343) q[2];
sx q[2];
rz(-1.8089393) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86262396) q[1];
sx q[1];
rz(-1.875307) q[1];
sx q[1];
rz(0.83359756) q[1];
rz(-1.5300203) q[3];
sx q[3];
rz(-0.47104657) q[3];
sx q[3];
rz(0.99179964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0703766) q[2];
sx q[2];
rz(-1.8698317) q[2];
sx q[2];
rz(1.1191204) q[2];
rz(1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041466) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(1.5995837) q[0];
rz(-1.0150602) q[1];
sx q[1];
rz(-1.5856182) q[1];
sx q[1];
rz(-0.17280811) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5046523) q[0];
sx q[0];
rz(-1.0811792) q[0];
sx q[0];
rz(2.2235653) q[0];
x q[1];
rz(-1.9394933) q[2];
sx q[2];
rz(-2.2208344) q[2];
sx q[2];
rz(1.7391313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0513251) q[1];
sx q[1];
rz(-2.6410612) q[1];
sx q[1];
rz(2.1085018) q[1];
rz(-1.0181581) q[3];
sx q[3];
rz(-2.2595836) q[3];
sx q[3];
rz(-1.3017201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66199866) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(2.5665723) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(-1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8136895) q[0];
sx q[0];
rz(-2.8459025) q[0];
sx q[0];
rz(3.1258702) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(1.9421008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48391438) q[0];
sx q[0];
rz(-1.8059314) q[0];
sx q[0];
rz(-2.0245309) q[0];
rz(-pi) q[1];
rz(-1.6686317) q[2];
sx q[2];
rz(-0.49287686) q[2];
sx q[2];
rz(0.6334444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6024692) q[1];
sx q[1];
rz(-1.1014465) q[1];
sx q[1];
rz(0.55901171) q[1];
x q[2];
rz(-2.1425072) q[3];
sx q[3];
rz(-1.8322819) q[3];
sx q[3];
rz(0.7657683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1477995) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(1.3851059) q[2];
rz(-0.6238474) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(-1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856336) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(-0.69325915) q[0];
rz(-0.26501003) q[1];
sx q[1];
rz(-2.3131504) q[1];
sx q[1];
rz(1.6185435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361621) q[0];
sx q[0];
rz(-2.3443065) q[0];
sx q[0];
rz(1.5553586) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9829282) q[2];
sx q[2];
rz(-1.0318021) q[2];
sx q[2];
rz(-1.7807775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0761145) q[1];
sx q[1];
rz(-1.5947184) q[1];
sx q[1];
rz(-0.33965276) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6783488) q[3];
sx q[3];
rz(-2.5578024) q[3];
sx q[3];
rz(-1.3046169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(-1.4576853) q[2];
rz(0.95528209) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.38462287) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(0.48267522) q[0];
rz(-0.92974281) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(0.39628705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1209377) q[0];
sx q[0];
rz(-1.9299986) q[0];
sx q[0];
rz(0.44393702) q[0];
rz(-pi) q[1];
rz(2.9626289) q[2];
sx q[2];
rz(-1.4216627) q[2];
sx q[2];
rz(-2.3863132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28010269) q[1];
sx q[1];
rz(-0.74833732) q[1];
sx q[1];
rz(1.8381717) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72209218) q[3];
sx q[3];
rz(-0.51755899) q[3];
sx q[3];
rz(0.90921569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3489909) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(0.3271884) q[2];
rz(0.78091019) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(-1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1031716) q[0];
sx q[0];
rz(-0.99089834) q[0];
sx q[0];
rz(0.25767576) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(1.5255899) q[2];
sx q[2];
rz(-2.2635985) q[2];
sx q[2];
rz(0.23512693) q[2];
rz(-0.38402186) q[3];
sx q[3];
rz(-2.1441318) q[3];
sx q[3];
rz(-1.287788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
