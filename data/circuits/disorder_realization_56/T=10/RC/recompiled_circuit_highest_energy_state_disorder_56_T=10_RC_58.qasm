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
rz(2.7409878) q[0];
sx q[0];
rz(-0.40979835) q[0];
sx q[0];
rz(1.0709437) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(-2.5735811) q[1];
sx q[1];
rz(-0.84363371) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5045568) q[0];
sx q[0];
rz(-1.3063653) q[0];
sx q[0];
rz(0.030585551) q[0];
x q[1];
rz(0.29542342) q[2];
sx q[2];
rz(-1.43338) q[2];
sx q[2];
rz(0.58770056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8275717) q[1];
sx q[1];
rz(-2.2306475) q[1];
sx q[1];
rz(-0.41484264) q[1];
rz(-pi) q[2];
rz(2.0655965) q[3];
sx q[3];
rz(-1.1462948) q[3];
sx q[3];
rz(-2.692846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7817276) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(-3.1245933) q[2];
rz(-1.7403691) q[3];
sx q[3];
rz(-0.56178105) q[3];
sx q[3];
rz(-1.9248272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725937) q[0];
sx q[0];
rz(-1.3324791) q[0];
sx q[0];
rz(2.3542985) q[0];
rz(-0.17838082) q[1];
sx q[1];
rz(-1.2063824) q[1];
sx q[1];
rz(-1.23752) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3792619) q[0];
sx q[0];
rz(-1.6325765) q[0];
sx q[0];
rz(0.5835377) q[0];
rz(-pi) q[1];
rz(2.1844615) q[2];
sx q[2];
rz(-1.8120013) q[2];
sx q[2];
rz(-1.7231736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2678821) q[1];
sx q[1];
rz(-1.137966) q[1];
sx q[1];
rz(-0.3496561) q[1];
rz(2.7010553) q[3];
sx q[3];
rz(-1.8228662) q[3];
sx q[3];
rz(-1.019875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7201207) q[2];
sx q[2];
rz(-1.4482435) q[2];
sx q[2];
rz(-0.063701542) q[2];
rz(-1.2740159) q[3];
sx q[3];
rz(-0.84309045) q[3];
sx q[3];
rz(1.564285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0434697) q[0];
sx q[0];
rz(-1.948057) q[0];
sx q[0];
rz(-2.3980339) q[0];
rz(2.6877563) q[1];
sx q[1];
rz(-1.3498787) q[1];
sx q[1];
rz(0.41890621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6565257) q[0];
sx q[0];
rz(-0.48975268) q[0];
sx q[0];
rz(-1.144125) q[0];
rz(-pi) q[1];
rz(-2.3055196) q[2];
sx q[2];
rz(-1.8030858) q[2];
sx q[2];
rz(0.79233263) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53330227) q[1];
sx q[1];
rz(-0.69606298) q[1];
sx q[1];
rz(2.0230369) q[1];
rz(-0.13678582) q[3];
sx q[3];
rz(-2.2453551) q[3];
sx q[3];
rz(-3.0709895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.087674) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(1.2308925) q[2];
rz(0.538921) q[3];
sx q[3];
rz(-1.475622) q[3];
sx q[3];
rz(1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26381275) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(-0.60212773) q[0];
rz(1.6944132) q[1];
sx q[1];
rz(-0.68232957) q[1];
sx q[1];
rz(-0.86300659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569538) q[0];
sx q[0];
rz(-0.85345972) q[0];
sx q[0];
rz(1.4199281) q[0];
rz(-pi) q[1];
rz(2.237488) q[2];
sx q[2];
rz(-0.46510425) q[2];
sx q[2];
rz(-2.1423774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80490404) q[1];
sx q[1];
rz(-2.4987578) q[1];
sx q[1];
rz(2.0175319) q[1];
rz(-pi) q[2];
rz(-0.3163655) q[3];
sx q[3];
rz(-2.2834407) q[3];
sx q[3];
rz(-2.0848839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0756695) q[2];
sx q[2];
rz(-1.9663591) q[2];
sx q[2];
rz(0.90265957) q[2];
rz(-2.3265808) q[3];
sx q[3];
rz(-2.5383526) q[3];
sx q[3];
rz(-2.7284315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.0937423) q[0];
sx q[0];
rz(-3.0939026) q[0];
sx q[0];
rz(0.93511859) q[0];
rz(-1.6085666) q[1];
sx q[1];
rz(-2.8159499) q[1];
sx q[1];
rz(0.26587048) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64447999) q[0];
sx q[0];
rz(-1.2605259) q[0];
sx q[0];
rz(3.029986) q[0];
rz(-pi) q[1];
rz(-1.8498672) q[2];
sx q[2];
rz(-0.84013591) q[2];
sx q[2];
rz(-2.1389291) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8740377) q[1];
sx q[1];
rz(-1.7312164) q[1];
sx q[1];
rz(0.35355132) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3060027) q[3];
sx q[3];
rz(-1.6700498) q[3];
sx q[3];
rz(-0.81391993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.461146) q[2];
sx q[2];
rz(-1.6411883) q[2];
sx q[2];
rz(-0.45073304) q[2];
rz(-1.1905253) q[3];
sx q[3];
rz(-2.4400986) q[3];
sx q[3];
rz(-2.7019971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9615237) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(-3.0971089) q[0];
rz(2.8728409) q[1];
sx q[1];
rz(-2.2046397) q[1];
sx q[1];
rz(0.43089795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92545729) q[0];
sx q[0];
rz(-1.8325984) q[0];
sx q[0];
rz(-2.5039704) q[0];
rz(-pi) q[1];
x q[1];
rz(0.459983) q[2];
sx q[2];
rz(-0.95427931) q[2];
sx q[2];
rz(3.1107855) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0468586) q[1];
sx q[1];
rz(-0.92093912) q[1];
sx q[1];
rz(0.51194329) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5765801) q[3];
sx q[3];
rz(-0.91051379) q[3];
sx q[3];
rz(0.74951142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1162794) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(-0.6244134) q[2];
rz(-2.0404909) q[3];
sx q[3];
rz(-2.3291984) q[3];
sx q[3];
rz(-0.98744923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.2295912) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(-2.4578995) q[0];
rz(-1.6992441) q[1];
sx q[1];
rz(-0.78386274) q[1];
sx q[1];
rz(-2.9396465) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037887427) q[0];
sx q[0];
rz(-1.572084) q[0];
sx q[0];
rz(-1.2426069) q[0];
x q[1];
rz(1.5187289) q[2];
sx q[2];
rz(-2.3935742) q[2];
sx q[2];
rz(-1.6287273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8117042) q[1];
sx q[1];
rz(-2.0913908) q[1];
sx q[1];
rz(1.4938526) q[1];
rz(0.087183909) q[3];
sx q[3];
rz(-0.4592866) q[3];
sx q[3];
rz(0.82827866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0553637) q[2];
sx q[2];
rz(-1.4198885) q[2];
sx q[2];
rz(-1.8602547) q[2];
rz(-0.66644871) q[3];
sx q[3];
rz(-2.0446348) q[3];
sx q[3];
rz(-2.5950529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619038) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(-2.6276278) q[0];
rz(2.1886096) q[1];
sx q[1];
rz(-0.62478462) q[1];
sx q[1];
rz(1.1187339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5394997) q[0];
sx q[0];
rz(-1.148466) q[0];
sx q[0];
rz(-1.6602449) q[0];
rz(-0.79134361) q[2];
sx q[2];
rz(-1.2553167) q[2];
sx q[2];
rz(2.3545654) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1744078) q[1];
sx q[1];
rz(-0.46665114) q[1];
sx q[1];
rz(-2.0626759) q[1];
rz(-pi) q[2];
rz(1.9268613) q[3];
sx q[3];
rz(-2.0907932) q[3];
sx q[3];
rz(0.69737747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.755456) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(2.238671) q[2];
rz(1.4593982) q[3];
sx q[3];
rz(-1.7995588) q[3];
sx q[3];
rz(-0.82227796) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67895401) q[0];
sx q[0];
rz(-0.32005388) q[0];
sx q[0];
rz(1.1902887) q[0];
rz(2.8207488) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(1.920059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0319043) q[0];
sx q[0];
rz(-1.3656989) q[0];
sx q[0];
rz(1.7102276) q[0];
x q[1];
rz(1.2668816) q[2];
sx q[2];
rz(-0.7579782) q[2];
sx q[2];
rz(-0.2371108) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62432837) q[1];
sx q[1];
rz(-0.78189865) q[1];
sx q[1];
rz(-2.8502591) q[1];
x q[2];
rz(-1.3876347) q[3];
sx q[3];
rz(-1.3792999) q[3];
sx q[3];
rz(-2.2091256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(-0.3332738) q[2];
rz(2.1244369) q[3];
sx q[3];
rz(-0.95484304) q[3];
sx q[3];
rz(-2.8216383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3453813) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(-0.91878015) q[0];
rz(-0.19110075) q[1];
sx q[1];
rz(-2.7156576) q[1];
sx q[1];
rz(0.79413116) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0988172) q[0];
sx q[0];
rz(-2.2198027) q[0];
sx q[0];
rz(-2.2521583) q[0];
x q[1];
rz(1.7816069) q[2];
sx q[2];
rz(-1.9242058) q[2];
sx q[2];
rz(-1.4251874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9136466) q[1];
sx q[1];
rz(-1.3034879) q[1];
sx q[1];
rz(2.2993907) q[1];
rz(-pi) q[2];
rz(1.3598241) q[3];
sx q[3];
rz(-1.2559685) q[3];
sx q[3];
rz(-1.6877768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0309151) q[2];
sx q[2];
rz(-2.5258749) q[2];
sx q[2];
rz(-0.75882971) q[2];
rz(0.87016726) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-2.0738535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0017241521) q[0];
sx q[0];
rz(-1.682946) q[0];
sx q[0];
rz(-1.2276822) q[0];
rz(-2.7036746) q[1];
sx q[1];
rz(-1.4358078) q[1];
sx q[1];
rz(-1.5461071) q[1];
rz(-1.9867867) q[2];
sx q[2];
rz(-2.0218973) q[2];
sx q[2];
rz(0.99662957) q[2];
rz(-2.1049166) q[3];
sx q[3];
rz(-0.86363367) q[3];
sx q[3];
rz(0.96994079) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
