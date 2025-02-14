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
rz(2.7317943) q[0];
sx q[0];
rz(8.3538342) q[0];
rz(-2.5228956) q[1];
sx q[1];
rz(-0.56801152) q[1];
sx q[1];
rz(0.84363371) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742357) q[0];
sx q[0];
rz(-1.5412742) q[0];
sx q[0];
rz(1.3062472) q[0];
x q[1];
rz(1.4272408) q[2];
sx q[2];
rz(-1.8633522) q[2];
sx q[2];
rz(-2.1168328) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1484969) q[1];
sx q[1];
rz(-1.246713) q[1];
sx q[1];
rz(-2.2739972) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0759962) q[3];
sx q[3];
rz(-1.9952979) q[3];
sx q[3];
rz(2.692846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3598651) q[2];
sx q[2];
rz(-2.3024776) q[2];
sx q[2];
rz(-3.1245933) q[2];
rz(-1.4012236) q[3];
sx q[3];
rz(-2.5798116) q[3];
sx q[3];
rz(1.2167654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.568999) q[0];
sx q[0];
rz(-1.3324791) q[0];
sx q[0];
rz(-2.3542985) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(-1.23752) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3792619) q[0];
sx q[0];
rz(-1.6325765) q[0];
sx q[0];
rz(2.558055) q[0];
x q[1];
rz(-2.8493153) q[2];
sx q[2];
rz(-0.9773796) q[2];
sx q[2];
rz(-3.127272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54531714) q[1];
sx q[1];
rz(-1.2545689) q[1];
sx q[1];
rz(1.113722) q[1];
x q[2];
rz(-2.7010553) q[3];
sx q[3];
rz(-1.3187265) q[3];
sx q[3];
rz(-1.019875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42147192) q[2];
sx q[2];
rz(-1.6933491) q[2];
sx q[2];
rz(-0.063701542) q[2];
rz(1.8675768) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(1.5773076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0434697) q[0];
sx q[0];
rz(-1.1935357) q[0];
sx q[0];
rz(-2.3980339) q[0];
rz(0.45383635) q[1];
sx q[1];
rz(-1.791714) q[1];
sx q[1];
rz(0.41890621) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1322238) q[0];
sx q[0];
rz(-1.1282578) q[0];
sx q[0];
rz(2.924463) q[0];
rz(-pi) q[1];
rz(-2.3055196) q[2];
sx q[2];
rz(-1.8030858) q[2];
sx q[2];
rz(-2.34926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.53330227) q[1];
sx q[1];
rz(-0.69606298) q[1];
sx q[1];
rz(-1.1185557) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13678582) q[3];
sx q[3];
rz(-0.89623755) q[3];
sx q[3];
rz(3.0709895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0539187) q[2];
sx q[2];
rz(-0.81834617) q[2];
sx q[2];
rz(-1.2308925) q[2];
rz(2.6026717) q[3];
sx q[3];
rz(-1.475622) q[3];
sx q[3];
rz(-1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8777799) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(-0.60212773) q[0];
rz(-1.4471794) q[1];
sx q[1];
rz(-0.68232957) q[1];
sx q[1];
rz(-0.86300659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857638) q[0];
sx q[0];
rz(-1.6842972) q[0];
sx q[0];
rz(-0.72300006) q[0];
rz(-pi) q[1];
rz(2.8406937) q[2];
sx q[2];
rz(-1.2105807) q[2];
sx q[2];
rz(0.27733251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80490404) q[1];
sx q[1];
rz(-2.4987578) q[1];
sx q[1];
rz(2.0175319) q[1];
rz(-pi) q[2];
rz(-2.3087048) q[3];
sx q[3];
rz(-1.3331659) q[3];
sx q[3];
rz(-0.72494635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.065923125) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(2.2389331) q[2];
rz(-2.3265808) q[3];
sx q[3];
rz(-2.5383526) q[3];
sx q[3];
rz(-2.7284315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0478504) q[0];
sx q[0];
rz(-3.0939026) q[0];
sx q[0];
rz(2.2064741) q[0];
rz(1.6085666) q[1];
sx q[1];
rz(-0.32564274) q[1];
sx q[1];
rz(0.26587048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1810721) q[0];
sx q[0];
rz(-1.4645394) q[0];
sx q[0];
rz(-1.8828859) q[0];
rz(0.29823096) q[2];
sx q[2];
rz(-0.77285337) q[2];
sx q[2];
rz(1.4082343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.362124) q[1];
sx q[1];
rz(-1.2219795) q[1];
sx q[1];
rz(-1.3999983) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2072125) q[3];
sx q[3];
rz(-0.28237469) q[3];
sx q[3];
rz(2.7350712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.461146) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(2.6908596) q[2];
rz(-1.9510673) q[3];
sx q[3];
rz(-2.4400986) q[3];
sx q[3];
rz(2.7019971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9615237) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(0.044483749) q[0];
rz(-2.8728409) q[1];
sx q[1];
rz(-2.2046397) q[1];
sx q[1];
rz(2.7106947) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6857032) q[0];
sx q[0];
rz(-2.1834032) q[0];
sx q[0];
rz(-1.2489178) q[0];
rz(0.90164945) q[2];
sx q[2];
rz(-1.9414273) q[2];
sx q[2];
rz(-1.2610051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1483254) q[1];
sx q[1];
rz(-1.9714515) q[1];
sx q[1];
rz(-2.2877778) q[1];
rz(-pi) q[2];
rz(-0.66029064) q[3];
sx q[3];
rz(-1.5753645) q[3];
sx q[3];
rz(-2.3167603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0253133) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(0.6244134) q[2];
rz(2.0404909) q[3];
sx q[3];
rz(-2.3291984) q[3];
sx q[3];
rz(-2.1541434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2295912) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(0.68369317) q[0];
rz(1.6992441) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(-2.9396465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.604902) q[0];
sx q[0];
rz(-2.8134007) q[0];
sx q[0];
rz(-1.5668014) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0933385) q[2];
sx q[2];
rz(-0.82403467) q[2];
sx q[2];
rz(1.6997018) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20257945) q[1];
sx q[1];
rz(-1.6375305) q[1];
sx q[1];
rz(-0.52187397) q[1];
rz(-pi) q[2];
rz(-0.45777623) q[3];
sx q[3];
rz(-1.5321863) q[3];
sx q[3];
rz(0.82070551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0553637) q[2];
sx q[2];
rz(-1.7217041) q[2];
sx q[2];
rz(1.8602547) q[2];
rz(-2.4751439) q[3];
sx q[3];
rz(-2.0446348) q[3];
sx q[3];
rz(2.5950529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619038) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(0.51396489) q[0];
rz(0.95298302) q[1];
sx q[1];
rz(-0.62478462) q[1];
sx q[1];
rz(2.0228588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361439) q[0];
sx q[0];
rz(-1.4892254) q[0];
sx q[0];
rz(2.7177627) q[0];
x q[1];
rz(-2.7114026) q[2];
sx q[2];
rz(-2.3025844) q[2];
sx q[2];
rz(-2.6553287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7756164) q[1];
sx q[1];
rz(-1.9785641) q[1];
sx q[1];
rz(2.9080176) q[1];
rz(-pi) q[2];
rz(-0.548377) q[3];
sx q[3];
rz(-1.8781239) q[3];
sx q[3];
rz(-0.69068324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.755456) q[2];
sx q[2];
rz(-2.7894207) q[2];
sx q[2];
rz(-2.238671) q[2];
rz(-1.6821945) q[3];
sx q[3];
rz(-1.3420339) q[3];
sx q[3];
rz(-2.3193147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67895401) q[0];
sx q[0];
rz(-0.32005388) q[0];
sx q[0];
rz(-1.1902887) q[0];
rz(2.8207488) q[1];
sx q[1];
rz(-1.6879993) q[1];
sx q[1];
rz(-1.920059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71305481) q[0];
sx q[0];
rz(-2.894141) q[0];
sx q[0];
rz(-2.5525981) q[0];
x q[1];
rz(2.8655445) q[2];
sx q[2];
rz(-0.85535565) q[2];
sx q[2];
rz(-0.1705585) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1561191) q[1];
sx q[1];
rz(-1.7745943) q[1];
sx q[1];
rz(-0.76038313) q[1];
rz(0.19467312) q[3];
sx q[3];
rz(-1.750573) q[3];
sx q[3];
rz(-0.67357066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90423501) q[2];
sx q[2];
rz(-2.9866437) q[2];
sx q[2];
rz(0.3332738) q[2];
rz(-1.0171558) q[3];
sx q[3];
rz(-0.95484304) q[3];
sx q[3];
rz(0.3199544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79621133) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(2.2228125) q[0];
rz(2.9504919) q[1];
sx q[1];
rz(-2.7156576) q[1];
sx q[1];
rz(-2.3474615) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0988172) q[0];
sx q[0];
rz(-0.92178994) q[0];
sx q[0];
rz(2.2521583) q[0];
rz(-pi) q[1];
rz(-0.36074071) q[2];
sx q[2];
rz(-1.7684002) q[2];
sx q[2];
rz(0.2195356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0865143) q[1];
sx q[1];
rz(-2.3740413) q[1];
sx q[1];
rz(1.961019) q[1];
rz(2.5701282) q[3];
sx q[3];
rz(-0.37701615) q[3];
sx q[3];
rz(-2.0588889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0309151) q[2];
sx q[2];
rz(-0.6157178) q[2];
sx q[2];
rz(-0.75882971) q[2];
rz(-0.87016726) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1398685) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(0.43791804) q[1];
sx q[1];
rz(-1.4358078) q[1];
sx q[1];
rz(-1.5461071) q[1];
rz(1.9867867) q[2];
sx q[2];
rz(-1.1196954) q[2];
sx q[2];
rz(-2.1449631) q[2];
rz(1.036676) q[3];
sx q[3];
rz(-0.86363367) q[3];
sx q[3];
rz(0.96994079) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
