OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.714158) q[0];
sx q[0];
rz(-2.7058869) q[0];
sx q[0];
rz(-0.92619196) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46967888) q[0];
sx q[0];
rz(-0.9979453) q[0];
sx q[0];
rz(-1.3290622) q[0];
rz(-pi) q[1];
rz(-2.9973642) q[2];
sx q[2];
rz(-0.63604522) q[2];
sx q[2];
rz(-0.47338212) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.86245) q[1];
sx q[1];
rz(-1.7596054) q[1];
sx q[1];
rz(-0.11629176) q[1];
x q[2];
rz(-2.0472237) q[3];
sx q[3];
rz(-2.8989887) q[3];
sx q[3];
rz(-2.0403595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7131876) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(0.92450809) q[2];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(-1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532042) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(0.19451441) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-0.054873437) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7354436) q[0];
sx q[0];
rz(-2.9950954) q[0];
sx q[0];
rz(2.2551401) q[0];
rz(-pi) q[1];
rz(-0.61972159) q[2];
sx q[2];
rz(-2.0267068) q[2];
sx q[2];
rz(-2.1330619) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84150746) q[1];
sx q[1];
rz(-1.9133139) q[1];
sx q[1];
rz(-3.1373346) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67316405) q[3];
sx q[3];
rz(-1.4041956) q[3];
sx q[3];
rz(1.8002321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43869552) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(1.4820209) q[2];
rz(2.1510018) q[3];
sx q[3];
rz(-1.7329268) q[3];
sx q[3];
rz(-1.3668485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(-0.24770501) q[1];
sx q[1];
rz(-1.1876371) q[1];
sx q[1];
rz(1.8050271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35555392) q[0];
sx q[0];
rz(-2.9125014) q[0];
sx q[0];
rz(-1.3130207) q[0];
rz(-pi) q[1];
rz(0.58972085) q[2];
sx q[2];
rz(-1.1955185) q[2];
sx q[2];
rz(0.42082618) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3414351) q[1];
sx q[1];
rz(-2.2812124) q[1];
sx q[1];
rz(-0.40409778) q[1];
rz(-1.9732479) q[3];
sx q[3];
rz(-0.84512049) q[3];
sx q[3];
rz(2.890051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0388564) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(-1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(1.6259441) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17722002) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(1.7376602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0899635) q[0];
sx q[0];
rz(-1.3768059) q[0];
sx q[0];
rz(0.35332638) q[0];
rz(-pi) q[1];
rz(2.3061182) q[2];
sx q[2];
rz(-2.2611141) q[2];
sx q[2];
rz(-2.8719605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85642203) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(1.1986033) q[1];
rz(-pi) q[2];
rz(-0.87938829) q[3];
sx q[3];
rz(-1.6367568) q[3];
sx q[3];
rz(-1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(-0.02031859) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(2.7385353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-2.0251944) q[0];
sx q[0];
rz(2.1257341) q[0];
rz(1.5083183) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(-3.1075081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.236892) q[0];
sx q[0];
rz(-1.1780945) q[0];
sx q[0];
rz(-1.7444201) q[0];
rz(1.3047406) q[2];
sx q[2];
rz(-1.012946) q[2];
sx q[2];
rz(1.3706052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7923813) q[1];
sx q[1];
rz(-1.5639018) q[1];
sx q[1];
rz(-1.8551808) q[1];
rz(-pi) q[2];
rz(-2.6328153) q[3];
sx q[3];
rz(-1.0049337) q[3];
sx q[3];
rz(2.1600427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(1.1711228) q[2];
rz(1.3327538) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(-0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4532582) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(2.7057498) q[0];
rz(1.1054989) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(-1.9357392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642731) q[0];
sx q[0];
rz(-0.84007971) q[0];
sx q[0];
rz(-2.8217836) q[0];
rz(3.1255683) q[2];
sx q[2];
rz(-1.6659701) q[2];
sx q[2];
rz(-1.3118088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0059433) q[1];
sx q[1];
rz(-1.597153) q[1];
sx q[1];
rz(-0.51075682) q[1];
rz(-pi) q[2];
x q[2];
rz(0.033127012) q[3];
sx q[3];
rz(-1.9702385) q[3];
sx q[3];
rz(-1.2101733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-3.0899866) q[2];
rz(-0.40766019) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(-0.78397059) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(-2.6046682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85978956) q[0];
sx q[0];
rz(-1.4872695) q[0];
sx q[0];
rz(-3.0755755) q[0];
rz(-pi) q[1];
rz(1.827048) q[2];
sx q[2];
rz(-1.429024) q[2];
sx q[2];
rz(1.8585376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.19837025) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(-0.13065773) q[1];
rz(-pi) q[2];
rz(-0.22646871) q[3];
sx q[3];
rz(-0.60301757) q[3];
sx q[3];
rz(-0.23803593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15988222) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(2.8602709) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(-2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0328338) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-2.7539745) q[0];
rz(3.0265813) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.176929) q[0];
sx q[0];
rz(-1.7690072) q[0];
sx q[0];
rz(2.6508509) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6667716) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(2.1542187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3927314) q[1];
sx q[1];
rz(-1.0034605) q[1];
sx q[1];
rz(-2.9859403) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5070595) q[3];
sx q[3];
rz(-2.3252441) q[3];
sx q[3];
rz(-0.07894978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9541624) q[2];
sx q[2];
rz(-2.607441) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(0.77504843) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.9545814) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(-2.9206081) q[0];
rz(2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(-1.0029213) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.040398) q[0];
sx q[0];
rz(-1.4324491) q[0];
sx q[0];
rz(0.25427108) q[0];
rz(-0.9858176) q[2];
sx q[2];
rz(-2.0009811) q[2];
sx q[2];
rz(-0.9790203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79473125) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(0.55333432) q[1];
rz(-1.8995729) q[3];
sx q[3];
rz(-2.4283096) q[3];
sx q[3];
rz(-0.79771358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(-2.6878099) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(2.2973072) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(2.9504543) q[0];
rz(0.29516164) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(0.89231649) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494209) q[0];
sx q[0];
rz(-0.76602174) q[0];
sx q[0];
rz(-1.8269405) q[0];
rz(-pi) q[1];
rz(-0.92992444) q[2];
sx q[2];
rz(-1.495549) q[2];
sx q[2];
rz(1.6844695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.031242) q[1];
sx q[1];
rz(-2.894245) q[1];
sx q[1];
rz(-2.8691643) q[1];
rz(-pi) q[2];
rz(-0.43182208) q[3];
sx q[3];
rz(-1.8027657) q[3];
sx q[3];
rz(2.1524129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(0.85956335) q[2];
rz(1.2236979) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201465) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(1.6620811) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.66172) q[2];
sx q[2];
rz(-2.8488013) q[2];
sx q[2];
rz(-1.6691096) q[2];
rz(1.0449833) q[3];
sx q[3];
rz(-0.079418728) q[3];
sx q[3];
rz(3.1374745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
