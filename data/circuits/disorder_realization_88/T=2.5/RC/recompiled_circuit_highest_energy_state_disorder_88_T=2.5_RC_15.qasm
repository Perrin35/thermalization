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
rz(-3.1336194) q[0];
sx q[0];
rz(-2.70533) q[0];
sx q[0];
rz(1.0101969) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(-2.3354524) q[1];
sx q[1];
rz(-1.0239209) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181241) q[0];
sx q[0];
rz(-1.082927) q[0];
sx q[0];
rz(0.0077631891) q[0];
rz(-pi) q[1];
rz(1.802581) q[2];
sx q[2];
rz(-1.1884153) q[2];
sx q[2];
rz(-0.2222375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4657106) q[1];
sx q[1];
rz(-0.16168586) q[1];
sx q[1];
rz(1.449701) q[1];
x q[2];
rz(2.8918582) q[3];
sx q[3];
rz(-1.6563376) q[3];
sx q[3];
rz(-1.5126603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94850928) q[2];
sx q[2];
rz(-2.812959) q[2];
sx q[2];
rz(-0.23552093) q[2];
rz(-2.113302) q[3];
sx q[3];
rz(-1.8833269) q[3];
sx q[3];
rz(1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4368206) q[0];
sx q[0];
rz(-1.7828159) q[0];
sx q[0];
rz(-0.92192465) q[0];
rz(-3.0251265) q[1];
sx q[1];
rz(-1.7200229) q[1];
sx q[1];
rz(-0.88752735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5796611) q[0];
sx q[0];
rz(-1.4401739) q[0];
sx q[0];
rz(-1.3829253) q[0];
rz(-0.15917425) q[2];
sx q[2];
rz(-0.35458699) q[2];
sx q[2];
rz(-1.2693506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.973399) q[1];
sx q[1];
rz(-2.138294) q[1];
sx q[1];
rz(-1.4367815) q[1];
x q[2];
rz(-3.001962) q[3];
sx q[3];
rz(-0.88362304) q[3];
sx q[3];
rz(2.4202514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55887115) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(-2.9026046) q[2];
rz(2.41921) q[3];
sx q[3];
rz(-0.84010774) q[3];
sx q[3];
rz(-1.1649789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-2.8257985) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(-2.353299) q[0];
rz(-1.7514508) q[1];
sx q[1];
rz(-2.7056521) q[1];
sx q[1];
rz(-1.699126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9787131) q[0];
sx q[0];
rz(-1.7121234) q[0];
sx q[0];
rz(-0.070959758) q[0];
x q[1];
rz(-2.2295775) q[2];
sx q[2];
rz(-0.46123966) q[2];
sx q[2];
rz(-1.3818503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6953814) q[1];
sx q[1];
rz(-1.3331116) q[1];
sx q[1];
rz(2.2878245) q[1];
rz(-pi) q[2];
rz(-1.0679108) q[3];
sx q[3];
rz(-1.7278132) q[3];
sx q[3];
rz(-2.5105421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4060789) q[2];
sx q[2];
rz(-1.1512076) q[2];
sx q[2];
rz(-2.8847412) q[2];
rz(2.2779321) q[3];
sx q[3];
rz(-2.05859) q[3];
sx q[3];
rz(0.5736205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470806) q[0];
sx q[0];
rz(-0.032945078) q[0];
sx q[0];
rz(1.6324014) q[0];
rz(-0.31653658) q[1];
sx q[1];
rz(-1.5959975) q[1];
sx q[1];
rz(-2.1260156) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94293649) q[0];
sx q[0];
rz(-1.1186677) q[0];
sx q[0];
rz(-2.6200635) q[0];
rz(-2.8695753) q[2];
sx q[2];
rz(-2.3292856) q[2];
sx q[2];
rz(-2.7886632) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9719661) q[1];
sx q[1];
rz(-1.3463273) q[1];
sx q[1];
rz(1.9588934) q[1];
rz(-pi) q[2];
rz(1.7523132) q[3];
sx q[3];
rz(-2.3861775) q[3];
sx q[3];
rz(-2.8528573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29869002) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(1.451937) q[2];
rz(-1.2339833) q[3];
sx q[3];
rz(-1.0943639) q[3];
sx q[3];
rz(2.2598677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3877617) q[0];
sx q[0];
rz(-1.8330638) q[0];
sx q[0];
rz(2.0552788) q[0];
rz(1.7408675) q[1];
sx q[1];
rz(-2.3298405) q[1];
sx q[1];
rz(-3.1382255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28890507) q[0];
sx q[0];
rz(-2.5470535) q[0];
sx q[0];
rz(-3.1311228) q[0];
x q[1];
rz(0.98262991) q[2];
sx q[2];
rz(-1.8525043) q[2];
sx q[2];
rz(2.1717193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5382522) q[1];
sx q[1];
rz(-1.0639186) q[1];
sx q[1];
rz(1.5065306) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60008757) q[3];
sx q[3];
rz(-1.5494124) q[3];
sx q[3];
rz(-2.1789973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1668261) q[2];
sx q[2];
rz(-0.70009118) q[2];
sx q[2];
rz(-0.33494803) q[2];
rz(2.8454928) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(3.0193442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.629338) q[0];
sx q[0];
rz(-0.4991931) q[0];
sx q[0];
rz(0.41803023) q[0];
rz(-0.66608518) q[1];
sx q[1];
rz(-1.8981551) q[1];
sx q[1];
rz(-2.8935249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2193766) q[0];
sx q[0];
rz(-1.2398949) q[0];
sx q[0];
rz(-1.0537345) q[0];
rz(-1.6841959) q[2];
sx q[2];
rz(-2.6893986) q[2];
sx q[2];
rz(-0.8314641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29344968) q[1];
sx q[1];
rz(-2.0781233) q[1];
sx q[1];
rz(0.44815843) q[1];
rz(2.2002688) q[3];
sx q[3];
rz(-1.3622614) q[3];
sx q[3];
rz(2.1311614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4416113) q[2];
sx q[2];
rz(-1.2629843) q[2];
sx q[2];
rz(1.9786394) q[2];
rz(0.59398389) q[3];
sx q[3];
rz(-1.6006399) q[3];
sx q[3];
rz(2.1514413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6043337) q[0];
sx q[0];
rz(-0.27892497) q[0];
sx q[0];
rz(3.0766686) q[0];
rz(-2.7760432) q[1];
sx q[1];
rz(-1.431501) q[1];
sx q[1];
rz(-0.72986594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.176031) q[0];
sx q[0];
rz(-0.65098) q[0];
sx q[0];
rz(-1.3887547) q[0];
rz(2.3785602) q[2];
sx q[2];
rz(-1.4307466) q[2];
sx q[2];
rz(-2.829388) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5090885) q[1];
sx q[1];
rz(-1.6659686) q[1];
sx q[1];
rz(-1.0126963) q[1];
rz(-pi) q[2];
rz(-0.49298005) q[3];
sx q[3];
rz(-0.74876174) q[3];
sx q[3];
rz(1.6822588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87021747) q[2];
sx q[2];
rz(-1.7729746) q[2];
sx q[2];
rz(-2.0591056) q[2];
rz(-1.6535951) q[3];
sx q[3];
rz(-2.7797785) q[3];
sx q[3];
rz(-0.35058072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38700405) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(-0.3311232) q[0];
rz(1.4777199) q[1];
sx q[1];
rz(-0.96550566) q[1];
sx q[1];
rz(-0.32041916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3715558) q[0];
sx q[0];
rz(-1.5450243) q[0];
sx q[0];
rz(0.098866247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4244979) q[2];
sx q[2];
rz(-1.4182036) q[2];
sx q[2];
rz(-0.998978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81605212) q[1];
sx q[1];
rz(-0.42652247) q[1];
sx q[1];
rz(0.037530516) q[1];
rz(-pi) q[2];
rz(2.255104) q[3];
sx q[3];
rz(-1.6152665) q[3];
sx q[3];
rz(-1.1680166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7122345) q[2];
sx q[2];
rz(-1.7244491) q[2];
sx q[2];
rz(-0.36163914) q[2];
rz(0.9137736) q[3];
sx q[3];
rz(-0.29699609) q[3];
sx q[3];
rz(0.44338068) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906032) q[0];
sx q[0];
rz(-2.3539703) q[0];
sx q[0];
rz(0.11601624) q[0];
rz(-1.8138255) q[1];
sx q[1];
rz(-1.1632183) q[1];
sx q[1];
rz(1.9414925) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93107779) q[0];
sx q[0];
rz(-1.3793945) q[0];
sx q[0];
rz(-1.5899508) q[0];
rz(0.16559439) q[2];
sx q[2];
rz(-0.8919009) q[2];
sx q[2];
rz(2.0904581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1397651) q[1];
sx q[1];
rz(-1.5854302) q[1];
sx q[1];
rz(-2.57994) q[1];
rz(-pi) q[2];
rz(-0.96941822) q[3];
sx q[3];
rz(-0.28801258) q[3];
sx q[3];
rz(1.7565816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3001331) q[2];
sx q[2];
rz(-1.7068784) q[2];
sx q[2];
rz(2.8099698) q[2];
rz(-0.2837818) q[3];
sx q[3];
rz(-1.1016568) q[3];
sx q[3];
rz(0.42154977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83061853) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(0.16657883) q[0];
rz(-0.84959787) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(-0.073624484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4723785) q[0];
sx q[0];
rz(-1.3166691) q[0];
sx q[0];
rz(0.73331244) q[0];
rz(0.87009116) q[2];
sx q[2];
rz(-2.4656762) q[2];
sx q[2];
rz(-1.5004339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5391748) q[1];
sx q[1];
rz(-0.6197558) q[1];
sx q[1];
rz(-2.5386993) q[1];
rz(-pi) q[2];
rz(2.5154153) q[3];
sx q[3];
rz(-1.1707843) q[3];
sx q[3];
rz(1.7020066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.9061884) q[2];
sx q[2];
rz(-0.81934682) q[2];
sx q[2];
rz(2.4833615) q[2];
rz(1.341358) q[3];
sx q[3];
rz(-0.96786371) q[3];
sx q[3];
rz(-0.85097504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8258719) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(-0.52147621) q[1];
sx q[1];
rz(-1.5670525) q[1];
sx q[1];
rz(-1.570931) q[1];
rz(1.0339708) q[2];
sx q[2];
rz(-2.2848717) q[2];
sx q[2];
rz(0.90938766) q[2];
rz(-0.6220606) q[3];
sx q[3];
rz(-1.4348772) q[3];
sx q[3];
rz(0.025195599) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
