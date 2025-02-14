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
rz(-1.7614814) q[0];
sx q[0];
rz(7.8742134) q[0];
sx q[0];
rz(5.9013517) q[0];
rz(-0.39185697) q[1];
sx q[1];
rz(4.691603) q[1];
sx q[1];
rz(7.2075972) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5474553) q[0];
sx q[0];
rz(-2.1004803) q[0];
sx q[0];
rz(0.87181566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7713791) q[2];
sx q[2];
rz(-1.5892803) q[2];
sx q[2];
rz(-0.12272515) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83340911) q[1];
sx q[1];
rz(-1.7004968) q[1];
sx q[1];
rz(-0.42713487) q[1];
x q[2];
rz(2.1788039) q[3];
sx q[3];
rz(-2.5217524) q[3];
sx q[3];
rz(2.7328797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(1.5748242) q[2];
rz(-1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-0.80802149) q[0];
sx q[0];
rz(-1.1760733) q[0];
rz(0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(-0.31879058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87321157) q[0];
sx q[0];
rz(-0.86641524) q[0];
sx q[0];
rz(1.9834788) q[0];
rz(2.712557) q[2];
sx q[2];
rz(-1.7861767) q[2];
sx q[2];
rz(-1.693415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76004878) q[1];
sx q[1];
rz(-0.69763073) q[1];
sx q[1];
rz(-0.95894869) q[1];
rz(-pi) q[2];
rz(1.1652496) q[3];
sx q[3];
rz(-0.31003788) q[3];
sx q[3];
rz(1.3998717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0626283) q[2];
sx q[2];
rz(-2.2603409) q[2];
sx q[2];
rz(-0.47743615) q[2];
rz(1.7834974) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(-0.0099946578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7497712) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(2.0023477) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(-2.4635945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8157522) q[0];
sx q[0];
rz(-1.6597972) q[0];
sx q[0];
rz(-0.98317105) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35786104) q[2];
sx q[2];
rz(-2.1659019) q[2];
sx q[2];
rz(-2.8107363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19723141) q[1];
sx q[1];
rz(-1.1121096) q[1];
sx q[1];
rz(0.58802559) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0150681) q[3];
sx q[3];
rz(-0.58150269) q[3];
sx q[3];
rz(2.8648928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71538007) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-1.0769843) q[2];
rz(-1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(0.25585678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.857321) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(2.381109) q[0];
rz(2.7898232) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(2.9913091) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9780804) q[0];
sx q[0];
rz(-0.0038359782) q[0];
sx q[0];
rz(-1.0413961) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1614805) q[2];
sx q[2];
rz(-2.7285353) q[2];
sx q[2];
rz(-1.6020136) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9517512) q[1];
sx q[1];
rz(-2.449936) q[1];
sx q[1];
rz(1.0361978) q[1];
x q[2];
rz(0.45771367) q[3];
sx q[3];
rz(-1.8031075) q[3];
sx q[3];
rz(-1.7328229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1111697) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(-2.5566768) q[2];
rz(0.065718204) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(2.0871302) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86972648) q[0];
sx q[0];
rz(-1.9166742) q[0];
sx q[0];
rz(-0.72750339) q[0];
rz(-1.2959405) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(-0.71472275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22435221) q[0];
sx q[0];
rz(-0.2324129) q[0];
sx q[0];
rz(-0.31347855) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6863912) q[2];
sx q[2];
rz(-1.2718448) q[2];
sx q[2];
rz(1.40708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59698518) q[1];
sx q[1];
rz(-1.7555573) q[1];
sx q[1];
rz(1.1853335) q[1];
x q[2];
rz(1.0017702) q[3];
sx q[3];
rz(-2.4845893) q[3];
sx q[3];
rz(-3.0083619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(-1.145251) q[2];
rz(-0.20259419) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957434) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(2.9737245) q[0];
rz(0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(-1.3289183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3053554) q[0];
sx q[0];
rz(-1.0389881) q[0];
sx q[0];
rz(-2.1270069) q[0];
rz(-pi) q[1];
rz(-3.0014056) q[2];
sx q[2];
rz(-1.6194238) q[2];
sx q[2];
rz(1.0076154) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5649453) q[1];
sx q[1];
rz(-0.35035808) q[1];
sx q[1];
rz(-2.3957361) q[1];
rz(-0.82002016) q[3];
sx q[3];
rz(-0.32147543) q[3];
sx q[3];
rz(-3.0177651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7557175) q[2];
sx q[2];
rz(-1.8942602) q[2];
sx q[2];
rz(0.34783777) q[2];
rz(-0.31473413) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(-2.0033964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9921853) q[0];
sx q[0];
rz(-1.6353761) q[0];
sx q[0];
rz(2.1904679) q[0];
rz(-0.29257193) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(2.1975885) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2068588) q[0];
sx q[0];
rz(-2.8331596) q[0];
sx q[0];
rz(1.6215723) q[0];
rz(-pi) q[1];
rz(-2.0538141) q[2];
sx q[2];
rz(-1.1406787) q[2];
sx q[2];
rz(1.6295691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4170226) q[1];
sx q[1];
rz(-0.17462433) q[1];
sx q[1];
rz(-1.440446) q[1];
rz(-pi) q[2];
rz(1.0349501) q[3];
sx q[3];
rz(-2.6831045) q[3];
sx q[3];
rz(-1.637984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49147478) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(-1.4415119) q[2];
rz(3.1169543) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(0.98954454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330419) q[0];
sx q[0];
rz(-1.7311743) q[0];
sx q[0];
rz(-0.1380052) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(-2.2053351) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9662921) q[0];
sx q[0];
rz(-1.9387687) q[0];
sx q[0];
rz(-1.247768) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1571884) q[2];
sx q[2];
rz(-2.3086562) q[2];
sx q[2];
rz(-0.10228233) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4605751) q[1];
sx q[1];
rz(-1.5854757) q[1];
sx q[1];
rz(-1.3880066) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8209723) q[3];
sx q[3];
rz(-0.71133876) q[3];
sx q[3];
rz(-2.5742755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1428947) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(1.3188837) q[2];
rz(-0.12254347) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97813022) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(2.7769856) q[0];
rz(-1.3847903) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(0.91420954) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6624266) q[0];
sx q[0];
rz(-1.7044596) q[0];
sx q[0];
rz(-2.0019309) q[0];
x q[1];
rz(1.6169102) q[2];
sx q[2];
rz(-2.3048525) q[2];
sx q[2];
rz(3.0457599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85687145) q[1];
sx q[1];
rz(-2.4865013) q[1];
sx q[1];
rz(0.22677179) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8960885) q[3];
sx q[3];
rz(-2.5872018) q[3];
sx q[3];
rz(1.7911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30283516) q[2];
sx q[2];
rz(-2.5573533) q[2];
sx q[2];
rz(-1.6287104) q[2];
rz(-2.5527939) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(2.485386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47422472) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(-2.8571416) q[0];
rz(-2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(0.42025748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81440825) q[0];
sx q[0];
rz(-1.1564009) q[0];
sx q[0];
rz(-2.9073858) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0345694) q[2];
sx q[2];
rz(-0.83597413) q[2];
sx q[2];
rz(-0.42586621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3455167) q[1];
sx q[1];
rz(-0.14847595) q[1];
sx q[1];
rz(0.92904894) q[1];
x q[2];
rz(2.612667) q[3];
sx q[3];
rz(-0.63415895) q[3];
sx q[3];
rz(-2.8932539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1493211) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(2.1743656) q[2];
rz(2.1238756) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348542) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(-3.0178487) q[1];
sx q[1];
rz(-2.1950304) q[1];
sx q[1];
rz(-1.310941) q[1];
rz(2.7563949) q[2];
sx q[2];
rz(-2.7976947) q[2];
sx q[2];
rz(-1.7277389) q[2];
rz(-0.31207257) q[3];
sx q[3];
rz(-1.5371116) q[3];
sx q[3];
rz(-0.1826382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
