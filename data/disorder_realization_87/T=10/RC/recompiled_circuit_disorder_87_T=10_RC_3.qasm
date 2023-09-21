OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(5.5570931) q[0];
sx q[0];
rz(9.2232016) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(2.6013241) q[1];
sx q[1];
rz(7.2202914) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0159338) q[0];
sx q[0];
rz(-2.4056245) q[0];
sx q[0];
rz(0.65471162) q[0];
x q[1];
rz(-1.0785157) q[2];
sx q[2];
rz(-2.1247851) q[2];
sx q[2];
rz(2.2433777) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1260682) q[1];
sx q[1];
rz(-1.4926732) q[1];
sx q[1];
rz(1.3874345) q[1];
rz(-pi) q[2];
rz(-0.097221656) q[3];
sx q[3];
rz(-1.8386278) q[3];
sx q[3];
rz(-0.83084805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(-0.086159555) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57698292) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(1.8857229) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(2.870141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2651228) q[0];
sx q[0];
rz(-2.070283) q[0];
sx q[0];
rz(0.45544099) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51606744) q[2];
sx q[2];
rz(-1.8415673) q[2];
sx q[2];
rz(-0.64112907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27122341) q[1];
sx q[1];
rz(-1.9454114) q[1];
sx q[1];
rz(-2.3323374) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1984962) q[3];
sx q[3];
rz(-2.3447403) q[3];
sx q[3];
rz(-2.4380395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0469971) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(-2.2180637) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-2.999021) q[0];
rz(-1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-2.9325063) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3269743) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(0.66280611) q[0];
x q[1];
rz(2.4019037) q[2];
sx q[2];
rz(-1.2581173) q[2];
sx q[2];
rz(2.3198421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0774539) q[1];
sx q[1];
rz(-1.1541379) q[1];
sx q[1];
rz(2.0616848) q[1];
x q[2];
rz(1.3554058) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(-2.4733558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(2.1787815) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-3.1406291) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9199333) q[0];
sx q[0];
rz(-2.1533305) q[0];
sx q[0];
rz(-0.68908738) q[0];
rz(-pi) q[1];
rz(-3.1374627) q[2];
sx q[2];
rz(-3.0175856) q[2];
sx q[2];
rz(-2.4069402) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5597792) q[1];
sx q[1];
rz(-2.1372791) q[1];
sx q[1];
rz(2.986746) q[1];
rz(-pi) q[2];
rz(3.0022995) q[3];
sx q[3];
rz(-0.5738429) q[3];
sx q[3];
rz(-1.3822671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0659539) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(2.6351392) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(-0.26062632) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31291744) q[0];
sx q[0];
rz(-2.9936491) q[0];
sx q[0];
rz(-2.8539128) q[0];
rz(1.3172651) q[2];
sx q[2];
rz(-1.3273444) q[2];
sx q[2];
rz(-1.4075116) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.645748) q[1];
sx q[1];
rz(-0.9965047) q[1];
sx q[1];
rz(1.5498284) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3652606) q[3];
sx q[3];
rz(-0.20271248) q[3];
sx q[3];
rz(2.4622038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9115209) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(-0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(-2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3689573) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(-1.649958) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.7274436) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0522239) q[0];
sx q[0];
rz(-1.5234689) q[0];
sx q[0];
rz(-2.8663859) q[0];
rz(-1.3041777) q[2];
sx q[2];
rz(-1.068371) q[2];
sx q[2];
rz(2.7555639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27281877) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(0.52656071) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7870951) q[3];
sx q[3];
rz(-1.6445451) q[3];
sx q[3];
rz(-2.4019965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1386537) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(0.49267832) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(0.12167715) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(-2.7391403) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2936195) q[0];
sx q[0];
rz(-1.3012039) q[0];
sx q[0];
rz(0.86975354) q[0];
rz(-pi) q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(2.7260821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23558815) q[1];
sx q[1];
rz(-1.2885805) q[1];
sx q[1];
rz(-1.7547592) q[1];
rz(1.4943069) q[3];
sx q[3];
rz(-2.7923931) q[3];
sx q[3];
rz(-2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(-0.86223117) q[2];
rz(-2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(-0.73295897) q[0];
rz(-0.14006242) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(1.0345116) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0089793423) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(2.5402252) q[2];
sx q[2];
rz(-0.12430087) q[2];
sx q[2];
rz(2.3483495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1840399) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(-2.1368105) q[1];
rz(-1.4626059) q[3];
sx q[3];
rz(-1.1547609) q[3];
sx q[3];
rz(3.0707404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(0.72426978) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(-1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6999321) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-0.016816703) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-2.3628078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.461207) q[0];
sx q[0];
rz(-1.2347504) q[0];
sx q[0];
rz(-1.8340322) q[0];
rz(-pi) q[1];
rz(2.07431) q[2];
sx q[2];
rz(-0.52769606) q[2];
sx q[2];
rz(-0.50349456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15021819) q[1];
sx q[1];
rz(-2.7017936) q[1];
sx q[1];
rz(-2.1437953) q[1];
rz(-pi) q[2];
rz(-1.2905144) q[3];
sx q[3];
rz(-0.22961578) q[3];
sx q[3];
rz(-1.5748101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(-1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(0.15144908) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(-2.418628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13847362) q[0];
sx q[0];
rz(-1.4644633) q[0];
sx q[0];
rz(-1.7220201) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9777708) q[2];
sx q[2];
rz(-0.40728912) q[2];
sx q[2];
rz(-2.2182857) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0235698) q[1];
sx q[1];
rz(-1.5155063) q[1];
sx q[1];
rz(-1.0667332) q[1];
rz(-pi) q[2];
rz(-0.39977269) q[3];
sx q[3];
rz(-0.82953605) q[3];
sx q[3];
rz(1.5012036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(-0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(-2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(2.2407871) q[2];
sx q[2];
rz(-1.1849891) q[2];
sx q[2];
rz(1.6580788) q[2];
rz(-0.59393926) q[3];
sx q[3];
rz(-0.38729061) q[3];
sx q[3];
rz(-0.95872986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];