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
rz(0.12953144) q[0];
sx q[0];
rz(-0.13106267) q[0];
sx q[0];
rz(2.5370606) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(-0.25783917) q[1];
sx q[1];
rz(-1.2300307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2211232) q[0];
sx q[0];
rz(-1.9029494) q[0];
sx q[0];
rz(-1.3060547) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31754874) q[2];
sx q[2];
rz(-0.63830909) q[2];
sx q[2];
rz(0.93250634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22495077) q[1];
sx q[1];
rz(-1.0826546) q[1];
sx q[1];
rz(-0.78190885) q[1];
rz(-pi) q[2];
rz(1.8038407) q[3];
sx q[3];
rz(-2.6342158) q[3];
sx q[3];
rz(-0.22358596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0543542) q[2];
sx q[2];
rz(-1.9638991) q[2];
sx q[2];
rz(-3.0053906) q[2];
rz(2.6916091) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(-0.78342485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1012652) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(-0.26327565) q[0];
rz(1.0385665) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(2.3668049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.446774) q[0];
sx q[0];
rz(-1.7763486) q[0];
sx q[0];
rz(-0.62384042) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2293533) q[2];
sx q[2];
rz(-0.99102913) q[2];
sx q[2];
rz(1.5383681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0475743) q[1];
sx q[1];
rz(-1.1915922) q[1];
sx q[1];
rz(0.039888558) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0618013) q[3];
sx q[3];
rz(-0.92872075) q[3];
sx q[3];
rz(-1.1428558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7210377) q[2];
sx q[2];
rz(-1.6246395) q[2];
sx q[2];
rz(-0.59845412) q[2];
rz(1.7789486) q[3];
sx q[3];
rz(-2.2612031) q[3];
sx q[3];
rz(-2.8708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28629985) q[0];
sx q[0];
rz(-2.6787651) q[0];
sx q[0];
rz(1.5077952) q[0];
rz(-2.9871121) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(-0.78937626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1868527) q[0];
sx q[0];
rz(-0.53760872) q[0];
sx q[0];
rz(-2.9718999) q[0];
x q[1];
rz(-1.71307) q[2];
sx q[2];
rz(-1.2923489) q[2];
sx q[2];
rz(1.5253893) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7273318) q[1];
sx q[1];
rz(-1.8751029) q[1];
sx q[1];
rz(-2.2060966) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79929914) q[3];
sx q[3];
rz(-1.2033312) q[3];
sx q[3];
rz(-2.7225843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7827451) q[2];
sx q[2];
rz(-0.91155702) q[2];
sx q[2];
rz(1.6645128) q[2];
rz(1.9278256) q[3];
sx q[3];
rz(-1.341308) q[3];
sx q[3];
rz(-1.7267797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724991) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(-0.61071998) q[0];
rz(-2.4702813) q[1];
sx q[1];
rz(-1.6339615) q[1];
sx q[1];
rz(1.3549365) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1917065) q[0];
sx q[0];
rz(-2.0512132) q[0];
sx q[0];
rz(1.0597764) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1995537) q[2];
sx q[2];
rz(-1.6880195) q[2];
sx q[2];
rz(3.0407112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6411533) q[1];
sx q[1];
rz(-1.3257003) q[1];
sx q[1];
rz(-0.064747253) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1576304) q[3];
sx q[3];
rz(-1.4915183) q[3];
sx q[3];
rz(-2.6266499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2063724) q[2];
sx q[2];
rz(-0.78038961) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(0.74784652) q[3];
sx q[3];
rz(-2.3220389) q[3];
sx q[3];
rz(-1.4486754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.185323) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(2.4638033) q[0];
rz(2.9970844) q[1];
sx q[1];
rz(-0.78737193) q[1];
sx q[1];
rz(-1.6544624) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3412171) q[0];
sx q[0];
rz(-1.7066597) q[0];
sx q[0];
rz(3.1151943) q[0];
rz(-pi) q[1];
rz(-3.0838548) q[2];
sx q[2];
rz(-1.1694093) q[2];
sx q[2];
rz(0.58338469) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8890052) q[1];
sx q[1];
rz(-2.41373) q[1];
sx q[1];
rz(0.83830203) q[1];
x q[2];
rz(-3.0934382) q[3];
sx q[3];
rz(-1.1954525) q[3];
sx q[3];
rz(-2.1233692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80374485) q[2];
sx q[2];
rz(-0.2636815) q[2];
sx q[2];
rz(-0.35550508) q[2];
rz(-1.0454987) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(2.579328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1657555) q[0];
sx q[0];
rz(-0.58335692) q[0];
sx q[0];
rz(-3.052886) q[0];
rz(1.5345796) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(-0.075693695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3189118) q[0];
sx q[0];
rz(-1.3782129) q[0];
sx q[0];
rz(0.67355021) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6275109) q[2];
sx q[2];
rz(-2.4091588) q[2];
sx q[2];
rz(-0.54679856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5867501) q[1];
sx q[1];
rz(-1.4289218) q[1];
sx q[1];
rz(-1.0125404) q[1];
x q[2];
rz(2.663732) q[3];
sx q[3];
rz(-2.5091672) q[3];
sx q[3];
rz(-0.53148182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3169516) q[2];
sx q[2];
rz(-1.757916) q[2];
sx q[2];
rz(-1.0392044) q[2];
rz(-0.21229395) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742663) q[0];
sx q[0];
rz(-0.75263158) q[0];
sx q[0];
rz(-1.0443895) q[0];
rz(2.3692865) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(-2.4078802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45424592) q[0];
sx q[0];
rz(-2.2139619) q[0];
sx q[0];
rz(-2.1976794) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1999424) q[2];
sx q[2];
rz(-2.7463253) q[2];
sx q[2];
rz(-2.7140009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1024688) q[1];
sx q[1];
rz(-2.8274635) q[1];
sx q[1];
rz(0.38508319) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83603199) q[3];
sx q[3];
rz(-1.6510909) q[3];
sx q[3];
rz(1.3463595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5250728) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(-2.5329242) q[2];
rz(2.0742553) q[3];
sx q[3];
rz(-0.87456861) q[3];
sx q[3];
rz(2.5909891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50255018) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(-1.4962366) q[0];
rz(1.3642338) q[1];
sx q[1];
rz(-2.5271006) q[1];
sx q[1];
rz(-0.38633698) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2992512) q[0];
sx q[0];
rz(-2.4066397) q[0];
sx q[0];
rz(0.39836653) q[0];
rz(1.9278717) q[2];
sx q[2];
rz(-0.197796) q[2];
sx q[2];
rz(-1.7864428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.92703351) q[1];
sx q[1];
rz(-2.7874569) q[1];
sx q[1];
rz(-2.781809) q[1];
x q[2];
rz(-0.058606996) q[3];
sx q[3];
rz(-0.59133619) q[3];
sx q[3];
rz(3.1268596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51314696) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(0.15303843) q[2];
rz(2.2506574) q[3];
sx q[3];
rz(-0.95518437) q[3];
sx q[3];
rz(-2.4002767) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9855758) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(2.7942221) q[0];
rz(0.037847606) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(-0.52245021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.396401) q[0];
sx q[0];
rz(-0.42310444) q[0];
sx q[0];
rz(-1.7897096) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36549296) q[2];
sx q[2];
rz(-1.5704944) q[2];
sx q[2];
rz(-0.049916849) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3047239) q[1];
sx q[1];
rz(-0.68075276) q[1];
sx q[1];
rz(0.70226837) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96790989) q[3];
sx q[3];
rz(-2.4914352) q[3];
sx q[3];
rz(0.12891836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5106421) q[2];
sx q[2];
rz(-2.6783671) q[2];
sx q[2];
rz(-1.9412712) q[2];
rz(0.55780324) q[3];
sx q[3];
rz(-1.6226945) q[3];
sx q[3];
rz(-0.38448486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8484304) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(1.4137319) q[0];
rz(-3.0005455) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(0.65473762) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075874627) q[0];
sx q[0];
rz(-1.7786586) q[0];
sx q[0];
rz(-1.3285358) q[0];
rz(-pi) q[1];
rz(2.5467196) q[2];
sx q[2];
rz(-1.1385749) q[2];
sx q[2];
rz(3.1405666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.84147553) q[1];
sx q[1];
rz(-1.8130596) q[1];
sx q[1];
rz(-2.9650142) q[1];
rz(-pi) q[2];
rz(-0.68924381) q[3];
sx q[3];
rz(-0.74604496) q[3];
sx q[3];
rz(0.41710284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(0.72511017) q[2];
rz(-2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(0.58722535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643322) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(2.4779392) q[2];
sx q[2];
rz(-1.2812231) q[2];
sx q[2];
rz(-0.44322586) q[2];
rz(-1.7305456) q[3];
sx q[3];
rz(-0.86426576) q[3];
sx q[3];
rz(-0.15440253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
