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
rz(-2.969279) q[0];
sx q[0];
rz(-0.20792374) q[0];
sx q[0];
rz(1.8508258) q[0];
rz(0.15200226) q[1];
sx q[1];
rz(-2.4988889) q[1];
sx q[1];
rz(0.42833498) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5458459) q[0];
sx q[0];
rz(-1.6753249) q[0];
sx q[0];
rz(-1.680611) q[0];
x q[1];
rz(-1.8646127) q[2];
sx q[2];
rz(-0.91962469) q[2];
sx q[2];
rz(3.0319954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1312771) q[1];
sx q[1];
rz(-2.4020264) q[1];
sx q[1];
rz(3.0775978) q[1];
x q[2];
rz(0.21875225) q[3];
sx q[3];
rz(-1.7667337) q[3];
sx q[3];
rz(2.4234555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20994818) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(-1.9317365) q[2];
rz(3.0620388) q[3];
sx q[3];
rz(-2.7456561) q[3];
sx q[3];
rz(-2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2570268) q[0];
sx q[0];
rz(-2.631812) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(-0.89029038) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(-2.8382137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6174406) q[0];
sx q[0];
rz(-0.57311237) q[0];
sx q[0];
rz(1.2138496) q[0];
x q[1];
rz(-1.0531002) q[2];
sx q[2];
rz(-1.3441022) q[2];
sx q[2];
rz(-1.1013307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1282004) q[1];
sx q[1];
rz(-1.6329995) q[1];
sx q[1];
rz(-1.4799467) q[1];
x q[2];
rz(0.56191617) q[3];
sx q[3];
rz(-0.73560916) q[3];
sx q[3];
rz(-0.058678415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38637912) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(0.73491043) q[2];
rz(-0.84878659) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(-0.36062226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3291149) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(-0.078711674) q[0];
rz(-2.3178237) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(1.8471921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5536907) q[0];
sx q[0];
rz(-1.5103834) q[0];
sx q[0];
rz(-1.509045) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50991518) q[2];
sx q[2];
rz(-1.0140061) q[2];
sx q[2];
rz(0.82600466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0261198) q[1];
sx q[1];
rz(-1.2900262) q[1];
sx q[1];
rz(-1.6927648) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2749316) q[3];
sx q[3];
rz(-1.7539644) q[3];
sx q[3];
rz(-0.66732348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89982975) q[2];
sx q[2];
rz(-2.7624942) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(1.5166616) q[3];
sx q[3];
rz(-0.95976019) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4888332) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(-0.56831992) q[0];
rz(1.142451) q[1];
sx q[1];
rz(-1.7894952) q[1];
sx q[1];
rz(1.4521339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947493) q[0];
sx q[0];
rz(-0.94943014) q[0];
sx q[0];
rz(1.2970379) q[0];
x q[1];
rz(0.4770437) q[2];
sx q[2];
rz(-2.9397268) q[2];
sx q[2];
rz(2.9143726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4751079) q[1];
sx q[1];
rz(-2.1748073) q[1];
sx q[1];
rz(-0.38352769) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0978974) q[3];
sx q[3];
rz(-2.7644025) q[3];
sx q[3];
rz(-2.8998714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51231724) q[2];
sx q[2];
rz(-0.53048152) q[2];
sx q[2];
rz(0.076889195) q[2];
rz(0.39938375) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870616) q[0];
sx q[0];
rz(-2.7847544) q[0];
sx q[0];
rz(-0.29990184) q[0];
rz(-1.1497644) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(-2.6531175) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91575275) q[0];
sx q[0];
rz(-1.3855278) q[0];
sx q[0];
rz(0.033323296) q[0];
rz(-pi) q[1];
rz(2.5484094) q[2];
sx q[2];
rz(-2.0016243) q[2];
sx q[2];
rz(-2.5938428) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8676973) q[1];
sx q[1];
rz(-0.58389837) q[1];
sx q[1];
rz(-0.03429596) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8397055) q[3];
sx q[3];
rz(-0.9465512) q[3];
sx q[3];
rz(1.7549124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38146314) q[2];
sx q[2];
rz(-0.13801408) q[2];
sx q[2];
rz(1.4786973) q[2];
rz(1.1100769) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(-2.4983675) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4019302) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(-0.31841835) q[0];
rz(-3.1069801) q[1];
sx q[1];
rz(-0.56762677) q[1];
sx q[1];
rz(2.6908223) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0552669) q[0];
sx q[0];
rz(-1.0604211) q[0];
sx q[0];
rz(-3.1253405) q[0];
x q[1];
rz(2.4813406) q[2];
sx q[2];
rz(-2.1532032) q[2];
sx q[2];
rz(0.16124111) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3035982) q[1];
sx q[1];
rz(-2.384409) q[1];
sx q[1];
rz(0.96214575) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1799906) q[3];
sx q[3];
rz(-0.79422659) q[3];
sx q[3];
rz(-3.1193747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.047711756) q[2];
sx q[2];
rz(-2.9517951) q[2];
sx q[2];
rz(-3.0799358) q[2];
rz(0.098585248) q[3];
sx q[3];
rz(-2.3969789) q[3];
sx q[3];
rz(-0.70257598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.75673574) q[0];
sx q[0];
rz(-1.1133794) q[0];
sx q[0];
rz(3.1016438) q[0];
rz(0.7705676) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(-0.22824731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21439274) q[0];
sx q[0];
rz(-1.6480371) q[0];
sx q[0];
rz(1.6595675) q[0];
x q[1];
rz(1.8769774) q[2];
sx q[2];
rz(-1.5723229) q[2];
sx q[2];
rz(0.36097872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3769987) q[1];
sx q[1];
rz(-1.5176763) q[1];
sx q[1];
rz(0.082100987) q[1];
x q[2];
rz(-0.41577783) q[3];
sx q[3];
rz(-0.97857394) q[3];
sx q[3];
rz(-1.7433188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.063529) q[2];
sx q[2];
rz(-1.0588131) q[2];
sx q[2];
rz(-2.883319) q[2];
rz(2.9410948) q[3];
sx q[3];
rz(-0.79810464) q[3];
sx q[3];
rz(2.958278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16266009) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(0.3072511) q[0];
rz(1.2112674) q[1];
sx q[1];
rz(-0.4937506) q[1];
sx q[1];
rz(-0.58120751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2659822) q[0];
sx q[0];
rz(-1.6020244) q[0];
sx q[0];
rz(-1.9755548) q[0];
rz(-pi) q[1];
rz(0.18980726) q[2];
sx q[2];
rz(-0.71868374) q[2];
sx q[2];
rz(-0.41658066) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8381184) q[1];
sx q[1];
rz(-1.3708769) q[1];
sx q[1];
rz(-0.78472991) q[1];
rz(-pi) q[2];
rz(-1.2138661) q[3];
sx q[3];
rz(-2.078131) q[3];
sx q[3];
rz(1.1905673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4622978) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(-0.031166859) q[2];
rz(-0.22552414) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533326) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(-0.70892507) q[0];
rz(-2.9290579) q[1];
sx q[1];
rz(-1.0147076) q[1];
sx q[1];
rz(-2.2841891) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564726) q[0];
sx q[0];
rz(-1.5758638) q[0];
sx q[0];
rz(1.4410254) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6078244) q[2];
sx q[2];
rz(-1.8532751) q[2];
sx q[2];
rz(-0.9953645) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.556223) q[1];
sx q[1];
rz(-1.959728) q[1];
sx q[1];
rz(-2.5635368) q[1];
rz(-pi) q[2];
rz(-0.69714947) q[3];
sx q[3];
rz(-1.7350041) q[3];
sx q[3];
rz(2.177161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7015486) q[2];
sx q[2];
rz(-0.35228071) q[2];
sx q[2];
rz(2.3360543) q[2];
rz(0.37197477) q[3];
sx q[3];
rz(-1.6476846) q[3];
sx q[3];
rz(0.39723799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.032967903) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(0.97880542) q[0];
rz(0.72273123) q[1];
sx q[1];
rz(-2.0131854) q[1];
sx q[1];
rz(0.61789787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4198087) q[0];
sx q[0];
rz(-0.81994826) q[0];
sx q[0];
rz(2.7616114) q[0];
rz(-2.0892136) q[2];
sx q[2];
rz(-0.978865) q[2];
sx q[2];
rz(-3.0822138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5717575) q[1];
sx q[1];
rz(-1.7392842) q[1];
sx q[1];
rz(1.7525826) q[1];
rz(-0.79094751) q[3];
sx q[3];
rz(-3.0229048) q[3];
sx q[3];
rz(-3.0968248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.939398) q[2];
sx q[2];
rz(-1.7859744) q[2];
sx q[2];
rz(-3.1414269) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(0.59529006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547739) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(-1.3407002) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(0.98967057) q[2];
sx q[2];
rz(-2.4151617) q[2];
sx q[2];
rz(-2.5834609) q[2];
rz(0.67047337) q[3];
sx q[3];
rz(-2.5994876) q[3];
sx q[3];
rz(2.2710298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
