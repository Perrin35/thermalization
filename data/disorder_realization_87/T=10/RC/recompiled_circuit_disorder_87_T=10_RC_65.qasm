OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(-2.9400163) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21298458) q[0];
sx q[0];
rz(-2.1323418) q[0];
sx q[0];
rz(-1.0667849) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4891698) q[2];
sx q[2];
rz(-0.7235652) q[2];
sx q[2];
rz(0.10318081) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1569251) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(-1.9763293) q[1];
rz(1.2308146) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(1.1839379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(-3.0554331) q[2];
rz(2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57698292) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(1.3264867) q[0];
rz(-1.8857229) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(0.27145162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4685681) q[0];
sx q[0];
rz(-0.66267555) q[0];
sx q[0];
rz(-0.89232348) q[0];
x q[1];
rz(0.51241264) q[2];
sx q[2];
rz(-2.5645442) q[2];
sx q[2];
rz(-1.3702099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2084864) q[1];
sx q[1];
rz(-0.83175627) q[1];
sx q[1];
rz(2.0887124) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7852887) q[3];
sx q[3];
rz(-2.2998527) q[3];
sx q[3];
rz(-0.19405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0945956) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-0.14257167) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3554879) q[0];
sx q[0];
rz(-0.84083637) q[0];
sx q[0];
rz(-0.79746042) q[0];
rz(0.44720165) q[2];
sx q[2];
rz(-0.7913835) q[2];
sx q[2];
rz(-0.42391047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8612954) q[1];
sx q[1];
rz(-2.0164844) q[1];
sx q[1];
rz(2.6764826) q[1];
x q[2];
rz(-2.530982) q[3];
sx q[3];
rz(-0.36452499) q[3];
sx q[3];
rz(-0.025346905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(2.7374632) q[2];
rz(1.8557619) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(-0.96281111) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2179746) q[0];
sx q[0];
rz(-1.0110564) q[0];
sx q[0];
rz(2.2773507) q[0];
x q[1];
rz(-1.5702815) q[2];
sx q[2];
rz(-1.4467903) q[2];
sx q[2];
rz(-2.4027783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8648659) q[1];
sx q[1];
rz(-0.58502561) q[1];
sx q[1];
rz(-1.8086955) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13929316) q[3];
sx q[3];
rz(-2.5677498) q[3];
sx q[3];
rz(1.7593256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(3.0299419) q[2];
rz(0.81104898) q[3];
sx q[3];
rz(-0.45447293) q[3];
sx q[3];
rz(0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(2.7850889) q[0];
rz(-0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(2.8809663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31291744) q[0];
sx q[0];
rz(-0.14794359) q[0];
sx q[0];
rz(-2.8539128) q[0];
rz(-pi) q[1];
rz(-0.79029681) q[2];
sx q[2];
rz(-2.7919263) q[2];
sx q[2];
rz(2.5555573) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.534429) q[1];
sx q[1];
rz(-0.57463127) q[1];
sx q[1];
rz(-3.1092005) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9959833) q[3];
sx q[3];
rz(-1.7123316) q[3];
sx q[3];
rz(-1.4841929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2300718) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(-1.4916346) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(1.7274436) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3154253) q[0];
sx q[0];
rz(-2.8624479) q[0];
sx q[0];
rz(2.9690353) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6238407) q[2];
sx q[2];
rz(-1.8038097) q[2];
sx q[2];
rz(-2.0875967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1097088) q[1];
sx q[1];
rz(-2.0307396) q[1];
sx q[1];
rz(-2.1214532) q[1];
rz(-pi) q[2];
x q[2];
rz(0.075501637) q[3];
sx q[3];
rz(-1.3550948) q[3];
sx q[3];
rz(0.84738934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.002939) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(-3.1138528) q[2];
rz(2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.3940575) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7572927) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(3.0199155) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(0.40245232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84797317) q[0];
sx q[0];
rz(-1.8403887) q[0];
sx q[0];
rz(0.86975354) q[0];
x q[1];
rz(3.0955663) q[2];
sx q[2];
rz(-1.5791025) q[2];
sx q[2];
rz(1.9866895) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35343364) q[1];
sx q[1];
rz(-0.33553365) q[1];
sx q[1];
rz(0.56281705) q[1];
rz(-pi) q[2];
rz(1.6472858) q[3];
sx q[3];
rz(-0.34919958) q[3];
sx q[3];
rz(-2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.014331269) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(0.47752738) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(0.91807085) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0089793423) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(-pi) q[1];
rz(3.0389298) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(1.7662802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1840399) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(2.1368105) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7233852) q[3];
sx q[3];
rz(-1.6697262) q[3];
sx q[3];
rz(-1.4560771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(2.365716) q[2];
rz(-0.72426978) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(2.3628078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22395615) q[0];
sx q[0];
rz(-2.7178239) q[0];
sx q[0];
rz(2.5013322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0987894) q[2];
sx q[2];
rz(-1.325377) q[2];
sx q[2];
rz(0.62308842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3719337) q[1];
sx q[1];
rz(-1.2049335) q[1];
sx q[1];
rz(2.8918173) q[1];
x q[2];
rz(3.0770244) q[3];
sx q[3];
rz(-1.7913006) q[3];
sx q[3];
rz(1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.3042973) q[2];
sx q[2];
rz(-0.71643913) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(-1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0338106) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(2.418628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3176214) q[0];
sx q[0];
rz(-2.9569607) q[0];
sx q[0];
rz(0.95438192) q[0];
x q[1];
rz(-1.9777708) q[2];
sx q[2];
rz(-0.40728912) q[2];
sx q[2];
rz(-2.2182857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55262676) q[1];
sx q[1];
rz(-2.6347661) q[1];
sx q[1];
rz(-1.456702) q[1];
x q[2];
rz(-0.78852699) q[3];
sx q[3];
rz(-1.8619814) q[3];
sx q[3];
rz(-0.2083208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(0.52552137) q[2];
rz(-2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239607) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(2.0422968) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(0.99156689) q[2];
sx q[2];
rz(-2.3835923) q[2];
sx q[2];
rz(0.53072416) q[2];
rz(1.3463734) q[3];
sx q[3];
rz(-1.8891469) q[3];
sx q[3];
rz(-0.3286152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];