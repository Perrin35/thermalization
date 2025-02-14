OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(-0.23393272) q[0];
rz(-2.2995931) q[1];
sx q[1];
rz(-1.41398) q[1];
sx q[1];
rz(-1.6230621) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18478014) q[0];
sx q[0];
rz(-1.0597469) q[0];
sx q[0];
rz(1.2132116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0064924) q[2];
sx q[2];
rz(-2.6681136) q[2];
sx q[2];
rz(1.4896637) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3668574) q[1];
sx q[1];
rz(-2.1090057) q[1];
sx q[1];
rz(0.66267207) q[1];
rz(-pi) q[2];
rz(-0.47187658) q[3];
sx q[3];
rz(-2.1778244) q[3];
sx q[3];
rz(-0.31479731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6003549) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(-2.3360628) q[2];
rz(-0.39189288) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(-2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58981744) q[0];
sx q[0];
rz(-2.4601695) q[0];
sx q[0];
rz(2.7031194) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(-0.10800392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7399044) q[0];
sx q[0];
rz(-1.5459014) q[0];
sx q[0];
rz(-1.6817001) q[0];
rz(-pi) q[1];
rz(0.82307016) q[2];
sx q[2];
rz(-2.0688754) q[2];
sx q[2];
rz(0.66253875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9671772) q[1];
sx q[1];
rz(-1.6545873) q[1];
sx q[1];
rz(-1.3166974) q[1];
rz(-pi) q[2];
rz(-1.9181973) q[3];
sx q[3];
rz(-2.6250771) q[3];
sx q[3];
rz(-1.7361189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7763623) q[2];
sx q[2];
rz(-0.47351101) q[2];
sx q[2];
rz(-3.0253809) q[2];
rz(0.41415563) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6666343) q[0];
sx q[0];
rz(-1.8284429) q[0];
sx q[0];
rz(0.83589244) q[0];
rz(1.5180786) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(1.9035043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19003632) q[0];
sx q[0];
rz(-1.0685896) q[0];
sx q[0];
rz(-0.42597187) q[0];
x q[1];
rz(-0.61254259) q[2];
sx q[2];
rz(-2.7296737) q[2];
sx q[2];
rz(-1.4985794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1819396) q[1];
sx q[1];
rz(-1.6283855) q[1];
sx q[1];
rz(-1.6883541) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8594871) q[3];
sx q[3];
rz(-1.8477401) q[3];
sx q[3];
rz(-1.1462584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.28430024) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(-2.4227552) q[2];
rz(2.5607064) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.523664) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(1.1267598) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(-2.0984971) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1820871) q[0];
sx q[0];
rz(-2.8694177) q[0];
sx q[0];
rz(1.4578117) q[0];
x q[1];
rz(-1.0342233) q[2];
sx q[2];
rz(-1.787428) q[2];
sx q[2];
rz(2.72675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14341893) q[1];
sx q[1];
rz(-0.34356782) q[1];
sx q[1];
rz(-1.0982537) q[1];
rz(-pi) q[2];
rz(0.11123379) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(1.8452132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8635233) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(2.1045904) q[2];
rz(0.95748025) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(-2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1400414) q[0];
sx q[0];
rz(-2.934444) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(1.8278488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2976332) q[0];
sx q[0];
rz(-2.5974413) q[0];
sx q[0];
rz(2.0354969) q[0];
rz(0.037363923) q[2];
sx q[2];
rz(-1.4722344) q[2];
sx q[2];
rz(-2.7600206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.524595) q[1];
sx q[1];
rz(-1.6596848) q[1];
sx q[1];
rz(1.7434381) q[1];
x q[2];
rz(-2.143712) q[3];
sx q[3];
rz(-0.30227236) q[3];
sx q[3];
rz(0.59461601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4801415) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(-0.56337774) q[2];
rz(1.2207458) q[3];
sx q[3];
rz(-1.3910339) q[3];
sx q[3];
rz(-2.5105072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(0.011938183) q[0];
rz(-2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(0.24519244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5545776) q[0];
sx q[0];
rz(-2.3759807) q[0];
sx q[0];
rz(-2.2157726) q[0];
rz(-2.4293173) q[2];
sx q[2];
rz(-0.33604188) q[2];
sx q[2];
rz(-3.0672376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4575544) q[1];
sx q[1];
rz(-2.0054617) q[1];
sx q[1];
rz(-2.5358389) q[1];
rz(-pi) q[2];
rz(0.40254503) q[3];
sx q[3];
rz(-1.8410826) q[3];
sx q[3];
rz(-2.9077934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37990722) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(0.53064972) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(-1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.1928007) q[0];
sx q[0];
rz(-1.8185607) q[0];
sx q[0];
rz(-0.2970933) q[0];
rz(1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(0.43103257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6490899) q[0];
sx q[0];
rz(-1.8592872) q[0];
sx q[0];
rz(1.7899075) q[0];
x q[1];
rz(2.9210429) q[2];
sx q[2];
rz(-2.9052581) q[2];
sx q[2];
rz(0.55921171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37087155) q[1];
sx q[1];
rz(-1.2464735) q[1];
sx q[1];
rz(-1.0350219) q[1];
rz(-pi) q[2];
rz(0.2934612) q[3];
sx q[3];
rz(-1.0612504) q[3];
sx q[3];
rz(-2.5603385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3057574) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(-2.8744899) q[2];
rz(-3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(-0.10572461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06374643) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(0.63999501) q[0];
rz(0.85482875) q[1];
sx q[1];
rz(-1.6565485) q[1];
sx q[1];
rz(-1.7291501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28977355) q[0];
sx q[0];
rz(-2.4443194) q[0];
sx q[0];
rz(0.66582219) q[0];
rz(-pi) q[1];
rz(-2.0630025) q[2];
sx q[2];
rz(-1.6160496) q[2];
sx q[2];
rz(2.947888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0113653) q[1];
sx q[1];
rz(-1.9107358) q[1];
sx q[1];
rz(-1.4814427) q[1];
rz(-1.1453923) q[3];
sx q[3];
rz(-1.2266876) q[3];
sx q[3];
rz(-2.1895125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5414446) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(0.42204648) q[2];
rz(2.6680434) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(-2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710829) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(1.3516082) q[0];
rz(-2.8260258) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(1.9884761) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0457927) q[0];
sx q[0];
rz(-1.519916) q[0];
sx q[0];
rz(1.5183543) q[0];
rz(-pi) q[1];
rz(0.081772371) q[2];
sx q[2];
rz(-0.53682971) q[2];
sx q[2];
rz(2.5935136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35495423) q[1];
sx q[1];
rz(-1.6980722) q[1];
sx q[1];
rz(2.5580225) q[1];
rz(-pi) q[2];
rz(1.9824636) q[3];
sx q[3];
rz(-1.2198866) q[3];
sx q[3];
rz(1.0699492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45712581) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(0.37377629) q[2];
rz(1.9155546) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.8208338) q[0];
rz(-2.2044115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(-2.6722867) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26312882) q[0];
sx q[0];
rz(-1.6471905) q[0];
sx q[0];
rz(1.6883649) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1524544) q[2];
sx q[2];
rz(-1.4114221) q[2];
sx q[2];
rz(-1.5053269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1782068) q[1];
sx q[1];
rz(-1.4904516) q[1];
sx q[1];
rz(2.9909913) q[1];
rz(1.3239884) q[3];
sx q[3];
rz(-1.2780407) q[3];
sx q[3];
rz(0.61652641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4361973) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(1.9311284) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(-2.1657522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198915) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(-2.2484491) q[1];
sx q[1];
rz(-1.8291263) q[1];
sx q[1];
rz(-1.6222454) q[1];
rz(1.860113) q[2];
sx q[2];
rz(-1.3500967) q[2];
sx q[2];
rz(-1.4934749) q[2];
rz(0.94598573) q[3];
sx q[3];
rz(-2.47249) q[3];
sx q[3];
rz(-2.5853018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
