OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(4.1424799) q[0];
sx q[0];
rz(9.6371798) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(-1.2815055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6706657) q[0];
sx q[0];
rz(-2.7424194) q[0];
sx q[0];
rz(-2.1098718) q[0];
x q[1];
rz(1.6002866) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(0.32683795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68583381) q[1];
sx q[1];
rz(-1.7384006) q[1];
sx q[1];
rz(-2.0007677) q[1];
rz(-pi) q[2];
rz(-0.063135677) q[3];
sx q[3];
rz(-1.4423443) q[3];
sx q[3];
rz(-2.0383143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(-2.4860399) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(-0.15629388) q[0];
rz(0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(-1.4888391) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50870897) q[0];
sx q[0];
rz(-1.7299621) q[0];
sx q[0];
rz(-0.55523086) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.312617) q[2];
sx q[2];
rz(-1.0114705) q[2];
sx q[2];
rz(1.0499357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3540346) q[1];
sx q[1];
rz(-2.2984142) q[1];
sx q[1];
rz(-2.3393199) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.00499) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(-0.68340106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(-2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(1.9077574) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082224) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(0.5439533) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.9525607) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9246763) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-0.49305537) q[0];
rz(-pi) q[1];
rz(0.24818111) q[2];
sx q[2];
rz(-1.6096874) q[2];
sx q[2];
rz(1.3523906) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5570453) q[1];
sx q[1];
rz(-1.0439596) q[1];
sx q[1];
rz(0.79750632) q[1];
rz(-2.4431908) q[3];
sx q[3];
rz(-1.2050556) q[3];
sx q[3];
rz(1.48667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(0.35513487) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(0.60607934) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(1.9569424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9979447) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(2.4585637) q[0];
rz(-0.46063936) q[2];
sx q[2];
rz(-1.1038291) q[2];
sx q[2];
rz(-2.4796782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8682755) q[1];
sx q[1];
rz(-1.6874716) q[1];
sx q[1];
rz(0.97881808) q[1];
x q[2];
rz(-0.70373669) q[3];
sx q[3];
rz(-0.53547137) q[3];
sx q[3];
rz(-2.717201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(-0.8574287) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(0.52039352) q[0];
rz(-0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(-0.72881126) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4686574) q[0];
sx q[0];
rz(-2.9142083) q[0];
sx q[0];
rz(2.6926281) q[0];
x q[1];
rz(-0.8466709) q[2];
sx q[2];
rz(-1.6933105) q[2];
sx q[2];
rz(0.12144897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5853873) q[1];
sx q[1];
rz(-0.3857179) q[1];
sx q[1];
rz(-1.0005887) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6008984) q[3];
sx q[3];
rz(-1.027642) q[3];
sx q[3];
rz(-2.3054996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(-0.63306159) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(-0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(0.9206413) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716498) q[0];
sx q[0];
rz(-2.2075704) q[0];
sx q[0];
rz(0.79983478) q[0];
x q[1];
rz(-1.6502041) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(0.90027819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27053988) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(-0.096710042) q[1];
x q[2];
rz(0.88455172) q[3];
sx q[3];
rz(-0.93730799) q[3];
sx q[3];
rz(-1.8201049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-0.075604288) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(-0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(-1.9627337) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(1.1694318) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021381) q[0];
sx q[0];
rz(-1.6374145) q[0];
sx q[0];
rz(-1.735795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.85497) q[2];
sx q[2];
rz(-1.9467762) q[2];
sx q[2];
rz(0.77501955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84050256) q[1];
sx q[1];
rz(-1.2355348) q[1];
sx q[1];
rz(-0.5718949) q[1];
rz(-1.5285989) q[3];
sx q[3];
rz(-1.183126) q[3];
sx q[3];
rz(-1.8488415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.823267) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(1.7100122) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1201046) q[0];
sx q[0];
rz(-1.095626) q[0];
sx q[0];
rz(3.0475463) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2051815) q[2];
sx q[2];
rz(-0.80447703) q[2];
sx q[2];
rz(0.22125439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8843958) q[1];
sx q[1];
rz(-2.307502) q[1];
sx q[1];
rz(-0.019330545) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1222955) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(-2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-0.28253728) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(-2.4422586) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(-1.2618077) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5489053) q[0];
sx q[0];
rz(-1.840001) q[0];
sx q[0];
rz(0.44434987) q[0];
x q[1];
rz(-2.8724573) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(0.035949635) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.251235) q[1];
sx q[1];
rz(-2.2489378) q[1];
sx q[1];
rz(1.2441586) q[1];
x q[2];
rz(0.52957876) q[3];
sx q[3];
rz(-1.7486608) q[3];
sx q[3];
rz(1.9999268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(0.66649377) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(-0.054140422) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150072) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(1.8146145) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80600593) q[2];
sx q[2];
rz(-1.5984048) q[2];
sx q[2];
rz(-2.1105786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66026238) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(-2.7923613) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.530738) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(0.50869298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3728309) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(-2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6282745) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(0.8387768) q[2];
sx q[2];
rz(-1.5427187) q[2];
sx q[2];
rz(-1.6533921) q[2];
rz(1.0735687) q[3];
sx q[3];
rz(-2.0484925) q[3];
sx q[3];
rz(-1.9852553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
