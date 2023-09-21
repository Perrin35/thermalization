OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(2.9291908) q[0];
rz(0.71495932) q[1];
sx q[1];
rz(-2.3541048) q[1];
sx q[1];
rz(1.2815055) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40385383) q[0];
sx q[0];
rz(-1.7716584) q[0];
sx q[0];
rz(-1.9181197) q[0];
x q[1];
rz(-0.69252695) q[2];
sx q[2];
rz(-1.5934957) q[2];
sx q[2];
rz(1.9164617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.233853) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(-1.9563667) q[1];
rz(-pi) q[2];
rz(1.1164066) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(1.56173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(1.6248576) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1600906) q[0];
sx q[0];
rz(-2.1182051) q[0];
sx q[0];
rz(1.7574969) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4376051) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(0.068254452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3540346) q[1];
sx q[1];
rz(-0.84317849) q[1];
sx q[1];
rz(2.3393199) q[1];
rz(-0.1366027) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(-2.4581916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
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
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(1.9525607) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30071242) q[0];
sx q[0];
rz(-1.2749201) q[0];
sx q[0];
rz(-1.4062675) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6109153) q[2];
sx q[2];
rz(-1.3228068) q[2];
sx q[2];
rz(0.2085533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51057286) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(-0.87639798) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69840188) q[3];
sx q[3];
rz(-1.2050556) q[3];
sx q[3];
rz(-1.6549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90551463) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(-1.2515602) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(0.60607934) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(-1.1846503) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23683324) q[0];
sx q[0];
rz(-1.3968811) q[0];
sx q[0];
rz(-0.21763344) q[0];
rz(2.0834288) q[2];
sx q[2];
rz(-1.1626273) q[2];
sx q[2];
rz(-0.689091) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46890098) q[1];
sx q[1];
rz(-2.5395782) q[1];
sx q[1];
rz(-1.7778346) q[1];
x q[2];
rz(-0.70373669) q[3];
sx q[3];
rz(-0.53547137) q[3];
sx q[3];
rz(-2.717201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(0.63956368) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-0.48979428) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(0.52039352) q[0];
rz(-3.0386472) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6729352) q[0];
sx q[0];
rz(-0.22738439) q[0];
sx q[0];
rz(-0.44896455) q[0];
x q[1];
rz(-0.16291933) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(-1.7999072) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5853873) q[1];
sx q[1];
rz(-0.3857179) q[1];
sx q[1];
rz(-1.0005887) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95727386) q[3];
sx q[3];
rz(-2.0271218) q[3];
sx q[3];
rz(-0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(-0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(0.19764915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55006856) q[0];
sx q[0];
rz(-2.1854489) q[0];
sx q[0];
rz(-0.7556677) q[0];
rz(-pi) q[1];
rz(2.501776) q[2];
sx q[2];
rz(-0.13249995) q[2];
sx q[2];
rz(-2.8853531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3870961) q[1];
sx q[1];
rz(-1.6134141) q[1];
sx q[1];
rz(2.6864762) q[1];
rz(-0.88455172) q[3];
sx q[3];
rz(-0.93730799) q[3];
sx q[3];
rz(-1.3214878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-0.075604288) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(-2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(0.042073123) q[0];
rz(1.9627337) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.1694318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5116611) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(1.9566262) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61755007) q[2];
sx q[2];
rz(-2.6744161) q[2];
sx q[2];
rz(-3.0385366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9387503) q[1];
sx q[1];
rz(-0.65333594) q[1];
sx q[1];
rz(-0.57196879) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10294442) q[3];
sx q[3];
rz(-0.38984459) q[3];
sx q[3];
rz(-1.9600705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(2.5792504) q[0];
rz(1.7100122) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021488) q[0];
sx q[0];
rz(-1.095626) q[0];
sx q[0];
rz(0.094046353) q[0];
rz(-pi) q[1];
rz(-2.2675603) q[2];
sx q[2];
rz(-2.0119785) q[2];
sx q[2];
rz(-2.2639084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22842912) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(1.5921028) q[1];
x q[2];
rz(-1.0443586) q[3];
sx q[3];
rz(-1.8705771) q[3];
sx q[3];
rz(1.6696904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(2.8590554) q[2];
rz(2.1841168) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.7392993) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9937447) q[0];
sx q[0];
rz(-1.9980668) q[0];
sx q[0];
rz(-1.2742313) q[0];
x q[1];
rz(1.3627522) q[2];
sx q[2];
rz(-0.92804747) q[2];
sx q[2];
rz(-0.3026697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1101493) q[1];
sx q[1];
rz(-1.3182536) q[1];
sx q[1];
rz(2.4367711) q[1];
rz(-pi) q[2];
rz(-1.7761566) q[3];
sx q[3];
rz(-2.0911651) q[3];
sx q[3];
rz(-2.6092649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.6837439) q[2];
rz(0.66649377) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(1.9816459) q[0];
rz(3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265855) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(-1.3269781) q[0];
rz(-pi) q[1];
rz(-1.530933) q[2];
sx q[2];
rz(-2.3764052) q[2];
sx q[2];
rz(-2.6305692) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.479446) q[1];
sx q[1];
rz(-1.3227191) q[1];
sx q[1];
rz(-0.76920385) q[1];
rz(0.6108547) q[3];
sx q[3];
rz(-0.29963747) q[3];
sx q[3];
rz(-0.50869298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(-1.6147511) q[2];
rz(-2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-2.4894864) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(-2.3028159) q[2];
sx q[2];
rz(-1.5427187) q[2];
sx q[2];
rz(-1.6533921) q[2];
rz(0.53229971) q[3];
sx q[3];
rz(-2.0082062) q[3];
sx q[3];
rz(2.4826241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
