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
rz(1.8600872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0948439) q[0];
sx q[0];
rz(-1.2307414) q[0];
sx q[0];
rz(-0.21324555) q[0];
rz(-1.5413061) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(-2.8147547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.233853) q[1];
sx q[1];
rz(-2.6820161) q[1];
sx q[1];
rz(1.9563667) q[1];
x q[2];
rz(-1.1164066) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(-1.56173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32221258) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(1.5167351) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(-2.9852988) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(-1.6527536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50870897) q[0];
sx q[0];
rz(-1.4116305) q[0];
sx q[0];
rz(0.55523086) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82350125) q[2];
sx q[2];
rz(-0.89580065) q[2];
sx q[2];
rz(-2.0958971) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3540346) q[1];
sx q[1];
rz(-0.84317849) q[1];
sx q[1];
rz(2.3393199) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1366027) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(2.4581916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-0.29671159) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.5082224) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(1.189032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8231359) q[0];
sx q[0];
rz(-1.4134777) q[0];
sx q[0];
rz(-0.29968981) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5306773) q[2];
sx q[2];
rz(-1.8187858) q[2];
sx q[2];
rz(2.9330394) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51057286) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(-0.87639798) q[1];
x q[2];
rz(-0.69840188) q[3];
sx q[3];
rz(-1.9365371) q[3];
sx q[3];
rz(-1.6549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.2515602) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.8006515) q[0];
rz(2.5355133) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(1.1846503) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458762) q[0];
sx q[0];
rz(-1.7850951) q[0];
sx q[0];
rz(1.3927668) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6809533) q[2];
sx q[2];
rz(-1.1038291) q[2];
sx q[2];
rz(2.4796782) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6726917) q[1];
sx q[1];
rz(-0.60201445) q[1];
sx q[1];
rz(1.3637581) q[1];
rz(-pi) q[2];
rz(0.70373669) q[3];
sx q[3];
rz(-2.6061213) q[3];
sx q[3];
rz(-2.717201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-2.502029) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(3.0386472) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(2.4127814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6729352) q[0];
sx q[0];
rz(-0.22738439) q[0];
sx q[0];
rz(-0.44896455) q[0];
rz(-2.9786733) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(1.7999072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0923857) q[1];
sx q[1];
rz(-1.2485463) q[1];
sx q[1];
rz(-2.9258123) q[1];
x q[2];
rz(0.8647996) q[3];
sx q[3];
rz(-0.74665657) q[3];
sx q[3];
rz(1.4454696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(0.63306159) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.0387886) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(0.9206413) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(2.9439435) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5915241) q[0];
sx q[0];
rz(-2.1854489) q[0];
sx q[0];
rz(0.7556677) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0350787) q[2];
sx q[2];
rz(-1.649756) q[2];
sx q[2];
rz(0.67895141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16285322) q[1];
sx q[1];
rz(-2.0254685) q[1];
sx q[1];
rz(1.5233558) q[1];
x q[2];
rz(2.382155) q[3];
sx q[3];
rz(-2.1067838) q[3];
sx q[3];
rz(-2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7224715) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(-1.1694318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021381) q[0];
sx q[0];
rz(-1.5041782) q[0];
sx q[0];
rz(-1.4057977) q[0];
rz(-pi) q[1];
rz(0.39016907) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(0.90261501) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52160727) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(1.1779838) q[1];
rz(0.10294442) q[3];
sx q[3];
rz(-2.7517481) q[3];
sx q[3];
rz(-1.9600705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(3.1271093) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(-2.6838578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6477752) q[0];
sx q[0];
rz(-1.6543979) q[0];
sx q[0];
rz(-2.0477707) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5896795) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(2.1053095) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25719682) q[1];
sx q[1];
rz(-2.307502) q[1];
sx q[1];
rz(0.019330545) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3433414) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(-2.8727369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3056425) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(-1.4022934) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.8797849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5489053) q[0];
sx q[0];
rz(-1.3015916) q[0];
sx q[0];
rz(-0.44434987) q[0];
rz(-pi) q[1];
rz(-2.8724573) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(-0.035949635) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7463747) q[1];
sx q[1];
rz(-0.74133855) q[1];
sx q[1];
rz(0.37903255) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7997167) q[3];
sx q[3];
rz(-0.55594) q[3];
sx q[3];
rz(3.0059909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-1.0704401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72296952) q[0];
sx q[0];
rz(-1.4040935) q[0];
sx q[0];
rz(-2.3136501) q[0];
x q[1];
rz(3.1033377) q[2];
sx q[2];
rz(-2.3352211) q[2];
sx q[2];
rz(-2.5753266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1420319) q[1];
sx q[1];
rz(-0.83082098) q[1];
sx q[1];
rz(1.2318516) q[1];
rz(-2.8937267) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(-0.47249139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(1.6147511) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(1.0340446) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-1.5288011) q[2];
sx q[2];
rz(-2.4091346) q[2];
sx q[2];
rz(-0.11382881) q[2];
rz(-1.0735687) q[3];
sx q[3];
rz(-1.0931001) q[3];
sx q[3];
rz(1.1563374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
