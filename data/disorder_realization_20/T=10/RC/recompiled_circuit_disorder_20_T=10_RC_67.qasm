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
rz(3.9290805) q[1];
sx q[1];
rz(10.706283) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0948439) q[0];
sx q[0];
rz(-1.2307414) q[0];
sx q[0];
rz(2.9283471) q[0];
rz(1.6002866) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(0.32683795) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80860955) q[1];
sx q[1];
rz(-1.994351) q[1];
sx q[1];
rz(-2.9575648) q[1];
x q[2];
rz(-3.078457) q[3];
sx q[3];
rz(-1.4423443) q[3];
sx q[3];
rz(2.0383143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8193801) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(1.5167351) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-2.4860399) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(2.9852988) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(-1.6527536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.981502) q[0];
sx q[0];
rz(-1.0233876) q[0];
sx q[0];
rz(-1.3840958) q[0];
x q[1];
rz(0.70398753) q[2];
sx q[2];
rz(-0.96103243) q[2];
sx q[2];
rz(-0.068254452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21178791) q[1];
sx q[1];
rz(-1.025052) q[1];
sx q[1];
rz(-2.2499229) q[1];
x q[2];
rz(-2.923521) q[3];
sx q[3];
rz(-3.0017188) q[3];
sx q[3];
rz(-1.1034031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8504101) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082224) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(-0.53031522) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(-1.9525607) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30071242) q[0];
sx q[0];
rz(-1.2749201) q[0];
sx q[0];
rz(1.7353252) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8934115) q[2];
sx q[2];
rz(-1.6096874) q[2];
sx q[2];
rz(-1.789202) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6310198) q[1];
sx q[1];
rz(-2.2377308) q[1];
sx q[1];
rz(-2.2651947) q[1];
rz(-2.6044106) q[3];
sx q[3];
rz(-0.77386412) q[3];
sx q[3];
rz(-2.8230132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90551463) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.8900324) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(-0.35513487) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.8006515) q[0];
rz(0.60607934) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(1.9569424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.143648) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(-0.68302897) q[0];
x q[1];
rz(2.2934154) q[2];
sx q[2];
rz(-0.64372534) q[2];
sx q[2];
rz(-1.6456457) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6726917) q[1];
sx q[1];
rz(-2.5395782) q[1];
sx q[1];
rz(1.7778346) q[1];
x q[2];
rz(-2.7167927) q[3];
sx q[3];
rz(-1.9072755) q[3];
sx q[3];
rz(-0.51589033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23345315) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(-2.6211991) q[0];
rz(3.0386472) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-2.4127814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6729352) q[0];
sx q[0];
rz(-0.22738439) q[0];
sx q[0];
rz(-0.44896455) q[0];
rz(0.8466709) q[2];
sx q[2];
rz(-1.4482822) q[2];
sx q[2];
rz(-3.0201437) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0923857) q[1];
sx q[1];
rz(-1.8930463) q[1];
sx q[1];
rz(-2.9258123) q[1];
rz(-pi) q[2];
rz(0.54069424) q[3];
sx q[3];
rz(-1.027642) q[3];
sx q[3];
rz(2.3054996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(-2.4318802) q[2];
rz(0.63306159) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028041) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(2.9439435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5699428) q[0];
sx q[0];
rz(-2.2075704) q[0];
sx q[0];
rz(-2.3417579) q[0];
rz(-pi) q[1];
rz(-0.10651393) q[2];
sx q[2];
rz(-1.649756) q[2];
sx q[2];
rz(0.67895141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9787394) q[1];
sx q[1];
rz(-2.0254685) q[1];
sx q[1];
rz(1.6182369) q[1];
rz(0.71182735) q[3];
sx q[3];
rz(-0.89755745) q[3];
sx q[3];
rz(-0.87513808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.027823042) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.9627337) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(1.1694318) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4394546) q[0];
sx q[0];
rz(-1.6374145) q[0];
sx q[0];
rz(-1.4057977) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7514236) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(-2.2389776) q[2];
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
rz(-1.9636088) q[1];
x q[2];
rz(0.10294442) q[3];
sx q[3];
rz(-2.7517481) q[3];
sx q[3];
rz(-1.9600705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(-0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91677529) q[0];
sx q[0];
rz(-2.6579034) q[0];
sx q[0];
rz(-1.7513357) q[0];
x q[1];
rz(2.5896795) q[2];
sx q[2];
rz(-0.95171026) q[2];
sx q[2];
rz(2.1053095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.815005) q[1];
sx q[1];
rz(-1.5564789) q[1];
sx q[1];
rz(-0.83399764) q[1];
rz(-0.3433414) q[3];
sx q[3];
rz(-2.0715054) q[3];
sx q[3];
rz(-2.8727369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83595014) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(-2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.8797849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6540419) q[0];
sx q[0];
rz(-0.51484171) q[0];
sx q[0];
rz(-2.5709855) q[0];
x q[1];
rz(-0.26913531) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(3.105643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.251235) q[1];
sx q[1];
rz(-0.8926549) q[1];
sx q[1];
rz(-1.2441586) q[1];
x q[2];
rz(-1.3654361) q[3];
sx q[3];
rz(-1.0504275) q[3];
sx q[3];
rz(-2.6092649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4385779) q[2];
sx q[2];
rz(-1.6377178) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.9816459) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72296952) q[0];
sx q[0];
rz(-1.4040935) q[0];
sx q[0];
rz(2.3136501) q[0];
rz(-pi) q[1];
rz(-1.530933) q[2];
sx q[2];
rz(-2.3764052) q[2];
sx q[2];
rz(0.51102343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.66026238) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(-2.7923613) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24786592) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(0.47249139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(-2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-0.8387768) q[2];
sx q[2];
rz(-1.598874) q[2];
sx q[2];
rz(1.4882006) q[2];
rz(-2.6092929) q[3];
sx q[3];
rz(-2.0082062) q[3];
sx q[3];
rz(2.4826241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
