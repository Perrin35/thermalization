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
rz(-2.7377388) q[0];
sx q[0];
rz(-1.7716584) q[0];
sx q[0];
rz(1.9181197) q[0];
rz(1.6002866) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(-2.8147547) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.233853) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(1.185226) q[1];
rz(-pi) q[2];
rz(0.063135677) q[3];
sx q[3];
rz(-1.6992484) q[3];
sx q[3];
rz(-2.0383143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4441929) q[0];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8119368) q[0];
sx q[0];
rz(-0.57528472) q[0];
sx q[0];
rz(-2.8459957) q[0];
rz(-pi) q[1];
rz(-0.70398753) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(3.0733382) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7557775) q[1];
sx q[1];
rz(-1.0040682) q[1];
sx q[1];
rz(-2.4789026) q[1];
rz(0.1366027) q[3];
sx q[3];
rz(-1.6009637) q[3];
sx q[3];
rz(0.68340106) q[3];
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
rz(2.8091649) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082224) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(-2.6112774) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(-1.189032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9246763) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-2.6485373) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9844935) q[2];
sx q[2];
rz(-0.25114775) q[2];
sx q[2];
rz(-2.7709393) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6310198) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(2.2651947) q[1];
rz(-pi) q[2];
rz(2.0344814) q[3];
sx q[3];
rz(-2.2148805) q[3];
sx q[3];
rz(-0.37582276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-2.7864578) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.1846503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957164) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(1.3927668) q[0];
rz(-0.84817727) q[2];
sx q[2];
rz(-2.4978673) q[2];
sx q[2];
rz(-1.4959469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2733172) q[1];
sx q[1];
rz(-1.6874716) q[1];
sx q[1];
rz(0.97881808) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9373478) q[3];
sx q[3];
rz(-1.1712211) q[3];
sx q[3];
rz(-1.2031581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(2.502029) q[2];
rz(-2.284164) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(-0.52039352) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8048808) q[0];
sx q[0];
rz(-1.4727955) q[0];
sx q[0];
rz(0.20551198) q[0];
rz(-1.38703) q[2];
sx q[2];
rz(-0.73256058) q[2];
sx q[2];
rz(-1.5549321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5908969) q[1];
sx q[1];
rz(-1.7753074) q[1];
sx q[1];
rz(1.900161) q[1];
rz(0.8647996) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(1.696123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1028041) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(-2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-2.9439435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716498) q[0];
sx q[0];
rz(-0.93402223) q[0];
sx q[0];
rz(-2.3417579) q[0];
x q[1];
rz(-1.6502041) q[2];
sx q[2];
rz(-1.4646155) q[2];
sx q[2];
rz(-0.90027819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8710528) q[1];
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
rz(-pi) q[1];
rz(-3.1137696) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(3.0659884) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.1694318) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6299316) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(1.1849665) q[0];
x q[1];
rz(2.5240426) q[2];
sx q[2];
rz(-2.6744161) q[2];
sx q[2];
rz(3.0385366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.84050256) q[1];
sx q[1];
rz(-1.9060578) q[1];
sx q[1];
rz(-2.5696978) q[1];
rz(-0.10294442) q[3];
sx q[3];
rz(-0.38984459) q[3];
sx q[3];
rz(-1.9600705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(3.1271093) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-2.5884132) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977215) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(2.5792504) q[0];
rz(1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(-0.45773488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021488) q[0];
sx q[0];
rz(-1.095626) q[0];
sx q[0];
rz(-0.094046353) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2051815) q[2];
sx q[2];
rz(-2.3371156) q[2];
sx q[2];
rz(2.9203383) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3265877) q[1];
sx q[1];
rz(-1.5851138) q[1];
sx q[1];
rz(2.307595) q[1];
rz(-pi) q[2];
rz(-2.0972341) q[3];
sx q[3];
rz(-1.8705771) q[3];
sx q[3];
rz(-1.6696904) q[3];
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
rz(-0.28253728) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95865059) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48755074) q[0];
sx q[0];
rz(-0.51484171) q[0];
sx q[0];
rz(-0.57060711) q[0];
rz(-pi) q[1];
rz(1.7788404) q[2];
sx q[2];
rz(-2.2135452) q[2];
sx q[2];
rz(-0.3026697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1101493) q[1];
sx q[1];
rz(-1.3182536) q[1];
sx q[1];
rz(-2.4367711) q[1];
rz(0.34187596) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(0.13560175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4385779) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(-2.4750989) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-1.0704401) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4186231) q[0];
sx q[0];
rz(-1.4040935) q[0];
sx q[0];
rz(-2.3136501) q[0];
rz(-pi) q[1];
rz(-0.038254914) q[2];
sx q[2];
rz(-0.80637156) q[2];
sx q[2];
rz(2.5753266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9995607) q[1];
sx q[1];
rz(-2.3107717) q[1];
sx q[1];
rz(-1.2318516) q[1];
x q[2];
rz(-2.530738) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(-2.6328997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(-1.6147511) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6282745) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(-1.5288011) q[2];
sx q[2];
rz(-2.4091346) q[2];
sx q[2];
rz(-0.11382881) q[2];
rz(0.74450775) q[3];
sx q[3];
rz(-2.4662938) q[3];
sx q[3];
rz(0.288356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
