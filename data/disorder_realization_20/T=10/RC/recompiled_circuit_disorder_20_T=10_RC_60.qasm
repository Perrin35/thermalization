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
rz(-2.1407054) q[0];
sx q[0];
rz(0.21240182) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40385383) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(1.9181197) q[0];
rz(0.69252695) q[2];
sx q[2];
rz(-1.5480969) q[2];
sx q[2];
rz(1.9164617) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68583381) q[1];
sx q[1];
rz(-1.403192) q[1];
sx q[1];
rz(-2.0007677) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.699502) q[3];
sx q[3];
rz(-1.5081815) q[3];
sx q[3];
rz(-2.6821729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(-2.8939698) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(-0.65555278) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(2.9852988) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(-1.4888391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1600906) q[0];
sx q[0];
rz(-1.0233876) q[0];
sx q[0];
rz(-1.7574969) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82350125) q[2];
sx q[2];
rz(-0.89580065) q[2];
sx q[2];
rz(1.0456955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3858151) q[1];
sx q[1];
rz(-2.1375244) q[1];
sx q[1];
rz(-0.66269006) q[1];
rz(0.1366027) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(2.4581916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(-2.8091649) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(-1.1214316) q[1];
sx q[1];
rz(1.189032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9246763) q[0];
sx q[0];
rz(-0.33736704) q[0];
sx q[0];
rz(0.49305537) q[0];
rz(-1.5306773) q[2];
sx q[2];
rz(-1.8187858) q[2];
sx q[2];
rz(0.2085533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6310198) q[1];
sx q[1];
rz(-2.2377308) q[1];
sx q[1];
rz(0.87639798) q[1];
rz(-0.69840188) q[3];
sx q[3];
rz(-1.2050556) q[3];
sx q[3];
rz(-1.48667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90551463) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(-1.8006515) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(-1.1846503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23683324) q[0];
sx q[0];
rz(-1.7447115) q[0];
sx q[0];
rz(-2.9239592) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0834288) q[2];
sx q[2];
rz(-1.9789654) q[2];
sx q[2];
rz(0.689091) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6726917) q[1];
sx q[1];
rz(-0.60201445) q[1];
sx q[1];
rz(1.7778346) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9373478) q[3];
sx q[3];
rz(-1.9703715) q[3];
sx q[3];
rz(-1.2031581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(0.63956368) q[2];
rz(-2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(0.52039352) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(-2.4127814) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9279014) q[0];
sx q[0];
rz(-1.7753082) q[0];
sx q[0];
rz(1.4707028) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16291933) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(1.3416854) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.049206991) q[1];
sx q[1];
rz(-1.8930463) q[1];
sx q[1];
rz(0.21578034) q[1];
rz(-pi) q[2];
rz(-0.8647996) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(-1.696123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52299243) q[0];
sx q[0];
rz(-0.97609659) q[0];
sx q[0];
rz(-2.3408875) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6502041) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(0.90027819) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8710528) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(-3.0448826) q[1];
x q[2];
rz(-2.4297653) q[3];
sx q[3];
rz(-0.89755745) q[3];
sx q[3];
rz(-0.87513808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(0.075604288) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.1694318) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4394546) q[0];
sx q[0];
rz(-1.5041782) q[0];
sx q[0];
rz(-1.735795) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5240426) q[2];
sx q[2];
rz(-2.6744161) q[2];
sx q[2];
rz(3.0385366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6199854) q[1];
sx q[1];
rz(-2.1072525) q[1];
sx q[1];
rz(1.1779838) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0386482) q[3];
sx q[3];
rz(-0.38984459) q[3];
sx q[3];
rz(-1.1815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(-0.45773488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4938175) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(-2.0477707) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5896795) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(1.0362831) q[2];
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
rz(-2.7982513) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(-2.8727369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.7392993) q[0];
rz(-2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.2618077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9937447) q[0];
sx q[0];
rz(-1.1435259) q[0];
sx q[0];
rz(-1.2742313) q[0];
rz(2.8724573) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(-3.105643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0314434) q[1];
sx q[1];
rz(-1.823339) q[1];
sx q[1];
rz(-2.4367711) q[1];
x q[2];
rz(1.3654361) q[3];
sx q[3];
rz(-2.0911651) q[3];
sx q[3];
rz(-2.6092649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.6377178) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2465729) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(1.0704401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1150072) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(-1.8146145) q[0];
rz(1.530933) q[2];
sx q[2];
rz(-0.76518744) q[2];
sx q[2];
rz(0.51102343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66026238) q[1];
sx q[1];
rz(-0.80033014) q[1];
sx q[1];
rz(-2.7923613) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.530738) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(-2.6328997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(-0.57957831) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(1.5288011) q[2];
sx q[2];
rz(-0.73245807) q[2];
sx q[2];
rz(3.0277638) q[2];
rz(2.0680239) q[3];
sx q[3];
rz(-1.0931001) q[3];
sx q[3];
rz(1.1563374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
