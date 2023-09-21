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
rz(-0.21240182) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7377388) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(-1.9181197) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5413061) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(2.8147547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68583381) q[1];
sx q[1];
rz(-1.7384006) q[1];
sx q[1];
rz(2.0007677) q[1];
rz(-pi) q[2];
rz(-0.063135677) q[3];
sx q[3];
rz(-1.4423443) q[3];
sx q[3];
rz(1.1032784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(-2.9852988) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(-1.4888391) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3296559) q[0];
sx q[0];
rz(-2.5663079) q[0];
sx q[0];
rz(0.29559691) q[0];
rz(2.312617) q[2];
sx q[2];
rz(-2.1301221) q[2];
sx q[2];
rz(1.0499357) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3858151) q[1];
sx q[1];
rz(-1.0040682) q[1];
sx q[1];
rz(0.66269006) q[1];
rz(-pi) q[2];
rz(0.21807166) q[3];
sx q[3];
rz(-0.13987386) q[3];
sx q[3];
rz(1.1034031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(0.29671159) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(-2.6112774) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.189032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9246763) q[0];
sx q[0];
rz(-0.33736704) q[0];
sx q[0];
rz(2.6485373) q[0];
rz(-pi) q[1];
rz(1.5306773) q[2];
sx q[2];
rz(-1.3228068) q[2];
sx q[2];
rz(0.2085533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6991068) q[1];
sx q[1];
rz(-0.92256303) q[1];
sx q[1];
rz(0.68251619) q[1];
rz(0.53718209) q[3];
sx q[3];
rz(-0.77386412) q[3];
sx q[3];
rz(0.31857946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90551463) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.2515602) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.8006515) q[0];
rz(2.5355133) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(-1.9569424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047594) q[0];
sx q[0];
rz(-1.3968811) q[0];
sx q[0];
rz(-0.21763344) q[0];
rz(-pi) q[1];
rz(-2.6809533) q[2];
sx q[2];
rz(-1.1038291) q[2];
sx q[2];
rz(2.4796782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6726917) q[1];
sx q[1];
rz(-0.60201445) q[1];
sx q[1];
rz(1.7778346) q[1];
rz(-pi) q[2];
rz(-1.9373478) q[3];
sx q[3];
rz(-1.1712211) q[3];
sx q[3];
rz(1.9384345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
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
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63067591) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(2.4127814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136912) q[0];
sx q[0];
rz(-1.3662845) q[0];
sx q[0];
rz(-1.4707028) q[0];
rz(-pi) q[1];
rz(0.8466709) q[2];
sx q[2];
rz(-1.4482822) q[2];
sx q[2];
rz(0.12144897) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.049206991) q[1];
sx q[1];
rz(-1.2485463) q[1];
sx q[1];
rz(0.21578034) q[1];
x q[2];
rz(2.1843188) q[3];
sx q[3];
rz(-2.0271218) q[3];
sx q[3];
rz(0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.0795152) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6186002) q[0];
sx q[0];
rz(-0.97609659) q[0];
sx q[0];
rz(-0.80070514) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6502041) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(-2.2413145) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27053988) q[1];
sx q[1];
rz(-0.45696837) q[1];
sx q[1];
rz(0.096710042) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2570409) q[3];
sx q[3];
rz(-2.2042847) q[3];
sx q[3];
rz(-1.3214878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(2.3125125) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(3.0995195) q[0];
rz(1.9627337) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.9721608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1202576) q[0];
sx q[0];
rz(-1.406167) q[0];
sx q[0];
rz(0.067532587) q[0];
x q[1];
rz(-1.2866227) q[2];
sx q[2];
rz(-1.1948164) q[2];
sx q[2];
rz(-2.3665731) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6199854) q[1];
sx q[1];
rz(-2.1072525) q[1];
sx q[1];
rz(-1.9636088) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10294442) q[3];
sx q[3];
rz(-0.38984459) q[3];
sx q[3];
rz(1.1815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977215) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(2.6838578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4938175) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(2.0477707) q[0];
x q[1];
rz(-2.5896795) q[2];
sx q[2];
rz(-0.95171026) q[2];
sx q[2];
rz(1.0362831) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22842912) q[1];
sx q[1];
rz(-0.73691165) q[1];
sx q[1];
rz(1.5921028) q[1];
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
x q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(2.8590554) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(-1.2618077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14784797) q[0];
sx q[0];
rz(-1.1435259) q[0];
sx q[0];
rz(1.8673613) q[0];
rz(-pi) q[1];
rz(-1.7788404) q[2];
sx q[2];
rz(-2.2135452) q[2];
sx q[2];
rz(-2.838923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8903576) q[1];
sx q[1];
rz(-0.8926549) q[1];
sx q[1];
rz(1.8974341) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7997167) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(-0.13560175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70301473) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2465729) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.9816459) q[0];
rz(3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(2.0711526) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265855) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(1.3269781) q[0];
rz(-pi) q[1];
rz(3.1033377) q[2];
sx q[2];
rz(-2.3352211) q[2];
sx q[2];
rz(-2.5753266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66026238) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(-2.7923613) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24786592) q[3];
sx q[3];
rz(-1.4006747) q[3];
sx q[3];
rz(2.6691013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3728309) q[2];
sx q[2];
rz(-2.5925345) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51331818) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(0.037739567) q[2];
sx q[2];
rz(-0.83913091) q[2];
sx q[2];
rz(-0.057374949) q[2];
rz(2.6092929) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];