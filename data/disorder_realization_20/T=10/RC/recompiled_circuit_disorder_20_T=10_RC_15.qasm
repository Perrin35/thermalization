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
rz(0.71495932) q[1];
sx q[1];
rz(-2.3541048) q[1];
sx q[1];
rz(1.2815055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6706657) q[0];
sx q[0];
rz(-0.39917329) q[0];
sx q[0];
rz(2.1098718) q[0];
rz(-pi) q[1];
rz(-1.5413061) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(-0.32683795) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.1164066) q[3];
sx q[3];
rz(-0.14305275) q[3];
sx q[3];
rz(-1.56173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(-2.8939698) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(2.4860399) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(-0.15629388) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(-1.4888391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3296559) q[0];
sx q[0];
rz(-2.5663079) q[0];
sx q[0];
rz(2.8459957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.312617) q[2];
sx q[2];
rz(-2.1301221) q[2];
sx q[2];
rz(-2.0916569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3858151) q[1];
sx q[1];
rz(-2.1375244) q[1];
sx q[1];
rz(0.66269006) q[1];
rz(0.1366027) q[3];
sx q[3];
rz(-1.6009637) q[3];
sx q[3];
rz(-2.4581916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(-0.53031522) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(1.9525607) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3184568) q[0];
sx q[0];
rz(-1.4134777) q[0];
sx q[0];
rz(-2.8419028) q[0];
rz(-pi) q[1];
rz(0.15709917) q[2];
sx q[2];
rz(-2.8904449) q[2];
sx q[2];
rz(2.7709393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5570453) q[1];
sx q[1];
rz(-2.097633) q[1];
sx q[1];
rz(-2.3440863) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69840188) q[3];
sx q[3];
rz(-1.9365371) q[3];
sx q[3];
rz(-1.48667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(-2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77804756) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(2.5355133) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.1846503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458762) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(1.7488259) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6809533) q[2];
sx q[2];
rz(-2.0377635) q[2];
sx q[2];
rz(-0.66191445) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8682755) q[1];
sx q[1];
rz(-1.454121) q[1];
sx q[1];
rz(-0.97881808) q[1];
rz(-2.437856) q[3];
sx q[3];
rz(-2.6061213) q[3];
sx q[3];
rz(0.4243917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(-0.8574287) q[3];
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
rz(-pi) q[2];
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
rz(-2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(-0.52039352) q[0];
rz(-3.0386472) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-0.72881126) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33671185) q[0];
sx q[0];
rz(-1.4727955) q[0];
sx q[0];
rz(2.9360807) q[0];
rz(-pi) q[1];
rz(2.9786733) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(1.7999072) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5506958) q[1];
sx q[1];
rz(-1.7753074) q[1];
sx q[1];
rz(1.900161) q[1];
rz(-pi) q[2];
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
sx q[1];
rz(pi/2) q[1];
rz(-2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(-2.5085311) q[3];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(0.88371712) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(0.19764915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52299243) q[0];
sx q[0];
rz(-2.1654961) q[0];
sx q[0];
rz(-2.3408875) q[0];
rz(-1.4913885) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(-0.90027819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27053988) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(0.096710042) q[1];
x q[2];
rz(-0.71182735) q[3];
sx q[3];
rz(-0.89755745) q[3];
sx q[3];
rz(0.87513808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(-0.82908019) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
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
rz(-1.1584287) q[1];
sx q[1];
rz(1.9721608) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4394546) q[0];
sx q[0];
rz(-1.6374145) q[0];
sx q[0];
rz(1.735795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2866227) q[2];
sx q[2];
rz(-1.1948164) q[2];
sx q[2];
rz(-2.3665731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.52160727) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(1.1779838) q[1];
x q[2];
rz(-0.38798214) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(-2.8475873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(3.1271093) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-2.5792504) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(-0.45773488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4938175) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(-2.0477707) q[0];
x q[1];
rz(-0.93641113) q[2];
sx q[2];
rz(-0.80447703) q[2];
sx q[2];
rz(2.9203383) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9131635) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(-1.5494898) q[1];
rz(-pi) q[2];
rz(-0.3433414) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(-0.2688558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(0.28253728) q[2];
rz(0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(1.7392993) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(1.8797849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48755074) q[0];
sx q[0];
rz(-0.51484171) q[0];
sx q[0];
rz(-2.5709855) q[0];
x q[1];
rz(-0.65323921) q[2];
sx q[2];
rz(-1.7368894) q[2];
sx q[2];
rz(1.3939898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1101493) q[1];
sx q[1];
rz(-1.3182536) q[1];
sx q[1];
rz(2.4367711) q[1];
rz(2.6120139) q[3];
sx q[3];
rz(-1.7486608) q[3];
sx q[3];
rz(-1.9999268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2465729) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.9816459) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(-2.0711526) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1150072) q[0];
sx q[0];
rz(-0.75782776) q[0];
sx q[0];
rz(-1.3269781) q[0];
x q[1];
rz(-1.6106597) q[2];
sx q[2];
rz(-0.76518744) q[2];
sx q[2];
rz(-2.6305692) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66214661) q[1];
sx q[1];
rz(-1.3227191) q[1];
sx q[1];
rz(-0.76920385) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8937267) q[3];
sx q[3];
rz(-1.4006747) q[3];
sx q[3];
rz(-2.6691013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(0.57957831) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(0.65210623) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(1.6127916) q[2];
sx q[2];
rz(-2.4091346) q[2];
sx q[2];
rz(-0.11382881) q[2];
rz(-0.74450775) q[3];
sx q[3];
rz(-0.67529882) q[3];
sx q[3];
rz(-2.8532367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
