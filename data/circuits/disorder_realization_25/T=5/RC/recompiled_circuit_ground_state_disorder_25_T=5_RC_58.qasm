OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8018262) q[0];
sx q[0];
rz(-2.5751994) q[0];
sx q[0];
rz(-1.1993778) q[0];
rz(-0.83130032) q[1];
sx q[1];
rz(4.755862) q[1];
sx q[1];
rz(8.9250467) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3008376) q[0];
sx q[0];
rz(-0.84512701) q[0];
sx q[0];
rz(0.27780224) q[0];
rz(1.209916) q[2];
sx q[2];
rz(-1.7583876) q[2];
sx q[2];
rz(-2.9071992) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.079165212) q[1];
sx q[1];
rz(-1.2794331) q[1];
sx q[1];
rz(-2.3645762) q[1];
rz(-1.6263032) q[3];
sx q[3];
rz(-2.422547) q[3];
sx q[3];
rz(-2.2628502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1685593) q[2];
sx q[2];
rz(-1.3848105) q[2];
sx q[2];
rz(-0.92156571) q[2];
rz(1.5714931) q[3];
sx q[3];
rz(-2.470033) q[3];
sx q[3];
rz(0.69935548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8916931) q[0];
sx q[0];
rz(-2.1136916) q[0];
sx q[0];
rz(1.6722884) q[0];
rz(1.6647313) q[1];
sx q[1];
rz(-1.3427443) q[1];
sx q[1];
rz(0.5161759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1219541) q[0];
sx q[0];
rz(-0.40651822) q[0];
sx q[0];
rz(-1.818114) q[0];
rz(-2.8968229) q[2];
sx q[2];
rz(-1.2483856) q[2];
sx q[2];
rz(0.66157326) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6440638) q[1];
sx q[1];
rz(-1.6402771) q[1];
sx q[1];
rz(-2.890088) q[1];
rz(-2.1352876) q[3];
sx q[3];
rz(-0.7745452) q[3];
sx q[3];
rz(1.2513127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76253647) q[2];
sx q[2];
rz(-1.9409337) q[2];
sx q[2];
rz(2.8442247) q[2];
rz(2.6250046) q[3];
sx q[3];
rz(-1.9818431) q[3];
sx q[3];
rz(0.62492257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49279889) q[0];
sx q[0];
rz(-0.53557098) q[0];
sx q[0];
rz(0.18185644) q[0];
rz(0.44405538) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(3.0603538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59653283) q[0];
sx q[0];
rz(-1.1756011) q[0];
sx q[0];
rz(1.6020653) q[0];
x q[1];
rz(-1.9047348) q[2];
sx q[2];
rz(-2.4231964) q[2];
sx q[2];
rz(0.12018724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.013147203) q[1];
sx q[1];
rz(-2.0390292) q[1];
sx q[1];
rz(0.20856671) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3879561) q[3];
sx q[3];
rz(-2.5732079) q[3];
sx q[3];
rz(0.050384132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9265201) q[2];
sx q[2];
rz(-2.7762065) q[2];
sx q[2];
rz(0.60058769) q[2];
rz(1.8961204) q[3];
sx q[3];
rz(-2.5966817) q[3];
sx q[3];
rz(0.039552461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0092225) q[0];
sx q[0];
rz(-0.70766574) q[0];
sx q[0];
rz(-0.02455499) q[0];
rz(0.37697667) q[1];
sx q[1];
rz(-1.8507277) q[1];
sx q[1];
rz(2.7797508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.637476) q[0];
sx q[0];
rz(-1.8808805) q[0];
sx q[0];
rz(1.732336) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.078388647) q[2];
sx q[2];
rz(-1.9447717) q[2];
sx q[2];
rz(-2.0736935) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1805318) q[1];
sx q[1];
rz(-2.1994563) q[1];
sx q[1];
rz(-1.3500477) q[1];
rz(2.0482158) q[3];
sx q[3];
rz(-1.1934115) q[3];
sx q[3];
rz(-1.4879102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5234066) q[2];
sx q[2];
rz(-1.8768825) q[2];
sx q[2];
rz(-2.3183909) q[2];
rz(-0.19436714) q[3];
sx q[3];
rz(-2.7394962) q[3];
sx q[3];
rz(0.91485867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0219367) q[0];
sx q[0];
rz(-1.5384262) q[0];
sx q[0];
rz(0.62171474) q[0];
rz(0.74553982) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(3.0527557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5058891) q[0];
sx q[0];
rz(-2.997809) q[0];
sx q[0];
rz(-2.6583688) q[0];
rz(-pi) q[1];
rz(0.88314573) q[2];
sx q[2];
rz(-2.7252203) q[2];
sx q[2];
rz(-1.0036381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6703265) q[1];
sx q[1];
rz(-2.2818292) q[1];
sx q[1];
rz(1.9416351) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1301378) q[3];
sx q[3];
rz(-0.36962907) q[3];
sx q[3];
rz(-0.91418788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9125774) q[2];
sx q[2];
rz(-1.503016) q[2];
sx q[2];
rz(1.5360606) q[2];
rz(0.80794263) q[3];
sx q[3];
rz(-0.84351051) q[3];
sx q[3];
rz(1.1725175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.262893) q[0];
sx q[0];
rz(-1.7182588) q[0];
sx q[0];
rz(0.95243564) q[0];
rz(1.7772504) q[1];
sx q[1];
rz(-1.9352501) q[1];
sx q[1];
rz(-0.84699026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5085707) q[0];
sx q[0];
rz(-1.7135597) q[0];
sx q[0];
rz(-1.4748013) q[0];
rz(-pi) q[1];
rz(2.3479727) q[2];
sx q[2];
rz(-1.334903) q[2];
sx q[2];
rz(-1.867804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1141194) q[1];
sx q[1];
rz(-2.4714258) q[1];
sx q[1];
rz(-0.7188188) q[1];
rz(-pi) q[2];
rz(-1.5853508) q[3];
sx q[3];
rz(-1.6739917) q[3];
sx q[3];
rz(-1.9510912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18818894) q[2];
sx q[2];
rz(-1.0651411) q[2];
sx q[2];
rz(2.7597001) q[2];
rz(2.0415107) q[3];
sx q[3];
rz(-2.3298658) q[3];
sx q[3];
rz(-2.7303117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57322684) q[0];
sx q[0];
rz(-1.3245642) q[0];
sx q[0];
rz(0.49760094) q[0];
rz(2.7865903) q[1];
sx q[1];
rz(-1.4811131) q[1];
sx q[1];
rz(-0.76638806) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7805745) q[0];
sx q[0];
rz(-1.0583335) q[0];
sx q[0];
rz(1.6105525) q[0];
x q[1];
rz(-2.1649579) q[2];
sx q[2];
rz(-2.4554002) q[2];
sx q[2];
rz(-0.90168574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.97624) q[1];
sx q[1];
rz(-1.2818977) q[1];
sx q[1];
rz(1.6331178) q[1];
x q[2];
rz(-1.822352) q[3];
sx q[3];
rz(-0.90853229) q[3];
sx q[3];
rz(1.4317625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4088722) q[2];
sx q[2];
rz(-1.5343821) q[2];
sx q[2];
rz(-1.4917779) q[2];
rz(2.0626227) q[3];
sx q[3];
rz(-1.8337199) q[3];
sx q[3];
rz(-2.5274247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33477467) q[0];
sx q[0];
rz(-0.65009999) q[0];
sx q[0];
rz(-0.41710576) q[0];
rz(-0.053622309) q[1];
sx q[1];
rz(-1.6157179) q[1];
sx q[1];
rz(-2.0448304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039711397) q[0];
sx q[0];
rz(-2.5371612) q[0];
sx q[0];
rz(1.7834316) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.057889537) q[2];
sx q[2];
rz(-2.0308924) q[2];
sx q[2];
rz(-1.5050846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1138823) q[1];
sx q[1];
rz(-1.7358297) q[1];
sx q[1];
rz(-1.5636958) q[1];
x q[2];
rz(-1.4194059) q[3];
sx q[3];
rz(-0.89821363) q[3];
sx q[3];
rz(-2.0167493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2199478) q[2];
sx q[2];
rz(-1.767445) q[2];
sx q[2];
rz(-2.5448223) q[2];
rz(2.2611332) q[3];
sx q[3];
rz(-1.7170649) q[3];
sx q[3];
rz(2.5318291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45744687) q[0];
sx q[0];
rz(-2.5992751) q[0];
sx q[0];
rz(2.82161) q[0];
rz(-0.34842247) q[1];
sx q[1];
rz(-0.44487822) q[1];
sx q[1];
rz(0.50331032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43673462) q[0];
sx q[0];
rz(-2.9665369) q[0];
sx q[0];
rz(1.0147834) q[0];
rz(-2.3186683) q[2];
sx q[2];
rz(-1.9562998) q[2];
sx q[2];
rz(-0.67367919) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.383068) q[1];
sx q[1];
rz(-0.43979859) q[1];
sx q[1];
rz(-1.7563933) q[1];
rz(-pi) q[2];
rz(-1.4473404) q[3];
sx q[3];
rz(-1.0528101) q[3];
sx q[3];
rz(-0.65604612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8866426) q[2];
sx q[2];
rz(-0.4526259) q[2];
sx q[2];
rz(-0.63549271) q[2];
rz(-1.6057711) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(0.087285727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4601634) q[0];
sx q[0];
rz(-0.59578139) q[0];
sx q[0];
rz(1.7728565) q[0];
rz(-0.5303371) q[1];
sx q[1];
rz(-1.4884721) q[1];
sx q[1];
rz(-2.368685) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094495234) q[0];
sx q[0];
rz(-0.90029923) q[0];
sx q[0];
rz(2.0740202) q[0];
rz(-pi) q[1];
rz(-1.673416) q[2];
sx q[2];
rz(-0.78067242) q[2];
sx q[2];
rz(0.80678015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49131888) q[1];
sx q[1];
rz(-1.4161619) q[1];
sx q[1];
rz(2.0238961) q[1];
rz(-pi) q[2];
rz(-2.8498296) q[3];
sx q[3];
rz(-2.004753) q[3];
sx q[3];
rz(-1.4572373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7698001) q[2];
sx q[2];
rz(-0.95855203) q[2];
sx q[2];
rz(-2.8690763) q[2];
rz(2.6767327) q[3];
sx q[3];
rz(-1.9404989) q[3];
sx q[3];
rz(-2.4242937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40806017) q[0];
sx q[0];
rz(-1.6049186) q[0];
sx q[0];
rz(-1.4419755) q[0];
rz(1.6170665) q[1];
sx q[1];
rz(-2.1433612) q[1];
sx q[1];
rz(-2.8467766) q[1];
rz(-2.8256399) q[2];
sx q[2];
rz(-1.3262333) q[2];
sx q[2];
rz(-1.4953351) q[2];
rz(2.9197027) q[3];
sx q[3];
rz(-1.1380914) q[3];
sx q[3];
rz(-1.9590845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
