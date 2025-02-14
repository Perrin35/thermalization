OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.003886) q[0];
sx q[0];
rz(-2.6258111) q[0];
sx q[0];
rz(-0.04096026) q[0];
rz(0.75086683) q[1];
sx q[1];
rz(-0.48315307) q[1];
sx q[1];
rz(-2.9274489) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0182604) q[0];
sx q[0];
rz(-1.7361479) q[0];
sx q[0];
rz(1.6613217) q[0];
rz(3.0532367) q[2];
sx q[2];
rz(-1.7680829) q[2];
sx q[2];
rz(-3.0947229) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63068608) q[1];
sx q[1];
rz(-1.4171466) q[1];
sx q[1];
rz(-2.2935778) q[1];
rz(-pi) q[2];
rz(2.6539702) q[3];
sx q[3];
rz(-2.0242526) q[3];
sx q[3];
rz(2.3359156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9302049) q[2];
sx q[2];
rz(-1.7519506) q[2];
sx q[2];
rz(-0.99838057) q[2];
rz(-1.3409279) q[3];
sx q[3];
rz(-2.0710335) q[3];
sx q[3];
rz(0.65795952) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9407161) q[0];
sx q[0];
rz(-2.2287892) q[0];
sx q[0];
rz(-0.87769133) q[0];
rz(-2.8857723) q[1];
sx q[1];
rz(-0.95953774) q[1];
sx q[1];
rz(-2.6928601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2923778) q[0];
sx q[0];
rz(-1.5149982) q[0];
sx q[0];
rz(-0.32955881) q[0];
rz(-pi) q[1];
rz(1.3569758) q[2];
sx q[2];
rz(-1.4656642) q[2];
sx q[2];
rz(2.2315764) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4325368) q[1];
sx q[1];
rz(-0.87608713) q[1];
sx q[1];
rz(2.0164151) q[1];
rz(-pi) q[2];
rz(-0.93495448) q[3];
sx q[3];
rz(-1.2024013) q[3];
sx q[3];
rz(1.660543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9933219) q[2];
sx q[2];
rz(-1.4485056) q[2];
sx q[2];
rz(-1.5987965) q[2];
rz(-2.350542) q[3];
sx q[3];
rz(-1.461921) q[3];
sx q[3];
rz(-2.3492474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032884447) q[0];
sx q[0];
rz(-2.5486163) q[0];
sx q[0];
rz(-1.8610883) q[0];
rz(-2.9325824) q[1];
sx q[1];
rz(-1.1384532) q[1];
sx q[1];
rz(-1.4271522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9379332) q[0];
sx q[0];
rz(-2.1224995) q[0];
sx q[0];
rz(2.0996086) q[0];
rz(-0.99781783) q[2];
sx q[2];
rz(-2.9429475) q[2];
sx q[2];
rz(-1.3494267) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5722981) q[1];
sx q[1];
rz(-1.8406017) q[1];
sx q[1];
rz(2.5266179) q[1];
x q[2];
rz(2.9455037) q[3];
sx q[3];
rz(-2.5485989) q[3];
sx q[3];
rz(0.65341572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3543432) q[2];
sx q[2];
rz(-2.3086583) q[2];
sx q[2];
rz(-0.99861097) q[2];
rz(-0.35501114) q[3];
sx q[3];
rz(-1.7576926) q[3];
sx q[3];
rz(-1.2242873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76096475) q[0];
sx q[0];
rz(-2.9686718) q[0];
sx q[0];
rz(-2.2291016) q[0];
rz(-1.6418705) q[1];
sx q[1];
rz(-1.8955756) q[1];
sx q[1];
rz(1.6778827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6448224) q[0];
sx q[0];
rz(-2.8319785) q[0];
sx q[0];
rz(2.7965318) q[0];
x q[1];
rz(-2.5374012) q[2];
sx q[2];
rz(-1.6407654) q[2];
sx q[2];
rz(1.7330012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6341175) q[1];
sx q[1];
rz(-1.1833131) q[1];
sx q[1];
rz(1.3204367) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8837508) q[3];
sx q[3];
rz(-2.3944562) q[3];
sx q[3];
rz(1.136387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0616167) q[2];
sx q[2];
rz(-1.4392263) q[2];
sx q[2];
rz(-1.7495135) q[2];
rz(-1.1213087) q[3];
sx q[3];
rz(-1.0871355) q[3];
sx q[3];
rz(-0.054718941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.012080972) q[0];
sx q[0];
rz(-1.5750146) q[0];
sx q[0];
rz(-2.6636301) q[0];
rz(-1.6732008) q[1];
sx q[1];
rz(-1.5719527) q[1];
sx q[1];
rz(-2.0782616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2713094) q[0];
sx q[0];
rz(-1.0510322) q[0];
sx q[0];
rz(-2.5867046) q[0];
x q[1];
rz(0.14866406) q[2];
sx q[2];
rz(-1.0188661) q[2];
sx q[2];
rz(2.8665598) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2995475) q[1];
sx q[1];
rz(-2.2631524) q[1];
sx q[1];
rz(1.0177708) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78467031) q[3];
sx q[3];
rz(-0.41483799) q[3];
sx q[3];
rz(-1.1718599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.31074506) q[2];
sx q[2];
rz(-1.9644535) q[2];
sx q[2];
rz(-1.2486521) q[2];
rz(2.2511075) q[3];
sx q[3];
rz(-1.9200446) q[3];
sx q[3];
rz(0.42731592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.86522) q[0];
sx q[0];
rz(-2.7557912) q[0];
sx q[0];
rz(1.2060097) q[0];
rz(1.6845866) q[1];
sx q[1];
rz(-0.99452368) q[1];
sx q[1];
rz(-2.1280033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26631334) q[0];
sx q[0];
rz(-1.5277098) q[0];
sx q[0];
rz(-0.049751374) q[0];
rz(-pi) q[1];
rz(0.59047575) q[2];
sx q[2];
rz(-1.9370859) q[2];
sx q[2];
rz(2.4136638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9070753) q[1];
sx q[1];
rz(-0.6014809) q[1];
sx q[1];
rz(0.3882577) q[1];
rz(-pi) q[2];
rz(-1.9928855) q[3];
sx q[3];
rz(-1.5076901) q[3];
sx q[3];
rz(-1.0168187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69063416) q[2];
sx q[2];
rz(-2.3834507) q[2];
sx q[2];
rz(-0.071619384) q[2];
rz(3.0028263) q[3];
sx q[3];
rz(-1.3774739) q[3];
sx q[3];
rz(0.57851401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0305158) q[0];
sx q[0];
rz(-2.0589622) q[0];
sx q[0];
rz(-1.6336596) q[0];
rz(-1.1208447) q[1];
sx q[1];
rz(-1.0034674) q[1];
sx q[1];
rz(0.20787636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250666) q[0];
sx q[0];
rz(-1.2057852) q[0];
sx q[0];
rz(2.9004296) q[0];
rz(-pi) q[1];
rz(-2.8817564) q[2];
sx q[2];
rz(-1.4186315) q[2];
sx q[2];
rz(0.92607385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1285591) q[1];
sx q[1];
rz(-0.39617523) q[1];
sx q[1];
rz(1.0349367) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7925345) q[3];
sx q[3];
rz(-2.9051637) q[3];
sx q[3];
rz(0.010836212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8230744) q[2];
sx q[2];
rz(-1.069331) q[2];
sx q[2];
rz(11/(7*pi)) q[2];
rz(2.9253166) q[3];
sx q[3];
rz(-1.6612771) q[3];
sx q[3];
rz(-1.1094619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08865393) q[0];
sx q[0];
rz(-1.5733938) q[0];
sx q[0];
rz(-1.1608359) q[0];
rz(2.898518) q[1];
sx q[1];
rz(-1.4832152) q[1];
sx q[1];
rz(-1.9879139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6223654) q[0];
sx q[0];
rz(-1.3524242) q[0];
sx q[0];
rz(1.1535991) q[0];
rz(1.3166381) q[2];
sx q[2];
rz(-1.1593737) q[2];
sx q[2];
rz(-2.9947077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8582871) q[1];
sx q[1];
rz(-1.6231307) q[1];
sx q[1];
rz(3.1238265) q[1];
x q[2];
rz(-2.1716464) q[3];
sx q[3];
rz(-1.613253) q[3];
sx q[3];
rz(0.24985931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.941075) q[2];
sx q[2];
rz(-2.2633471) q[2];
sx q[2];
rz(-0.78594691) q[2];
rz(3.1089697) q[3];
sx q[3];
rz(-2.2226108) q[3];
sx q[3];
rz(-2.6526764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0068479) q[0];
sx q[0];
rz(-2.5866046) q[0];
sx q[0];
rz(0.81163374) q[0];
rz(0.45400485) q[1];
sx q[1];
rz(-1.3521399) q[1];
sx q[1];
rz(-2.6754191) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9627169) q[0];
sx q[0];
rz(-2.6486925) q[0];
sx q[0];
rz(1.3702223) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2978372) q[2];
sx q[2];
rz(-2.5409219) q[2];
sx q[2];
rz(-2.7417083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5890903) q[1];
sx q[1];
rz(-1.2118166) q[1];
sx q[1];
rz(0.41429452) q[1];
rz(-pi) q[2];
rz(-2.2884263) q[3];
sx q[3];
rz(-0.93255842) q[3];
sx q[3];
rz(0.555942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0780645) q[2];
sx q[2];
rz(-1.177265) q[2];
sx q[2];
rz(-1.9130116) q[2];
rz(-2.581253) q[3];
sx q[3];
rz(-1.2816343) q[3];
sx q[3];
rz(1.3250215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90801564) q[0];
sx q[0];
rz(-0.19291872) q[0];
sx q[0];
rz(-2.8237901) q[0];
rz(-2.4494412) q[1];
sx q[1];
rz(-2.0798123) q[1];
sx q[1];
rz(-0.94921509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.229521) q[0];
sx q[0];
rz(-1.5186274) q[0];
sx q[0];
rz(2.0020141) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91429488) q[2];
sx q[2];
rz(-1.339448) q[2];
sx q[2];
rz(2.0670239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80133841) q[1];
sx q[1];
rz(-2.0726554) q[1];
sx q[1];
rz(1.5141634) q[1];
x q[2];
rz(-1.0383706) q[3];
sx q[3];
rz(-0.11603131) q[3];
sx q[3];
rz(0.86933245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.074241) q[2];
sx q[2];
rz(-2.1965252) q[2];
sx q[2];
rz(2.2828339) q[2];
rz(-0.37352118) q[3];
sx q[3];
rz(-1.7645323) q[3];
sx q[3];
rz(-0.080373272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141948) q[0];
sx q[0];
rz(-0.62794958) q[0];
sx q[0];
rz(-3.0218883) q[0];
rz(1.9627199) q[1];
sx q[1];
rz(-2.1736455) q[1];
sx q[1];
rz(-1.4812352) q[1];
rz(2.1356077) q[2];
sx q[2];
rz(-2.0567675) q[2];
sx q[2];
rz(-1.3153402) q[2];
rz(1.6681485) q[3];
sx q[3];
rz(-0.90383263) q[3];
sx q[3];
rz(-2.2624349) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
