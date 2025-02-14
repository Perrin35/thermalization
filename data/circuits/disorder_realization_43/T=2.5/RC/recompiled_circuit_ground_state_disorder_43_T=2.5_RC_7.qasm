OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3620152) q[0];
sx q[0];
rz(-2.9147122) q[0];
sx q[0];
rz(1.3751295) q[0];
rz(3.0472164) q[1];
sx q[1];
rz(-0.91369319) q[1];
sx q[1];
rz(-2.8606666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048427933) q[0];
sx q[0];
rz(-0.68138382) q[0];
sx q[0];
rz(1.8244034) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9380577) q[2];
sx q[2];
rz(-2.8312771) q[2];
sx q[2];
rz(-0.65904891) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1350159) q[1];
sx q[1];
rz(-2.6408051) q[1];
sx q[1];
rz(2.1441475) q[1];
rz(-pi) q[2];
rz(2.5689718) q[3];
sx q[3];
rz(-1.2702281) q[3];
sx q[3];
rz(-0.68651783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0278339) q[2];
sx q[2];
rz(-1.7366624) q[2];
sx q[2];
rz(3.0244381) q[2];
rz(-2.8095918) q[3];
sx q[3];
rz(-2.3947075) q[3];
sx q[3];
rz(-2.3565256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.4431045) q[0];
sx q[0];
rz(-2.3790058) q[0];
sx q[0];
rz(-0.37345988) q[0];
rz(2.7442878) q[1];
sx q[1];
rz(-1.1445069) q[1];
sx q[1];
rz(-2.4041046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339379) q[0];
sx q[0];
rz(-0.93385252) q[0];
sx q[0];
rz(-2.4600814) q[0];
rz(-pi) q[1];
rz(2.5869484) q[2];
sx q[2];
rz(-1.1078385) q[2];
sx q[2];
rz(0.92051586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7389679) q[1];
sx q[1];
rz(-2.1161656) q[1];
sx q[1];
rz(-1.7508372) q[1];
rz(-1.2190388) q[3];
sx q[3];
rz(-2.2969807) q[3];
sx q[3];
rz(-1.6399469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.071659) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(0.37754479) q[2];
rz(-2.8318882) q[3];
sx q[3];
rz(-1.0891424) q[3];
sx q[3];
rz(1.496605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51652235) q[0];
sx q[0];
rz(-0.83928883) q[0];
sx q[0];
rz(1.2155493) q[0];
rz(-1.0141605) q[1];
sx q[1];
rz(-1.7182257) q[1];
sx q[1];
rz(-0.97022143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033701115) q[0];
sx q[0];
rz(-1.6112856) q[0];
sx q[0];
rz(1.0711566) q[0];
rz(-pi) q[1];
rz(0.90486352) q[2];
sx q[2];
rz(-2.3653125) q[2];
sx q[2];
rz(1.640445) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.045779) q[1];
sx q[1];
rz(-1.9111553) q[1];
sx q[1];
rz(-2.9540121) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35200624) q[3];
sx q[3];
rz(-1.4624634) q[3];
sx q[3];
rz(-1.7217404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0479451) q[2];
sx q[2];
rz(-2.2990172) q[2];
sx q[2];
rz(-2.688431) q[2];
rz(1.2293182) q[3];
sx q[3];
rz(-1.8639576) q[3];
sx q[3];
rz(0.17190988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969084) q[0];
sx q[0];
rz(-0.6854282) q[0];
sx q[0];
rz(2.7759283) q[0];
rz(2.2293495) q[1];
sx q[1];
rz(-1.8625926) q[1];
sx q[1];
rz(-2.7361187) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.643754) q[0];
sx q[0];
rz(-2.9037243) q[0];
sx q[0];
rz(-2.0980623) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3581966) q[2];
sx q[2];
rz(-2.3153164) q[2];
sx q[2];
rz(0.43897334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7995389) q[1];
sx q[1];
rz(-0.56302445) q[1];
sx q[1];
rz(-2.8232226) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2527189) q[3];
sx q[3];
rz(-2.1858099) q[3];
sx q[3];
rz(-2.033704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.071986467) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(-1.5988505) q[2];
rz(0.37374464) q[3];
sx q[3];
rz(-1.475324) q[3];
sx q[3];
rz(-1.4603978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.54295802) q[0];
sx q[0];
rz(-0.77474189) q[0];
sx q[0];
rz(-2.4454818) q[0];
rz(1.8822582) q[1];
sx q[1];
rz(-1.101661) q[1];
sx q[1];
rz(-2.7764244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8697978) q[0];
sx q[0];
rz(-0.68469098) q[0];
sx q[0];
rz(0.71918804) q[0];
x q[1];
rz(-1.1579726) q[2];
sx q[2];
rz(-1.8637805) q[2];
sx q[2];
rz(2.3940046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7988749) q[1];
sx q[1];
rz(-1.5694071) q[1];
sx q[1];
rz(-0.54552127) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9837512) q[3];
sx q[3];
rz(-0.98414183) q[3];
sx q[3];
rz(1.1053305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68088561) q[2];
sx q[2];
rz(-1.978771) q[2];
sx q[2];
rz(0.17670259) q[2];
rz(-1.7548615) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(2.7669014) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6143167) q[0];
sx q[0];
rz(-2.2820331) q[0];
sx q[0];
rz(-1.5981307) q[0];
rz(1.3044926) q[1];
sx q[1];
rz(-1.0544798) q[1];
sx q[1];
rz(-0.80610448) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5209608) q[0];
sx q[0];
rz(-1.8072819) q[0];
sx q[0];
rz(-2.3804401) q[0];
x q[1];
rz(0.97053846) q[2];
sx q[2];
rz(-2.0036864) q[2];
sx q[2];
rz(-2.8899756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7215868) q[1];
sx q[1];
rz(-1.3665821) q[1];
sx q[1];
rz(0.050724357) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5664653) q[3];
sx q[3];
rz(-2.025423) q[3];
sx q[3];
rz(-1.2808587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1690037) q[2];
sx q[2];
rz(-0.4946332) q[2];
sx q[2];
rz(3.0779823) q[2];
rz(0.58498597) q[3];
sx q[3];
rz(-1.7413185) q[3];
sx q[3];
rz(1.6271648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.1503898) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(0.72917953) q[0];
rz(1.179262) q[1];
sx q[1];
rz(-0.74960342) q[1];
sx q[1];
rz(2.8470305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9091932) q[0];
sx q[0];
rz(-2.5721484) q[0];
sx q[0];
rz(0.48869407) q[0];
rz(-3.0315517) q[2];
sx q[2];
rz(-1.1400643) q[2];
sx q[2];
rz(0.34039341) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8940436) q[1];
sx q[1];
rz(-1.4534543) q[1];
sx q[1];
rz(-2.3168922) q[1];
rz(2.0414646) q[3];
sx q[3];
rz(-0.72689712) q[3];
sx q[3];
rz(1.6290381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6509167) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(-0.73224625) q[2];
rz(-2.2722774) q[3];
sx q[3];
rz(-1.1704051) q[3];
sx q[3];
rz(-1.7349582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.231584) q[0];
sx q[0];
rz(-1.5861479) q[0];
sx q[0];
rz(2.3439132) q[0];
rz(-0.32294598) q[1];
sx q[1];
rz(-2.1452417) q[1];
sx q[1];
rz(-1.4083883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7777611) q[0];
sx q[0];
rz(-2.1182502) q[0];
sx q[0];
rz(-1.3307894) q[0];
rz(-2.0297208) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(-1.9314507) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3288142) q[1];
sx q[1];
rz(-1.8056889) q[1];
sx q[1];
rz(2.3702904) q[1];
rz(-pi) q[2];
rz(-0.50813913) q[3];
sx q[3];
rz(-2.7781825) q[3];
sx q[3];
rz(-1.5606708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8352167) q[2];
sx q[2];
rz(-2.0312467) q[2];
sx q[2];
rz(2.9442673) q[2];
rz(2.3399682) q[3];
sx q[3];
rz(-0.97951952) q[3];
sx q[3];
rz(-2.7200123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9824958) q[0];
sx q[0];
rz(-1.4747341) q[0];
sx q[0];
rz(2.2267447) q[0];
rz(0.21028701) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(2.4519144) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4431157) q[0];
sx q[0];
rz(-1.7962828) q[0];
sx q[0];
rz(1.4754292) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6946227) q[2];
sx q[2];
rz(-3.0231907) q[2];
sx q[2];
rz(-1.3294544) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27305973) q[1];
sx q[1];
rz(-2.0439721) q[1];
sx q[1];
rz(2.563963) q[1];
rz(-pi) q[2];
rz(2.0965936) q[3];
sx q[3];
rz(-1.8711578) q[3];
sx q[3];
rz(-1.5142358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0476394) q[2];
sx q[2];
rz(-2.2932105) q[2];
sx q[2];
rz(0.32996714) q[2];
rz(1.0432358) q[3];
sx q[3];
rz(-1.4179351) q[3];
sx q[3];
rz(0.51658336) q[3];
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
rz(2.1107776) q[0];
sx q[0];
rz(-1.8106221) q[0];
sx q[0];
rz(2.0844841) q[0];
rz(2.8874176) q[1];
sx q[1];
rz(-1.1421685) q[1];
sx q[1];
rz(-2.4370297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0443665) q[0];
sx q[0];
rz(-1.2597879) q[0];
sx q[0];
rz(-1.0754271) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3866587) q[2];
sx q[2];
rz(-1.5103075) q[2];
sx q[2];
rz(-2.9764701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3290289) q[1];
sx q[1];
rz(-2.6575251) q[1];
sx q[1];
rz(1.1776393) q[1];
rz(-2.4289565) q[3];
sx q[3];
rz(-0.27240005) q[3];
sx q[3];
rz(2.2505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5083984) q[2];
sx q[2];
rz(-1.1722379) q[2];
sx q[2];
rz(2.634826) q[2];
rz(0.1861598) q[3];
sx q[3];
rz(-2.9045744) q[3];
sx q[3];
rz(-2.6059634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3770461) q[0];
sx q[0];
rz(-1.4519539) q[0];
sx q[0];
rz(-2.0013381) q[0];
rz(2.2346732) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(-0.89217734) q[2];
sx q[2];
rz(-1.1363251) q[2];
sx q[2];
rz(-2.688495) q[2];
rz(2.818454) q[3];
sx q[3];
rz(-1.825014) q[3];
sx q[3];
rz(-2.0925044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
