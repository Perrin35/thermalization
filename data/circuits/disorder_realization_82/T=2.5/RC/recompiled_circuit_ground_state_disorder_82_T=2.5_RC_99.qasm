OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3759484) q[0];
sx q[0];
rz(-2.7274237) q[0];
sx q[0];
rz(0.8015269) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(3.0473696) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.191358) q[0];
sx q[0];
rz(-1.1947898) q[0];
sx q[0];
rz(-2.3467499) q[0];
rz(-pi) q[1];
rz(-1.6802189) q[2];
sx q[2];
rz(-0.27421194) q[2];
sx q[2];
rz(0.69137979) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-0.5597502) q[1];
sx q[1];
rz(1.00882) q[1];
rz(-pi) q[2];
rz(-0.014299794) q[3];
sx q[3];
rz(-2.9889782) q[3];
sx q[3];
rz(-0.26210704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(-1.630265) q[2];
rz(-0.18621914) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(-0.65096861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179825) q[0];
sx q[0];
rz(-1.9600493) q[0];
sx q[0];
rz(-0.13667662) q[0];
rz(0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(0.76900855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5154607) q[0];
sx q[0];
rz(-2.7410586) q[0];
sx q[0];
rz(2.787091) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42119512) q[2];
sx q[2];
rz(-0.73891089) q[2];
sx q[2];
rz(0.74904672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5040759) q[1];
sx q[1];
rz(-2.3055162) q[1];
sx q[1];
rz(-1.3512263) q[1];
x q[2];
rz(-1.3202928) q[3];
sx q[3];
rz(-0.13938306) q[3];
sx q[3];
rz(-1.936862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39917055) q[2];
sx q[2];
rz(-0.48226446) q[2];
sx q[2];
rz(-1.4313618) q[2];
rz(2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2630149) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(1.0947134) q[0];
rz(-1.8393983) q[1];
sx q[1];
rz(-2.1523988) q[1];
sx q[1];
rz(3.1375695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1510295) q[0];
sx q[0];
rz(-2.4044665) q[0];
sx q[0];
rz(1.766979) q[0];
rz(-pi) q[1];
rz(2.268774) q[2];
sx q[2];
rz(-0.96257639) q[2];
sx q[2];
rz(-1.5814511) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6371256) q[1];
sx q[1];
rz(-1.5404623) q[1];
sx q[1];
rz(-2.8843969) q[1];
rz(-pi) q[2];
x q[2];
rz(1.324292) q[3];
sx q[3];
rz(-0.87216264) q[3];
sx q[3];
rz(-2.7072631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5912938) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(-2.6934521) q[2];
rz(0.48629931) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(1.1264616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(-2.4718156) q[0];
rz(-1.9122596) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(0.545048) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16825039) q[0];
sx q[0];
rz(-1.1118044) q[0];
sx q[0];
rz(-2.5516872) q[0];
rz(0.82370211) q[2];
sx q[2];
rz(-2.3266226) q[2];
sx q[2];
rz(-2.4832249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2983568) q[1];
sx q[1];
rz(-1.7853754) q[1];
sx q[1];
rz(-3.0797348) q[1];
rz(2.0470951) q[3];
sx q[3];
rz(-2.3102488) q[3];
sx q[3];
rz(-0.0090713105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58686078) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(1.2083017) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5912882) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(-3.1377129) q[0];
rz(-2.8344391) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(-2.602813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85941797) q[0];
sx q[0];
rz(-0.27056405) q[0];
sx q[0];
rz(2.5859588) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7096552) q[2];
sx q[2];
rz(-1.9966239) q[2];
sx q[2];
rz(0.52639293) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9591892) q[1];
sx q[1];
rz(-2.1418874) q[1];
sx q[1];
rz(-1.9923927) q[1];
rz(-pi) q[2];
rz(0.49541028) q[3];
sx q[3];
rz(-0.79548478) q[3];
sx q[3];
rz(2.0919679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4469845) q[2];
sx q[2];
rz(-2.0422523) q[2];
sx q[2];
rz(-1.2896607) q[2];
rz(-1.6192216) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(1.6600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.30037844) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(0.73219055) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(2.5644459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4704551) q[0];
sx q[0];
rz(-1.5671434) q[0];
sx q[0];
rz(1.5421583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24546282) q[2];
sx q[2];
rz(-0.6010913) q[2];
sx q[2];
rz(2.9562151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53458339) q[1];
sx q[1];
rz(-1.4614236) q[1];
sx q[1];
rz(-1.6499728) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7558877) q[3];
sx q[3];
rz(-1.3284747) q[3];
sx q[3];
rz(1.7923812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.131549) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(1.5232167) q[2];
rz(-3.0736382) q[3];
sx q[3];
rz(-0.32502919) q[3];
sx q[3];
rz(2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286569) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(2.3701684) q[0];
rz(2.6029288) q[1];
sx q[1];
rz(-1.3464876) q[1];
sx q[1];
rz(-1.0661941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26771256) q[0];
sx q[0];
rz(-1.4940049) q[0];
sx q[0];
rz(-1.8241054) q[0];
x q[1];
rz(-2.0838924) q[2];
sx q[2];
rz(-1.3670245) q[2];
sx q[2];
rz(2.9625826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33246189) q[1];
sx q[1];
rz(-1.5576502) q[1];
sx q[1];
rz(-2.9244208) q[1];
rz(-pi) q[2];
rz(-1.3049502) q[3];
sx q[3];
rz(-1.8943333) q[3];
sx q[3];
rz(0.27184799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44925877) q[2];
sx q[2];
rz(-0.80654311) q[2];
sx q[2];
rz(0.44245693) q[2];
rz(-2.9509406) q[3];
sx q[3];
rz(-2.4176044) q[3];
sx q[3];
rz(-2.382522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8964748) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.8316487) q[0];
rz(2.7438121) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(-1.3936874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01263857) q[0];
sx q[0];
rz(-2.0067747) q[0];
sx q[0];
rz(2.1672744) q[0];
rz(-pi) q[1];
rz(-0.23598052) q[2];
sx q[2];
rz(-1.8644635) q[2];
sx q[2];
rz(0.27494173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85227784) q[1];
sx q[1];
rz(-2.5767269) q[1];
sx q[1];
rz(0.58128618) q[1];
rz(-pi) q[2];
x q[2];
rz(2.05768) q[3];
sx q[3];
rz(-0.22398914) q[3];
sx q[3];
rz(-1.4518224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1132249) q[2];
sx q[2];
rz(-2.4312879) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(2.4469589) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-2.1139483) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46591127) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(0.95941108) q[0];
rz(2.361946) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(-0.27228212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.773943) q[0];
sx q[0];
rz(-0.01423562) q[0];
sx q[0];
rz(1.9444501) q[0];
rz(1.2113038) q[2];
sx q[2];
rz(-2.4036416) q[2];
sx q[2];
rz(2.8518845) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0148264) q[1];
sx q[1];
rz(-2.0527168) q[1];
sx q[1];
rz(0.46238203) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5078154) q[3];
sx q[3];
rz(-2.2867775) q[3];
sx q[3];
rz(-0.26308003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73725629) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(-0.20757248) q[2];
rz(-3.0370144) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(-1.6149717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.680147) q[0];
sx q[0];
rz(-1.9141645) q[0];
sx q[0];
rz(2.0848059) q[0];
rz(-0.54569221) q[1];
sx q[1];
rz(-2.1634407) q[1];
sx q[1];
rz(-2.812885) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0775111) q[0];
sx q[0];
rz(-1.2405378) q[0];
sx q[0];
rz(-2.0547377) q[0];
rz(-0.07241197) q[2];
sx q[2];
rz(-1.5161523) q[2];
sx q[2];
rz(-1.8397728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2519319) q[1];
sx q[1];
rz(-1.228782) q[1];
sx q[1];
rz(3.1238848) q[1];
rz(-2.3700033) q[3];
sx q[3];
rz(-1.3719012) q[3];
sx q[3];
rz(-2.1344913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1591961) q[2];
sx q[2];
rz(-3.0089162) q[2];
sx q[2];
rz(1.1245493) q[2];
rz(3.0647965) q[3];
sx q[3];
rz(-2.7087637) q[3];
sx q[3];
rz(-2.2176149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9818253) q[0];
sx q[0];
rz(-1.8851017) q[0];
sx q[0];
rz(-1.5633352) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-2.0682206) q[2];
sx q[2];
rz(-2.2523027) q[2];
sx q[2];
rz(-0.61423617) q[2];
rz(-0.40533752) q[3];
sx q[3];
rz(-1.574285) q[3];
sx q[3];
rz(-2.3226572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
