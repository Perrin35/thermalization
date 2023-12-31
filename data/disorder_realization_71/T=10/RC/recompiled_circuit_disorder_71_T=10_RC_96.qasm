OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60183817) q[0];
sx q[0];
rz(-1.1798501) q[0];
sx q[0];
rz(1.8200309) q[0];
rz(-0.9853739) q[2];
sx q[2];
rz(-1.7314163) q[2];
sx q[2];
rz(-2.0656245) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4235916) q[1];
sx q[1];
rz(-0.33034409) q[1];
sx q[1];
rz(1.6172536) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7238486) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(-0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-2.0200502) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(0.87444011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16337285) q[0];
sx q[0];
rz(-1.997943) q[0];
sx q[0];
rz(-0.33421974) q[0];
rz(2.1463257) q[2];
sx q[2];
rz(-1.7003635) q[2];
sx q[2];
rz(2.1339983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0622105) q[1];
sx q[1];
rz(-2.1746238) q[1];
sx q[1];
rz(-2.126333) q[1];
rz(-pi) q[2];
rz(-3.087895) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(-2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4008537) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-0.43593105) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(-2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(2.7456465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229189) q[0];
sx q[0];
rz(-2.4791105) q[0];
sx q[0];
rz(0.89612095) q[0];
rz(-pi) q[1];
rz(1.5244353) q[2];
sx q[2];
rz(-0.29435396) q[2];
sx q[2];
rz(-2.4183395) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6763941) q[1];
sx q[1];
rz(-2.4583543) q[1];
sx q[1];
rz(3.1264683) q[1];
rz(0.08926908) q[3];
sx q[3];
rz(-0.25651989) q[3];
sx q[3];
rz(-1.741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6039156) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72162119) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(0.23342361) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.061302) q[0];
sx q[0];
rz(-1.4740605) q[0];
sx q[0];
rz(-0.060782766) q[0];
rz(1.8065622) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(2.6218888) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.213821) q[1];
sx q[1];
rz(-1.2328316) q[1];
sx q[1];
rz(-0.73422276) q[1];
x q[2];
rz(0.38846429) q[3];
sx q[3];
rz(-2.3187175) q[3];
sx q[3];
rz(-2.8188761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0507811) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(-2.3941669) q[2];
rz(0.22339544) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-3.0773556) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-2.3805526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8781389) q[0];
sx q[0];
rz(-1.8315151) q[0];
sx q[0];
rz(-3.1104452) q[0];
rz(-pi) q[1];
rz(0.20247395) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(1.7970049) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0423454) q[1];
sx q[1];
rz(-1.5464916) q[1];
sx q[1];
rz(-2.821032) q[1];
x q[2];
rz(2.8588572) q[3];
sx q[3];
rz(-1.8394107) q[3];
sx q[3];
rz(-2.3408567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80660194) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(-2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(0.32593265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3098785) q[0];
sx q[0];
rz(-1.9448115) q[0];
sx q[0];
rz(2.5783587) q[0];
rz(0.58165254) q[2];
sx q[2];
rz(-2.622421) q[2];
sx q[2];
rz(0.66228629) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1833916) q[1];
sx q[1];
rz(-0.32918731) q[1];
sx q[1];
rz(0.42894657) q[1];
rz(-2.9767354) q[3];
sx q[3];
rz(-3.0391209) q[3];
sx q[3];
rz(-0.41302478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(1.4298965) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(2.4196999) q[0];
rz(1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(-0.11925764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3463466) q[0];
sx q[0];
rz(-2.3756785) q[0];
sx q[0];
rz(0.38608293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7330405) q[2];
sx q[2];
rz(-1.679323) q[2];
sx q[2];
rz(-2.721399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81469401) q[1];
sx q[1];
rz(-2.6242995) q[1];
sx q[1];
rz(-1.3002404) q[1];
rz(-pi) q[2];
rz(0.70763208) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(-1.6361145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(-1.0634364) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0793593) q[0];
sx q[0];
rz(-1.0842807) q[0];
sx q[0];
rz(-1.833672) q[0];
rz(3.1259414) q[2];
sx q[2];
rz(-2.1502697) q[2];
sx q[2];
rz(1.5755115) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13247989) q[1];
sx q[1];
rz(-1.1785893) q[1];
sx q[1];
rz(1.023804) q[1];
x q[2];
rz(-2.3732244) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(-1.8661631) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(2.0057604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4275883) q[0];
sx q[0];
rz(-1.7363318) q[0];
sx q[0];
rz(1.0323314) q[0];
rz(-pi) q[1];
rz(-1.5485974) q[2];
sx q[2];
rz(-1.1361085) q[2];
sx q[2];
rz(-1.4093902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.51582) q[1];
sx q[1];
rz(-1.1449709) q[1];
sx q[1];
rz(1.5172525) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1080092) q[3];
sx q[3];
rz(-2.4164003) q[3];
sx q[3];
rz(-0.68554002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(2.2686968) q[2];
rz(0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(1.8632442) q[0];
rz(-2.1168013) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(-1.9445673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.409317) q[0];
sx q[0];
rz(-2.0079552) q[0];
sx q[0];
rz(-0.93426312) q[0];
rz(-pi) q[1];
rz(1.808666) q[2];
sx q[2];
rz(-1.6878205) q[2];
sx q[2];
rz(2.7931917) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7750818) q[1];
sx q[1];
rz(-2.0644036) q[1];
sx q[1];
rz(2.2713186) q[1];
x q[2];
rz(0.50935575) q[3];
sx q[3];
rz(-1.2378113) q[3];
sx q[3];
rz(1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3060351) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(0.79997921) q[2];
rz(-1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-0.75688731) q[2];
sx q[2];
rz(-0.49513985) q[2];
sx q[2];
rz(-2.0078299) q[2];
rz(-2.342631) q[3];
sx q[3];
rz(-2.3286455) q[3];
sx q[3];
rz(0.13959985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
