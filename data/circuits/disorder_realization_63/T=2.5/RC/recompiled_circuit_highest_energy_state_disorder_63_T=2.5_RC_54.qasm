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
rz(1.0032049) q[0];
sx q[0];
rz(1.9411074) q[0];
sx q[0];
rz(8.3984126) q[0];
rz(0.29844555) q[1];
sx q[1];
rz(-1.6331853) q[1];
sx q[1];
rz(1.0025947) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74513197) q[0];
sx q[0];
rz(-2.9365109) q[0];
sx q[0];
rz(1.6830326) q[0];
rz(-pi) q[1];
rz(-3.0709362) q[2];
sx q[2];
rz(-2.4774654) q[2];
sx q[2];
rz(0.16985591) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2015675) q[1];
sx q[1];
rz(-2.7301412) q[1];
sx q[1];
rz(0.33951851) q[1];
rz(2.1777654) q[3];
sx q[3];
rz(-2.1133528) q[3];
sx q[3];
rz(2.3574717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5528494) q[2];
sx q[2];
rz(-2.1407949) q[2];
sx q[2];
rz(0.0351077) q[2];
rz(1.2037753) q[3];
sx q[3];
rz(-1.8193918) q[3];
sx q[3];
rz(-0.28425851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069000706) q[0];
sx q[0];
rz(-2.6995316) q[0];
sx q[0];
rz(1.0042071) q[0];
rz(-3.0191811) q[1];
sx q[1];
rz(-1.4675843) q[1];
sx q[1];
rz(2.4845128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5733684) q[0];
sx q[0];
rz(-1.4889297) q[0];
sx q[0];
rz(-0.38827814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2934712) q[2];
sx q[2];
rz(-1.4155626) q[2];
sx q[2];
rz(-2.5350132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7065694) q[1];
sx q[1];
rz(-2.2945255) q[1];
sx q[1];
rz(0.91591452) q[1];
x q[2];
rz(-2.8960896) q[3];
sx q[3];
rz(-2.3424753) q[3];
sx q[3];
rz(-1.4202858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.55679) q[2];
sx q[2];
rz(-2.5636797) q[2];
sx q[2];
rz(0.59188262) q[2];
rz(1.4764192) q[3];
sx q[3];
rz(-1.6074601) q[3];
sx q[3];
rz(2.2679451) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0037435) q[0];
sx q[0];
rz(-2.9326404) q[0];
sx q[0];
rz(-1.9269706) q[0];
rz(-0.95642033) q[1];
sx q[1];
rz(-1.950187) q[1];
sx q[1];
rz(-1.1716243) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7192243) q[0];
sx q[0];
rz(-2.5107493) q[0];
sx q[0];
rz(-0.8623841) q[0];
rz(-pi) q[1];
rz(-1.8128305) q[2];
sx q[2];
rz(-0.71513745) q[2];
sx q[2];
rz(2.8492343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4332648) q[1];
sx q[1];
rz(-1.2339795) q[1];
sx q[1];
rz(-0.59971209) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0517393) q[3];
sx q[3];
rz(-1.5260586) q[3];
sx q[3];
rz(-2.5344332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9413635) q[2];
sx q[2];
rz(-2.0912632) q[2];
sx q[2];
rz(1.7639147) q[2];
rz(2.7348943) q[3];
sx q[3];
rz(-0.64928693) q[3];
sx q[3];
rz(2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.7895301) q[0];
sx q[0];
rz(-1.8822414) q[0];
sx q[0];
rz(1.2404741) q[0];
rz(0.71818304) q[1];
sx q[1];
rz(-1.8156464) q[1];
sx q[1];
rz(-0.2389508) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093042308) q[0];
sx q[0];
rz(-1.4497546) q[0];
sx q[0];
rz(-2.8516475) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3117606) q[2];
sx q[2];
rz(-0.90832635) q[2];
sx q[2];
rz(1.5892513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5217704) q[1];
sx q[1];
rz(-1.8293293) q[1];
sx q[1];
rz(0.14743989) q[1];
rz(1.4364868) q[3];
sx q[3];
rz(-2.1013936) q[3];
sx q[3];
rz(-2.1710896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6255528) q[2];
sx q[2];
rz(-1.7480787) q[2];
sx q[2];
rz(-1.7986521) q[2];
rz(2.2602153) q[3];
sx q[3];
rz(-0.16888976) q[3];
sx q[3];
rz(2.3605997) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7278904) q[0];
sx q[0];
rz(-2.8975633) q[0];
sx q[0];
rz(1.4843041) q[0];
rz(-0.65385747) q[1];
sx q[1];
rz(-1.5212719) q[1];
sx q[1];
rz(-3.003655) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333899) q[0];
sx q[0];
rz(-1.7195367) q[0];
sx q[0];
rz(1.9277907) q[0];
rz(-0.65278585) q[2];
sx q[2];
rz(-1.7695656) q[2];
sx q[2];
rz(-0.38236375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9508507) q[1];
sx q[1];
rz(-0.44176451) q[1];
sx q[1];
rz(-1.8015693) q[1];
x q[2];
rz(1.3119427) q[3];
sx q[3];
rz(-1.5167243) q[3];
sx q[3];
rz(-1.8340221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5518034) q[2];
sx q[2];
rz(-1.1729596) q[2];
sx q[2];
rz(-3.0544082) q[2];
rz(-3.0269571) q[3];
sx q[3];
rz(-0.81511027) q[3];
sx q[3];
rz(2.7242928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36599416) q[0];
sx q[0];
rz(-1.6021148) q[0];
sx q[0];
rz(-2.6566246) q[0];
rz(-2.1980749) q[1];
sx q[1];
rz(-0.8129932) q[1];
sx q[1];
rz(1.9224723) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2526357) q[0];
sx q[0];
rz(-0.85272861) q[0];
sx q[0];
rz(3.0120993) q[0];
rz(3.1252507) q[2];
sx q[2];
rz(-1.8927843) q[2];
sx q[2];
rz(2.0193375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.74138) q[1];
sx q[1];
rz(-1.4948436) q[1];
sx q[1];
rz(1.6981359) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7597105) q[3];
sx q[3];
rz(-1.1897161) q[3];
sx q[3];
rz(0.71807968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9014827) q[2];
sx q[2];
rz(-2.0811452) q[2];
sx q[2];
rz(2.2002952) q[2];
rz(-0.24602041) q[3];
sx q[3];
rz(-1.6451719) q[3];
sx q[3];
rz(1.7814319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0435903) q[0];
sx q[0];
rz(-2.0397546) q[0];
sx q[0];
rz(-1.164042) q[0];
rz(1.7835435) q[1];
sx q[1];
rz(-0.86654228) q[1];
sx q[1];
rz(-1.39303) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355125) q[0];
sx q[0];
rz(-2.7580259) q[0];
sx q[0];
rz(-1.3719029) q[0];
rz(0.54945182) q[2];
sx q[2];
rz(-1.3260815) q[2];
sx q[2];
rz(-3.0915054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5389312) q[1];
sx q[1];
rz(-1.6762513) q[1];
sx q[1];
rz(2.6698673) q[1];
rz(-pi) q[2];
rz(-0.89207666) q[3];
sx q[3];
rz(-1.3952655) q[3];
sx q[3];
rz(-0.20993671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67036575) q[2];
sx q[2];
rz(-1.7134075) q[2];
sx q[2];
rz(-0.29895374) q[2];
rz(0.40545884) q[3];
sx q[3];
rz(-2.7751228) q[3];
sx q[3];
rz(2.5772742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4295171) q[0];
sx q[0];
rz(-2.8686664) q[0];
sx q[0];
rz(2.3522229) q[0];
rz(-2.7366267) q[1];
sx q[1];
rz(-1.3507495) q[1];
sx q[1];
rz(2.6466323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69823658) q[0];
sx q[0];
rz(-0.77455097) q[0];
sx q[0];
rz(0.86408918) q[0];
rz(-2.1407897) q[2];
sx q[2];
rz(-1.8256639) q[2];
sx q[2];
rz(-2.3895181) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9443173) q[1];
sx q[1];
rz(-2.3305571) q[1];
sx q[1];
rz(-1.2326101) q[1];
rz(0.52431732) q[3];
sx q[3];
rz(-0.90745196) q[3];
sx q[3];
rz(1.2330221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1627545) q[2];
sx q[2];
rz(-1.9163722) q[2];
sx q[2];
rz(0.2552574) q[2];
rz(-0.034959547) q[3];
sx q[3];
rz(-1.6182263) q[3];
sx q[3];
rz(2.8953654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73664767) q[0];
sx q[0];
rz(-1.0628137) q[0];
sx q[0];
rz(-2.431562) q[0];
rz(1.0011477) q[1];
sx q[1];
rz(-2.2444057) q[1];
sx q[1];
rz(2.5206916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5516736) q[0];
sx q[0];
rz(-2.7749333) q[0];
sx q[0];
rz(-1.7398341) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65695564) q[2];
sx q[2];
rz(-1.1816506) q[2];
sx q[2];
rz(-0.78788131) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58644303) q[1];
sx q[1];
rz(-1.4634919) q[1];
sx q[1];
rz(0.34894033) q[1];
rz(0.7767949) q[3];
sx q[3];
rz(-1.3703386) q[3];
sx q[3];
rz(1.7944698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8093449) q[2];
sx q[2];
rz(-1.6184018) q[2];
sx q[2];
rz(0.28918949) q[2];
rz(-1.1122164) q[3];
sx q[3];
rz(-0.29505348) q[3];
sx q[3];
rz(1.0203863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4777098) q[0];
sx q[0];
rz(-1.3417256) q[0];
sx q[0];
rz(2.2148602) q[0];
rz(-2.0345188) q[1];
sx q[1];
rz(-1.5137545) q[1];
sx q[1];
rz(0.56799299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04161373) q[0];
sx q[0];
rz(-1.8868514) q[0];
sx q[0];
rz(0.63412447) q[0];
x q[1];
rz(-0.44804087) q[2];
sx q[2];
rz(-0.5574286) q[2];
sx q[2];
rz(0.008226062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2462594) q[1];
sx q[1];
rz(-2.3998108) q[1];
sx q[1];
rz(-1.0146902) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0509012) q[3];
sx q[3];
rz(-2.1445159) q[3];
sx q[3];
rz(0.75999505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0893112) q[2];
sx q[2];
rz(-2.3312882) q[2];
sx q[2];
rz(-1.7029765) q[2];
rz(0.29664052) q[3];
sx q[3];
rz(-2.7999925) q[3];
sx q[3];
rz(0.44704416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.457837) q[0];
sx q[0];
rz(-1.0975657) q[0];
sx q[0];
rz(0.56957635) q[0];
rz(-0.33869047) q[1];
sx q[1];
rz(-1.6117922) q[1];
sx q[1];
rz(1.502996) q[1];
rz(1.7768301) q[2];
sx q[2];
rz(-1.3246957) q[2];
sx q[2];
rz(-0.93304721) q[2];
rz(2.972907) q[3];
sx q[3];
rz(-1.8416234) q[3];
sx q[3];
rz(2.4467322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
