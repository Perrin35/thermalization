OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(1.3258452) q[0];
rz(1.3658547) q[1];
sx q[1];
rz(-0.39443016) q[1];
sx q[1];
rz(-2.3529513) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45730725) q[0];
sx q[0];
rz(-0.72170242) q[0];
sx q[0];
rz(2.0474322) q[0];
rz(0.85429384) q[2];
sx q[2];
rz(-2.2546356) q[2];
sx q[2];
rz(1.8707866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0683191) q[1];
sx q[1];
rz(-1.2879268) q[1];
sx q[1];
rz(1.7712996) q[1];
rz(-pi) q[2];
rz(1.3943421) q[3];
sx q[3];
rz(-0.6354699) q[3];
sx q[3];
rz(-0.42144767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49393645) q[2];
sx q[2];
rz(-1.4467354) q[2];
sx q[2];
rz(0.20230618) q[2];
rz(0.10937396) q[3];
sx q[3];
rz(-0.60905639) q[3];
sx q[3];
rz(-0.085453184) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0789455) q[0];
sx q[0];
rz(-2.2181692) q[0];
sx q[0];
rz(1.9664636) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-1.0025832) q[1];
sx q[1];
rz(-0.35951231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8749411) q[0];
sx q[0];
rz(-1.734899) q[0];
sx q[0];
rz(3.0709188) q[0];
rz(-pi) q[1];
rz(-2.5841332) q[2];
sx q[2];
rz(-1.895708) q[2];
sx q[2];
rz(0.2211472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2613758) q[1];
sx q[1];
rz(-1.2296074) q[1];
sx q[1];
rz(0.64067322) q[1];
rz(-pi) q[2];
rz(0.97378181) q[3];
sx q[3];
rz(-2.0691853) q[3];
sx q[3];
rz(0.9998876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29490617) q[2];
sx q[2];
rz(-2.0945695) q[2];
sx q[2];
rz(0.27493757) q[2];
rz(-1.8168195) q[3];
sx q[3];
rz(-1.6144269) q[3];
sx q[3];
rz(-1.8435318) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777968) q[0];
sx q[0];
rz(-2.9834788) q[0];
sx q[0];
rz(-0.39704821) q[0];
rz(0.60931698) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(1.5416001) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5865094) q[0];
sx q[0];
rz(-1.4662881) q[0];
sx q[0];
rz(1.2782607) q[0];
rz(-pi) q[1];
rz(-1.6360388) q[2];
sx q[2];
rz(-1.5833304) q[2];
sx q[2];
rz(1.6103075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1101602) q[1];
sx q[1];
rz(-0.84430391) q[1];
sx q[1];
rz(-3.1189819) q[1];
rz(-pi) q[2];
rz(-2.2619369) q[3];
sx q[3];
rz(-2.1113951) q[3];
sx q[3];
rz(-0.72002711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93727532) q[2];
sx q[2];
rz(-1.5466362) q[2];
sx q[2];
rz(0.23359648) q[2];
rz(-1.2967845) q[3];
sx q[3];
rz(-1.1917043) q[3];
sx q[3];
rz(-0.72062033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5066756) q[0];
sx q[0];
rz(-1.9838061) q[0];
sx q[0];
rz(-1.3651715) q[0];
rz(-1.2719951) q[1];
sx q[1];
rz(-1.1993661) q[1];
sx q[1];
rz(1.8796399) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.711447) q[0];
sx q[0];
rz(-1.0017348) q[0];
sx q[0];
rz(1.3422768) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8795525) q[2];
sx q[2];
rz(-1.5456887) q[2];
sx q[2];
rz(-1.0756191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8820928) q[1];
sx q[1];
rz(-1.5061404) q[1];
sx q[1];
rz(-2.9586709) q[1];
rz(-pi) q[2];
rz(2.3172195) q[3];
sx q[3];
rz(-1.715987) q[3];
sx q[3];
rz(2.8100325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23123732) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(-2.3528698) q[2];
rz(1.8606868) q[3];
sx q[3];
rz(-1.3642045) q[3];
sx q[3];
rz(-1.0219215) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3874409) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(-0.66106853) q[0];
rz(-2.0263653) q[1];
sx q[1];
rz(-1.3122357) q[1];
sx q[1];
rz(2.1818395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5710053) q[0];
sx q[0];
rz(-2.5474605) q[0];
sx q[0];
rz(-2.6460365) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9767545) q[2];
sx q[2];
rz(-1.0292813) q[2];
sx q[2];
rz(-2.7680226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57454039) q[1];
sx q[1];
rz(-0.31992074) q[1];
sx q[1];
rz(-2.0053276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3101391) q[3];
sx q[3];
rz(-2.4803526) q[3];
sx q[3];
rz(2.6186297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5614718) q[2];
sx q[2];
rz(-2.6802907) q[2];
sx q[2];
rz(2.0816154) q[2];
rz(-2.3914242) q[3];
sx q[3];
rz(-1.7629905) q[3];
sx q[3];
rz(1.2257956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8691413) q[0];
sx q[0];
rz(-1.5831818) q[0];
sx q[0];
rz(-2.9172752) q[0];
rz(-1.51651) q[1];
sx q[1];
rz(-0.61619157) q[1];
sx q[1];
rz(3.065965) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13091732) q[0];
sx q[0];
rz(-1.0300127) q[0];
sx q[0];
rz(0.45387876) q[0];
rz(1.8706029) q[2];
sx q[2];
rz(-2.8245016) q[2];
sx q[2];
rz(-0.89978774) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5995248) q[1];
sx q[1];
rz(-2.0318527) q[1];
sx q[1];
rz(0.49054029) q[1];
x q[2];
rz(2.7073955) q[3];
sx q[3];
rz(-2.2480474) q[3];
sx q[3];
rz(-1.3788669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6077967) q[2];
sx q[2];
rz(-2.3888612) q[2];
sx q[2];
rz(3.094574) q[2];
rz(-0.67695824) q[3];
sx q[3];
rz(-1.8432063) q[3];
sx q[3];
rz(-0.025378749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.15636477) q[0];
sx q[0];
rz(-1.4893091) q[0];
sx q[0];
rz(1.697668) q[0];
rz(-2.2185183) q[1];
sx q[1];
rz(-2.1363027) q[1];
sx q[1];
rz(-0.18327555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3175439) q[0];
sx q[0];
rz(-2.5185761) q[0];
sx q[0];
rz(-2.7035294) q[0];
x q[1];
rz(-0.44281339) q[2];
sx q[2];
rz(-1.9573136) q[2];
sx q[2];
rz(2.7903976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1122656) q[1];
sx q[1];
rz(-1.8958175) q[1];
sx q[1];
rz(-1.0461665) q[1];
rz(-1.9952946) q[3];
sx q[3];
rz(-2.5190926) q[3];
sx q[3];
rz(1.5869035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.034059374) q[2];
sx q[2];
rz(-1.4564161) q[2];
sx q[2];
rz(-3.122186) q[2];
rz(2.9803045) q[3];
sx q[3];
rz(-2.1262157) q[3];
sx q[3];
rz(-1.9136782) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1041383) q[0];
sx q[0];
rz(-2.7327974) q[0];
sx q[0];
rz(-3.0492875) q[0];
rz(-1.9654988) q[1];
sx q[1];
rz(-0.25830019) q[1];
sx q[1];
rz(1.4091122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8831439) q[0];
sx q[0];
rz(-1.5165268) q[0];
sx q[0];
rz(-1.4656214) q[0];
x q[1];
rz(0.46940501) q[2];
sx q[2];
rz(-2.1300211) q[2];
sx q[2];
rz(-2.6162868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17615023) q[1];
sx q[1];
rz(-1.1854082) q[1];
sx q[1];
rz(-1.8606645) q[1];
rz(-pi) q[2];
rz(0.61586611) q[3];
sx q[3];
rz(-1.0857673) q[3];
sx q[3];
rz(0.46014338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1000503) q[2];
sx q[2];
rz(-1.5644194) q[2];
sx q[2];
rz(1.9680295) q[2];
rz(0.38343492) q[3];
sx q[3];
rz(-1.8337367) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0122871) q[0];
sx q[0];
rz(-1.4947991) q[0];
sx q[0];
rz(0.20558414) q[0];
rz(2.3221817) q[1];
sx q[1];
rz(-1.2148379) q[1];
sx q[1];
rz(-2.6507071) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48614472) q[0];
sx q[0];
rz(-1.1974338) q[0];
sx q[0];
rz(0.70151897) q[0];
rz(-pi) q[1];
rz(1.2621636) q[2];
sx q[2];
rz(-0.96519816) q[2];
sx q[2];
rz(3.1250033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5832688) q[1];
sx q[1];
rz(-2.138639) q[1];
sx q[1];
rz(-0.056599157) q[1];
rz(-pi) q[2];
rz(-0.098693178) q[3];
sx q[3];
rz(-2.0798363) q[3];
sx q[3];
rz(2.4157897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72863355) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(-1.884985) q[2];
rz(2.7058153) q[3];
sx q[3];
rz(-1.1083009) q[3];
sx q[3];
rz(2.2573439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7213781) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(2.9086928) q[0];
rz(0.13889343) q[1];
sx q[1];
rz(-1.9878309) q[1];
sx q[1];
rz(-0.5084261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8247037) q[0];
sx q[0];
rz(-2.1068067) q[0];
sx q[0];
rz(-2.8404923) q[0];
rz(1.0632238) q[2];
sx q[2];
rz(-1.3295577) q[2];
sx q[2];
rz(-0.2080179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1436746) q[1];
sx q[1];
rz(-1.7047722) q[1];
sx q[1];
rz(-3.133598) q[1];
rz(-pi) q[2];
rz(0.53896972) q[3];
sx q[3];
rz(-2.1930755) q[3];
sx q[3];
rz(1.301018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4109219) q[2];
sx q[2];
rz(-1.4620917) q[2];
sx q[2];
rz(-2.3250568) q[2];
rz(2.7517892) q[3];
sx q[3];
rz(-1.0023508) q[3];
sx q[3];
rz(1.7635112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6941866) q[0];
sx q[0];
rz(-2.2185855) q[0];
sx q[0];
rz(2.9119281) q[0];
rz(1.6387088) q[1];
sx q[1];
rz(-1.3353744) q[1];
sx q[1];
rz(-2.2257805) q[1];
rz(-1.3982795) q[2];
sx q[2];
rz(-1.1926706) q[2];
sx q[2];
rz(-0.43422912) q[2];
rz(-1.3744204) q[3];
sx q[3];
rz(-0.71153258) q[3];
sx q[3];
rz(1.1179433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
