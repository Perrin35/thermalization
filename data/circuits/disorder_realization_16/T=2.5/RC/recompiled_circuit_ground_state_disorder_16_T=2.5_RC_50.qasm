OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.10915169) q[0];
sx q[0];
rz(4.2946058) q[0];
sx q[0];
rz(9.499318) q[0];
rz(1.6115161) q[1];
sx q[1];
rz(-3.0822152) q[1];
sx q[1];
rz(-0.52409726) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68428129) q[0];
sx q[0];
rz(-2.1403098) q[0];
sx q[0];
rz(-2.7049541) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33676001) q[2];
sx q[2];
rz(-1.3408061) q[2];
sx q[2];
rz(2.5652792) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84108564) q[1];
sx q[1];
rz(-2.1504554) q[1];
sx q[1];
rz(-2.4546844) q[1];
x q[2];
rz(1.662781) q[3];
sx q[3];
rz(-1.5301203) q[3];
sx q[3];
rz(2.3856861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3268299) q[2];
sx q[2];
rz(-2.6430898) q[2];
sx q[2];
rz(0.89547431) q[2];
rz(-0.26252663) q[3];
sx q[3];
rz(-0.34590507) q[3];
sx q[3];
rz(-0.088570647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.217591) q[0];
sx q[0];
rz(-1.9460678) q[0];
sx q[0];
rz(3.0254645) q[0];
rz(-0.63236347) q[1];
sx q[1];
rz(-2.875681) q[1];
sx q[1];
rz(-1.4858474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0362186) q[0];
sx q[0];
rz(-1.5040795) q[0];
sx q[0];
rz(-1.5936038) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9048278) q[2];
sx q[2];
rz(-2.1000266) q[2];
sx q[2];
rz(2.6108612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33148674) q[1];
sx q[1];
rz(-1.4432194) q[1];
sx q[1];
rz(1.4639616) q[1];
rz(3.0108475) q[3];
sx q[3];
rz(-2.2669889) q[3];
sx q[3];
rz(-1.2519022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31132013) q[2];
sx q[2];
rz(-0.98036426) q[2];
sx q[2];
rz(-0.92973989) q[2];
rz(-2.7021507) q[3];
sx q[3];
rz(-1.6012871) q[3];
sx q[3];
rz(-0.75696993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065652549) q[0];
sx q[0];
rz(-2.3891698) q[0];
sx q[0];
rz(0.886985) q[0];
rz(1.2607964) q[1];
sx q[1];
rz(-1.1075243) q[1];
sx q[1];
rz(-1.6541803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9943574) q[0];
sx q[0];
rz(-1.571313) q[0];
sx q[0];
rz(-1.6883255) q[0];
x q[1];
rz(-1.4753129) q[2];
sx q[2];
rz(-1.9621358) q[2];
sx q[2];
rz(3.0961159) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2339685) q[1];
sx q[1];
rz(-2.8115936) q[1];
sx q[1];
rz(3.0872295) q[1];
x q[2];
rz(1.2441176) q[3];
sx q[3];
rz(-1.0138113) q[3];
sx q[3];
rz(-2.9426334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6416574) q[2];
sx q[2];
rz(-0.96330088) q[2];
sx q[2];
rz(-0.43528834) q[2];
rz(-0.38517243) q[3];
sx q[3];
rz(-1.2110854) q[3];
sx q[3];
rz(2.1729573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4061072) q[0];
sx q[0];
rz(-0.084097363) q[0];
sx q[0];
rz(0.054585833) q[0];
rz(-1.6361884) q[1];
sx q[1];
rz(-1.9278229) q[1];
sx q[1];
rz(-0.41753599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62359257) q[0];
sx q[0];
rz(-1.5125649) q[0];
sx q[0];
rz(2.7043155) q[0];
x q[1];
rz(-2.7052077) q[2];
sx q[2];
rz(-0.67989319) q[2];
sx q[2];
rz(1.5117548) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8167953) q[1];
sx q[1];
rz(-0.15843219) q[1];
sx q[1];
rz(-3.0976803) q[1];
rz(-pi) q[2];
rz(1.8154816) q[3];
sx q[3];
rz(-2.8716794) q[3];
sx q[3];
rz(3.1343366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2893082) q[2];
sx q[2];
rz(-0.70312971) q[2];
sx q[2];
rz(0.0011778041) q[2];
rz(-1.7581455) q[3];
sx q[3];
rz(-2.8774084) q[3];
sx q[3];
rz(1.0435102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99172878) q[0];
sx q[0];
rz(-0.1606476) q[0];
sx q[0];
rz(-0.37586656) q[0];
rz(2.1916892) q[1];
sx q[1];
rz(-1.4774731) q[1];
sx q[1];
rz(0.72521597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071142) q[0];
sx q[0];
rz(-1.6867945) q[0];
sx q[0];
rz(-3.0222361) q[0];
rz(-pi) q[1];
rz(2.2967417) q[2];
sx q[2];
rz(-0.84502673) q[2];
sx q[2];
rz(-1.8595075) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23622295) q[1];
sx q[1];
rz(-0.097644383) q[1];
sx q[1];
rz(3.0821441) q[1];
rz(-1.3106724) q[3];
sx q[3];
rz(-0.88444607) q[3];
sx q[3];
rz(-2.7174866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3788562) q[2];
sx q[2];
rz(-1.6894268) q[2];
sx q[2];
rz(0.47259304) q[2];
rz(-3.0224814) q[3];
sx q[3];
rz(-2.5178858) q[3];
sx q[3];
rz(-2.3314893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3177719) q[0];
sx q[0];
rz(-0.20109421) q[0];
sx q[0];
rz(0.066545181) q[0];
rz(-0.19509527) q[1];
sx q[1];
rz(-2.0654443) q[1];
sx q[1];
rz(1.9871575) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6267363) q[0];
sx q[0];
rz(-1.3382698) q[0];
sx q[0];
rz(1.0794845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45656352) q[2];
sx q[2];
rz(-1.6491659) q[2];
sx q[2];
rz(-2.5541277) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9137568) q[1];
sx q[1];
rz(-1.933473) q[1];
sx q[1];
rz(0.11432027) q[1];
rz(0.0129778) q[3];
sx q[3];
rz(-2.1514031) q[3];
sx q[3];
rz(-1.5185771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2164312) q[2];
sx q[2];
rz(-1.5387646) q[2];
sx q[2];
rz(-0.72539854) q[2];
rz(1.7321436) q[3];
sx q[3];
rz(-1.1812482) q[3];
sx q[3];
rz(-0.60666549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055939097) q[0];
sx q[0];
rz(-2.5655209) q[0];
sx q[0];
rz(-2.2357909) q[0];
rz(-2.5869351) q[1];
sx q[1];
rz(-0.83811086) q[1];
sx q[1];
rz(0.14403266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2650288) q[0];
sx q[0];
rz(-1.3110037) q[0];
sx q[0];
rz(-2.0869528) q[0];
x q[1];
rz(-0.085113581) q[2];
sx q[2];
rz(-2.0817698) q[2];
sx q[2];
rz(0.86513317) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8160496) q[1];
sx q[1];
rz(-1.1874949) q[1];
sx q[1];
rz(0.30168962) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8215173) q[3];
sx q[3];
rz(-1.7169239) q[3];
sx q[3];
rz(3.0878029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6683228) q[2];
sx q[2];
rz(-0.39445764) q[2];
sx q[2];
rz(-2.1486166) q[2];
rz(2.9571577) q[3];
sx q[3];
rz(-0.95576972) q[3];
sx q[3];
rz(0.9930281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011878012) q[0];
sx q[0];
rz(-2.5373902) q[0];
sx q[0];
rz(2.7469444) q[0];
rz(2.6036085) q[1];
sx q[1];
rz(-2.4514276) q[1];
sx q[1];
rz(1.949955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25090835) q[0];
sx q[0];
rz(-0.60870752) q[0];
sx q[0];
rz(-2.6842791) q[0];
rz(-pi) q[1];
x q[1];
rz(2.678474) q[2];
sx q[2];
rz(-2.7292433) q[2];
sx q[2];
rz(0.3316484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0318089) q[1];
sx q[1];
rz(-1.1369103) q[1];
sx q[1];
rz(0.96200301) q[1];
rz(-pi) q[2];
rz(-0.58624506) q[3];
sx q[3];
rz(-0.72544155) q[3];
sx q[3];
rz(-0.44437757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.789232) q[2];
sx q[2];
rz(-1.5287986) q[2];
sx q[2];
rz(0.94190502) q[2];
rz(-0.76147979) q[3];
sx q[3];
rz(-1.1250291) q[3];
sx q[3];
rz(-3.1010845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2725459) q[0];
sx q[0];
rz(-1.7750374) q[0];
sx q[0];
rz(-2.0696562) q[0];
rz(2.5275224) q[1];
sx q[1];
rz(-2.2449988) q[1];
sx q[1];
rz(0.44642064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75335129) q[0];
sx q[0];
rz(-2.2470625) q[0];
sx q[0];
rz(-0.36305289) q[0];
rz(-2.5183886) q[2];
sx q[2];
rz(-0.7936306) q[2];
sx q[2];
rz(2.6422184) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.804396) q[1];
sx q[1];
rz(-2.2257518) q[1];
sx q[1];
rz(-1.4319929) q[1];
rz(0.29576755) q[3];
sx q[3];
rz(-0.52433521) q[3];
sx q[3];
rz(1.7400896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8322231) q[2];
sx q[2];
rz(-2.1840405) q[2];
sx q[2];
rz(1.8831801) q[2];
rz(-1.9084515) q[3];
sx q[3];
rz(-2.5150053) q[3];
sx q[3];
rz(-1.3174177) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719139) q[0];
sx q[0];
rz(-0.8154251) q[0];
sx q[0];
rz(-1.423214) q[0];
rz(2.8990959) q[1];
sx q[1];
rz(-0.47912326) q[1];
sx q[1];
rz(1.8877782) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44556955) q[0];
sx q[0];
rz(-1.7544839) q[0];
sx q[0];
rz(1.1653698) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0331633) q[2];
sx q[2];
rz(-0.067010894) q[2];
sx q[2];
rz(-2.5052414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7775443) q[1];
sx q[1];
rz(-1.8692888) q[1];
sx q[1];
rz(2.1095898) q[1];
rz(-pi) q[2];
x q[2];
rz(2.092179) q[3];
sx q[3];
rz(-0.42179042) q[3];
sx q[3];
rz(-2.3613195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8157876) q[2];
sx q[2];
rz(-2.3141404) q[2];
sx q[2];
rz(3.117756) q[2];
rz(1.2787974) q[3];
sx q[3];
rz(-2.8102504) q[3];
sx q[3];
rz(0.7767902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9780289) q[0];
sx q[0];
rz(-1.704287) q[0];
sx q[0];
rz(-1.1821672) q[0];
rz(-0.72721807) q[1];
sx q[1];
rz(-1.6738418) q[1];
sx q[1];
rz(2.3152836) q[1];
rz(-1.6979065) q[2];
sx q[2];
rz(-2.9893176) q[2];
sx q[2];
rz(1.9592374) q[2];
rz(1.4489849) q[3];
sx q[3];
rz(-1.5234608) q[3];
sx q[3];
rz(0.85352637) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
