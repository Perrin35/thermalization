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
rz(1.2443378) q[0];
sx q[0];
rz(2.3490348) q[0];
sx q[0];
rz(9.6777182) q[0];
rz(2.3674372) q[1];
sx q[1];
rz(-2.7634662) q[1];
sx q[1];
rz(2.7546496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9772676) q[0];
sx q[0];
rz(-1.2480134) q[0];
sx q[0];
rz(-0.27077814) q[0];
x q[1];
rz(-3.0091506) q[2];
sx q[2];
rz(-0.58183607) q[2];
sx q[2];
rz(-2.066032) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2110212) q[1];
sx q[1];
rz(-1.4912349) q[1];
sx q[1];
rz(0.08427517) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0798595) q[3];
sx q[3];
rz(-1.0359952) q[3];
sx q[3];
rz(-2.8859069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.046254961) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(1.4026027) q[2];
rz(1.4143573) q[3];
sx q[3];
rz(-2.4284913) q[3];
sx q[3];
rz(-2.6426017) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801179) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(0.71037355) q[0];
rz(1.0164227) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(0.76510915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9020664) q[0];
sx q[0];
rz(-1.2529116) q[0];
sx q[0];
rz(2.5796298) q[0];
x q[1];
rz(-2.3813407) q[2];
sx q[2];
rz(-0.85101247) q[2];
sx q[2];
rz(-0.30530294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2110846) q[1];
sx q[1];
rz(-2.395906) q[1];
sx q[1];
rz(-1.0812283) q[1];
rz(-pi) q[2];
rz(-1.9809057) q[3];
sx q[3];
rz(-2.7795459) q[3];
sx q[3];
rz(-2.1459879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1171099) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(2.4369241) q[2];
rz(-2.9674528) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(0.38159889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0524549) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-0.88743368) q[0];
rz(-1.2499836) q[1];
sx q[1];
rz(-1.2108112) q[1];
sx q[1];
rz(-1.3608305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75653532) q[0];
sx q[0];
rz(-2.4891315) q[0];
sx q[0];
rz(-2.4908309) q[0];
x q[1];
rz(1.4329696) q[2];
sx q[2];
rz(-1.4915474) q[2];
sx q[2];
rz(-0.85286959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4579297) q[1];
sx q[1];
rz(-0.81240801) q[1];
sx q[1];
rz(-0.13813604) q[1];
rz(-pi) q[2];
rz(1.4668457) q[3];
sx q[3];
rz(-0.39109215) q[3];
sx q[3];
rz(-1.1268738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36797324) q[2];
sx q[2];
rz(-0.36586389) q[2];
sx q[2];
rz(1.492929) q[2];
rz(-2.6922373) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0538977) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.4087403) q[0];
rz(-2.159481) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(1.7353479) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9179233) q[0];
sx q[0];
rz(-1.6817998) q[0];
sx q[0];
rz(1.3551606) q[0];
rz(-pi) q[1];
rz(-3.1020145) q[2];
sx q[2];
rz(-1.5180032) q[2];
sx q[2];
rz(2.9556208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2212053) q[1];
sx q[1];
rz(-1.3246598) q[1];
sx q[1];
rz(2.1046776) q[1];
rz(-pi) q[2];
rz(2.8729183) q[3];
sx q[3];
rz(-2.2019535) q[3];
sx q[3];
rz(2.1674066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1534319) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(0.15596685) q[2];
rz(0.65034741) q[3];
sx q[3];
rz(-1.9675083) q[3];
sx q[3];
rz(-0.11421886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0735737) q[0];
sx q[0];
rz(-1.8719712) q[0];
sx q[0];
rz(-0.20075783) q[0];
rz(1.9643895) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(-1.1600561) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667127) q[0];
sx q[0];
rz(-1.794941) q[0];
sx q[0];
rz(0.17205162) q[0];
rz(-1.6408553) q[2];
sx q[2];
rz(-1.4558917) q[2];
sx q[2];
rz(-2.9821378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1694035) q[1];
sx q[1];
rz(-0.47220017) q[1];
sx q[1];
rz(0.71581033) q[1];
x q[2];
rz(-0.20108624) q[3];
sx q[3];
rz(-1.0498091) q[3];
sx q[3];
rz(-2.9231284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7572299) q[2];
sx q[2];
rz(-0.94567662) q[2];
sx q[2];
rz(2.6833656) q[2];
rz(2.5107757) q[3];
sx q[3];
rz(-1.2354555) q[3];
sx q[3];
rz(0.49921504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45796564) q[0];
sx q[0];
rz(-0.85245913) q[0];
sx q[0];
rz(3.1347347) q[0];
rz(2.9099756) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(2.0793656) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7201014) q[0];
sx q[0];
rz(-1.4473697) q[0];
sx q[0];
rz(0.31436679) q[0];
rz(-pi) q[1];
rz(-0.3756282) q[2];
sx q[2];
rz(-0.53002073) q[2];
sx q[2];
rz(-0.14800581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4783027) q[1];
sx q[1];
rz(-1.8243102) q[1];
sx q[1];
rz(0.67914825) q[1];
rz(-pi) q[2];
rz(-2.2799479) q[3];
sx q[3];
rz(-1.6966736) q[3];
sx q[3];
rz(-1.6050086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0007533) q[2];
sx q[2];
rz(-1.8893628) q[2];
sx q[2];
rz(1.8737277) q[2];
rz(2.2800692) q[3];
sx q[3];
rz(-2.0778766) q[3];
sx q[3];
rz(-1.4388194) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93290257) q[0];
sx q[0];
rz(-2.8080495) q[0];
sx q[0];
rz(0.21555899) q[0];
rz(2.6453099) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(2.9821679) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27860363) q[0];
sx q[0];
rz(-1.5734324) q[0];
sx q[0];
rz(0.9359261) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3770476) q[2];
sx q[2];
rz(-1.1074578) q[2];
sx q[2];
rz(-3.0510356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9558952) q[1];
sx q[1];
rz(-1.5296827) q[1];
sx q[1];
rz(1.8061045) q[1];
rz(-pi) q[2];
rz(-0.46298085) q[3];
sx q[3];
rz(-2.0714875) q[3];
sx q[3];
rz(-0.41741308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7364007) q[2];
sx q[2];
rz(-1.1242194) q[2];
sx q[2];
rz(-2.6250725) q[2];
rz(-1.9308331) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(1.3794544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(2.6013689) q[0];
rz(-1.5243439) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(1.2670955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62477797) q[0];
sx q[0];
rz(-2.5886726) q[0];
sx q[0];
rz(0.42527683) q[0];
rz(-pi) q[1];
rz(-0.71155602) q[2];
sx q[2];
rz(-2.2268953) q[2];
sx q[2];
rz(3.1392821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6798576) q[1];
sx q[1];
rz(-1.9377794) q[1];
sx q[1];
rz(-1.7253897) q[1];
rz(-pi) q[2];
rz(2.0319875) q[3];
sx q[3];
rz(-2.4465585) q[3];
sx q[3];
rz(-1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8871062) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(-2.9642588) q[2];
rz(-1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(-2.9593318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3607445) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(0.10072197) q[0];
rz(1.1489457) q[1];
sx q[1];
rz(-1.1958242) q[1];
sx q[1];
rz(-1.1740059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.188499) q[0];
sx q[0];
rz(-1.3818912) q[0];
sx q[0];
rz(-0.61236463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4743382) q[2];
sx q[2];
rz(-2.3347179) q[2];
sx q[2];
rz(-0.45713666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8619949) q[1];
sx q[1];
rz(-1.5902385) q[1];
sx q[1];
rz(-1.5330381) q[1];
x q[2];
rz(0.87971808) q[3];
sx q[3];
rz(-1.6501325) q[3];
sx q[3];
rz(-2.9226231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9775057) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(-2.738415) q[2];
rz(-2.4781503) q[3];
sx q[3];
rz(-2.5622029) q[3];
sx q[3];
rz(0.0042393953) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73096257) q[0];
sx q[0];
rz(-2.0616489) q[0];
sx q[0];
rz(0.078911111) q[0];
rz(2.2321189) q[1];
sx q[1];
rz(-1.0458922) q[1];
sx q[1];
rz(2.7452819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27061227) q[0];
sx q[0];
rz(-1.5617234) q[0];
sx q[0];
rz(-1.5864282) q[0];
rz(1.5780007) q[2];
sx q[2];
rz(-0.77709891) q[2];
sx q[2];
rz(2.8030628) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24416284) q[1];
sx q[1];
rz(-1.3363991) q[1];
sx q[1];
rz(0.48903709) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2811398) q[3];
sx q[3];
rz(-2.3858968) q[3];
sx q[3];
rz(-2.6967654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69184715) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(2.3425102) q[2];
rz(-3.0633022) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.4770724) q[0];
sx q[0];
rz(-1.399853) q[0];
sx q[0];
rz(0.54367263) q[0];
rz(1.9480582) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(0.1706201) q[2];
sx q[2];
rz(-2.0224051) q[2];
sx q[2];
rz(-0.97013459) q[2];
rz(-1.0738437) q[3];
sx q[3];
rz(-2.2324149) q[3];
sx q[3];
rz(0.80692337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
