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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(2.7106078) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(0.42981848) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0551222) q[0];
sx q[0];
rz(-1.4807379) q[0];
sx q[0];
rz(0.49205972) q[0];
rz(-pi) q[1];
x q[1];
rz(1.138428) q[2];
sx q[2];
rz(-1.585344) q[2];
sx q[2];
rz(0.61847875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7215772) q[1];
sx q[1];
rz(-0.41999451) q[1];
sx q[1];
rz(0.53424044) q[1];
x q[2];
rz(-0.63685959) q[3];
sx q[3];
rz(-2.5489293) q[3];
sx q[3];
rz(-0.67837472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(2.1303614) q[2];
rz(2.1003335) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(-2.7019555) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083199) q[0];
sx q[0];
rz(-2.1450873) q[0];
sx q[0];
rz(0.74747768) q[0];
rz(0.29391089) q[1];
sx q[1];
rz(-2.0103318) q[1];
sx q[1];
rz(0.25516137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145484) q[0];
sx q[0];
rz(-2.5393195) q[0];
sx q[0];
rz(1.8667579) q[0];
x q[1];
rz(-2.8361144) q[2];
sx q[2];
rz(-1.7490088) q[2];
sx q[2];
rz(-0.62159789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1595711) q[1];
sx q[1];
rz(-0.82948331) q[1];
sx q[1];
rz(-0.94653603) q[1];
x q[2];
rz(2.9547666) q[3];
sx q[3];
rz(-2.1869724) q[3];
sx q[3];
rz(-2.0957859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1110288) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(-2.2500706) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(-2.7454624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4075277) q[0];
sx q[0];
rz(-1.1363109) q[0];
sx q[0];
rz(-0.044128142) q[0];
rz(1.8924425) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(0.485802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51903782) q[0];
sx q[0];
rz(-1.3455032) q[0];
sx q[0];
rz(1.5507429) q[0];
rz(-pi) q[1];
rz(1.75778) q[2];
sx q[2];
rz(-1.4348363) q[2];
sx q[2];
rz(-0.41278186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1959744) q[1];
sx q[1];
rz(-0.83688078) q[1];
sx q[1];
rz(-0.7463781) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3084675) q[3];
sx q[3];
rz(-2.5841048) q[3];
sx q[3];
rz(-1.1380029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(0.69474727) q[2];
rz(-2.5658549) q[3];
sx q[3];
rz(-1.271194) q[3];
sx q[3];
rz(1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1908252) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(2.8488391) q[0];
rz(2.5726908) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(-2.2946766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900565) q[0];
sx q[0];
rz(-2.7256294) q[0];
sx q[0];
rz(-1.3777758) q[0];
rz(2.2402873) q[2];
sx q[2];
rz(-1.6274823) q[2];
sx q[2];
rz(-2.0079119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0140527) q[1];
sx q[1];
rz(-1.3779281) q[1];
sx q[1];
rz(2.8991153) q[1];
rz(-3.0740739) q[3];
sx q[3];
rz(-1.3775) q[3];
sx q[3];
rz(2.388978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7920821) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(2.9328031) q[2];
rz(2.4236692) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(2.2010402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013041) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(-0.33863246) q[0];
rz(0.70308095) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(-2.3032761) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378483) q[0];
sx q[0];
rz(-2.6347343) q[0];
sx q[0];
rz(-0.54308191) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19836004) q[2];
sx q[2];
rz(-0.83399978) q[2];
sx q[2];
rz(-1.552358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.053239659) q[1];
sx q[1];
rz(-2.2704098) q[1];
sx q[1];
rz(-3.0979473) q[1];
rz(-pi) q[2];
rz(-2.5346967) q[3];
sx q[3];
rz(-1.9268774) q[3];
sx q[3];
rz(-2.0739561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4579953) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(1.0028769) q[2];
rz(0.14497997) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(3.0785479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.018464) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(-0.96328324) q[0];
rz(-0.18114289) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(-1.9021665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604804) q[0];
sx q[0];
rz(-1.943893) q[0];
sx q[0];
rz(2.4831071) q[0];
rz(-pi) q[1];
rz(2.6341637) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(-1.5786174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2378578) q[1];
sx q[1];
rz(-1.1958201) q[1];
sx q[1];
rz(-1.9720248) q[1];
x q[2];
rz(-1.419098) q[3];
sx q[3];
rz(-1.4432419) q[3];
sx q[3];
rz(2.2383245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88508254) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(0.93290848) q[2];
rz(-1.7289303) q[3];
sx q[3];
rz(-1.7191073) q[3];
sx q[3];
rz(-2.9322114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949718) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(-2.3387261) q[0];
rz(-0.74202263) q[1];
sx q[1];
rz(-1.6525533) q[1];
sx q[1];
rz(-2.7313357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0502346) q[0];
sx q[0];
rz(-2.0977533) q[0];
sx q[0];
rz(1.1731253) q[0];
rz(-pi) q[1];
rz(-1.4531936) q[2];
sx q[2];
rz(-2.3315993) q[2];
sx q[2];
rz(-1.2076898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8852119) q[1];
sx q[1];
rz(-1.7325914) q[1];
sx q[1];
rz(-0.075116861) q[1];
rz(-pi) q[2];
rz(-0.65703742) q[3];
sx q[3];
rz(-1.9578303) q[3];
sx q[3];
rz(-1.5157733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86847574) q[2];
sx q[2];
rz(-1.5831524) q[2];
sx q[2];
rz(2.9203019) q[2];
rz(0.41915974) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(2.6656849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9047423) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(1.5297484) q[0];
rz(1.8086241) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.6882247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4299592) q[0];
sx q[0];
rz(-0.89444133) q[0];
sx q[0];
rz(-0.74672378) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4508993) q[2];
sx q[2];
rz(-1.7603619) q[2];
sx q[2];
rz(-1.8360863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21199521) q[1];
sx q[1];
rz(-2.4182165) q[1];
sx q[1];
rz(-2.9138447) q[1];
rz(-1.0681719) q[3];
sx q[3];
rz(-2.8262071) q[3];
sx q[3];
rz(-0.50704765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(-0.010738372) q[2];
rz(-2.8209316) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(-0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43309942) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(2.6668715) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-2.0359998) q[1];
sx q[1];
rz(-1.6620103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06874456) q[0];
sx q[0];
rz(-1.5340051) q[0];
sx q[0];
rz(-1.5307013) q[0];
rz(1.3298673) q[2];
sx q[2];
rz(-0.90902599) q[2];
sx q[2];
rz(2.2994201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2728647) q[1];
sx q[1];
rz(-1.1050779) q[1];
sx q[1];
rz(1.8293423) q[1];
rz(-0.18085705) q[3];
sx q[3];
rz(-1.81362) q[3];
sx q[3];
rz(-2.9202785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.097229615) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(2.7105159) q[2];
rz(-0.19032446) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(-0.88461191) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97656074) q[0];
sx q[0];
rz(-0.29104069) q[0];
sx q[0];
rz(-0.2188368) q[0];
rz(-2.1826375) q[1];
sx q[1];
rz(-0.7364277) q[1];
sx q[1];
rz(-1.7494019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0420211) q[0];
sx q[0];
rz(-0.19992689) q[0];
sx q[0];
rz(-0.26148343) q[0];
rz(-1.4264246) q[2];
sx q[2];
rz(-1.1544268) q[2];
sx q[2];
rz(-1.8954111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.47962418) q[1];
sx q[1];
rz(-1.2733409) q[1];
sx q[1];
rz(0.65315078) q[1];
x q[2];
rz(-0.99361657) q[3];
sx q[3];
rz(-2.4419067) q[3];
sx q[3];
rz(0.71551871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2739111) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(-1.8217746) q[2];
rz(0.53226081) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(1.5490612) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821447) q[0];
sx q[0];
rz(-1.5503217) q[0];
sx q[0];
rz(0.4009608) q[0];
rz(2.9136912) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(-2.7264222) q[2];
sx q[2];
rz(-1.3317458) q[2];
sx q[2];
rz(1.7785704) q[2];
rz(-2.9798672) q[3];
sx q[3];
rz(-0.99746708) q[3];
sx q[3];
rz(0.4494638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
