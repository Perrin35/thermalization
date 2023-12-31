OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(-2.3606665) q[0];
sx q[0];
rz(0.20679064) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(-2.9614255) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5313523) q[0];
sx q[0];
rz(-0.58824476) q[0];
sx q[0];
rz(-2.504185) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51282672) q[2];
sx q[2];
rz(-1.6851808) q[2];
sx q[2];
rz(2.0778542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.441592) q[1];
sx q[1];
rz(-0.68725005) q[1];
sx q[1];
rz(-1.6874466) q[1];
x q[2];
rz(0.68182919) q[3];
sx q[3];
rz(-1.7141984) q[3];
sx q[3];
rz(-1.486206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(-1.8475378) q[2];
rz(-0.40575746) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(-0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92293537) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(-3.0157715) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(-1.3719826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94700891) q[0];
sx q[0];
rz(-0.8964552) q[0];
sx q[0];
rz(1.6499004) q[0];
x q[1];
rz(-1.4470909) q[2];
sx q[2];
rz(-2.1134085) q[2];
sx q[2];
rz(1.5109085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0224277) q[1];
sx q[1];
rz(-0.25401527) q[1];
sx q[1];
rz(2.959842) q[1];
x q[2];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.7892924) q[3];
sx q[3];
rz(-0.0056497638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25386086) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(1.0590142) q[0];
rz(1.9937218) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36944593) q[0];
sx q[0];
rz(-1.9659974) q[0];
sx q[0];
rz(2.4482083) q[0];
x q[1];
rz(1.9469444) q[2];
sx q[2];
rz(-0.40824879) q[2];
sx q[2];
rz(-2.2779235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4050582) q[1];
sx q[1];
rz(-1.5824123) q[1];
sx q[1];
rz(-1.8809153) q[1];
rz(-2.2952609) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(-2.4785329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(-1.1188544) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.3942962) q[0];
rz(2.1030203) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(-2.0746453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2571714) q[0];
sx q[0];
rz(-1.0413678) q[0];
sx q[0];
rz(0.04495312) q[0];
rz(-1.489584) q[2];
sx q[2];
rz(-1.3021384) q[2];
sx q[2];
rz(3.0693698) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2414788) q[1];
sx q[1];
rz(-2.334324) q[1];
sx q[1];
rz(0.57358731) q[1];
rz(0.22575836) q[3];
sx q[3];
rz(-2.0421713) q[3];
sx q[3];
rz(-2.2798722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(-2.6211522) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(2.382544) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(-0.83475137) q[0];
rz(1.2443776) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(1.429819) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058250931) q[0];
sx q[0];
rz(-1.5081076) q[0];
sx q[0];
rz(1.5682674) q[0];
x q[1];
rz(-2.8565035) q[2];
sx q[2];
rz(-1.8008917) q[2];
sx q[2];
rz(-0.5321815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5606219) q[1];
sx q[1];
rz(-2.3498658) q[1];
sx q[1];
rz(3.1091299) q[1];
x q[2];
rz(-0.17794869) q[3];
sx q[3];
rz(-2.2214409) q[3];
sx q[3];
rz(1.734317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66118583) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(2.8018518) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(-2.5710411) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3756838) q[0];
sx q[0];
rz(-1.8609957) q[0];
sx q[0];
rz(-1.1955839) q[0];
rz(0.678755) q[2];
sx q[2];
rz(-2.5828913) q[2];
sx q[2];
rz(-0.15685454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30299444) q[1];
sx q[1];
rz(-2.0717151) q[1];
sx q[1];
rz(-0.8143199) q[1];
rz(0.029383226) q[3];
sx q[3];
rz(-0.64850649) q[3];
sx q[3];
rz(0.23119584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-2.3449507) q[2];
rz(-2.8213275) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.8626574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657848) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(-2.7835223) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.8310865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.755237) q[0];
sx q[0];
rz(-2.1023395) q[0];
sx q[0];
rz(1.7887572) q[0];
rz(2.481776) q[2];
sx q[2];
rz(-2.3265127) q[2];
sx q[2];
rz(1.3224524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9999034) q[1];
sx q[1];
rz(-2.018376) q[1];
sx q[1];
rz(0.56772851) q[1];
x q[2];
rz(2.2745423) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(1.5501319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33401176) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(-2.5618662) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(-1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26577935) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(1.9050725) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(-1.1901201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8461788) q[0];
sx q[0];
rz(-2.9018887) q[0];
sx q[0];
rz(-1.833605) q[0];
x q[1];
rz(1.4899848) q[2];
sx q[2];
rz(-1.6390071) q[2];
sx q[2];
rz(-2.8516172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78672532) q[1];
sx q[1];
rz(-2.7495972) q[1];
sx q[1];
rz(2.3945432) q[1];
x q[2];
rz(1.0134775) q[3];
sx q[3];
rz(-1.9625469) q[3];
sx q[3];
rz(1.980892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3999346) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-0.41440543) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(2.8167021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906616) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.855809) q[0];
sx q[0];
rz(-1.7109509) q[0];
sx q[0];
rz(-2.6967718) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29813913) q[2];
sx q[2];
rz(-0.65152822) q[2];
sx q[2];
rz(-0.27116129) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8543429) q[1];
sx q[1];
rz(-2.2209475) q[1];
sx q[1];
rz(-1.2681505) q[1];
rz(-0.98032326) q[3];
sx q[3];
rz(-1.983641) q[3];
sx q[3];
rz(-2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7342547) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(-0.43241832) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(-2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(0.7235136) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7454119) q[0];
sx q[0];
rz(-1.4416845) q[0];
sx q[0];
rz(-2.4095035) q[0];
rz(-pi) q[1];
rz(1.6610442) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(1.2573164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.77955058) q[1];
sx q[1];
rz(-2.3099646) q[1];
sx q[1];
rz(-1.9401624) q[1];
rz(-pi) q[2];
rz(-1.5756597) q[3];
sx q[3];
rz(-1.8412207) q[3];
sx q[3];
rz(1.5750386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92131203) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(2.5881361) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217459) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
rz(-0.28221054) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(2.4137905) q[2];
sx q[2];
rz(-0.62725485) q[2];
sx q[2];
rz(-2.2841452) q[2];
rz(-1.2890733) q[3];
sx q[3];
rz(-1.1829794) q[3];
sx q[3];
rz(1.1869528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
