OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6131634) q[0];
sx q[0];
rz(-2.0818721) q[0];
sx q[0];
rz(-0.73097316) q[0];
rz(1.641474) q[1];
sx q[1];
rz(5.2483622) q[1];
sx q[1];
rz(11.622826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65981728) q[0];
sx q[0];
rz(-0.075591139) q[0];
sx q[0];
rz(-0.48150058) q[0];
x q[1];
rz(-2.4906292) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-2.0219321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4409677) q[1];
sx q[1];
rz(-1.96083) q[1];
sx q[1];
rz(2.7729211) q[1];
rz(-pi) q[2];
rz(-0.54240534) q[3];
sx q[3];
rz(-0.35202682) q[3];
sx q[3];
rz(-0.56234081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(-1.4261774) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(-2.1616518) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(0.59894484) q[0];
rz(-1.8006181) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(-2.1751931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37106284) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(2.4735527) q[0];
rz(-1.6205377) q[2];
sx q[2];
rz(-2.192394) q[2];
sx q[2];
rz(-0.13812401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71827723) q[1];
sx q[1];
rz(-1.4916972) q[1];
sx q[1];
rz(3.0273816) q[1];
rz(-0.74921272) q[3];
sx q[3];
rz(-1.2470761) q[3];
sx q[3];
rz(-2.2183228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.085658375) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(-1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.1874636) q[0];
rz(1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(-1.3175861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353969) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(-2.8185185) q[0];
x q[1];
rz(-1.7682398) q[2];
sx q[2];
rz(-1.7806446) q[2];
sx q[2];
rz(0.91453493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8220362) q[1];
sx q[1];
rz(-1.469538) q[1];
sx q[1];
rz(0.42145573) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0242689) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(-0.26457149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(-0.2066361) q[2];
rz(0.7080428) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.9151691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7851631) q[0];
sx q[0];
rz(-0.7298846) q[0];
sx q[0];
rz(-1.4714144) q[0];
x q[1];
rz(-0.45986508) q[2];
sx q[2];
rz(-2.2042639) q[2];
sx q[2];
rz(-1.8749274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1434506) q[1];
sx q[1];
rz(-1.7805903) q[1];
sx q[1];
rz(1.7721121) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5049582) q[3];
sx q[3];
rz(-2.5599179) q[3];
sx q[3];
rz(3.0326774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(0.99299661) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(0.58247724) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156508) q[0];
sx q[0];
rz(-1.6424718) q[0];
sx q[0];
rz(-0.96836758) q[0];
x q[1];
rz(1.3864473) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(-2.9550936) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(0.71449844) q[1];
rz(-pi) q[2];
rz(1.5037687) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(-1.0073347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65486583) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(0.081710903) q[2];
rz(0.47406667) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896647) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-2.0369464) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4545126) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(1.7772654) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32128895) q[2];
sx q[2];
rz(-2.4761204) q[2];
sx q[2];
rz(1.2207536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0749952) q[1];
sx q[1];
rz(-1.0269594) q[1];
sx q[1];
rz(-2.4165513) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1052746) q[3];
sx q[3];
rz(-1.421531) q[3];
sx q[3];
rz(-2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(-2.5863623) q[2];
rz(0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5884488) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-2.0297091) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0985078) q[0];
sx q[0];
rz(-2.3102009) q[0];
sx q[0];
rz(2.1466473) q[0];
x q[1];
rz(-1.6484117) q[2];
sx q[2];
rz(-1.5570119) q[2];
sx q[2];
rz(0.26253653) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.249914) q[1];
sx q[1];
rz(-0.87712446) q[1];
sx q[1];
rz(1.2077043) q[1];
rz(1.3248596) q[3];
sx q[3];
rz(-1.2517559) q[3];
sx q[3];
rz(-2.6815573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-0.61541921) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62548816) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(-1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-0.63751784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2190327) q[0];
sx q[0];
rz(-2.1171283) q[0];
sx q[0];
rz(0.68740293) q[0];
x q[1];
rz(-1.9279187) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(1.9879607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70367614) q[1];
sx q[1];
rz(-0.99214593) q[1];
sx q[1];
rz(1.089536) q[1];
rz(-1.4333126) q[3];
sx q[3];
rz(-1.004389) q[3];
sx q[3];
rz(-2.8979104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(2.0054224) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4207941) q[0];
sx q[0];
rz(-1.7013229) q[0];
sx q[0];
rz(2.7230848) q[0];
rz(2.9416111) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(2.8560864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.059541313) q[1];
sx q[1];
rz(-2.6350628) q[1];
sx q[1];
rz(-0.7854714) q[1];
rz(-0.34992643) q[3];
sx q[3];
rz(-2.3438128) q[3];
sx q[3];
rz(-1.6328904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(0.96735111) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6185146) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(-0.05649795) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79849762) q[0];
sx q[0];
rz(-1.4416749) q[0];
sx q[0];
rz(-1.9341747) q[0];
rz(2.8374412) q[2];
sx q[2];
rz(-1.4545822) q[2];
sx q[2];
rz(-2.846037) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2056634) q[1];
sx q[1];
rz(-0.91055369) q[1];
sx q[1];
rz(0.6026938) q[1];
rz(-2.9910562) q[3];
sx q[3];
rz(-1.5621462) q[3];
sx q[3];
rz(-1.0516402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.848032) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(0.85110101) q[2];
sx q[2];
rz(-0.38817482) q[2];
sx q[2];
rz(2.9183802) q[2];
rz(-0.33070926) q[3];
sx q[3];
rz(-1.0996795) q[3];
sx q[3];
rz(-0.75547937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];