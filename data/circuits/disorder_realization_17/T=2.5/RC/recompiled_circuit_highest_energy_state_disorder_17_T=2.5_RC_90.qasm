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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(-2.5435574) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(-0.66507566) q[1];
sx q[1];
rz(-2.9887078) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4035002) q[0];
sx q[0];
rz(-3.0041143) q[0];
sx q[0];
rz(-2.1901778) q[0];
rz(-0.10993425) q[2];
sx q[2];
rz(-2.6510932) q[2];
sx q[2];
rz(-2.5562499) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28935948) q[1];
sx q[1];
rz(-2.8948862) q[1];
sx q[1];
rz(1.0323332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5737719) q[3];
sx q[3];
rz(-1.5248393) q[3];
sx q[3];
rz(-2.9480235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87024706) q[2];
sx q[2];
rz(-0.18834867) q[2];
sx q[2];
rz(-1.9625473) q[2];
rz(2.7315268) q[3];
sx q[3];
rz(-1.054801) q[3];
sx q[3];
rz(-0.89455354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9533185) q[0];
sx q[0];
rz(-2.4725547) q[0];
sx q[0];
rz(0.27729312) q[0];
rz(2.9412728) q[1];
sx q[1];
rz(-0.89626139) q[1];
sx q[1];
rz(3.0576113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1071716) q[0];
sx q[0];
rz(-1.1033397) q[0];
sx q[0];
rz(-0.22308992) q[0];
x q[1];
rz(2.4139187) q[2];
sx q[2];
rz(-1.6817589) q[2];
sx q[2];
rz(0.79906711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9586304) q[1];
sx q[1];
rz(-1.1144899) q[1];
sx q[1];
rz(2.0668304) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53392729) q[3];
sx q[3];
rz(-1.797343) q[3];
sx q[3];
rz(1.904303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9597441) q[2];
sx q[2];
rz(-0.66755787) q[2];
sx q[2];
rz(-2.3954771) q[2];
rz(-2.5203846) q[3];
sx q[3];
rz(-0.64380232) q[3];
sx q[3];
rz(-1.752689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32034945) q[0];
sx q[0];
rz(-0.72750434) q[0];
sx q[0];
rz(-1.9387091) q[0];
rz(-0.41162833) q[1];
sx q[1];
rz(-1.8609214) q[1];
sx q[1];
rz(1.113755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4899556) q[0];
sx q[0];
rz(-0.86511602) q[0];
sx q[0];
rz(-2.5791427) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95543785) q[2];
sx q[2];
rz(-1.2005383) q[2];
sx q[2];
rz(2.7317765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71761638) q[1];
sx q[1];
rz(-0.90002093) q[1];
sx q[1];
rz(-0.17016478) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90818543) q[3];
sx q[3];
rz(-0.81530967) q[3];
sx q[3];
rz(-2.1774928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0959629) q[2];
sx q[2];
rz(-0.7926422) q[2];
sx q[2];
rz(1.9795798) q[2];
rz(2.3368733) q[3];
sx q[3];
rz(-1.2730803) q[3];
sx q[3];
rz(-1.2812251) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74748892) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.5149186) q[0];
rz(2.6331242) q[1];
sx q[1];
rz(-0.6503121) q[1];
sx q[1];
rz(1.5079087) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8554022) q[0];
sx q[0];
rz(-1.5686745) q[0];
sx q[0];
rz(-1.4197664) q[0];
rz(0.99782439) q[2];
sx q[2];
rz(-0.25616872) q[2];
sx q[2];
rz(2.1081971) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3292309) q[1];
sx q[1];
rz(-2.8155102) q[1];
sx q[1];
rz(-1.9291758) q[1];
rz(-2.1056469) q[3];
sx q[3];
rz(-0.7300762) q[3];
sx q[3];
rz(1.4393782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7278829) q[2];
sx q[2];
rz(-1.1943613) q[2];
sx q[2];
rz(-0.74753648) q[2];
rz(-1.9109292) q[3];
sx q[3];
rz(-0.80941284) q[3];
sx q[3];
rz(2.6443853) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704945) q[0];
sx q[0];
rz(-2.0132988) q[0];
sx q[0];
rz(-1.6823912) q[0];
rz(-0.35405007) q[1];
sx q[1];
rz(-2.470128) q[1];
sx q[1];
rz(-2.5669602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0316236) q[0];
sx q[0];
rz(-1.4935483) q[0];
sx q[0];
rz(1.5573274) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8425521) q[2];
sx q[2];
rz(-2.3210706) q[2];
sx q[2];
rz(2.3472903) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4230391) q[1];
sx q[1];
rz(-1.7977814) q[1];
sx q[1];
rz(-1.6237698) q[1];
rz(-3.026612) q[3];
sx q[3];
rz(-1.774174) q[3];
sx q[3];
rz(0.16826828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5057126) q[2];
sx q[2];
rz(-1.1735703) q[2];
sx q[2];
rz(-2.9393348) q[2];
rz(-1.0328736) q[3];
sx q[3];
rz(-2.5346916) q[3];
sx q[3];
rz(3.0752944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1101087) q[0];
sx q[0];
rz(-0.38723543) q[0];
sx q[0];
rz(-2.7701344) q[0];
rz(2.8455041) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(-2.9782226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59054331) q[0];
sx q[0];
rz(-1.7544909) q[0];
sx q[0];
rz(-0.50409533) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1131606) q[2];
sx q[2];
rz(-2.073098) q[2];
sx q[2];
rz(0.1973803) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5886053) q[1];
sx q[1];
rz(-1.3901911) q[1];
sx q[1];
rz(-0.44579472) q[1];
rz(-pi) q[2];
rz(-2.9904731) q[3];
sx q[3];
rz(-0.88081473) q[3];
sx q[3];
rz(2.2595764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0614193) q[2];
sx q[2];
rz(-2.8714608) q[2];
sx q[2];
rz(2.488625) q[2];
rz(0.4717007) q[3];
sx q[3];
rz(-2.3931849) q[3];
sx q[3];
rz(2.4145943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93040526) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(1.0373254) q[0];
rz(-0.51517454) q[1];
sx q[1];
rz(-1.2778792) q[1];
sx q[1];
rz(-1.3264664) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65086397) q[0];
sx q[0];
rz(-2.0959637) q[0];
sx q[0];
rz(-0.45142169) q[0];
rz(-1.6112566) q[2];
sx q[2];
rz(-1.2190814) q[2];
sx q[2];
rz(2.2556698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80250128) q[1];
sx q[1];
rz(-2.3330124) q[1];
sx q[1];
rz(-1.436961) q[1];
rz(-pi) q[2];
rz(-3.0781704) q[3];
sx q[3];
rz(-1.8759955) q[3];
sx q[3];
rz(-2.2535107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6933763) q[2];
sx q[2];
rz(-0.46458149) q[2];
sx q[2];
rz(-2.8153937) q[2];
rz(0.97801963) q[3];
sx q[3];
rz(-1.8947442) q[3];
sx q[3];
rz(3.0514362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6982894) q[0];
sx q[0];
rz(-2.8841618) q[0];
sx q[0];
rz(0.28939104) q[0];
rz(-3.1293213) q[1];
sx q[1];
rz(-2.8616276) q[1];
sx q[1];
rz(1.6562921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78083437) q[0];
sx q[0];
rz(-0.88874431) q[0];
sx q[0];
rz(-2.2498796) q[0];
rz(-pi) q[1];
rz(2.9180235) q[2];
sx q[2];
rz(-1.6235329) q[2];
sx q[2];
rz(-1.2745672) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0721133) q[1];
sx q[1];
rz(-2.8918307) q[1];
sx q[1];
rz(0.62432557) q[1];
x q[2];
rz(2.8217153) q[3];
sx q[3];
rz(-1.1231224) q[3];
sx q[3];
rz(-3.0818528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5926008) q[2];
sx q[2];
rz(-1.5718549) q[2];
sx q[2];
rz(2.9686417) q[2];
rz(-1.3744099) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(1.4779199) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7462815) q[0];
sx q[0];
rz(-0.44773856) q[0];
sx q[0];
rz(-2.3315499) q[0];
rz(0.41656247) q[1];
sx q[1];
rz(-0.77605334) q[1];
sx q[1];
rz(0.3449482) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589011) q[0];
sx q[0];
rz(-1.693629) q[0];
sx q[0];
rz(-0.71101153) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8497963) q[2];
sx q[2];
rz(-1.4115745) q[2];
sx q[2];
rz(0.70406841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7564063) q[1];
sx q[1];
rz(-2.029772) q[1];
sx q[1];
rz(-2.7800757) q[1];
rz(-3.0250375) q[3];
sx q[3];
rz(-1.6892589) q[3];
sx q[3];
rz(-0.20878709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42744669) q[2];
sx q[2];
rz(-2.2553208) q[2];
sx q[2];
rz(-3.0575338) q[2];
rz(2.2345624) q[3];
sx q[3];
rz(-1.4182988) q[3];
sx q[3];
rz(2.1840054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3323988) q[0];
sx q[0];
rz(-3.0539303) q[0];
sx q[0];
rz(-2.8293389) q[0];
rz(-0.23278438) q[1];
sx q[1];
rz(-2.7504031) q[1];
sx q[1];
rz(1.8544474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904743) q[0];
sx q[0];
rz(-0.93997191) q[0];
sx q[0];
rz(1.1680956) q[0];
rz(-pi) q[1];
rz(0.63913362) q[2];
sx q[2];
rz(-1.5331563) q[2];
sx q[2];
rz(-2.2485178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.171265) q[1];
sx q[1];
rz(-1.1210151) q[1];
sx q[1];
rz(-0.27189092) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5165043) q[3];
sx q[3];
rz(-1.8755138) q[3];
sx q[3];
rz(-0.24686043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9707668) q[2];
sx q[2];
rz(-1.8584741) q[2];
sx q[2];
rz(1.6440561) q[2];
rz(-1.0981285) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(-0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1371786) q[0];
sx q[0];
rz(-1.7902086) q[0];
sx q[0];
rz(-0.93850346) q[0];
rz(1.8347523) q[1];
sx q[1];
rz(-0.88941457) q[1];
sx q[1];
rz(-0.79759146) q[1];
rz(1.7457444) q[2];
sx q[2];
rz(-0.69654225) q[2];
sx q[2];
rz(-2.2707743) q[2];
rz(-2.9994503) q[3];
sx q[3];
rz(-0.90928838) q[3];
sx q[3];
rz(1.7293775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
