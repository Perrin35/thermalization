OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8389799) q[0];
sx q[0];
rz(-1.7054727) q[0];
sx q[0];
rz(2.9247395) q[0];
rz(-0.4878374) q[1];
sx q[1];
rz(-0.87895972) q[1];
sx q[1];
rz(-0.15329696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3756806) q[0];
sx q[0];
rz(-1.8366251) q[0];
sx q[0];
rz(-1.2714809) q[0];
rz(1.054864) q[2];
sx q[2];
rz(-1.9867975) q[2];
sx q[2];
rz(1.6226559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5342321) q[1];
sx q[1];
rz(-0.56391729) q[1];
sx q[1];
rz(-0.47718559) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80634597) q[3];
sx q[3];
rz(-1.4571179) q[3];
sx q[3];
rz(-2.0557085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6903901) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(0.85980493) q[2];
rz(-0.23990038) q[3];
sx q[3];
rz(-2.4247215) q[3];
sx q[3];
rz(0.2828323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4848223) q[0];
sx q[0];
rz(-1.4161994) q[0];
sx q[0];
rz(2.2221478) q[0];
rz(0.8473618) q[1];
sx q[1];
rz(-2.5571926) q[1];
sx q[1];
rz(1.1712317) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15554504) q[0];
sx q[0];
rz(-1.520507) q[0];
sx q[0];
rz(1.3119447) q[0];
rz(-2.2861346) q[2];
sx q[2];
rz(-1.3163065) q[2];
sx q[2];
rz(1.364691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50368631) q[1];
sx q[1];
rz(-1.827938) q[1];
sx q[1];
rz(-1.9288344) q[1];
rz(-pi) q[2];
rz(0.53816157) q[3];
sx q[3];
rz(-1.5892913) q[3];
sx q[3];
rz(1.5876113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.072711572) q[2];
sx q[2];
rz(-2.2525807) q[2];
sx q[2];
rz(1.9286801) q[2];
rz(2.9291901) q[3];
sx q[3];
rz(-2.0272777) q[3];
sx q[3];
rz(0.67306486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251959) q[0];
sx q[0];
rz(-1.238751) q[0];
sx q[0];
rz(0.58309251) q[0];
rz(-1.8816226) q[1];
sx q[1];
rz(-2.5446353) q[1];
sx q[1];
rz(-1.8276385) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010600032) q[0];
sx q[0];
rz(-2.6525462) q[0];
sx q[0];
rz(2.2940192) q[0];
rz(-pi) q[1];
rz(-0.73380664) q[2];
sx q[2];
rz(-1.5948405) q[2];
sx q[2];
rz(2.1847385) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62210768) q[1];
sx q[1];
rz(-1.4478323) q[1];
sx q[1];
rz(-0.019799063) q[1];
rz(-pi) q[2];
rz(-1.9989955) q[3];
sx q[3];
rz(-1.7230265) q[3];
sx q[3];
rz(-2.5604932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49695697) q[2];
sx q[2];
rz(-2.3778215) q[2];
sx q[2];
rz(2.606707) q[2];
rz(0.48318091) q[3];
sx q[3];
rz(-1.5213685) q[3];
sx q[3];
rz(-0.11972891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508761) q[0];
sx q[0];
rz(-3.0829939) q[0];
sx q[0];
rz(-1.6706049) q[0];
rz(1.2141256) q[1];
sx q[1];
rz(-1.7045226) q[1];
sx q[1];
rz(-0.4462744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4497919) q[0];
sx q[0];
rz(-1.746382) q[0];
sx q[0];
rz(0.069616074) q[0];
rz(-pi) q[1];
rz(0.52881119) q[2];
sx q[2];
rz(-1.5545648) q[2];
sx q[2];
rz(2.0179786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7301) q[1];
sx q[1];
rz(-0.92918438) q[1];
sx q[1];
rz(1.8909251) q[1];
rz(1.5851444) q[3];
sx q[3];
rz(-1.8135462) q[3];
sx q[3];
rz(-1.2402463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92397592) q[2];
sx q[2];
rz(-2.6805704) q[2];
sx q[2];
rz(-0.32611845) q[2];
rz(2.7255132) q[3];
sx q[3];
rz(-1.5030428) q[3];
sx q[3];
rz(0.7777586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071082696) q[0];
sx q[0];
rz(-1.5776881) q[0];
sx q[0];
rz(0.34410205) q[0];
rz(0.35585078) q[1];
sx q[1];
rz(-2.3882723) q[1];
sx q[1];
rz(2.8867302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0094285) q[0];
sx q[0];
rz(-0.71064204) q[0];
sx q[0];
rz(-0.406213) q[0];
x q[1];
rz(-1.3717669) q[2];
sx q[2];
rz(-2.6712199) q[2];
sx q[2];
rz(-0.34188893) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30789852) q[1];
sx q[1];
rz(-2.1547518) q[1];
sx q[1];
rz(-2.6112982) q[1];
rz(-pi) q[2];
rz(-1.9790824) q[3];
sx q[3];
rz(-0.84813839) q[3];
sx q[3];
rz(1.4681787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43763375) q[2];
sx q[2];
rz(-0.86486977) q[2];
sx q[2];
rz(3.051009) q[2];
rz(0.80638742) q[3];
sx q[3];
rz(-1.839829) q[3];
sx q[3];
rz(-1.3069794) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0818455) q[0];
sx q[0];
rz(-1.2223926) q[0];
sx q[0];
rz(0.61862373) q[0];
rz(2.5289358) q[1];
sx q[1];
rz(-2.5109992) q[1];
sx q[1];
rz(0.15288606) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4893462) q[0];
sx q[0];
rz(-1.3382599) q[0];
sx q[0];
rz(1.812797) q[0];
rz(-pi) q[1];
rz(-1.6727078) q[2];
sx q[2];
rz(-2.1481124) q[2];
sx q[2];
rz(0.99330639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5093435) q[1];
sx q[1];
rz(-1.1763089) q[1];
sx q[1];
rz(-3.0837653) q[1];
x q[2];
rz(1.6083999) q[3];
sx q[3];
rz(-2.561224) q[3];
sx q[3];
rz(2.6036724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9360518) q[2];
sx q[2];
rz(-2.3369393) q[2];
sx q[2];
rz(-1.12977) q[2];
rz(-1.0063082) q[3];
sx q[3];
rz(-1.3200503) q[3];
sx q[3];
rz(1.6442851) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048700843) q[0];
sx q[0];
rz(-1.0660271) q[0];
sx q[0];
rz(-0.14376465) q[0];
rz(-2.9394506) q[1];
sx q[1];
rz(-0.527924) q[1];
sx q[1];
rz(-1.082083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4253628) q[0];
sx q[0];
rz(-0.83304616) q[0];
sx q[0];
rz(-2.0701755) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4074102) q[2];
sx q[2];
rz(-0.95962822) q[2];
sx q[2];
rz(-2.2847459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80716996) q[1];
sx q[1];
rz(-1.7826986) q[1];
sx q[1];
rz(2.2905877) q[1];
x q[2];
rz(0.73900004) q[3];
sx q[3];
rz(-2.0428052) q[3];
sx q[3];
rz(1.3907741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9988592) q[2];
sx q[2];
rz(-1.961668) q[2];
sx q[2];
rz(0.72706968) q[2];
rz(-1.3149698) q[3];
sx q[3];
rz(-1.246779) q[3];
sx q[3];
rz(1.6453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406782) q[0];
sx q[0];
rz(-1.0777363) q[0];
sx q[0];
rz(1.9986073) q[0];
rz(-2.9179528) q[1];
sx q[1];
rz(-2.5032005) q[1];
sx q[1];
rz(-1.2248096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1849104) q[0];
sx q[0];
rz(-1.3089797) q[0];
sx q[0];
rz(-1.8641406) q[0];
rz(-1.9838422) q[2];
sx q[2];
rz(-2.2487679) q[2];
sx q[2];
rz(-1.7630446) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8401007) q[1];
sx q[1];
rz(-2.9291541) q[1];
sx q[1];
rz(1.2919687) q[1];
rz(0.16815987) q[3];
sx q[3];
rz(-0.56795299) q[3];
sx q[3];
rz(1.1105383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5582922) q[2];
sx q[2];
rz(-1.6879098) q[2];
sx q[2];
rz(2.3395786) q[2];
rz(-1.2174886) q[3];
sx q[3];
rz(-3.0018482) q[3];
sx q[3];
rz(-1.2453992) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8121346) q[0];
sx q[0];
rz(-1.7725002) q[0];
sx q[0];
rz(-1.5611956) q[0];
rz(2.2691057) q[1];
sx q[1];
rz(-0.28631887) q[1];
sx q[1];
rz(0.40503851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0035343) q[0];
sx q[0];
rz(-1.0558319) q[0];
sx q[0];
rz(2.1206417) q[0];
x q[1];
rz(-0.15200673) q[2];
sx q[2];
rz(-0.59796732) q[2];
sx q[2];
rz(0.94546181) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2916221) q[1];
sx q[1];
rz(-2.1287781) q[1];
sx q[1];
rz(2.4172952) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5902577) q[3];
sx q[3];
rz(-1.7177337) q[3];
sx q[3];
rz(-0.16836873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8672455) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(0.61152968) q[2];
rz(1.3036171) q[3];
sx q[3];
rz(-0.36208624) q[3];
sx q[3];
rz(-1.1360137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63854727) q[0];
sx q[0];
rz(-1.847581) q[0];
sx q[0];
rz(0.053330388) q[0];
rz(2.5697925) q[1];
sx q[1];
rz(-1.5170393) q[1];
sx q[1];
rz(-1.2791876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8227575) q[0];
sx q[0];
rz(-1.3190206) q[0];
sx q[0];
rz(-0.29522268) q[0];
x q[1];
rz(1.8498184) q[2];
sx q[2];
rz(-0.49683647) q[2];
sx q[2];
rz(-1.9328839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0030583) q[1];
sx q[1];
rz(-1.6489677) q[1];
sx q[1];
rz(-0.98346904) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9235696) q[3];
sx q[3];
rz(-1.9338528) q[3];
sx q[3];
rz(-1.3040964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46644396) q[2];
sx q[2];
rz(-1.5326591) q[2];
sx q[2];
rz(2.4565728) q[2];
rz(-0.78392616) q[3];
sx q[3];
rz(-2.7214366) q[3];
sx q[3];
rz(0.73721957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6017799) q[0];
sx q[0];
rz(-1.89986) q[0];
sx q[0];
rz(1.4358406) q[0];
rz(-0.80541366) q[1];
sx q[1];
rz(-0.46331159) q[1];
sx q[1];
rz(0.70809271) q[1];
rz(-2.5067301) q[2];
sx q[2];
rz(-2.786953) q[2];
sx q[2];
rz(2.0155596) q[2];
rz(-2.7266034) q[3];
sx q[3];
rz(-1.4616953) q[3];
sx q[3];
rz(-2.2898522) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
