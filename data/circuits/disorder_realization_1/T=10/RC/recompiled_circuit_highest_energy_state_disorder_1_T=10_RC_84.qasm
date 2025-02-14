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
rz(1.7684608) q[0];
sx q[0];
rz(-1.7113577) q[0];
sx q[0];
rz(0.89682427) q[0];
rz(-1.7984017) q[1];
sx q[1];
rz(3.922037) q[1];
sx q[1];
rz(8.2889397) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97082393) q[0];
sx q[0];
rz(-1.5005732) q[0];
sx q[0];
rz(-1.4115566) q[0];
rz(0.48223038) q[2];
sx q[2];
rz(-0.94677529) q[2];
sx q[2];
rz(-2.4718747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.070277409) q[1];
sx q[1];
rz(-2.1143374) q[1];
sx q[1];
rz(0.75884968) q[1];
rz(-0.43026409) q[3];
sx q[3];
rz(-1.2183905) q[3];
sx q[3];
rz(-0.5510181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4844126) q[2];
sx q[2];
rz(-1.5766532) q[2];
sx q[2];
rz(0.29259345) q[2];
rz(1.6491133) q[3];
sx q[3];
rz(-1.7238659) q[3];
sx q[3];
rz(-2.9610146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0507616) q[0];
sx q[0];
rz(-1.2240336) q[0];
sx q[0];
rz(-2.5378788) q[0];
rz(-1.0366038) q[1];
sx q[1];
rz(-0.34244582) q[1];
sx q[1];
rz(-2.2110914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11545746) q[0];
sx q[0];
rz(-1.8107426) q[0];
sx q[0];
rz(-1.7895997) q[0];
x q[1];
rz(-1.7199688) q[2];
sx q[2];
rz(-2.1270555) q[2];
sx q[2];
rz(0.38146293) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82236921) q[1];
sx q[1];
rz(-1.7041053) q[1];
sx q[1];
rz(-2.4001166) q[1];
rz(-pi) q[2];
rz(2.1648079) q[3];
sx q[3];
rz(-1.8476764) q[3];
sx q[3];
rz(-2.8897485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36640627) q[2];
sx q[2];
rz(-1.2440871) q[2];
sx q[2];
rz(-1.6413956) q[2];
rz(-0.90534798) q[3];
sx q[3];
rz(-1.2572181) q[3];
sx q[3];
rz(2.1998028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6833078) q[0];
sx q[0];
rz(-1.4508805) q[0];
sx q[0];
rz(-0.84253755) q[0];
rz(1.9375577) q[1];
sx q[1];
rz(-1.2843818) q[1];
sx q[1];
rz(-1.1010928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6146832) q[0];
sx q[0];
rz(-1.2430265) q[0];
sx q[0];
rz(-3.0614168) q[0];
x q[1];
rz(-0.74275334) q[2];
sx q[2];
rz(-1.4223546) q[2];
sx q[2];
rz(-0.63723931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3345015) q[1];
sx q[1];
rz(-2.0223631) q[1];
sx q[1];
rz(1.0784763) q[1];
x q[2];
rz(2.1258611) q[3];
sx q[3];
rz(-0.69472488) q[3];
sx q[3];
rz(0.55391698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6847685) q[2];
sx q[2];
rz(-0.55906877) q[2];
sx q[2];
rz(-0.10747257) q[2];
rz(-1.0476073) q[3];
sx q[3];
rz(-1.8045629) q[3];
sx q[3];
rz(-1.5260772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90077129) q[0];
sx q[0];
rz(-2.695684) q[0];
sx q[0];
rz(1.4500424) q[0];
rz(-2.5598473) q[1];
sx q[1];
rz(-2.3090239) q[1];
sx q[1];
rz(2.0411101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9387526) q[0];
sx q[0];
rz(-1.0577518) q[0];
sx q[0];
rz(2.0111397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1911439) q[2];
sx q[2];
rz(-2.0850402) q[2];
sx q[2];
rz(2.4730027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2128513) q[1];
sx q[1];
rz(-0.87777661) q[1];
sx q[1];
rz(-1.568896) q[1];
rz(-2.8488345) q[3];
sx q[3];
rz(-2.2653711) q[3];
sx q[3];
rz(0.048489038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45740286) q[2];
sx q[2];
rz(-1.0418912) q[2];
sx q[2];
rz(0.602496) q[2];
rz(-1.6920793) q[3];
sx q[3];
rz(-1.4256198) q[3];
sx q[3];
rz(0.95363936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36798894) q[0];
sx q[0];
rz(-1.8743176) q[0];
sx q[0];
rz(1.3195272) q[0];
rz(-2.6124182) q[1];
sx q[1];
rz(-1.2472943) q[1];
sx q[1];
rz(-0.41437638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4236365) q[0];
sx q[0];
rz(-2.0265105) q[0];
sx q[0];
rz(0.44063963) q[0];
rz(0.82444382) q[2];
sx q[2];
rz(-1.0811812) q[2];
sx q[2];
rz(2.010422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16769174) q[1];
sx q[1];
rz(-2.0030336) q[1];
sx q[1];
rz(-2.8532989) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46196725) q[3];
sx q[3];
rz(-2.7865692) q[3];
sx q[3];
rz(-0.17661653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6759912) q[2];
sx q[2];
rz(-0.71054116) q[2];
sx q[2];
rz(-2.4856534) q[2];
rz(-1.9115619) q[3];
sx q[3];
rz(-2.0345119) q[3];
sx q[3];
rz(2.6225577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8116654) q[0];
sx q[0];
rz(-1.0141319) q[0];
sx q[0];
rz(2.2299679) q[0];
rz(1.3765593) q[1];
sx q[1];
rz(-0.39386097) q[1];
sx q[1];
rz(-2.6011661) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2776162) q[0];
sx q[0];
rz(-1.2860398) q[0];
sx q[0];
rz(-2.6984614) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0824605) q[2];
sx q[2];
rz(-2.793047) q[2];
sx q[2];
rz(-0.71427554) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47460654) q[1];
sx q[1];
rz(-1.1724533) q[1];
sx q[1];
rz(1.5035267) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94169887) q[3];
sx q[3];
rz(-1.3543918) q[3];
sx q[3];
rz(3.0890478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2889169) q[2];
sx q[2];
rz(-1.6923075) q[2];
sx q[2];
rz(1.1360315) q[2];
rz(-1.8063258) q[3];
sx q[3];
rz(-1.5583594) q[3];
sx q[3];
rz(0.11245888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5920608) q[0];
sx q[0];
rz(-2.9887152) q[0];
sx q[0];
rz(-2.1524647) q[0];
rz(-1.2024744) q[1];
sx q[1];
rz(-1.6705284) q[1];
sx q[1];
rz(2.8340526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06194845) q[0];
sx q[0];
rz(-0.35828081) q[0];
sx q[0];
rz(2.4419191) q[0];
x q[1];
rz(1.9498242) q[2];
sx q[2];
rz(-2.267365) q[2];
sx q[2];
rz(-1.3427757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9779098) q[1];
sx q[1];
rz(-1.3505757) q[1];
sx q[1];
rz(-1.5942595) q[1];
x q[2];
rz(-0.6749344) q[3];
sx q[3];
rz(-2.4111742) q[3];
sx q[3];
rz(2.6715913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5146553) q[2];
sx q[2];
rz(-0.81521002) q[2];
sx q[2];
rz(1.5283654) q[2];
rz(2.3441687) q[3];
sx q[3];
rz(-1.5484836) q[3];
sx q[3];
rz(3.048866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1180457) q[0];
sx q[0];
rz(-0.67414701) q[0];
sx q[0];
rz(2.4349037) q[0];
rz(0.30330172) q[1];
sx q[1];
rz(-1.773833) q[1];
sx q[1];
rz(1.2695262) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6578419) q[0];
sx q[0];
rz(-2.7829956) q[0];
sx q[0];
rz(2.7593385) q[0];
rz(-pi) q[1];
rz(-2.6466188) q[2];
sx q[2];
rz(-1.3522864) q[2];
sx q[2];
rz(1.72118) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4321338) q[1];
sx q[1];
rz(-1.7466326) q[1];
sx q[1];
rz(-1.9587112) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4702099) q[3];
sx q[3];
rz(-1.0142418) q[3];
sx q[3];
rz(0.33410698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0728486) q[2];
sx q[2];
rz(-0.35126433) q[2];
sx q[2];
rz(0.52034411) q[2];
rz(-1.2460234) q[3];
sx q[3];
rz(-1.7249853) q[3];
sx q[3];
rz(2.0791159) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.635427) q[0];
sx q[0];
rz(-1.6599673) q[0];
sx q[0];
rz(0.51314276) q[0];
rz(-2.6601833) q[1];
sx q[1];
rz(-2.3115999) q[1];
sx q[1];
rz(0.76024461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8177196) q[0];
sx q[0];
rz(-1.5580252) q[0];
sx q[0];
rz(-1.5991999) q[0];
rz(-0.77519007) q[2];
sx q[2];
rz(-1.7235867) q[2];
sx q[2];
rz(0.82681954) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.34955353) q[1];
sx q[1];
rz(-2.4124618) q[1];
sx q[1];
rz(3.0302073) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.097858345) q[3];
sx q[3];
rz(-1.9966085) q[3];
sx q[3];
rz(1.2657566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61797872) q[2];
sx q[2];
rz(-0.60680497) q[2];
sx q[2];
rz(2.4904909) q[2];
rz(-1.2111604) q[3];
sx q[3];
rz(-1.9821207) q[3];
sx q[3];
rz(-2.7291362) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3245658) q[0];
sx q[0];
rz(-2.4232219) q[0];
sx q[0];
rz(0.73085648) q[0];
rz(0.51365799) q[1];
sx q[1];
rz(-2.3274603) q[1];
sx q[1];
rz(-0.51630744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8532787) q[0];
sx q[0];
rz(-2.1876039) q[0];
sx q[0];
rz(-0.46648394) q[0];
rz(2.7514156) q[2];
sx q[2];
rz(-1.220229) q[2];
sx q[2];
rz(-1.0165018) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67077209) q[1];
sx q[1];
rz(-1.9665475) q[1];
sx q[1];
rz(2.938063) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4993419) q[3];
sx q[3];
rz(-2.7786479) q[3];
sx q[3];
rz(-1.7098984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98661304) q[2];
sx q[2];
rz(-1.2519138) q[2];
sx q[2];
rz(-0.66407859) q[2];
rz(-2.1443071) q[3];
sx q[3];
rz(-1.7136796) q[3];
sx q[3];
rz(2.6663781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8127139) q[0];
sx q[0];
rz(-1.7408149) q[0];
sx q[0];
rz(-0.029205532) q[0];
rz(3.008814) q[1];
sx q[1];
rz(-1.7280424) q[1];
sx q[1];
rz(-0.8676563) q[1];
rz(0.50162195) q[2];
sx q[2];
rz(-1.2615962) q[2];
sx q[2];
rz(-0.35485219) q[2];
rz(-2.0302782) q[3];
sx q[3];
rz(-1.9412771) q[3];
sx q[3];
rz(-2.8893445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
