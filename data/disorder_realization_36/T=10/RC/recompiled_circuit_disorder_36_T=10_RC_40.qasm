OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(2.5640092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89864697) q[0];
sx q[0];
rz(-0.67910128) q[0];
sx q[0];
rz(2.8773017) q[0];
rz(-pi) q[1];
rz(2.6063927) q[2];
sx q[2];
rz(-1.1241962) q[2];
sx q[2];
rz(-1.610178) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5756302) q[1];
sx q[1];
rz(-1.9735676) q[1];
sx q[1];
rz(-2.8494542) q[1];
rz(-pi) q[2];
rz(-0.64896119) q[3];
sx q[3];
rz(-1.0421841) q[3];
sx q[3];
rz(1.1543857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.2600391) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-0.84567436) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0613522) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(2.4387226) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77180736) q[2];
sx q[2];
rz(-2.9559921) q[2];
sx q[2];
rz(1.1112569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37296346) q[1];
sx q[1];
rz(-2.0634723) q[1];
sx q[1];
rz(1.3586587) q[1];
x q[2];
rz(-2.1917079) q[3];
sx q[3];
rz(-1.9904899) q[3];
sx q[3];
rz(1.8267531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(1.9539179) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(2.0551596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.438293) q[0];
sx q[0];
rz(-0.46143954) q[0];
sx q[0];
rz(-1.5745387) q[0];
rz(-pi) q[1];
rz(-3.0263607) q[2];
sx q[2];
rz(-0.92667246) q[2];
sx q[2];
rz(-2.6098721) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.31049) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(2.8527841) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1826035) q[3];
sx q[3];
rz(-0.43366323) q[3];
sx q[3];
rz(0.38503034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.4455459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1235385) q[0];
sx q[0];
rz(-0.92045438) q[0];
sx q[0];
rz(-1.3440078) q[0];
rz(-pi) q[1];
rz(-1.4807329) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(2.4243674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7062263) q[1];
sx q[1];
rz(-1.5063783) q[1];
sx q[1];
rz(0.21965841) q[1];
x q[2];
rz(1.0936071) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(1.5447865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(0.23322341) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695278) q[0];
sx q[0];
rz(-2.2846691) q[0];
sx q[0];
rz(0.010849997) q[0];
x q[1];
rz(-2.6536077) q[2];
sx q[2];
rz(-0.93163604) q[2];
sx q[2];
rz(-2.0069063) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2367868) q[1];
sx q[1];
rz(-2.2854837) q[1];
sx q[1];
rz(-2.3805815) q[1];
rz(-pi) q[2];
rz(2.3417926) q[3];
sx q[3];
rz(-2.4197257) q[3];
sx q[3];
rz(-1.1218027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(0.12621005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.431625) q[0];
sx q[0];
rz(-0.21490782) q[0];
sx q[0];
rz(2.9931195) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38527617) q[2];
sx q[2];
rz(-0.81309536) q[2];
sx q[2];
rz(-2.3713881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63657657) q[1];
sx q[1];
rz(-2.2762183) q[1];
sx q[1];
rz(0.24800639) q[1];
x q[2];
rz(-0.19212171) q[3];
sx q[3];
rz(-1.8936833) q[3];
sx q[3];
rz(2.9411112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90298992) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59259748) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(2.5820406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768809) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(-2.5401831) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5188811) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.7920997) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2110062) q[1];
sx q[1];
rz(-1.483327) q[1];
sx q[1];
rz(-1.4409815) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5114602) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(2.7255448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-2.7977978) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(-2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(-2.2699845) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14411892) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(0.075450443) q[0];
rz(2.7360104) q[2];
sx q[2];
rz(-1.7103346) q[2];
sx q[2];
rz(-1.1428733) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(-2.4396067) q[1];
rz(0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(-1.0761716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9386439) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(1.5244012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071951889) q[0];
sx q[0];
rz(-1.4884293) q[0];
sx q[0];
rz(-3.1125493) q[0];
rz(1.379307) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(0.42275235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7414654) q[1];
sx q[1];
rz(-1.8997846) q[1];
sx q[1];
rz(2.9743183) q[1];
x q[2];
rz(-2.5025326) q[3];
sx q[3];
rz(-1.3165054) q[3];
sx q[3];
rz(-1.1921079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(0.25340733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2449269) q[0];
sx q[0];
rz(-1.1765119) q[0];
sx q[0];
rz(-0.11163296) q[0];
rz(-2.017574) q[2];
sx q[2];
rz(-0.3728711) q[2];
sx q[2];
rz(0.98137059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39292654) q[1];
sx q[1];
rz(-1.3303419) q[1];
sx q[1];
rz(2.7886831) q[1];
rz(0.57755034) q[3];
sx q[3];
rz(-0.69637075) q[3];
sx q[3];
rz(-1.7601354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(0.5029451) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-1.8297557) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(-2.4651299) q[3];
sx q[3];
rz(-0.87920311) q[3];
sx q[3];
rz(1.9999947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];