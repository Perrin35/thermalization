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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(-0.21633202) q[0];
rz(0.95739111) q[1];
sx q[1];
rz(-2.4654145) q[1];
sx q[1];
rz(1.5114991) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0656242) q[0];
sx q[0];
rz(-2.127021) q[0];
sx q[0];
rz(-1.5821304) q[0];
rz(-pi) q[1];
rz(0.47343238) q[2];
sx q[2];
rz(-1.2780927) q[2];
sx q[2];
rz(-2.6930489) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18876702) q[1];
sx q[1];
rz(-2.5270542) q[1];
sx q[1];
rz(-0.60093083) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0262296) q[3];
sx q[3];
rz(-2.3854227) q[3];
sx q[3];
rz(0.84634534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.034885255) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(-2.0931639) q[2];
rz(-0.18167051) q[3];
sx q[3];
rz(-2.1766267) q[3];
sx q[3];
rz(-2.488193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6492017) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(2.6948068) q[0];
rz(-0.88042879) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(0.78278881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7715817) q[0];
sx q[0];
rz(-0.25213045) q[0];
sx q[0];
rz(1.3229516) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61666416) q[2];
sx q[2];
rz(-1.9266124) q[2];
sx q[2];
rz(-0.43105506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.782541) q[1];
sx q[1];
rz(-0.56582574) q[1];
sx q[1];
rz(-0.58360175) q[1];
x q[2];
rz(0.47329013) q[3];
sx q[3];
rz(-0.9447228) q[3];
sx q[3];
rz(1.2155217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.11144) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(0.97935575) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(-0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591831) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(0.78980494) q[0];
rz(2.9549331) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(-2.1479215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45897608) q[0];
sx q[0];
rz(-1.0280711) q[0];
sx q[0];
rz(-1.2786464) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8254059) q[2];
sx q[2];
rz(-0.63717604) q[2];
sx q[2];
rz(-0.030046163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6911654) q[1];
sx q[1];
rz(-1.1466999) q[1];
sx q[1];
rz(1.0662088) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8422115) q[3];
sx q[3];
rz(-2.2732353) q[3];
sx q[3];
rz(-2.0347119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8197202) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(-0.37128386) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6868941) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.7373079) q[0];
rz(-2.4644201) q[1];
sx q[1];
rz(-1.9858457) q[1];
sx q[1];
rz(-1.4926532) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20805222) q[0];
sx q[0];
rz(-1.786199) q[0];
sx q[0];
rz(-2.893704) q[0];
rz(-pi) q[1];
rz(2.1588232) q[2];
sx q[2];
rz(-0.50543907) q[2];
sx q[2];
rz(0.44014058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4863534) q[1];
sx q[1];
rz(-2.7905373) q[1];
sx q[1];
rz(-2.4454861) q[1];
rz(2.2954313) q[3];
sx q[3];
rz(-1.807888) q[3];
sx q[3];
rz(2.2711636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45381418) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(0.34823927) q[2];
rz(-1.6866775) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(-1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652047) q[0];
sx q[0];
rz(-0.56469733) q[0];
sx q[0];
rz(2.2494466) q[0];
rz(-2.6773894) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(1.8468599) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034445914) q[0];
sx q[0];
rz(-2.4305516) q[0];
sx q[0];
rz(-0.76816316) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6991529) q[2];
sx q[2];
rz(-2.1048628) q[2];
sx q[2];
rz(0.72557025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4720104) q[1];
sx q[1];
rz(-1.8947487) q[1];
sx q[1];
rz(-1.1060017) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2512359) q[3];
sx q[3];
rz(-1.0728175) q[3];
sx q[3];
rz(-2.339956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(2.5757705) q[2];
rz(0.081341751) q[3];
sx q[3];
rz(-0.9809202) q[3];
sx q[3];
rz(-2.5210023) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.464798) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(-1.5860522) q[0];
rz(-1.0606891) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(0.63180822) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2151129) q[0];
sx q[0];
rz(-0.2479015) q[0];
sx q[0];
rz(0.32899022) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6423694) q[2];
sx q[2];
rz(-2.2167335) q[2];
sx q[2];
rz(1.2690085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7673847) q[1];
sx q[1];
rz(-0.65113089) q[1];
sx q[1];
rz(-0.99094772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9618044) q[3];
sx q[3];
rz(-0.69151141) q[3];
sx q[3];
rz(-2.9464242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98171988) q[2];
sx q[2];
rz(-2.0231415) q[2];
sx q[2];
rz(-2.6427606) q[2];
rz(1.3129129) q[3];
sx q[3];
rz(-0.73176089) q[3];
sx q[3];
rz(-1.7074728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53443921) q[0];
sx q[0];
rz(-2.7767015) q[0];
sx q[0];
rz(0.69951192) q[0];
rz(-2.6761159) q[1];
sx q[1];
rz(-2.270348) q[1];
sx q[1];
rz(-2.2043998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584211) q[0];
sx q[0];
rz(-1.7390378) q[0];
sx q[0];
rz(-2.9267163) q[0];
rz(-pi) q[1];
rz(2.9019321) q[2];
sx q[2];
rz(-1.4693038) q[2];
sx q[2];
rz(3.0200151) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4129137) q[1];
sx q[1];
rz(-1.2866486) q[1];
sx q[1];
rz(-2.9212664) q[1];
x q[2];
rz(2.3714957) q[3];
sx q[3];
rz(-2.7316964) q[3];
sx q[3];
rz(2.485825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.297544) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(-2.76827) q[2];
rz(-1.9518055) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(-0.51026195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74476403) q[0];
sx q[0];
rz(-2.0592392) q[0];
sx q[0];
rz(-0.74791351) q[0];
rz(2.3751936) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(0.0029729923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52067962) q[0];
sx q[0];
rz(-1.7629884) q[0];
sx q[0];
rz(-0.8780794) q[0];
rz(1.4998798) q[2];
sx q[2];
rz(-1.5902963) q[2];
sx q[2];
rz(-1.8766581) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6312478) q[1];
sx q[1];
rz(-1.5273792) q[1];
sx q[1];
rz(-1.4166635) q[1];
rz(-1.7693172) q[3];
sx q[3];
rz(-1.3796419) q[3];
sx q[3];
rz(-1.6894345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.030674) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(1.0670916) q[2];
rz(0.083960697) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(-0.77897227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(3.0352266) q[0];
rz(0.97995177) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(0.76593691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04127114) q[0];
sx q[0];
rz(-1.0491956) q[0];
sx q[0];
rz(-2.0727511) q[0];
rz(-pi) q[1];
rz(1.6607051) q[2];
sx q[2];
rz(-1.769763) q[2];
sx q[2];
rz(-1.0294017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3801743) q[1];
sx q[1];
rz(-1.6055709) q[1];
sx q[1];
rz(-0.67616391) q[1];
rz(0.014102293) q[3];
sx q[3];
rz(-1.5670793) q[3];
sx q[3];
rz(0.93059082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6672259) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(0.70029798) q[2];
rz(2.4750366) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(-2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67548442) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(-1.745537) q[0];
rz(-1.4736157) q[1];
sx q[1];
rz(-1.2584078) q[1];
sx q[1];
rz(-2.4748763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6121126) q[0];
sx q[0];
rz(-2.1938938) q[0];
sx q[0];
rz(0.77147958) q[0];
rz(0.36138968) q[2];
sx q[2];
rz(-0.69212428) q[2];
sx q[2];
rz(3.1049984) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.279505) q[1];
sx q[1];
rz(-0.34395978) q[1];
sx q[1];
rz(2.7441447) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75717775) q[3];
sx q[3];
rz(-0.72523967) q[3];
sx q[3];
rz(0.66981572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0746158) q[2];
sx q[2];
rz(-2.4846027) q[2];
sx q[2];
rz(-0.89078772) q[2];
rz(-0.43205076) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.73410949) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(1.2474077) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(-1.9966077) q[2];
sx q[2];
rz(-2.7663284) q[2];
sx q[2];
rz(-0.13079499) q[2];
rz(1.6040989) q[3];
sx q[3];
rz(-2.0596444) q[3];
sx q[3];
rz(2.8465908) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
