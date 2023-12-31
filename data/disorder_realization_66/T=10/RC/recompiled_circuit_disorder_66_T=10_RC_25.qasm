OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(-1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(0.02286214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1296596) q[0];
sx q[0];
rz(-1.7190022) q[0];
sx q[0];
rz(-0.11797842) q[0];
rz(-pi) q[1];
rz(-0.51486751) q[2];
sx q[2];
rz(-1.9549184) q[2];
sx q[2];
rz(1.9799973) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80427985) q[1];
sx q[1];
rz(-2.7645281) q[1];
sx q[1];
rz(1.9878597) q[1];
rz(1.2223577) q[3];
sx q[3];
rz(-2.493353) q[3];
sx q[3];
rz(-3.0188308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239685) q[0];
sx q[0];
rz(-0.94857615) q[0];
sx q[0];
rz(-1.5505962) q[0];
rz(-pi) q[1];
rz(-2.000196) q[2];
sx q[2];
rz(-2.4944802) q[2];
sx q[2];
rz(0.61691689) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47463402) q[1];
sx q[1];
rz(-0.14161319) q[1];
sx q[1];
rz(2.2255564) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0410791) q[3];
sx q[3];
rz(-1.2840052) q[3];
sx q[3];
rz(1.6581397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.4651728) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-0.9056257) q[0];
rz(0.33272818) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4577643) q[0];
sx q[0];
rz(-1.0965075) q[0];
sx q[0];
rz(1.477924) q[0];
x q[1];
rz(-2.076782) q[2];
sx q[2];
rz(-1.2227321) q[2];
sx q[2];
rz(0.51007523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3789931) q[1];
sx q[1];
rz(-1.5931207) q[1];
sx q[1];
rz(-2.3861045) q[1];
x q[2];
rz(2.3178187) q[3];
sx q[3];
rz(-2.4274821) q[3];
sx q[3];
rz(0.87574524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5642693) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(1.3576515) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3992594) q[2];
sx q[2];
rz(-2.9651387) q[2];
sx q[2];
rz(1.2982969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9878848) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(0.99516408) q[1];
rz(0.025709318) q[3];
sx q[3];
rz(-2.6298012) q[3];
sx q[3];
rz(-0.6230841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(0.80835289) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-2.6079544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69913188) q[0];
sx q[0];
rz(-1.2958382) q[0];
sx q[0];
rz(-1.86637) q[0];
rz(-1.3445504) q[2];
sx q[2];
rz(-0.81586736) q[2];
sx q[2];
rz(-2.0362542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4644949) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(2.5527843) q[1];
rz(-pi) q[2];
rz(0.20547262) q[3];
sx q[3];
rz(-0.95147248) q[3];
sx q[3];
rz(1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9020033) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.9653962) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-1.0046545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83197901) q[0];
sx q[0];
rz(-1.9404611) q[0];
sx q[0];
rz(1.9395104) q[0];
rz(-pi) q[1];
rz(-0.3968233) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(2.4900988) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9690296) q[1];
sx q[1];
rz(-1.8168212) q[1];
sx q[1];
rz(2.5804156) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1296022) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(0.12366611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(2.4694494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37049609) q[0];
sx q[0];
rz(-1.5633977) q[0];
sx q[0];
rz(1.4047983) q[0];
rz(-pi) q[1];
rz(2.0124112) q[2];
sx q[2];
rz(-1.9600944) q[2];
sx q[2];
rz(2.0814975) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0828515) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(-0.65545603) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9716306) q[3];
sx q[3];
rz(-2.4331577) q[3];
sx q[3];
rz(-1.8426614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(-0.82459015) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(0.80668443) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8950302) q[0];
sx q[0];
rz(-0.44996214) q[0];
sx q[0];
rz(1.8438086) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(0.92668698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1572691) q[1];
sx q[1];
rz(-1.2829797) q[1];
sx q[1];
rz(2.4411574) q[1];
x q[2];
rz(1.8716378) q[3];
sx q[3];
rz(-1.6013718) q[3];
sx q[3];
rz(-1.2310864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-0.34067571) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239615) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(0.4171564) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9195243) q[2];
sx q[2];
rz(-2.1246315) q[2];
sx q[2];
rz(-2.1332707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3714911) q[1];
sx q[1];
rz(-2.1226468) q[1];
sx q[1];
rz(-0.21367195) q[1];
rz(-0.70298985) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(0.27858946) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(-3.0648807) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.840832) q[0];
sx q[0];
rz(-1.6941841) q[0];
sx q[0];
rz(-1.8873909) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1439728) q[2];
sx q[2];
rz(-1.7152889) q[2];
sx q[2];
rz(2.0720553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(1.5484023) q[1];
rz(-pi) q[2];
rz(-1.4823227) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(-0.74151553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(2.0942681) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(-0.69910819) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
