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
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4532783) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(-2.2384089) q[0];
rz(2.0055662) q[2];
sx q[2];
rz(-1.0966986) q[2];
sx q[2];
rz(-0.20027645) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3373128) q[1];
sx q[1];
rz(-0.37706456) q[1];
sx q[1];
rz(1.1537329) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25306563) q[3];
sx q[3];
rz(-0.96732891) q[3];
sx q[3];
rz(-0.30482182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-0.76618761) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413977) q[0];
sx q[0];
rz(-1.5543823) q[0];
sx q[0];
rz(-2.5192758) q[0];
rz(-pi) q[1];
rz(-2.8367963) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(0.09588974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18499204) q[1];
sx q[1];
rz(-1.6829832) q[1];
sx q[1];
rz(3.0549906) q[1];
rz(-1.1005136) q[3];
sx q[3];
rz(-1.2840052) q[3];
sx q[3];
rz(-1.4834529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860486) q[0];
sx q[0];
rz(-1.6533924) q[0];
sx q[0];
rz(2.6655469) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2134174) q[2];
sx q[2];
rz(-0.60545063) q[2];
sx q[2];
rz(-1.6124992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76259957) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(-0.75548817) q[1];
rz(-pi) q[2];
rz(0.53220766) q[3];
sx q[3];
rz(-1.0695219) q[3];
sx q[3];
rz(-1.3794848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(-1.3282233) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(-0.20733325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610791) q[0];
sx q[0];
rz(-1.5797537) q[0];
sx q[0];
rz(1.5293967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7423332) q[2];
sx q[2];
rz(-2.9651387) q[2];
sx q[2];
rz(1.2982969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85992766) q[1];
sx q[1];
rz(-1.4178226) q[1];
sx q[1];
rz(-1.8106736) q[1];
rz(1.585235) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(-0.80777848) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-2.6079544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14347178) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(-0.801416) q[0];
x q[1];
rz(2.9075378) q[2];
sx q[2];
rz(-2.3599527) q[2];
sx q[2];
rz(-0.78125886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.071760885) q[1];
sx q[1];
rz(-1.3511718) q[1];
sx q[1];
rz(-1.4211618) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2002033) q[3];
sx q[3];
rz(-1.7377059) q[3];
sx q[3];
rz(-3.133579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(0.56838244) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83197901) q[0];
sx q[0];
rz(-1.9404611) q[0];
sx q[0];
rz(-1.9395104) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9506358) q[2];
sx q[2];
rz(-2.7381884) q[2];
sx q[2];
rz(1.0952589) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.028588258) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(0.44087704) q[1];
x q[2];
rz(-3.1296022) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(-0.12366611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(2.693434) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1561688) q[0];
sx q[0];
rz(-2.9754313) q[0];
sx q[0];
rz(-1.5260494) q[0];
rz(-pi) q[1];
rz(-2.7156098) q[2];
sx q[2];
rz(-1.1642712) q[2];
sx q[2];
rz(0.68824088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0587412) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(2.4861366) q[1];
rz(1.4268731) q[3];
sx q[3];
rz(-2.2669499) q[3];
sx q[3];
rz(1.5211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(-0.82459015) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(0.80668443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8950302) q[0];
sx q[0];
rz(-2.6916305) q[0];
sx q[0];
rz(-1.297784) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47862349) q[2];
sx q[2];
rz(-2.9106986) q[2];
sx q[2];
rz(-0.17639562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9843236) q[1];
sx q[1];
rz(-1.2829797) q[1];
sx q[1];
rz(2.4411574) q[1];
rz(-pi) q[2];
x q[2];
rz(0.032012352) q[3];
sx q[3];
rz(-1.2700998) q[3];
sx q[3];
rz(-2.7923982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(2.1419443) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-0.28655562) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(2.8009169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4083456) q[0];
sx q[0];
rz(-2.354963) q[0];
sx q[0];
rz(2.0287081) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9195243) q[2];
sx q[2];
rz(-2.1246315) q[2];
sx q[2];
rz(1.008322) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3714911) q[1];
sx q[1];
rz(-2.1226468) q[1];
sx q[1];
rz(-0.21367195) q[1];
rz(-pi) q[2];
rz(-1.7394003) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-0.51914674) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72934812) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4519276) q[0];
sx q[0];
rz(-1.8849012) q[0];
sx q[0];
rz(3.0118224) q[0];
x q[1];
rz(-2.9701091) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-2.7330074) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0753239) q[1];
sx q[1];
rz(-1.2053718) q[1];
sx q[1];
rz(-1.5484023) q[1];
rz(-pi) q[2];
rz(0.17681392) q[3];
sx q[3];
rz(-1.483695) q[3];
sx q[3];
rz(0.84482312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(0.58615276) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(-1.0473245) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(2.4424845) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
