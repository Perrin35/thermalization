OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(0.56086993) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(1.9265494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62984798) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(-1.758979) q[0];
x q[1];
rz(-1.2826074) q[2];
sx q[2];
rz(-0.93027861) q[2];
sx q[2];
rz(-3.1079907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3455968) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(-2.7624346) q[1];
x q[2];
rz(-1.3052985) q[3];
sx q[3];
rz(-1.7146829) q[3];
sx q[3];
rz(3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-0.077117292) q[0];
rz(-2.8027957) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(1.5391301) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9248283) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(-3.0490962) q[0];
rz(1.7878754) q[2];
sx q[2];
rz(-0.84393822) q[2];
sx q[2];
rz(2.196687) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6290366) q[1];
sx q[1];
rz(-1.7450383) q[1];
sx q[1];
rz(-2.0577355) q[1];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.9841521) q[3];
sx q[3];
rz(1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(-2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-2.7084896) q[0];
rz(-1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(2.5862397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(0.81911266) q[0];
rz(2.136134) q[2];
sx q[2];
rz(-1.3034504) q[2];
sx q[2];
rz(-0.29758673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8086116) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.3036742) q[1];
rz(1.2461353) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(-3.0666921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(-1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(-0.25517685) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53084757) q[0];
sx q[0];
rz(-1.2753715) q[0];
sx q[0];
rz(2.7644964) q[0];
x q[1];
rz(-1.9011263) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(1.4532879) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0604917) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(2.9213419) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.051024036) q[3];
sx q[3];
rz(-1.050356) q[3];
sx q[3];
rz(0.40363064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1306886) q[0];
sx q[0];
rz(-0.70106693) q[0];
sx q[0];
rz(-2.5581193) q[0];
x q[1];
rz(2.0300843) q[2];
sx q[2];
rz(-1.2570724) q[2];
sx q[2];
rz(-2.8788061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1069378) q[1];
sx q[1];
rz(-1.5562623) q[1];
sx q[1];
rz(-2.8388139) q[1];
x q[2];
rz(1.7124743) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(0.39839881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1405979) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(0.9544968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1178841) q[0];
sx q[0];
rz(-2.5881898) q[0];
sx q[0];
rz(1.2721567) q[0];
x q[1];
rz(-1.5231832) q[2];
sx q[2];
rz(-2.5104668) q[2];
sx q[2];
rz(-2.7472251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33659014) q[1];
sx q[1];
rz(-2.046236) q[1];
sx q[1];
rz(-0.57979433) q[1];
rz(-pi) q[2];
rz(-1.2797221) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(0.94369704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-2.5315703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22247032) q[0];
sx q[0];
rz(-1.3195992) q[0];
sx q[0];
rz(-0.50719502) q[0];
x q[1];
rz(-0.1158175) q[2];
sx q[2];
rz(-2.2216185) q[2];
sx q[2];
rz(-1.7388369) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.94177946) q[1];
sx q[1];
rz(-0.66775125) q[1];
sx q[1];
rz(1.5303395) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(2.4050074) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5442473) q[0];
sx q[0];
rz(-1.5816507) q[0];
sx q[0];
rz(1.6008475) q[0];
x q[1];
rz(-3.0326764) q[2];
sx q[2];
rz(-1.8913336) q[2];
sx q[2];
rz(3.0963754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0524307) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(2.9176941) q[1];
x q[2];
rz(-2.8183476) q[3];
sx q[3];
rz(-1.1439699) q[3];
sx q[3];
rz(-0.95110287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(2.4576808) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-0.7448147) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9841524) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-2.0122583) q[0];
rz(-pi) q[1];
rz(-0.79297519) q[2];
sx q[2];
rz(-1.8599531) q[2];
sx q[2];
rz(-0.58794978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4312268) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(0.23780312) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3585988) q[3];
sx q[3];
rz(-1.2851614) q[3];
sx q[3];
rz(1.0234969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.1709447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0572646) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(1.5469993) q[0];
rz(0.77565907) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(-1.0722216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.4198235) q[1];
sx q[1];
rz(-0.024200736) q[1];
sx q[1];
rz(1.9355378) q[1];
rz(2.714614) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82407172) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(-0.099427632) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-2.4377433) q[2];
sx q[2];
rz(-0.8740295) q[2];
sx q[2];
rz(-0.340273) q[2];
rz(1.2373274) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
