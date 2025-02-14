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
rz(-3.1373625) q[0];
sx q[0];
rz(-2.5047996) q[0];
sx q[0];
rz(1.1871583) q[0];
rz(0.98183331) q[1];
sx q[1];
rz(-0.47458664) q[1];
sx q[1];
rz(2.2003953) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1008074) q[0];
sx q[0];
rz(-1.5691162) q[0];
sx q[0];
rz(-0.4953065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.773052) q[2];
sx q[2];
rz(-1.4448243) q[2];
sx q[2];
rz(-1.6304315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8601712) q[1];
sx q[1];
rz(-2.1263391) q[1];
sx q[1];
rz(-1.1475599) q[1];
rz(-pi) q[2];
rz(3.0676316) q[3];
sx q[3];
rz(-1.1152905) q[3];
sx q[3];
rz(0.8616254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59419751) q[2];
sx q[2];
rz(-1.5028468) q[2];
sx q[2];
rz(1.0198062) q[2];
rz(0.29113537) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(1.8118793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52784598) q[0];
sx q[0];
rz(-0.83714038) q[0];
sx q[0];
rz(1.5035195) q[0];
rz(-1.9095406) q[1];
sx q[1];
rz(-0.82545311) q[1];
sx q[1];
rz(-1.8881316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37663662) q[0];
sx q[0];
rz(-2.6692163) q[0];
sx q[0];
rz(1.0563904) q[0];
rz(-1.8623146) q[2];
sx q[2];
rz(-2.0644858) q[2];
sx q[2];
rz(1.3048628) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3121705) q[1];
sx q[1];
rz(-2.8924184) q[1];
sx q[1];
rz(0.62015688) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0157865) q[3];
sx q[3];
rz(-0.87712327) q[3];
sx q[3];
rz(1.7852269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85679179) q[2];
sx q[2];
rz(-2.5992726) q[2];
sx q[2];
rz(-2.8112603) q[2];
rz(0.10279837) q[3];
sx q[3];
rz(-1.2448575) q[3];
sx q[3];
rz(0.61906329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5264346) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(-2.0779628) q[0];
rz(0.10387736) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(-2.0361384) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.010546) q[0];
sx q[0];
rz(-1.5233114) q[0];
sx q[0];
rz(-3.1094883) q[0];
rz(-pi) q[1];
rz(2.8145829) q[2];
sx q[2];
rz(-0.78377027) q[2];
sx q[2];
rz(-1.4732337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7522831) q[1];
sx q[1];
rz(-2.0842617) q[1];
sx q[1];
rz(-3.0016243) q[1];
rz(-0.01171578) q[3];
sx q[3];
rz(-1.5289617) q[3];
sx q[3];
rz(-1.6138298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0217454) q[2];
sx q[2];
rz(-2.3357119) q[2];
sx q[2];
rz(-0.22551192) q[2];
rz(1.8146023) q[3];
sx q[3];
rz(-0.83695379) q[3];
sx q[3];
rz(2.5007611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91318146) q[0];
sx q[0];
rz(-1.5786194) q[0];
sx q[0];
rz(2.5033503) q[0];
rz(2.0300716) q[1];
sx q[1];
rz(-0.94739729) q[1];
sx q[1];
rz(2.9543455) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032124539) q[0];
sx q[0];
rz(-2.1817142) q[0];
sx q[0];
rz(-2.7892668) q[0];
rz(0.78088607) q[2];
sx q[2];
rz(-1.5263288) q[2];
sx q[2];
rz(1.5574297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9065215) q[1];
sx q[1];
rz(-1.491761) q[1];
sx q[1];
rz(2.6071878) q[1];
rz(-pi) q[2];
rz(-1.6583913) q[3];
sx q[3];
rz(-0.85932743) q[3];
sx q[3];
rz(-0.30184823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1049261) q[2];
sx q[2];
rz(-1.6380402) q[2];
sx q[2];
rz(0.46910134) q[2];
rz(-0.30096287) q[3];
sx q[3];
rz(-2.8807237) q[3];
sx q[3];
rz(-1.4891967) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8982573) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(-2.6014056) q[0];
rz(-0.46145269) q[1];
sx q[1];
rz(-1.5216454) q[1];
sx q[1];
rz(2.8823421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.741138) q[0];
sx q[0];
rz(-1.8807066) q[0];
sx q[0];
rz(-1.8551338) q[0];
rz(-pi) q[1];
rz(0.44633003) q[2];
sx q[2];
rz(-2.0392373) q[2];
sx q[2];
rz(2.6233752) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67422489) q[1];
sx q[1];
rz(-1.7265872) q[1];
sx q[1];
rz(0.99512659) q[1];
rz(-1.9018039) q[3];
sx q[3];
rz(-1.8192623) q[3];
sx q[3];
rz(1.0120165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79869142) q[2];
sx q[2];
rz(-2.4553802) q[2];
sx q[2];
rz(-2.2314824) q[2];
rz(2.5765007) q[3];
sx q[3];
rz(-1.7730954) q[3];
sx q[3];
rz(-0.71649396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703576) q[0];
sx q[0];
rz(-0.39327455) q[0];
sx q[0];
rz(-1.2275335) q[0];
rz(-2.6365989) q[1];
sx q[1];
rz(-2.125449) q[1];
sx q[1];
rz(1.203677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01194309) q[0];
sx q[0];
rz(-0.27515689) q[0];
sx q[0];
rz(0.45264523) q[0];
x q[1];
rz(-0.48690472) q[2];
sx q[2];
rz(-0.98003888) q[2];
sx q[2];
rz(2.7529773) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56702891) q[1];
sx q[1];
rz(-2.1334799) q[1];
sx q[1];
rz(-0.74498727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0685445) q[3];
sx q[3];
rz(-2.6274791) q[3];
sx q[3];
rz(-0.83465605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(-0.95899478) q[2];
rz(-0.2995019) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(-1.449077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2285948) q[0];
sx q[0];
rz(-1.2499502) q[0];
sx q[0];
rz(2.6477497) q[0];
rz(0.85810581) q[1];
sx q[1];
rz(-1.9789507) q[1];
sx q[1];
rz(1.7869305) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345024) q[0];
sx q[0];
rz(-2.1190152) q[0];
sx q[0];
rz(-2.8584418) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3184403) q[2];
sx q[2];
rz(-2.6359973) q[2];
sx q[2];
rz(2.8211968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87397447) q[1];
sx q[1];
rz(-2.3268891) q[1];
sx q[1];
rz(1.4741298) q[1];
rz(-2.2432547) q[3];
sx q[3];
rz(-1.9914989) q[3];
sx q[3];
rz(-1.6572233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4570423) q[2];
sx q[2];
rz(-0.6251308) q[2];
sx q[2];
rz(-2.5954512) q[2];
rz(0.1977194) q[3];
sx q[3];
rz(-1.0309912) q[3];
sx q[3];
rz(2.8945967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207183) q[0];
sx q[0];
rz(-1.1552759) q[0];
sx q[0];
rz(-2.0350631) q[0];
rz(0.58440343) q[1];
sx q[1];
rz(-1.9677275) q[1];
sx q[1];
rz(1.0027142) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4179392) q[0];
sx q[0];
rz(-2.0040211) q[0];
sx q[0];
rz(0.17670011) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24248388) q[2];
sx q[2];
rz(-1.2839181) q[2];
sx q[2];
rz(-0.74082021) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7948053) q[1];
sx q[1];
rz(-1.0587988) q[1];
sx q[1];
rz(-1.5293299) q[1];
x q[2];
rz(-0.61912304) q[3];
sx q[3];
rz(-1.8612218) q[3];
sx q[3];
rz(-0.30916902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5661085) q[2];
sx q[2];
rz(-2.2169952) q[2];
sx q[2];
rz(1.7522579) q[2];
rz(2.6050881) q[3];
sx q[3];
rz(-1.6335231) q[3];
sx q[3];
rz(0.84684053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.093512) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(-2.804948) q[0];
rz(-3.0632784) q[1];
sx q[1];
rz(-1.9322194) q[1];
sx q[1];
rz(1.9858817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1378095) q[0];
sx q[0];
rz(-2.3562618) q[0];
sx q[0];
rz(-1.0850052) q[0];
x q[1];
rz(0.066566932) q[2];
sx q[2];
rz(-1.9835501) q[2];
sx q[2];
rz(2.33605) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3382321) q[1];
sx q[1];
rz(-0.50932099) q[1];
sx q[1];
rz(1.8503016) q[1];
rz(-1.5404245) q[3];
sx q[3];
rz(-1.8874468) q[3];
sx q[3];
rz(0.79398549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0597421) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(2.3183863) q[2];
rz(-2.1246223) q[3];
sx q[3];
rz(-0.40139324) q[3];
sx q[3];
rz(0.038173525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391649) q[0];
sx q[0];
rz(-2.9186354) q[0];
sx q[0];
rz(-1.4270576) q[0];
rz(1.6629705) q[1];
sx q[1];
rz(-1.0717816) q[1];
sx q[1];
rz(-0.85085416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44247791) q[0];
sx q[0];
rz(-2.6529387) q[0];
sx q[0];
rz(1.3241934) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4434187) q[2];
sx q[2];
rz(-1.6568523) q[2];
sx q[2];
rz(-2.0620777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46063328) q[1];
sx q[1];
rz(-2.6151028) q[1];
sx q[1];
rz(-0.10641702) q[1];
x q[2];
rz(-3.0012742) q[3];
sx q[3];
rz(-1.5068973) q[3];
sx q[3];
rz(1.4177223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.458638) q[2];
sx q[2];
rz(-1.3502324) q[2];
sx q[2];
rz(-0.82009912) q[2];
rz(1.5470777) q[3];
sx q[3];
rz(-1.7087414) q[3];
sx q[3];
rz(1.9471751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1632669) q[0];
sx q[0];
rz(-1.5112725) q[0];
sx q[0];
rz(1.5164966) q[0];
rz(2.3691879) q[1];
sx q[1];
rz(-2.518387) q[1];
sx q[1];
rz(1.0060681) q[1];
rz(-2.8683581) q[2];
sx q[2];
rz(-0.71003503) q[2];
sx q[2];
rz(2.4207921) q[2];
rz(2.3631572) q[3];
sx q[3];
rz(-2.5166471) q[3];
sx q[3];
rz(-2.2770721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
