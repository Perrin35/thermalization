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
rz(-0.91312042) q[0];
sx q[0];
rz(-1.0364391) q[0];
sx q[0];
rz(-1.5358465) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52659863) q[0];
sx q[0];
rz(-2.8078571) q[0];
sx q[0];
rz(1.2646535) q[0];
rz(-pi) q[1];
rz(-0.65459697) q[2];
sx q[2];
rz(-2.4677548) q[2];
sx q[2];
rz(1.0985653) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84387842) q[1];
sx q[1];
rz(-1.529964) q[1];
sx q[1];
rz(-2.2740514) q[1];
x q[2];
rz(-1.3547475) q[3];
sx q[3];
rz(-0.69543823) q[3];
sx q[3];
rz(-1.9239294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6785757) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(-2.2185745) q[2];
rz(3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-3.104082) q[0];
rz(0.58049479) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(0.69211778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0626252) q[0];
sx q[0];
rz(-2.7699372) q[0];
sx q[0];
rz(1.070648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28022556) q[2];
sx q[2];
rz(-1.0839274) q[2];
sx q[2];
rz(3.1161199) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.925732) q[1];
sx q[1];
rz(-1.0714447) q[1];
sx q[1];
rz(-2.4768922) q[1];
x q[2];
rz(1.397527) q[3];
sx q[3];
rz(-2.3152346) q[3];
sx q[3];
rz(-2.2196774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.3956068) q[2];
rz(2.9546402) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66889399) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(0.5589267) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(-0.25159803) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0280625) q[0];
sx q[0];
rz(-2.0087272) q[0];
sx q[0];
rz(1.3597041) q[0];
rz(2.0992804) q[2];
sx q[2];
rz(-0.61868514) q[2];
sx q[2];
rz(1.9306192) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4797552) q[1];
sx q[1];
rz(-0.91714707) q[1];
sx q[1];
rz(2.8089671) q[1];
x q[2];
rz(-0.76520701) q[3];
sx q[3];
rz(-2.3022989) q[3];
sx q[3];
rz(-2.8347903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.805213) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(-1.7097998) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(-1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828204) q[0];
sx q[0];
rz(-0.50676218) q[0];
sx q[0];
rz(1.850542) q[0];
rz(-2.3123815) q[1];
sx q[1];
rz(-1.1081089) q[1];
sx q[1];
rz(-2.1270027) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982581) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(2.696585) q[0];
x q[1];
rz(2.5780771) q[2];
sx q[2];
rz(-1.2402099) q[2];
sx q[2];
rz(-2.3553768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3922799) q[1];
sx q[1];
rz(-1.5754414) q[1];
sx q[1];
rz(-2.1455812) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7144373) q[3];
sx q[3];
rz(-0.82089409) q[3];
sx q[3];
rz(-0.088975541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0097051) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(-0.76060549) q[2];
rz(2.6724114) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(0.73452264) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483265) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(2.560428) q[0];
rz(-1.5515074) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(-2.4466628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3669156) q[0];
sx q[0];
rz(-2.7645676) q[0];
sx q[0];
rz(0.89311262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.476711) q[2];
sx q[2];
rz(-1.861899) q[2];
sx q[2];
rz(-1.7895607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.123053) q[1];
sx q[1];
rz(-0.54529069) q[1];
sx q[1];
rz(-0.084597708) q[1];
x q[2];
rz(-1.3881237) q[3];
sx q[3];
rz(-1.338306) q[3];
sx q[3];
rz(-0.35216613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(2.8387873) q[2];
rz(-1.4263724) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(-2.5572131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645709) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(2.4466178) q[0];
rz(0.28680828) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(-1.274775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5999278) q[0];
sx q[0];
rz(-0.58690161) q[0];
sx q[0];
rz(2.15564) q[0];
rz(-pi) q[1];
rz(2.0389245) q[2];
sx q[2];
rz(-1.0777567) q[2];
sx q[2];
rz(1.4725034) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66836897) q[1];
sx q[1];
rz(-0.76938564) q[1];
sx q[1];
rz(-1.046087) q[1];
rz(-pi) q[2];
rz(1.8096089) q[3];
sx q[3];
rz(-2.1396355) q[3];
sx q[3];
rz(-0.59900008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0605165) q[2];
sx q[2];
rz(-2.923107) q[2];
sx q[2];
rz(0.9673574) q[2];
rz(0.71990144) q[3];
sx q[3];
rz(-2.1981809) q[3];
sx q[3];
rz(1.2758183) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8571781) q[0];
sx q[0];
rz(-0.030817742) q[0];
sx q[0];
rz(-1.6946633) q[0];
rz(0.46724304) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(-1.9460868) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0366096) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(-0.81654136) q[0];
rz(-pi) q[1];
rz(-1.2345957) q[2];
sx q[2];
rz(-0.78374353) q[2];
sx q[2];
rz(0.23434336) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1165645) q[1];
sx q[1];
rz(-1.0724853) q[1];
sx q[1];
rz(-0.070498585) q[1];
rz(1.9105114) q[3];
sx q[3];
rz(-3.1034443) q[3];
sx q[3];
rz(1.3919786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3168827) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(0.52428025) q[2];
rz(2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(-3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(1.9148781) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(-0.39047584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1487685) q[0];
sx q[0];
rz(-2.618194) q[0];
sx q[0];
rz(-2.5958909) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6975193) q[2];
sx q[2];
rz(-1.6055428) q[2];
sx q[2];
rz(-2.4028558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.272534) q[1];
sx q[1];
rz(-2.7400937) q[1];
sx q[1];
rz(2.1436006) q[1];
rz(-2.212575) q[3];
sx q[3];
rz(-1.1631771) q[3];
sx q[3];
rz(-0.81002012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32435027) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(1.3520757) q[2];
rz(2.9663441) q[3];
sx q[3];
rz(-1.3388355) q[3];
sx q[3];
rz(-1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7290203) q[0];
sx q[0];
rz(-0.46171284) q[0];
sx q[0];
rz(0.66226688) q[0];
rz(-0.63016713) q[1];
sx q[1];
rz(-1.2799193) q[1];
sx q[1];
rz(-2.9060649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5109638) q[0];
sx q[0];
rz(-2.8665344) q[0];
sx q[0];
rz(-2.4514626) q[0];
rz(-pi) q[1];
rz(1.477219) q[2];
sx q[2];
rz(-2.1380599) q[2];
sx q[2];
rz(2.7964489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3429361) q[1];
sx q[1];
rz(-1.6321811) q[1];
sx q[1];
rz(-0.023691698) q[1];
rz(-pi) q[2];
rz(-2.2673549) q[3];
sx q[3];
rz(-1.6739419) q[3];
sx q[3];
rz(-2.646605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0528637) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(2.5843184) q[2];
rz(-0.82734621) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(-1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92852229) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(-2.573977) q[0];
rz(2.7396743) q[1];
sx q[1];
rz(-0.64359507) q[1];
sx q[1];
rz(1.7793122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069893006) q[0];
sx q[0];
rz(-1.9687459) q[0];
sx q[0];
rz(-1.1342628) q[0];
rz(-2.882233) q[2];
sx q[2];
rz(-2.1742714) q[2];
sx q[2];
rz(3.0779787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.030008121) q[1];
sx q[1];
rz(-1.6492731) q[1];
sx q[1];
rz(0.93941883) q[1];
rz(-pi) q[2];
rz(-0.74257039) q[3];
sx q[3];
rz(-2.1705063) q[3];
sx q[3];
rz(0.35123435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(0.26025772) q[2];
rz(-2.4506954) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.6183841) q[3];
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
rz(2.8061168) q[0];
sx q[0];
rz(-1.5806883) q[0];
sx q[0];
rz(-1.6089532) q[0];
rz(0.43329049) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(-2.4337089) q[2];
sx q[2];
rz(-2.08692) q[2];
sx q[2];
rz(-0.89649123) q[2];
rz(2.3512202) q[3];
sx q[3];
rz(-1.6675259) q[3];
sx q[3];
rz(-0.70601757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
