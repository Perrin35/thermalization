OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.051209) q[0];
sx q[0];
rz(-0.78855711) q[0];
sx q[0];
rz(-0.30858827) q[0];
rz(0.42449549) q[1];
sx q[1];
rz(4.2257809) q[1];
sx q[1];
rz(9.0492166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0610866) q[0];
sx q[0];
rz(-1.8520141) q[0];
sx q[0];
rz(-0.078769509) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8396222) q[2];
sx q[2];
rz(-1.6014546) q[2];
sx q[2];
rz(1.7177402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5243317) q[1];
sx q[1];
rz(-1.4325805) q[1];
sx q[1];
rz(0.44869615) q[1];
rz(-pi) q[2];
rz(3.1377596) q[3];
sx q[3];
rz(-1.7143013) q[3];
sx q[3];
rz(-1.6832863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8033119) q[2];
sx q[2];
rz(-0.25140005) q[2];
sx q[2];
rz(0.54331642) q[2];
rz(1.733689) q[3];
sx q[3];
rz(-1.7270154) q[3];
sx q[3];
rz(-0.89731115) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5498085) q[0];
sx q[0];
rz(-1.9928638) q[0];
sx q[0];
rz(-2.5445004) q[0];
rz(1.0626556) q[1];
sx q[1];
rz(-2.5025044) q[1];
sx q[1];
rz(-2.7426681) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8455578) q[0];
sx q[0];
rz(-2.425023) q[0];
sx q[0];
rz(-0.96591732) q[0];
x q[1];
rz(-1.0137896) q[2];
sx q[2];
rz(-1.5913871) q[2];
sx q[2];
rz(2.2228732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0831757) q[1];
sx q[1];
rz(-2.1444369) q[1];
sx q[1];
rz(0.87223335) q[1];
x q[2];
rz(-0.90090294) q[3];
sx q[3];
rz(-0.76053166) q[3];
sx q[3];
rz(0.72872114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67542934) q[2];
sx q[2];
rz(-1.9025981) q[2];
sx q[2];
rz(-1.1413525) q[2];
rz(-2.2232248) q[3];
sx q[3];
rz(-0.71056241) q[3];
sx q[3];
rz(0.58951283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.454575) q[0];
sx q[0];
rz(-2.7641251) q[0];
sx q[0];
rz(-2.9662509) q[0];
rz(-1.31458) q[1];
sx q[1];
rz(-1.250123) q[1];
sx q[1];
rz(-2.4741727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31418598) q[0];
sx q[0];
rz(-1.6860322) q[0];
sx q[0];
rz(2.0699571) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4756617) q[2];
sx q[2];
rz(-2.4325822) q[2];
sx q[2];
rz(-0.60707742) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7857915) q[1];
sx q[1];
rz(-2.8257824) q[1];
sx q[1];
rz(0.0018382536) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61428689) q[3];
sx q[3];
rz(-1.3516055) q[3];
sx q[3];
rz(-0.012085513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1158925) q[2];
sx q[2];
rz(-0.95558715) q[2];
sx q[2];
rz(-2.058775) q[2];
rz(1.8358021) q[3];
sx q[3];
rz(-2.3256153) q[3];
sx q[3];
rz(2.344237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1378655) q[0];
sx q[0];
rz(-0.19345134) q[0];
sx q[0];
rz(-2.5306012) q[0];
rz(-2.5773279) q[1];
sx q[1];
rz(-1.0287501) q[1];
sx q[1];
rz(-2.4345523) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90908289) q[0];
sx q[0];
rz(-2.3444026) q[0];
sx q[0];
rz(-0.042122201) q[0];
x q[1];
rz(0.41019812) q[2];
sx q[2];
rz(-1.2648404) q[2];
sx q[2];
rz(-3.0434639) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1241856) q[1];
sx q[1];
rz(-2.0059791) q[1];
sx q[1];
rz(2.6098677) q[1];
x q[2];
rz(2.6534326) q[3];
sx q[3];
rz(-0.55503856) q[3];
sx q[3];
rz(2.6223573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45346144) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(-2.1235662) q[2];
rz(-0.38797837) q[3];
sx q[3];
rz(-1.6951025) q[3];
sx q[3];
rz(0.33897266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46868789) q[0];
sx q[0];
rz(-2.0141116) q[0];
sx q[0];
rz(-2.6197523) q[0];
rz(-2.8537967) q[1];
sx q[1];
rz(-2.5591873) q[1];
sx q[1];
rz(-2.5385181) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1408242) q[0];
sx q[0];
rz(-3.040977) q[0];
sx q[0];
rz(-0.510143) q[0];
rz(-2.4580946) q[2];
sx q[2];
rz(-1.1199691) q[2];
sx q[2];
rz(-1.3749387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9085616) q[1];
sx q[1];
rz(-2.7082293) q[1];
sx q[1];
rz(2.9410081) q[1];
rz(-pi) q[2];
rz(-0.92080812) q[3];
sx q[3];
rz(-0.42144708) q[3];
sx q[3];
rz(-1.309066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4294943) q[2];
sx q[2];
rz(-0.79271972) q[2];
sx q[2];
rz(-0.73954868) q[2];
rz(-1.0150917) q[3];
sx q[3];
rz(-0.52777094) q[3];
sx q[3];
rz(3.1225865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5571112) q[0];
sx q[0];
rz(-2.8040573) q[0];
sx q[0];
rz(-0.8648411) q[0];
rz(-2.5014014) q[1];
sx q[1];
rz(-2.4914111) q[1];
sx q[1];
rz(0.4943628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664864) q[0];
sx q[0];
rz(-0.49577489) q[0];
sx q[0];
rz(2.0074816) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.041448822) q[2];
sx q[2];
rz(-2.1728656) q[2];
sx q[2];
rz(2.1727501) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2701708) q[1];
sx q[1];
rz(-0.06690678) q[1];
sx q[1];
rz(-0.54065458) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0406144) q[3];
sx q[3];
rz(-1.1270071) q[3];
sx q[3];
rz(-0.63999301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.87259) q[2];
sx q[2];
rz(-0.63429093) q[2];
sx q[2];
rz(0.96160257) q[2];
rz(1.6327935) q[3];
sx q[3];
rz(-1.6274933) q[3];
sx q[3];
rz(0.076920286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85840571) q[0];
sx q[0];
rz(-2.7981693) q[0];
sx q[0];
rz(-1.0431694) q[0];
rz(2.4647602) q[1];
sx q[1];
rz(-0.64462858) q[1];
sx q[1];
rz(-0.72789311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5344369) q[0];
sx q[0];
rz(-2.712116) q[0];
sx q[0];
rz(1.8398102) q[0];
x q[1];
rz(-1.5743718) q[2];
sx q[2];
rz(-0.82873949) q[2];
sx q[2];
rz(-0.34012857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5965569) q[1];
sx q[1];
rz(-1.903773) q[1];
sx q[1];
rz(2.2117028) q[1];
rz(-0.48850208) q[3];
sx q[3];
rz(-2.4932749) q[3];
sx q[3];
rz(2.019995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4595043) q[2];
sx q[2];
rz(-2.3900718) q[2];
sx q[2];
rz(2.6085243) q[2];
rz(-1.8367977) q[3];
sx q[3];
rz(-0.47520906) q[3];
sx q[3];
rz(2.3615725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.099139) q[0];
sx q[0];
rz(-3.070153) q[0];
sx q[0];
rz(3.1169917) q[0];
rz(3.091231) q[1];
sx q[1];
rz(-2.5854526) q[1];
sx q[1];
rz(-0.76739001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25319812) q[0];
sx q[0];
rz(-2.1108187) q[0];
sx q[0];
rz(-1.5408976) q[0];
rz(-pi) q[1];
rz(-1.3774302) q[2];
sx q[2];
rz(-1.8948717) q[2];
sx q[2];
rz(2.0266768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9995867) q[1];
sx q[1];
rz(-0.97181994) q[1];
sx q[1];
rz(-1.0890278) q[1];
rz(-0.56782372) q[3];
sx q[3];
rz(-0.65560549) q[3];
sx q[3];
rz(-2.034831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6699827) q[2];
sx q[2];
rz(-0.47405425) q[2];
sx q[2];
rz(-1.6101884) q[2];
rz(1.0129741) q[3];
sx q[3];
rz(-1.4644724) q[3];
sx q[3];
rz(0.60232919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36039627) q[0];
sx q[0];
rz(-0.66806) q[0];
sx q[0];
rz(0.49041954) q[0];
rz(-0.40052739) q[1];
sx q[1];
rz(-0.35919765) q[1];
sx q[1];
rz(-1.5708956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2928481) q[0];
sx q[0];
rz(-1.9269173) q[0];
sx q[0];
rz(-0.57909052) q[0];
rz(-pi) q[1];
rz(-0.3272662) q[2];
sx q[2];
rz(-1.2246574) q[2];
sx q[2];
rz(-0.95771433) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3048858) q[1];
sx q[1];
rz(-1.9203102) q[1];
sx q[1];
rz(-3.034449) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79244198) q[3];
sx q[3];
rz(-1.0002478) q[3];
sx q[3];
rz(-1.5287409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5639497) q[2];
sx q[2];
rz(-0.22259139) q[2];
sx q[2];
rz(0.27601784) q[2];
rz(2.4769619) q[3];
sx q[3];
rz(-2.1641927) q[3];
sx q[3];
rz(-0.19653921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68263245) q[0];
sx q[0];
rz(-2.9660872) q[0];
sx q[0];
rz(-1.5631787) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-1.027532) q[1];
sx q[1];
rz(-3.077502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4558217) q[0];
sx q[0];
rz(-1.5470354) q[0];
sx q[0];
rz(-2.0288951) q[0];
x q[1];
rz(1.1214549) q[2];
sx q[2];
rz(-1.0021415) q[2];
sx q[2];
rz(1.519738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42228544) q[1];
sx q[1];
rz(-1.6537084) q[1];
sx q[1];
rz(1.0817429) q[1];
rz(2.1782297) q[3];
sx q[3];
rz(-2.2775536) q[3];
sx q[3];
rz(2.6092495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2854707) q[2];
sx q[2];
rz(-2.8557114) q[2];
sx q[2];
rz(2.2355283) q[2];
rz(1.2178577) q[3];
sx q[3];
rz(-1.2713615) q[3];
sx q[3];
rz(2.0876032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.8080407) q[0];
sx q[0];
rz(-1.6186436) q[0];
sx q[0];
rz(2.9343395) q[0];
rz(2.8605657) q[1];
sx q[1];
rz(-1.3916176) q[1];
sx q[1];
rz(2.4236046) q[1];
rz(1.6529374) q[2];
sx q[2];
rz(-1.1700656) q[2];
sx q[2];
rz(-1.1943371) q[2];
rz(2.6646176) q[3];
sx q[3];
rz(-1.5834482) q[3];
sx q[3];
rz(-0.96848828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
