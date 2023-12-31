OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(-2.8470319) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3573787) q[0];
sx q[0];
rz(-1.5445909) q[0];
sx q[0];
rz(1.635301) q[0];
rz(-pi) q[1];
rz(1.5316891) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(-1.4456911) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64695839) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(-1.3691982) q[1];
x q[2];
rz(2.9362039) q[3];
sx q[3];
rz(-0.60384149) q[3];
sx q[3];
rz(-2.3103034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.4650311) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1169491) q[0];
sx q[0];
rz(-1.2632217) q[0];
sx q[0];
rz(-1.7988388) q[0];
rz(-pi) q[1];
rz(0.53497603) q[2];
sx q[2];
rz(-1.8638532) q[2];
sx q[2];
rz(-1.5154293) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78039353) q[1];
sx q[1];
rz(-2.2269899) q[1];
sx q[1];
rz(-2.0304297) q[1];
rz(-pi) q[2];
rz(2.6133735) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(1.5147097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(-2.5487652) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-0.16361374) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2046976) q[0];
sx q[0];
rz(-2.097313) q[0];
sx q[0];
rz(-1.5419457) q[0];
rz(1.221237) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(-2.0958054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.219017) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(0.69505691) q[1];
rz(-0.22294873) q[3];
sx q[3];
rz(-1.795459) q[3];
sx q[3];
rz(-0.49062452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0107161) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(-0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-0.83909488) q[0];
rz(1.5377195) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71826868) q[0];
sx q[0];
rz(-1.6800796) q[0];
sx q[0];
rz(0.011358326) q[0];
rz(-pi) q[1];
rz(0.35595603) q[2];
sx q[2];
rz(-1.9320556) q[2];
sx q[2];
rz(-1.661983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3210541) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(1.1827521) q[1];
x q[2];
rz(2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.5031682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.443976) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-2.9042517) q[2];
rz(-2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4603235) q[0];
sx q[0];
rz(-1.8955505) q[0];
sx q[0];
rz(-0.89694174) q[0];
x q[1];
rz(-0.03582844) q[2];
sx q[2];
rz(-1.5140859) q[2];
sx q[2];
rz(-2.5113475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6254127) q[1];
sx q[1];
rz(-1.9429632) q[1];
sx q[1];
rz(-0.401293) q[1];
x q[2];
rz(2.2458514) q[3];
sx q[3];
rz(-2.3268019) q[3];
sx q[3];
rz(-0.18273396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-2.4772947) q[2];
rz(-0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(0.48318133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7532363) q[0];
sx q[0];
rz(-1.6325132) q[0];
sx q[0];
rz(-1.389099) q[0];
x q[1];
rz(1.6364355) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(2.133873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3997765) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(2.497655) q[1];
x q[2];
rz(1.9032848) q[3];
sx q[3];
rz(-1.6786492) q[3];
sx q[3];
rz(-0.40963848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3207021) q[0];
sx q[0];
rz(-1.989869) q[0];
sx q[0];
rz(-1.4120031) q[0];
x q[1];
rz(0.7438296) q[2];
sx q[2];
rz(-2.2093536) q[2];
sx q[2];
rz(-3.0628052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7455709) q[1];
sx q[1];
rz(-2.3704297) q[1];
sx q[1];
rz(-2.1019756) q[1];
rz(2.2456456) q[3];
sx q[3];
rz(-2.575867) q[3];
sx q[3];
rz(-1.5783527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(-1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(1.4138387) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.3716912) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28384128) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(-1.4197423) q[0];
rz(2.9210864) q[2];
sx q[2];
rz(-2.3683511) q[2];
sx q[2];
rz(-2.3436433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.129442) q[1];
sx q[1];
rz(-0.65629849) q[1];
sx q[1];
rz(2.6427569) q[1];
rz(-pi) q[2];
rz(2.7564704) q[3];
sx q[3];
rz(-0.38673863) q[3];
sx q[3];
rz(-0.56798565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9234377) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.6802457) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(-0.27059069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2345679) q[0];
sx q[0];
rz(-1.7922635) q[0];
sx q[0];
rz(1.7937167) q[0];
x q[1];
rz(2.0653535) q[2];
sx q[2];
rz(-0.38947546) q[2];
sx q[2];
rz(-1.7324093) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8671659) q[1];
sx q[1];
rz(-1.5741036) q[1];
sx q[1];
rz(-3.0043688) q[1];
rz(-3.0398265) q[3];
sx q[3];
rz(-2.0483077) q[3];
sx q[3];
rz(-2.4448642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(-1.6516997) q[2];
rz(2.4387032) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437317) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(2.9872966) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001227) q[0];
sx q[0];
rz(-0.98012692) q[0];
sx q[0];
rz(-2.3175879) q[0];
rz(-pi) q[1];
rz(0.43138357) q[2];
sx q[2];
rz(-1.0215534) q[2];
sx q[2];
rz(-1.0647578) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4780477) q[1];
sx q[1];
rz(-1.5819342) q[1];
sx q[1];
rz(0.34578172) q[1];
rz(-pi) q[2];
rz(-1.5417682) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(2.8673429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(-2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(0.32158357) q[2];
sx q[2];
rz(-2.2939773) q[2];
sx q[2];
rz(1.1800223) q[2];
rz(-0.7704173) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
