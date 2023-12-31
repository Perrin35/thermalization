OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(5.3806452) q[1];
sx q[1];
rz(7.5723958) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122983) q[0];
sx q[0];
rz(-2.0295143) q[0];
sx q[0];
rz(-0.32231583) q[0];
rz(-1.5640366) q[2];
sx q[2];
rz(-1.5951556) q[2];
sx q[2];
rz(1.3855795) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42613712) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(-0.21547683) q[1];
rz(-pi) q[2];
rz(1.3025769) q[3];
sx q[3];
rz(-1.5032288) q[3];
sx q[3];
rz(-2.4951631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5916799) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(2.8813598) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.82812) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27539) q[0];
sx q[0];
rz(-1.35944) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30828373) q[2];
sx q[2];
rz(-2.0364025) q[2];
sx q[2];
rz(-0.75454933) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5116509) q[1];
sx q[1];
rz(-2.959046) q[1];
sx q[1];
rz(-2.6444142) q[1];
rz(0.95150681) q[3];
sx q[3];
rz(-3.0118239) q[3];
sx q[3];
rz(-1.6956003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(0.95412811) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(0.69586786) q[0];
rz(-0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-2.9755039) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244708) q[0];
sx q[0];
rz(-2.8170601) q[0];
sx q[0];
rz(-2.0425509) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0735998) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(-1.3082248) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6925466) q[1];
sx q[1];
rz(-2.0448301) q[1];
sx q[1];
rz(2.0708592) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8852799) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(-0.29378763) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42472096) q[0];
sx q[0];
rz(-1.6443335) q[0];
sx q[0];
rz(1.5872692) q[0];
rz(-0.36802937) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(0.3993984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1244509) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(-0.084852858) q[1];
rz(-pi) q[2];
rz(2.9702441) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(-1.8360957) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(-1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(-2.9602125) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(0.10770527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2337423) q[0];
sx q[0];
rz(-1.478501) q[0];
sx q[0];
rz(1.924563) q[0];
rz(-pi) q[1];
rz(2.2568251) q[2];
sx q[2];
rz(-2.3713285) q[2];
sx q[2];
rz(-1.9212854) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6377423) q[1];
sx q[1];
rz(-2.3554194) q[1];
sx q[1];
rz(-2.1450858) q[1];
x q[2];
rz(-2.7961736) q[3];
sx q[3];
rz(-1.0343699) q[3];
sx q[3];
rz(2.9649343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(-3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058379563) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(-1.6311197) q[0];
rz(-1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-1.0169792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74422979) q[0];
sx q[0];
rz(-1.1083535) q[0];
sx q[0];
rz(-2.7683393) q[0];
x q[1];
rz(3.0252511) q[2];
sx q[2];
rz(-2.939072) q[2];
sx q[2];
rz(2.4570738) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0000227) q[1];
sx q[1];
rz(-1.1035898) q[1];
sx q[1];
rz(-1.0228553) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6050657) q[3];
sx q[3];
rz(-2.0201207) q[3];
sx q[3];
rz(-2.7262053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-0.44657782) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(0.84386688) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-0.73262501) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0680925) q[0];
sx q[0];
rz(-1.6627899) q[0];
sx q[0];
rz(1.6091225) q[0];
rz(-0.33118601) q[2];
sx q[2];
rz(-1.965431) q[2];
sx q[2];
rz(1.0489724) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0691094) q[1];
sx q[1];
rz(-1.3246037) q[1];
sx q[1];
rz(0.90087955) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3656093) q[3];
sx q[3];
rz(-2.725051) q[3];
sx q[3];
rz(0.52683559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(-0.037242446) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379631) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361355) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(-2.7917042) q[0];
x q[1];
rz(1.2572844) q[2];
sx q[2];
rz(-1.0672617) q[2];
sx q[2];
rz(1.7538278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85247773) q[1];
sx q[1];
rz(-2.5534938) q[1];
sx q[1];
rz(0.27681338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8141361) q[3];
sx q[3];
rz(-2.0753324) q[3];
sx q[3];
rz(2.6017021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-2.5869353) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(-1.6329637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976534) q[0];
sx q[0];
rz(-0.88730747) q[0];
sx q[0];
rz(2.9336714) q[0];
x q[1];
rz(0.091392322) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(-2.6932655) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0867882) q[1];
sx q[1];
rz(-1.0541704) q[1];
sx q[1];
rz(2.7727142) q[1];
rz(-pi) q[2];
rz(0.86767254) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(1.6004576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(2.7049086) q[2];
rz(-1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325539) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59795838) q[0];
sx q[0];
rz(-1.6126452) q[0];
sx q[0];
rz(0.50224395) q[0];
rz(-pi) q[1];
rz(-0.47423116) q[2];
sx q[2];
rz(-1.4719163) q[2];
sx q[2];
rz(-0.41254166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8379412) q[1];
sx q[1];
rz(-2.2397869) q[1];
sx q[1];
rz(-0.53955697) q[1];
x q[2];
rz(2.284427) q[3];
sx q[3];
rz(-1.790159) q[3];
sx q[3];
rz(-2.9227748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116466) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.272841) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(0.045885432) q[2];
sx q[2];
rz(-1.3165717) q[2];
sx q[2];
rz(-1.7236621) q[2];
rz(-1.4894555) q[3];
sx q[3];
rz(-2.079439) q[3];
sx q[3];
rz(0.54872201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
