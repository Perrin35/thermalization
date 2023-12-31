OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(-0.59564367) q[1];
sx q[1];
rz(-1.6593978) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971561) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(-0.42874254) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7035159) q[2];
sx q[2];
rz(-1.8258397) q[2];
sx q[2];
rz(1.0058574) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2130148) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(2.1133452) q[1];
rz(2.1336018) q[3];
sx q[3];
rz(-1.4853012) q[3];
sx q[3];
rz(-2.8443955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-0.86581725) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-2.1781133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022097691) q[0];
sx q[0];
rz(-0.95675981) q[0];
sx q[0];
rz(-3.1387781) q[0];
x q[1];
rz(1.0019148) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(-0.98904726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0986833) q[1];
sx q[1];
rz(-0.60815647) q[1];
sx q[1];
rz(0.98867464) q[1];
rz(-pi) q[2];
rz(1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-0.9427332) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117675) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(-1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(0.55999666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-0.43733894) q[0];
sx q[0];
rz(-1.0768946) q[0];
rz(-pi) q[1];
rz(1.5494924) q[2];
sx q[2];
rz(-1.8739803) q[2];
sx q[2];
rz(-1.6456749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9251688) q[1];
sx q[1];
rz(-0.41453002) q[1];
sx q[1];
rz(-2.6300936) q[1];
x q[2];
rz(0.47645724) q[3];
sx q[3];
rz(-1.1343079) q[3];
sx q[3];
rz(-0.95526327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(2.9690572) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(-2.8570783) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.2987312) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7553058) q[0];
sx q[0];
rz(-2.5884429) q[0];
sx q[0];
rz(2.0095216) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0012871731) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(-3.121701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1620996) q[1];
sx q[1];
rz(-2.5433308) q[1];
sx q[1];
rz(1.9059327) q[1];
rz(-pi) q[2];
rz(-1.8907359) q[3];
sx q[3];
rz(-0.3443998) q[3];
sx q[3];
rz(-2.875945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.6872905) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6032747) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(0.2968266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9223698) q[0];
sx q[0];
rz(-0.57225675) q[0];
sx q[0];
rz(-0.072364307) q[0];
rz(-pi) q[1];
rz(1.7237687) q[2];
sx q[2];
rz(-0.83380552) q[2];
sx q[2];
rz(-1.0066777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2424803) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(0.86109032) q[1];
x q[2];
rz(-1.4818707) q[3];
sx q[3];
rz(-0.85907912) q[3];
sx q[3];
rz(0.54642788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(0.60633916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7327001) q[0];
sx q[0];
rz(-3.0928287) q[0];
sx q[0];
rz(0.34838895) q[0];
rz(-0.89247993) q[2];
sx q[2];
rz(-1.9024897) q[2];
sx q[2];
rz(1.4075556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9669173) q[1];
sx q[1];
rz(-1.7511909) q[1];
sx q[1];
rz(1.0134407) q[1];
rz(-3.0665928) q[3];
sx q[3];
rz(-1.787775) q[3];
sx q[3];
rz(2.0892339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(-1.139337) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(-2.3513444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7471874) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(-1.2140843) q[0];
rz(-0.054025606) q[2];
sx q[2];
rz(-2.9276491) q[2];
sx q[2];
rz(-2.8097048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33752791) q[1];
sx q[1];
rz(-1.3009326) q[1];
sx q[1];
rz(-2.7359664) q[1];
rz(-pi) q[2];
rz(2.4942057) q[3];
sx q[3];
rz(-1.2643407) q[3];
sx q[3];
rz(-0.90741531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(-1.0472263) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(1.75288) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213288) q[0];
sx q[0];
rz(-1.5409894) q[0];
sx q[0];
rz(-1.0111615) q[0];
rz(-pi) q[1];
rz(2.6128204) q[2];
sx q[2];
rz(-1.9430338) q[2];
sx q[2];
rz(0.96226529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(-3.0180879) q[1];
x q[2];
rz(-0.29626131) q[3];
sx q[3];
rz(-0.98516609) q[3];
sx q[3];
rz(2.0010009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1887112) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(-1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700478) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(-0.55823278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687913) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(1.1084491) q[0];
rz(-pi) q[1];
rz(3.0874861) q[2];
sx q[2];
rz(-1.4634842) q[2];
sx q[2];
rz(-3.0215614) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25591125) q[1];
sx q[1];
rz(-2.3951888) q[1];
sx q[1];
rz(1.6649151) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3853119) q[3];
sx q[3];
rz(-1.0914601) q[3];
sx q[3];
rz(2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2925064) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.2020483) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6367209) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035318035) q[0];
sx q[0];
rz(-2.1786852) q[0];
sx q[0];
rz(1.4950698) q[0];
rz(-pi) q[1];
x q[1];
rz(2.142749) q[2];
sx q[2];
rz(-1.5339282) q[2];
sx q[2];
rz(1.2986623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20269468) q[1];
sx q[1];
rz(-1.2564066) q[1];
sx q[1];
rz(0.37757341) q[1];
x q[2];
rz(1.3396026) q[3];
sx q[3];
rz(-1.8095067) q[3];
sx q[3];
rz(-1.064144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(-0.26930299) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(1.833563) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(1.3310824) q[3];
sx q[3];
rz(-0.78835434) q[3];
sx q[3];
rz(2.2314856) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
