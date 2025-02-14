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
rz(1.9361629) q[0];
sx q[0];
rz(3.087145) q[0];
sx q[0];
rz(8.5589391) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(3.08334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39393878) q[0];
sx q[0];
rz(-0.17102111) q[0];
sx q[0];
rz(-1.8500516) q[0];
rz(3.0676663) q[2];
sx q[2];
rz(-1.0763028) q[2];
sx q[2];
rz(1.489515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21319751) q[1];
sx q[1];
rz(-1.0093912) q[1];
sx q[1];
rz(-2.2092845) q[1];
rz(-2.3184003) q[3];
sx q[3];
rz(-0.45837742) q[3];
sx q[3];
rz(2.2708054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3129348) q[2];
sx q[2];
rz(-2.2609495) q[2];
sx q[2];
rz(2.177296) q[2];
rz(-3.0621081) q[3];
sx q[3];
rz(-1.3026404) q[3];
sx q[3];
rz(0.77118072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1281779) q[0];
sx q[0];
rz(-0.34341136) q[0];
sx q[0];
rz(-1.6012023) q[0];
rz(-0.84114289) q[1];
sx q[1];
rz(-1.4079739) q[1];
sx q[1];
rz(2.2841563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6733137) q[0];
sx q[0];
rz(-2.7029917) q[0];
sx q[0];
rz(-1.5836141) q[0];
x q[1];
rz(2.2067864) q[2];
sx q[2];
rz(-1.9611214) q[2];
sx q[2];
rz(-1.1771415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.407436) q[1];
sx q[1];
rz(-2.4333471) q[1];
sx q[1];
rz(0.70869653) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6667834) q[3];
sx q[3];
rz(-1.848683) q[3];
sx q[3];
rz(-3.0888219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2555344) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(-1.9319755) q[2];
rz(-1.7602734) q[3];
sx q[3];
rz(-1.8807024) q[3];
sx q[3];
rz(0.59022006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34476122) q[0];
sx q[0];
rz(-1.3153356) q[0];
sx q[0];
rz(-1.8094081) q[0];
rz(-1.4765129) q[1];
sx q[1];
rz(-1.3308728) q[1];
sx q[1];
rz(0.61980334) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7999511) q[0];
sx q[0];
rz(-1.5864306) q[0];
sx q[0];
rz(0.077341103) q[0];
rz(-pi) q[1];
rz(-2.1479883) q[2];
sx q[2];
rz(-2.8393778) q[2];
sx q[2];
rz(2.4156918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18869723) q[1];
sx q[1];
rz(-2.8159375) q[1];
sx q[1];
rz(-0.22518994) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93094935) q[3];
sx q[3];
rz(-1.7849277) q[3];
sx q[3];
rz(-2.3420985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0009813112) q[2];
sx q[2];
rz(-0.32537127) q[2];
sx q[2];
rz(0.78835431) q[2];
rz(-1.2761448) q[3];
sx q[3];
rz(-1.9143462) q[3];
sx q[3];
rz(-0.86281323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1514423) q[0];
sx q[0];
rz(-1.3862415) q[0];
sx q[0];
rz(1.066712) q[0];
rz(1.3534631) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(1.1563168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8114093) q[0];
sx q[0];
rz(-1.5519987) q[0];
sx q[0];
rz(1.4292595) q[0];
x q[1];
rz(-0.026211003) q[2];
sx q[2];
rz(-1.6115341) q[2];
sx q[2];
rz(0.36529503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20322415) q[1];
sx q[1];
rz(-0.48048204) q[1];
sx q[1];
rz(2.234747) q[1];
rz(-1.5045936) q[3];
sx q[3];
rz(-2.2341223) q[3];
sx q[3];
rz(0.56831336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.452423) q[2];
sx q[2];
rz(-0.66340476) q[2];
sx q[2];
rz(-0.058825292) q[2];
rz(-2.5022653) q[3];
sx q[3];
rz(-1.2213734) q[3];
sx q[3];
rz(0.85734573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0462129) q[0];
sx q[0];
rz(-1.4154499) q[0];
sx q[0];
rz(3.131102) q[0];
rz(1.563021) q[1];
sx q[1];
rz(-0.69067162) q[1];
sx q[1];
rz(-0.48935997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90411579) q[0];
sx q[0];
rz(-1.7286219) q[0];
sx q[0];
rz(0.62741168) q[0];
rz(-1.0919369) q[2];
sx q[2];
rz(-0.74981028) q[2];
sx q[2];
rz(-2.9887226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6076679) q[1];
sx q[1];
rz(-1.2180389) q[1];
sx q[1];
rz(2.4466848) q[1];
rz(0.33133026) q[3];
sx q[3];
rz(-0.032535527) q[3];
sx q[3];
rz(1.0573385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8420777) q[2];
sx q[2];
rz(-1.0966417) q[2];
sx q[2];
rz(-0.00046029885) q[2];
rz(-2.4680303) q[3];
sx q[3];
rz(-2.2431777) q[3];
sx q[3];
rz(2.4323997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29086581) q[0];
sx q[0];
rz(-1.8274266) q[0];
sx q[0];
rz(-1.8379743) q[0];
rz(-0.29779008) q[1];
sx q[1];
rz(-2.2460263) q[1];
sx q[1];
rz(2.8477125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0654945) q[0];
sx q[0];
rz(-0.84653234) q[0];
sx q[0];
rz(0.80418555) q[0];
rz(0.43115012) q[2];
sx q[2];
rz(-2.2047055) q[2];
sx q[2];
rz(2.0946787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9781295) q[1];
sx q[1];
rz(-2.4678964) q[1];
sx q[1];
rz(-0.81836318) q[1];
x q[2];
rz(-2.8260293) q[3];
sx q[3];
rz(-1.3466918) q[3];
sx q[3];
rz(-1.4423646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6728354) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(-2.920816) q[2];
rz(-1.8686434) q[3];
sx q[3];
rz(-1.3947398) q[3];
sx q[3];
rz(1.1481736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4918168) q[0];
sx q[0];
rz(-2.5204372) q[0];
sx q[0];
rz(0.57449269) q[0];
rz(-2.5993787) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(0.39074674) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1339392) q[0];
sx q[0];
rz(-1.7152733) q[0];
sx q[0];
rz(2.5424315) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0894012) q[2];
sx q[2];
rz(-1.0269916) q[2];
sx q[2];
rz(0.79629818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83929944) q[1];
sx q[1];
rz(-1.5652735) q[1];
sx q[1];
rz(-0.41673613) q[1];
rz(-pi) q[2];
rz(3.0454759) q[3];
sx q[3];
rz(-2.7743154) q[3];
sx q[3];
rz(0.27378455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4603525) q[2];
sx q[2];
rz(-1.9443024) q[2];
sx q[2];
rz(-0.61275068) q[2];
rz(0.98226205) q[3];
sx q[3];
rz(-1.3145612) q[3];
sx q[3];
rz(2.6764892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9923582) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(-0.055334844) q[0];
rz(-1.5845567) q[1];
sx q[1];
rz(-1.0880071) q[1];
sx q[1];
rz(-0.43992821) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8986788) q[0];
sx q[0];
rz(-1.1132332) q[0];
sx q[0];
rz(-1.2957063) q[0];
x q[1];
rz(-0.49100809) q[2];
sx q[2];
rz(-1.3898808) q[2];
sx q[2];
rz(0.79800765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7288558) q[1];
sx q[1];
rz(-1.8860894) q[1];
sx q[1];
rz(-0.58087491) q[1];
x q[2];
rz(-1.1433268) q[3];
sx q[3];
rz(-0.60986589) q[3];
sx q[3];
rz(2.0855869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7679193) q[2];
sx q[2];
rz(-1.3481216) q[2];
sx q[2];
rz(0.31201735) q[2];
rz(0.9196552) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(-2.9254204) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31174082) q[0];
sx q[0];
rz(-0.98167247) q[0];
sx q[0];
rz(-3.1134636) q[0];
rz(1.4460538) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(2.6820954) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.979769) q[0];
sx q[0];
rz(-1.5422675) q[0];
sx q[0];
rz(-1.606444) q[0];
rz(-pi) q[1];
rz(-1.9295646) q[2];
sx q[2];
rz(-0.78132403) q[2];
sx q[2];
rz(1.5272905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7920835) q[1];
sx q[1];
rz(-2.4448423) q[1];
sx q[1];
rz(-1.3443483) q[1];
rz(-0.13925456) q[3];
sx q[3];
rz(-1.1920658) q[3];
sx q[3];
rz(-1.9851353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0405937) q[2];
sx q[2];
rz(-1.7969635) q[2];
sx q[2];
rz(2.8864268) q[2];
rz(0.98662871) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(-1.7184947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7607018) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(-1.799452) q[0];
rz(-3.1396719) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(2.0916746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0107897) q[0];
sx q[0];
rz(-1.783004) q[0];
sx q[0];
rz(0.15168587) q[0];
x q[1];
rz(-1.5639864) q[2];
sx q[2];
rz(-0.17596888) q[2];
sx q[2];
rz(2.3962452) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44651664) q[1];
sx q[1];
rz(-0.43276603) q[1];
sx q[1];
rz(-2.5997735) q[1];
x q[2];
rz(-2.3008651) q[3];
sx q[3];
rz(-1.7326712) q[3];
sx q[3];
rz(0.22541287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92901015) q[2];
sx q[2];
rz(-1.2139823) q[2];
sx q[2];
rz(2.6455961) q[2];
rz(1.4604733) q[3];
sx q[3];
rz(-1.7998327) q[3];
sx q[3];
rz(1.9546485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37954189) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(2.5024391) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(0.46001804) q[2];
sx q[2];
rz(-2.620853) q[2];
sx q[2];
rz(-3.1070409) q[2];
rz(-0.88981723) q[3];
sx q[3];
rz(-2.2636236) q[3];
sx q[3];
rz(-1.3286535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
