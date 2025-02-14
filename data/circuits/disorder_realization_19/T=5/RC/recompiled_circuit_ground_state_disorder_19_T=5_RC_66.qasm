OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8813397) q[0];
sx q[0];
rz(-0.94085675) q[0];
sx q[0];
rz(2.9139304) q[0];
rz(-2.8582299) q[1];
sx q[1];
rz(-0.41937399) q[1];
sx q[1];
rz(-1.8546606) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5628107) q[0];
sx q[0];
rz(-1.1039886) q[0];
sx q[0];
rz(-0.060725529) q[0];
rz(-pi) q[1];
rz(2.3782733) q[2];
sx q[2];
rz(-2.5322891) q[2];
sx q[2];
rz(-2.0610025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7007926) q[1];
sx q[1];
rz(-1.1733353) q[1];
sx q[1];
rz(-2.1093919) q[1];
rz(-0.41051045) q[3];
sx q[3];
rz(-0.68119739) q[3];
sx q[3];
rz(0.42689161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2490354) q[2];
sx q[2];
rz(-1.0113357) q[2];
sx q[2];
rz(-0.91157836) q[2];
rz(0.75561953) q[3];
sx q[3];
rz(-2.8204462) q[3];
sx q[3];
rz(0.49629456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.3590473) q[0];
sx q[0];
rz(-1.8308715) q[0];
sx q[0];
rz(2.0027335) q[0];
rz(-0.62659872) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(-2.0638594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2127578) q[0];
sx q[0];
rz(-2.1427615) q[0];
sx q[0];
rz(-1.1642745) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7863492) q[2];
sx q[2];
rz(-1.9890824) q[2];
sx q[2];
rz(2.929941) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.043916941) q[1];
sx q[1];
rz(-1.3070135) q[1];
sx q[1];
rz(1.5013807) q[1];
rz(0.84799453) q[3];
sx q[3];
rz(-1.902033) q[3];
sx q[3];
rz(3.1183463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6935912) q[2];
sx q[2];
rz(-0.75746626) q[2];
sx q[2];
rz(2.9288911) q[2];
rz(-1.5851783) q[3];
sx q[3];
rz(-2.3978265) q[3];
sx q[3];
rz(-1.0630382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14146516) q[0];
sx q[0];
rz(-2.8746222) q[0];
sx q[0];
rz(2.2947327) q[0];
rz(2.8894539) q[1];
sx q[1];
rz(-0.82770258) q[1];
sx q[1];
rz(-2.491378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0668117) q[0];
sx q[0];
rz(-1.909314) q[0];
sx q[0];
rz(-1.5241429) q[0];
rz(-2.3450801) q[2];
sx q[2];
rz(-2.079051) q[2];
sx q[2];
rz(0.69895335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47874988) q[1];
sx q[1];
rz(-1.2744941) q[1];
sx q[1];
rz(-3.0257312) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3350983) q[3];
sx q[3];
rz(-1.8915081) q[3];
sx q[3];
rz(0.93076462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3113159) q[2];
sx q[2];
rz(-2.4429784) q[2];
sx q[2];
rz(-2.1777731) q[2];
rz(1.3612932) q[3];
sx q[3];
rz(-2.6908974) q[3];
sx q[3];
rz(-3.0986339) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90218246) q[0];
sx q[0];
rz(-2.9189411) q[0];
sx q[0];
rz(-2.1233001) q[0];
rz(-0.88515431) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(-0.28908602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1219629) q[0];
sx q[0];
rz(-1.8719073) q[0];
sx q[0];
rz(2.798978) q[0];
rz(0.54479213) q[2];
sx q[2];
rz(-0.71342403) q[2];
sx q[2];
rz(-2.1868844) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70967453) q[1];
sx q[1];
rz(-0.54952114) q[1];
sx q[1];
rz(-1.2195682) q[1];
rz(-pi) q[2];
rz(0.52288309) q[3];
sx q[3];
rz(-0.25141806) q[3];
sx q[3];
rz(-2.4222684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79124147) q[2];
sx q[2];
rz(-1.505932) q[2];
sx q[2];
rz(1.852847) q[2];
rz(-0.57981235) q[3];
sx q[3];
rz(-0.47477397) q[3];
sx q[3];
rz(2.536186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8346005) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(0.32421625) q[0];
rz(0.038837198) q[1];
sx q[1];
rz(-0.7380929) q[1];
sx q[1];
rz(1.0968346) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841717) q[0];
sx q[0];
rz(-1.3270507) q[0];
sx q[0];
rz(-0.28864606) q[0];
rz(-1.4154424) q[2];
sx q[2];
rz(-0.58621472) q[2];
sx q[2];
rz(-1.4582576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8406163) q[1];
sx q[1];
rz(-0.30472091) q[1];
sx q[1];
rz(-1.4304377) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1504668) q[3];
sx q[3];
rz(-1.4864576) q[3];
sx q[3];
rz(1.727551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0754806) q[2];
sx q[2];
rz(-2.8754063) q[2];
sx q[2];
rz(-0.038662635) q[2];
rz(-1.4085116) q[3];
sx q[3];
rz(-1.6950636) q[3];
sx q[3];
rz(-2.8601638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51327813) q[0];
sx q[0];
rz(-0.81858855) q[0];
sx q[0];
rz(-1.0191089) q[0];
rz(-2.6844773) q[1];
sx q[1];
rz(-2.8725084) q[1];
sx q[1];
rz(1.7105182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5605049) q[0];
sx q[0];
rz(-1.3700587) q[0];
sx q[0];
rz(1.597658) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9408631) q[2];
sx q[2];
rz(-0.6218172) q[2];
sx q[2];
rz(1.3696826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0455605) q[1];
sx q[1];
rz(-1.9304196) q[1];
sx q[1];
rz(3.0477307) q[1];
x q[2];
rz(1.1248427) q[3];
sx q[3];
rz(-0.77174458) q[3];
sx q[3];
rz(1.2260557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21923253) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(-2.5617981) q[2];
rz(-2.842105) q[3];
sx q[3];
rz(-0.53286415) q[3];
sx q[3];
rz(-1.2286435) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806468) q[0];
sx q[0];
rz(-1.6609284) q[0];
sx q[0];
rz(-3.0378367) q[0];
rz(0.62458986) q[1];
sx q[1];
rz(-0.92085212) q[1];
sx q[1];
rz(1.3821028) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73653664) q[0];
sx q[0];
rz(-2.2782405) q[0];
sx q[0];
rz(0.28178431) q[0];
x q[1];
rz(0.34587282) q[2];
sx q[2];
rz(-2.8240339) q[2];
sx q[2];
rz(-1.4994743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8727726) q[1];
sx q[1];
rz(-1.3026857) q[1];
sx q[1];
rz(1.0467749) q[1];
x q[2];
rz(-1.8554888) q[3];
sx q[3];
rz(-2.2119388) q[3];
sx q[3];
rz(-2.5452328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3524807) q[2];
sx q[2];
rz(-1.4752957) q[2];
sx q[2];
rz(-0.86461198) q[2];
rz(-2.9027446) q[3];
sx q[3];
rz(-2.3819203) q[3];
sx q[3];
rz(0.71340942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9637941) q[0];
sx q[0];
rz(-2.5803784) q[0];
sx q[0];
rz(-2.906565) q[0];
rz(0.66559732) q[1];
sx q[1];
rz(-2.0033629) q[1];
sx q[1];
rz(-1.6682909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6971254) q[0];
sx q[0];
rz(-1.7791064) q[0];
sx q[0];
rz(2.134077) q[0];
rz(-pi) q[1];
rz(1.4334045) q[2];
sx q[2];
rz(-1.1072888) q[2];
sx q[2];
rz(2.3799294) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.044933407) q[1];
sx q[1];
rz(-2.9439264) q[1];
sx q[1];
rz(-0.37943073) q[1];
rz(-0.45110945) q[3];
sx q[3];
rz(-1.530297) q[3];
sx q[3];
rz(-2.9429352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16234806) q[2];
sx q[2];
rz(-0.20077106) q[2];
sx q[2];
rz(-2.6806504) q[2];
rz(-2.0567242) q[3];
sx q[3];
rz(-0.74551398) q[3];
sx q[3];
rz(-0.78456867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0775065) q[0];
sx q[0];
rz(-2.6106847) q[0];
sx q[0];
rz(-0.60894668) q[0];
rz(2.5755836) q[1];
sx q[1];
rz(-1.9809664) q[1];
sx q[1];
rz(1.9401248) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64513759) q[0];
sx q[0];
rz(-1.4352531) q[0];
sx q[0];
rz(3.1333424) q[0];
rz(1.1913774) q[2];
sx q[2];
rz(-0.51926368) q[2];
sx q[2];
rz(-1.0644703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2347327) q[1];
sx q[1];
rz(-2.3667049) q[1];
sx q[1];
rz(-3.0084684) q[1];
rz(-pi) q[2];
rz(1.8747599) q[3];
sx q[3];
rz(-1.1397151) q[3];
sx q[3];
rz(-1.4747185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10298771) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(-1.9496244) q[2];
rz(-2.1206756) q[3];
sx q[3];
rz(-2.9396785) q[3];
sx q[3];
rz(-1.5405704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182619) q[0];
sx q[0];
rz(-2.9409565) q[0];
sx q[0];
rz(0.9675135) q[0];
rz(-1.4406904) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(-0.38223019) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8865693) q[0];
sx q[0];
rz(-0.60211997) q[0];
sx q[0];
rz(1.2467136) q[0];
x q[1];
rz(0.43469825) q[2];
sx q[2];
rz(-1.9524196) q[2];
sx q[2];
rz(1.8392854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6401419) q[1];
sx q[1];
rz(-0.79265187) q[1];
sx q[1];
rz(0.16736302) q[1];
rz(-pi) q[2];
rz(2.5955673) q[3];
sx q[3];
rz(-1.5848397) q[3];
sx q[3];
rz(1.7469116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.050015673) q[2];
sx q[2];
rz(-2.2735368) q[2];
sx q[2];
rz(2.3559605) q[2];
rz(3.1146289) q[3];
sx q[3];
rz(-1.5188768) q[3];
sx q[3];
rz(0.54076076) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3383269) q[0];
sx q[0];
rz(-1.6852408) q[0];
sx q[0];
rz(2.2483873) q[0];
rz(2.4772353) q[1];
sx q[1];
rz(-1.2394445) q[1];
sx q[1];
rz(2.209421) q[1];
rz(2.4128466) q[2];
sx q[2];
rz(-1.797429) q[2];
sx q[2];
rz(0.28354473) q[2];
rz(0.6046927) q[3];
sx q[3];
rz(-1.562955) q[3];
sx q[3];
rz(2.4969586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
