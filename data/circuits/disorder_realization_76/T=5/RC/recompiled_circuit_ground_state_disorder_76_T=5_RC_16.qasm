OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(-1.0816242) q[0];
sx q[0];
rz(2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(1.8269202) q[1];
sx q[1];
rz(4.8575525) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8872966) q[0];
sx q[0];
rz(-2.097541) q[0];
sx q[0];
rz(-2.8077784) q[0];
x q[1];
rz(-0.47503586) q[2];
sx q[2];
rz(-2.5940707) q[2];
sx q[2];
rz(2.5753367) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4233746) q[1];
sx q[1];
rz(-1.4932695) q[1];
sx q[1];
rz(1.7192057) q[1];
x q[2];
rz(-1.6207148) q[3];
sx q[3];
rz(-2.422159) q[3];
sx q[3];
rz(-1.0843383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90929675) q[2];
sx q[2];
rz(-2.2984419) q[2];
sx q[2];
rz(-2.9006531) q[2];
rz(-0.029189261) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(-1.9624814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643243) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(-0.48467317) q[0];
rz(2.4618705) q[1];
sx q[1];
rz(-1.2815963) q[1];
sx q[1];
rz(-1.9814804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3368615) q[0];
sx q[0];
rz(-1.2668249) q[0];
sx q[0];
rz(1.5514403) q[0];
rz(-pi) q[1];
rz(2.4258575) q[2];
sx q[2];
rz(-1.01075) q[2];
sx q[2];
rz(-0.83362388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1966432) q[1];
sx q[1];
rz(-1.8183876) q[1];
sx q[1];
rz(2.5291614) q[1];
x q[2];
rz(2.6187702) q[3];
sx q[3];
rz(-1.2827875) q[3];
sx q[3];
rz(0.097974591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.844187) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(0.45787946) q[2];
rz(1.1229905) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(-0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8726525) q[0];
sx q[0];
rz(-0.70916969) q[0];
sx q[0];
rz(3.1100682) q[0];
rz(2.8541376) q[1];
sx q[1];
rz(-2.2645686) q[1];
sx q[1];
rz(1.9256176) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368204) q[0];
sx q[0];
rz(-0.88951123) q[0];
sx q[0];
rz(2.4347406) q[0];
rz(-pi) q[1];
rz(-0.38162614) q[2];
sx q[2];
rz(-0.51593057) q[2];
sx q[2];
rz(0.53781539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.029473631) q[1];
sx q[1];
rz(-0.39327565) q[1];
sx q[1];
rz(-2.359819) q[1];
x q[2];
rz(-0.9312882) q[3];
sx q[3];
rz(-1.0126601) q[3];
sx q[3];
rz(2.1344413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.03881255) q[2];
sx q[2];
rz(-1.5993885) q[2];
sx q[2];
rz(-2.3056324) q[2];
rz(1.7838259) q[3];
sx q[3];
rz(-0.72503763) q[3];
sx q[3];
rz(2.3939705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2826071) q[0];
sx q[0];
rz(-1.405412) q[0];
sx q[0];
rz(0.080168515) q[0];
rz(-1.0824341) q[1];
sx q[1];
rz(-0.23324649) q[1];
sx q[1];
rz(-1.3287883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6850615) q[0];
sx q[0];
rz(-2.079112) q[0];
sx q[0];
rz(1.1969181) q[0];
rz(-pi) q[1];
rz(-1.7432418) q[2];
sx q[2];
rz(-1.2524464) q[2];
sx q[2];
rz(2.6686252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7419974) q[1];
sx q[1];
rz(-0.4542225) q[1];
sx q[1];
rz(1.9059559) q[1];
rz(2.8058047) q[3];
sx q[3];
rz(-1.3779252) q[3];
sx q[3];
rz(-1.4953116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3749915) q[2];
sx q[2];
rz(-0.98321715) q[2];
sx q[2];
rz(-2.2341466) q[2];
rz(2.4713016) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(1.420174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(0.43148828) q[0];
rz(1.1314499) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(2.7010837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9517072) q[0];
sx q[0];
rz(-1.110297) q[0];
sx q[0];
rz(3.0104464) q[0];
x q[1];
rz(0.20655234) q[2];
sx q[2];
rz(-1.163854) q[2];
sx q[2];
rz(-0.47688866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9677093) q[1];
sx q[1];
rz(-1.5313799) q[1];
sx q[1];
rz(2.8110678) q[1];
rz(-pi) q[2];
rz(0.83310762) q[3];
sx q[3];
rz(-2.3596977) q[3];
sx q[3];
rz(1.1469008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12209192) q[2];
sx q[2];
rz(-0.92635218) q[2];
sx q[2];
rz(0.41395536) q[2];
rz(1.3881989) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(-3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269161) q[0];
sx q[0];
rz(-1.7358945) q[0];
sx q[0];
rz(0.93389121) q[0];
rz(2.4712708) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(1.3353039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437296) q[0];
sx q[0];
rz(-1.6707194) q[0];
sx q[0];
rz(2.1735682) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15768965) q[2];
sx q[2];
rz(-1.5458428) q[2];
sx q[2];
rz(-1.7316929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.851593) q[1];
sx q[1];
rz(-1.6834752) q[1];
sx q[1];
rz(-0.022701724) q[1];
rz(-pi) q[2];
rz(-0.63702668) q[3];
sx q[3];
rz(-1.5638132) q[3];
sx q[3];
rz(-2.1050326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2495217) q[2];
sx q[2];
rz(-0.27270174) q[2];
sx q[2];
rz(1.8708694) q[2];
rz(-1.3696085) q[3];
sx q[3];
rz(-1.2415875) q[3];
sx q[3];
rz(0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18913604) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(-0.15368803) q[0];
rz(-2.9391089) q[1];
sx q[1];
rz(-1.2362365) q[1];
sx q[1];
rz(-0.94857803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.089854) q[0];
sx q[0];
rz(-1.7719381) q[0];
sx q[0];
rz(-2.5178943) q[0];
rz(-pi) q[1];
rz(0.93164753) q[2];
sx q[2];
rz(-1.8338565) q[2];
sx q[2];
rz(1.6073038) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0041629) q[1];
sx q[1];
rz(-1.8297126) q[1];
sx q[1];
rz(0.24000411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7465215) q[3];
sx q[3];
rz(-0.65726377) q[3];
sx q[3];
rz(-2.6906309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.004868) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(0.0014121545) q[2];
rz(0.062156113) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(-2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75981265) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(1.2148452) q[0];
rz(0.75792056) q[1];
sx q[1];
rz(-2.6140723) q[1];
sx q[1];
rz(-0.55799276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4014977) q[0];
sx q[0];
rz(-1.3245653) q[0];
sx q[0];
rz(1.2241062) q[0];
x q[1];
rz(0.2336913) q[2];
sx q[2];
rz(-2.5776641) q[2];
sx q[2];
rz(-2.301931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87822014) q[1];
sx q[1];
rz(-2.0505295) q[1];
sx q[1];
rz(0.31444506) q[1];
rz(-pi) q[2];
rz(-1.4766944) q[3];
sx q[3];
rz(-0.27827874) q[3];
sx q[3];
rz(-0.32839963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(-0.26926678) q[2];
rz(-1.1632129) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729386) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(-2.0682251) q[0];
rz(2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(0.55666298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9087785) q[0];
sx q[0];
rz(-1.7126084) q[0];
sx q[0];
rz(1.2491666) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5684909) q[2];
sx q[2];
rz(-1.0759883) q[2];
sx q[2];
rz(3.0798517) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5334415) q[1];
sx q[1];
rz(-1.2131629) q[1];
sx q[1];
rz(0.74700345) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7438874) q[3];
sx q[3];
rz(-2.8880062) q[3];
sx q[3];
rz(1.2802779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90159455) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(1.1617917) q[2];
rz(-0.53317541) q[3];
sx q[3];
rz(-2.6170001) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6530782) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(2.5467806) q[0];
rz(2.5626903) q[1];
sx q[1];
rz(-1.6659104) q[1];
sx q[1];
rz(-2.4822809) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80583094) q[0];
sx q[0];
rz(-1.7044401) q[0];
sx q[0];
rz(0.62778309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5388902) q[2];
sx q[2];
rz(-1.9424159) q[2];
sx q[2];
rz(1.1204669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51344556) q[1];
sx q[1];
rz(-1.4288578) q[1];
sx q[1];
rz(-1.8177086) q[1];
rz(-pi) q[2];
rz(1.8439383) q[3];
sx q[3];
rz(-1.8255263) q[3];
sx q[3];
rz(1.8984853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.79002964) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(-0.053038049) q[2];
rz(-1.3430345) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(0.0742577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2863083) q[0];
sx q[0];
rz(-1.3636148) q[0];
sx q[0];
rz(-1.6028945) q[0];
rz(3.042649) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(1.5621875) q[2];
sx q[2];
rz(-1.1653524) q[2];
sx q[2];
rz(-0.63962519) q[2];
rz(-2.9704499) q[3];
sx q[3];
rz(-2.3928693) q[3];
sx q[3];
rz(0.00015043845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
