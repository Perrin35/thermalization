OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.054501) q[0];
sx q[0];
rz(-0.85252419) q[0];
sx q[0];
rz(-2.8791715) q[0];
rz(2.8334795) q[1];
sx q[1];
rz(-2.3051655) q[1];
sx q[1];
rz(-0.0094553789) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1605837) q[0];
sx q[0];
rz(-2.0045337) q[0];
sx q[0];
rz(-2.1061312) q[0];
x q[1];
rz(-2.1310102) q[2];
sx q[2];
rz(-3.0597914) q[2];
sx q[2];
rz(2.8103925) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7345924) q[1];
sx q[1];
rz(-1.0838036) q[1];
sx q[1];
rz(0.89501801) q[1];
x q[2];
rz(2.3289324) q[3];
sx q[3];
rz(-1.9403807) q[3];
sx q[3];
rz(-3.0280035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0414163) q[2];
sx q[2];
rz(-2.2569423) q[2];
sx q[2];
rz(2.6803988) q[2];
rz(2.6395116) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(-1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0111888) q[0];
sx q[0];
rz(-1.4606425) q[0];
sx q[0];
rz(-2.9050264) q[0];
rz(2.323281) q[1];
sx q[1];
rz(-1.9383483) q[1];
sx q[1];
rz(0.94258211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492185) q[0];
sx q[0];
rz(-1.5680917) q[0];
sx q[0];
rz(1.5399786) q[0];
rz(-pi) q[1];
rz(-1.4497527) q[2];
sx q[2];
rz(-1.5038052) q[2];
sx q[2];
rz(3.071272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4613104) q[1];
sx q[1];
rz(-2.4066917) q[1];
sx q[1];
rz(-0.0021541455) q[1];
rz(-2.0560741) q[3];
sx q[3];
rz(-1.266978) q[3];
sx q[3];
rz(-1.7485545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4054823) q[2];
sx q[2];
rz(-2.5030899) q[2];
sx q[2];
rz(-0.68518266) q[2];
rz(2.3403366) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(1.2509468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291229) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(-0.98010081) q[0];
rz(0.12067548) q[1];
sx q[1];
rz(-1.9735034) q[1];
sx q[1];
rz(-1.4792222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50608327) q[0];
sx q[0];
rz(-1.4091307) q[0];
sx q[0];
rz(-0.76753214) q[0];
x q[1];
rz(2.4083869) q[2];
sx q[2];
rz(-2.5016959) q[2];
sx q[2];
rz(0.24240968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66202098) q[1];
sx q[1];
rz(-1.5618854) q[1];
sx q[1];
rz(-1.5655976) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48101403) q[3];
sx q[3];
rz(-2.4949007) q[3];
sx q[3];
rz(2.8493136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7145308) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(2.96116) q[2];
rz(-2.7212972) q[3];
sx q[3];
rz(-2.1642978) q[3];
sx q[3];
rz(0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412398) q[0];
sx q[0];
rz(-1.7914597) q[0];
sx q[0];
rz(-3.0923162) q[0];
rz(-0.73206466) q[1];
sx q[1];
rz(-2.2923636) q[1];
sx q[1];
rz(-1.3759618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5696053) q[0];
sx q[0];
rz(-0.21776918) q[0];
sx q[0];
rz(-0.56977429) q[0];
rz(1.3591197) q[2];
sx q[2];
rz(-2.3251901) q[2];
sx q[2];
rz(-0.27118452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56616773) q[1];
sx q[1];
rz(-1.5034887) q[1];
sx q[1];
rz(-2.7507902) q[1];
rz(-2.2729418) q[3];
sx q[3];
rz(-1.0486866) q[3];
sx q[3];
rz(-2.3412395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0338318) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(2.4793009) q[2];
rz(-1.8108588) q[3];
sx q[3];
rz(-0.8225421) q[3];
sx q[3];
rz(1.2983769) q[3];
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
rz(pi/2) q[0];
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
rz(-1.0419643) q[0];
sx q[0];
rz(-0.30655107) q[0];
sx q[0];
rz(-2.1339259) q[0];
rz(3.0362466) q[1];
sx q[1];
rz(-0.65709972) q[1];
sx q[1];
rz(0.23060051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5499293) q[0];
sx q[0];
rz(-1.2386892) q[0];
sx q[0];
rz(-1.1989393) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9318387) q[2];
sx q[2];
rz(-1.3461543) q[2];
sx q[2];
rz(-1.063907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8457875) q[1];
sx q[1];
rz(-1.485386) q[1];
sx q[1];
rz(2.9763362) q[1];
rz(-pi) q[2];
rz(-2.9419961) q[3];
sx q[3];
rz(-0.65465876) q[3];
sx q[3];
rz(-2.9195291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0838123) q[2];
sx q[2];
rz(-2.1746706) q[2];
sx q[2];
rz(2.9528565) q[2];
rz(1.2356637) q[3];
sx q[3];
rz(-0.50714791) q[3];
sx q[3];
rz(1.88131) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1756303) q[0];
sx q[0];
rz(-1.2164793) q[0];
sx q[0];
rz(-2.8438582) q[0];
rz(-0.072602428) q[1];
sx q[1];
rz(-2.4563792) q[1];
sx q[1];
rz(-2.4593478) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52132359) q[0];
sx q[0];
rz(-2.1730039) q[0];
sx q[0];
rz(0.35767718) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93726341) q[2];
sx q[2];
rz(-1.450895) q[2];
sx q[2];
rz(-2.1351569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9798237) q[1];
sx q[1];
rz(-1.4228954) q[1];
sx q[1];
rz(-1.921341) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97203927) q[3];
sx q[3];
rz(-1.2597196) q[3];
sx q[3];
rz(-2.4233203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5785152) q[2];
sx q[2];
rz(-2.7882521) q[2];
sx q[2];
rz(-3.0541218) q[2];
rz(-2.2443917) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099365756) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(-1.9788096) q[0];
rz(0.21576628) q[1];
sx q[1];
rz(-0.81753221) q[1];
sx q[1];
rz(0.93528265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777841) q[0];
sx q[0];
rz(-2.3323614) q[0];
sx q[0];
rz(-2.8092395) q[0];
x q[1];
rz(0.14752702) q[2];
sx q[2];
rz(-1.8263683) q[2];
sx q[2];
rz(-2.7856261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21970785) q[1];
sx q[1];
rz(-0.7898191) q[1];
sx q[1];
rz(-3.0303257) q[1];
rz(2.8264826) q[3];
sx q[3];
rz(-1.8813881) q[3];
sx q[3];
rz(2.9950708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0420578) q[2];
sx q[2];
rz(-0.65379405) q[2];
sx q[2];
rz(-2.330244) q[2];
rz(2.8387496) q[3];
sx q[3];
rz(-1.1648014) q[3];
sx q[3];
rz(2.5109049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52043668) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(2.4503571) q[0];
rz(0.65903819) q[1];
sx q[1];
rz(-1.6233416) q[1];
sx q[1];
rz(0.59115994) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49159971) q[0];
sx q[0];
rz(-0.2831471) q[0];
sx q[0];
rz(-1.6121907) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0760087) q[2];
sx q[2];
rz(-2.587834) q[2];
sx q[2];
rz(-0.29386917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8670676) q[1];
sx q[1];
rz(-0.82885107) q[1];
sx q[1];
rz(-2.0193923) q[1];
x q[2];
rz(-1.2783613) q[3];
sx q[3];
rz(-1.5163444) q[3];
sx q[3];
rz(1.4002348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14472321) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(0.75374323) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(-0.47372216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562456) q[0];
sx q[0];
rz(-2.0401177) q[0];
sx q[0];
rz(0.22612485) q[0];
rz(-1.6926758) q[1];
sx q[1];
rz(-2.1782404) q[1];
sx q[1];
rz(-2.1400145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16765633) q[0];
sx q[0];
rz(-1.4261803) q[0];
sx q[0];
rz(-0.58982106) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6993352) q[2];
sx q[2];
rz(-1.8662037) q[2];
sx q[2];
rz(0.70198529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16538936) q[1];
sx q[1];
rz(-1.7586629) q[1];
sx q[1];
rz(-1.4632744) q[1];
rz(-pi) q[2];
rz(-1.8705192) q[3];
sx q[3];
rz(-0.90864881) q[3];
sx q[3];
rz(-1.3817092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61233026) q[2];
sx q[2];
rz(-2.3713106) q[2];
sx q[2];
rz(0.15753499) q[2];
rz(-2.9869288) q[3];
sx q[3];
rz(-1.9779343) q[3];
sx q[3];
rz(-2.7111588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(2.4100819) q[0];
sx q[0];
rz(-1.9022576) q[0];
sx q[0];
rz(2.4391158) q[0];
rz(-0.53600535) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(-1.747267) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6781021) q[0];
sx q[0];
rz(-2.9759088) q[0];
sx q[0];
rz(-1.413762) q[0];
x q[1];
rz(-1.461566) q[2];
sx q[2];
rz(-1.4773721) q[2];
sx q[2];
rz(2.1364853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5665219) q[1];
sx q[1];
rz(-2.5417323) q[1];
sx q[1];
rz(2.0886648) q[1];
x q[2];
rz(-1.5709747) q[3];
sx q[3];
rz(-1.1552253) q[3];
sx q[3];
rz(-3.0316169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(1.6949867) q[2];
rz(0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(-1.235777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21988729) q[0];
sx q[0];
rz(-2.4867262) q[0];
sx q[0];
rz(2.126271) q[0];
rz(-2.0657397) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(-0.37694945) q[2];
sx q[2];
rz(-1.6739346) q[2];
sx q[2];
rz(-2.0290537) q[2];
rz(2.6834834) q[3];
sx q[3];
rz(-2.3242713) q[3];
sx q[3];
rz(-0.8911688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
