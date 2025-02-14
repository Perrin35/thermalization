OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(-2.2890685) q[0];
sx q[0];
rz(2.8791715) q[0];
rz(-3.4497058) q[1];
sx q[1];
rz(0.83642712) q[1];
sx q[1];
rz(15.717419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0266761) q[0];
sx q[0];
rz(-2.466188) q[0];
sx q[0];
rz(0.83362718) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0105825) q[2];
sx q[2];
rz(-0.081801266) q[2];
sx q[2];
rz(2.8103925) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7345924) q[1];
sx q[1];
rz(-1.0838036) q[1];
sx q[1];
rz(-0.89501801) q[1];
rz(-pi) q[2];
rz(2.6515048) q[3];
sx q[3];
rz(-2.2667829) q[3];
sx q[3];
rz(2.0139607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0414163) q[2];
sx q[2];
rz(-2.2569423) q[2];
sx q[2];
rz(0.46119383) q[2];
rz(2.6395116) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(-1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13040386) q[0];
sx q[0];
rz(-1.4606425) q[0];
sx q[0];
rz(-2.9050264) q[0];
rz(-0.81831167) q[1];
sx q[1];
rz(-1.2032443) q[1];
sx q[1];
rz(2.1990105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078694669) q[0];
sx q[0];
rz(-1.601614) q[0];
sx q[0];
rz(0.0027058733) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6918399) q[2];
sx q[2];
rz(-1.5038052) q[2];
sx q[2];
rz(0.070320694) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89211226) q[1];
sx q[1];
rz(-1.5693519) q[1];
sx q[1];
rz(-0.73489983) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97901042) q[3];
sx q[3];
rz(-0.56603449) q[3];
sx q[3];
rz(0.69379679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4054823) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(2.45641) q[2];
rz(0.80125609) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41246978) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(0.98010081) q[0];
rz(-3.0209172) q[1];
sx q[1];
rz(-1.9735034) q[1];
sx q[1];
rz(-1.4792222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2421766) q[0];
sx q[0];
rz(-0.78095312) q[0];
sx q[0];
rz(-2.910898) q[0];
x q[1];
rz(-0.73320575) q[2];
sx q[2];
rz(-0.63989675) q[2];
sx q[2];
rz(2.899183) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.66202098) q[1];
sx q[1];
rz(-1.5618854) q[1];
sx q[1];
rz(-1.575995) q[1];
rz(-2.5517518) q[3];
sx q[3];
rz(-1.2882659) q[3];
sx q[3];
rz(0.88385201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4270619) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(2.96116) q[2];
rz(0.42029542) q[3];
sx q[3];
rz(-0.97729483) q[3];
sx q[3];
rz(2.596627) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20035289) q[0];
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
rz(1.1526718) q[0];
sx q[0];
rz(-1.3878569) q[0];
sx q[0];
rz(-1.4519917) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21993262) q[2];
sx q[2];
rz(-0.77785197) q[2];
sx q[2];
rz(0.57524112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56616773) q[1];
sx q[1];
rz(-1.638104) q[1];
sx q[1];
rz(2.7507902) q[1];
rz(-pi) q[2];
rz(-2.4957711) q[3];
sx q[3];
rz(-0.97664112) q[3];
sx q[3];
rz(2.770383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0338318) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(-2.4793009) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0419643) q[0];
sx q[0];
rz(-0.30655107) q[0];
sx q[0];
rz(-2.1339259) q[0];
rz(-3.0362466) q[1];
sx q[1];
rz(-0.65709972) q[1];
sx q[1];
rz(-0.23060051) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5916633) q[0];
sx q[0];
rz(-1.2386892) q[0];
sx q[0];
rz(1.1989393) q[0];
rz(-pi) q[1];
rz(-2.9020374) q[2];
sx q[2];
rz(-1.2192246) q[2];
sx q[2];
rz(-0.4229751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9437792) q[1];
sx q[1];
rz(-2.9557485) q[1];
sx q[1];
rz(2.6617104) q[1];
rz(-pi) q[2];
x q[2];
rz(2.496594) q[3];
sx q[3];
rz(-1.4497744) q[3];
sx q[3];
rz(-1.6337486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0577804) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(2.9528565) q[2];
rz(-1.2356637) q[3];
sx q[3];
rz(-2.6344447) q[3];
sx q[3];
rz(1.88131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96596232) q[0];
sx q[0];
rz(-1.2164793) q[0];
sx q[0];
rz(-0.29773444) q[0];
rz(-0.072602428) q[1];
sx q[1];
rz(-2.4563792) q[1];
sx q[1];
rz(0.6822449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6202691) q[0];
sx q[0];
rz(-2.1730039) q[0];
sx q[0];
rz(-0.35767718) q[0];
x q[1];
rz(2.9932026) q[2];
sx q[2];
rz(-2.1990657) q[2];
sx q[2];
rz(-2.4895957) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1156716) q[1];
sx q[1];
rz(-0.37927765) q[1];
sx q[1];
rz(-1.161452) q[1];
rz(-pi) q[2];
rz(-2.7704084) q[3];
sx q[3];
rz(-2.1371578) q[3];
sx q[3];
rz(2.0831828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5785152) q[2];
sx q[2];
rz(-2.7882521) q[2];
sx q[2];
rz(0.087470857) q[2];
rz(-2.2443917) q[3];
sx q[3];
rz(-1.2465435) q[3];
sx q[3];
rz(0.38952601) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099365756) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(-1.1627831) q[0];
rz(2.9258264) q[1];
sx q[1];
rz(-0.81753221) q[1];
sx q[1];
rz(-0.93528265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29998818) q[0];
sx q[0];
rz(-2.3242352) q[0];
sx q[0];
rz(-1.9005083) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14752702) q[2];
sx q[2];
rz(-1.8263683) q[2];
sx q[2];
rz(0.35596656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7120119) q[1];
sx q[1];
rz(-1.4918527) q[1];
sx q[1];
rz(-0.78671766) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31511001) q[3];
sx q[3];
rz(-1.2602046) q[3];
sx q[3];
rz(-0.14652182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0995348) q[2];
sx q[2];
rz(-0.65379405) q[2];
sx q[2];
rz(2.330244) q[2];
rz(-0.30284303) q[3];
sx q[3];
rz(-1.9767913) q[3];
sx q[3];
rz(0.6306878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.621156) q[0];
sx q[0];
rz(-0.32670894) q[0];
sx q[0];
rz(-0.69123554) q[0];
rz(2.4825545) q[1];
sx q[1];
rz(-1.6233416) q[1];
sx q[1];
rz(-0.59115994) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189445) q[0];
sx q[0];
rz(-1.582358) q[0];
sx q[0];
rz(-1.8537136) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.065584) q[2];
sx q[2];
rz(-0.55375868) q[2];
sx q[2];
rz(2.8477235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2745251) q[1];
sx q[1];
rz(-2.3127416) q[1];
sx q[1];
rz(-2.0193923) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0847315) q[3];
sx q[3];
rz(-1.2788075) q[3];
sx q[3];
rz(0.1869456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14472321) q[2];
sx q[2];
rz(-2.150712) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(-0.75374323) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(-2.6678705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.4562456) q[0];
sx q[0];
rz(-1.101475) q[0];
sx q[0];
rz(2.9154678) q[0];
rz(-1.4489168) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(1.0015782) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3069835) q[0];
sx q[0];
rz(-2.153647) q[0];
sx q[0];
rz(1.7442763) q[0];
rz(-1.1923157) q[2];
sx q[2];
rz(-2.2342302) q[2];
sx q[2];
rz(2.0326234) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4521617) q[1];
sx q[1];
rz(-0.2161444) q[1];
sx q[1];
rz(-0.5139022) q[1];
rz(-pi) q[2];
rz(1.8705192) q[3];
sx q[3];
rz(-2.2329438) q[3];
sx q[3];
rz(-1.3817092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5292624) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100819) q[0];
sx q[0];
rz(-1.2393351) q[0];
sx q[0];
rz(0.70247689) q[0];
rz(-2.6055873) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(1.747267) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4634906) q[0];
sx q[0];
rz(-2.9759088) q[0];
sx q[0];
rz(1.7278306) q[0];
rz(-pi) q[1];
rz(1.461566) q[2];
sx q[2];
rz(-1.6642206) q[2];
sx q[2];
rz(-1.0051073) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9623757) q[1];
sx q[1];
rz(-2.0834647) q[1];
sx q[1];
rz(2.8151399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0004041632) q[3];
sx q[3];
rz(-2.7260216) q[3];
sx q[3];
rz(-3.0311751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1348306) q[2];
sx q[2];
rz(-0.68209058) q[2];
sx q[2];
rz(-1.6949867) q[2];
rz(-2.8805978) q[3];
sx q[3];
rz(-1.010453) q[3];
sx q[3];
rz(-1.9058156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9217054) q[0];
sx q[0];
rz(-0.6548665) q[0];
sx q[0];
rz(-1.0153216) q[0];
rz(1.075853) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(1.6816611) q[2];
sx q[2];
rz(-1.9456429) q[2];
sx q[2];
rz(2.6425998) q[2];
rz(2.3786529) q[3];
sx q[3];
rz(-1.8992103) q[3];
sx q[3];
rz(-2.7873399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
